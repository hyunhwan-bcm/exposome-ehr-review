"""Tests for data-availability schema fields + focused LLM scan."""
import json
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from summarizer.schema import ManuscriptChecklist


# ── schema: new fields ───────────────────────────────────────────────────────


def _rec(**over):
    base = {
        "pmcid": "PMC1", "title": "t", "year": "2020",
        "ehr_used": False, "ehr_evidence": "n/a", "summary": "s",
    }
    base.update(over)
    return base


def test_data_availability_defaults_to_not_stated():
    c = ManuscriptChecklist(**_rec())
    assert c.data_availability == "not-stated"
    assert c.data_accession_links == []
    assert c.data_availability_statement == ""


def test_data_availability_enum_accepted():
    c = ManuscriptChecklist(**_rec(data_availability="public-repository",
                                   data_accession_links=["dbGaP phs000123", "https://github.com/x/y"]))
    assert c.data_availability == "public-repository"
    assert len(c.data_accession_links) == 2


def test_data_availability_normalizes_spaces_and_underscores():
    c = ManuscriptChecklist(**_rec(data_availability="available upon request"))
    assert c.data_availability == "available-upon-request"
    c2 = ManuscriptChecklist(**_rec(data_availability="public_repository"))
    assert c2.data_availability == "public-repository"


def test_data_availability_unknown_falls_to_not_stated():
    c = ManuscriptChecklist(**_rec(data_availability="maybe"))
    assert c.data_availability == "not-stated"


def test_accession_links_coerces_single_string():
    c = ManuscriptChecklist(**_rec(data_accession_links="GSE12345"))
    assert c.data_accession_links == ["GSE12345"]


def test_accession_links_drops_empties():
    c = ManuscriptChecklist(**_rec(data_accession_links=["GSE12345", "", "  "]))
    assert c.data_accession_links == ["GSE12345"]


def test_existing_summaries_without_field_still_validate():
    """Old per-paper JSONs (no data_availability) must still load after the schema change."""
    old = _rec()
    c = ManuscriptChecklist.model_validate(old)
    assert c.data_availability == "not-stated"


# ── scan script: updates per-paper JSON preserving existing fields ───────────


def test_scan_preserves_existing_fields(tmp_path: Path):
    """scan_one should merge the new data-availability fields onto an existing record."""
    import scan_data_availability as sda

    # existing per-paper JSON with full fields
    existing = _rec(pmcid="PMC42", title="Real Title", ehr_used=True,
                    ehr_evidence="HES codes", summary="real summary",
                    pathologies_diseases=["asthma"])
    out = tmp_path / "PMC42.json"
    out.write_text(json.dumps(existing, indent=2))

    # fake full-text file
    txt = tmp_path / "PMC42.txt"; txt.write_text("paper text")

    scan_result = {
        "data_availability": "public-repository",
        "data_accession_links": ["dbGaP phs000789"],
        "data_availability_statement": "Data available in dbGaP.",
    }
    with patch("scan_data_availability.extract", return_value=("text", "text")), \
         patch("scan_data_availability.ask_llm_data_availability", return_value=scan_result):
        result = sda.scan_one(txt, summary_dir=tmp_path)

    assert result["data_availability"] == "public-repository"
    # existing fields preserved
    assert result["title"] == "Real Title"
    assert result["ehr_used"] is True
    assert result["pathologies_diseases"] == ["asthma"]
    # written to disk + valid
    on_disk = json.loads(out.read_text())
    assert on_disk["data_availability"] == "public-repository"
    ManuscriptChecklist.model_validate(on_disk)  # round-trips


def test_scan_handles_no_existing_summary(tmp_path: Path):
    """If no per-paper JSON exists yet, scan_one creates one from metadata + scan."""
    import scan_data_availability as sda

    txt = tmp_path / "PMC99.txt"; txt.write_text("paper text")
    scan_result = {"data_availability": "in-house", "data_accession_links": [],
                   "data_availability_statement": "In-house cohort."}
    with patch("scan_data_availability.extract", return_value=("text", "text")), \
         patch("scan_data_availability.ask_llm_data_availability", return_value=scan_result):
        result = sda.scan_one(txt, summary_dir=tmp_path, meta={"PMC99": {"title": "T", "year": "2021"}})

    assert result["data_availability"] == "in-house"
    assert result["title"] == "T"
    on_disk = json.loads((tmp_path / "PMC99.json").read_text())
    assert on_disk["data_availability"] == "in-house"
    # required analytical fields absent -> written as partial (full summarize completes later)
    assert "ehr_used" not in on_disk


def test_ask_llm_returns_none_on_failure(tmp_path: Path):
    """ask_llm_data_availability returns None (graceful) instead of raising."""
    import scan_data_availability as sda
    mock_client = MagicMock()
    mock_client.chat.completions.create.side_effect = RuntimeError("boom")
    out = sda.ask_llm_data_availability(client=mock_client, model="m", text="x",
                                        pmcid="PMC1", title="t", year="2020")
    assert out is None


def test_ask_llm_recovers_markdown_field_response():
    """Gemma sometimes ignores JSON instructions and emits markdown fields."""
    import scan_data_availability as sda

    raw = """
* `data_availability`: "public-repository"
* `data_accession_links`:
  * ENA: PRJEB86099
  * Zenodo: https://doi.org/10.5281/zenodo.14913431
* `data_availability_statement`:
  * "Sequencing data have been deposited at the European Nucleotide Archive under project ENA: PRJEB86099."
"""
    mock_client = MagicMock()
    mock_client.chat.completions.create.return_value.choices = [
        MagicMock(message=MagicMock(content=raw))
    ]

    out = sda.ask_llm_data_availability(
        client=mock_client,
        model="m",
        text="The paper says data details are described in the availability section.",
        pmcid="PMC1",
        title="t",
        year="2020",
    )

    assert out["data_availability"] == "public-repository"
    assert "PRJEB86099" in out["data_accession_links"]
    assert "https://doi.org/10.5281/zenodo.14913431" in out["data_accession_links"]
    assert "Sequencing data have been deposited" in out["data_availability_statement"]


def test_accession_link_repairs_pdf_line_wrap():
    """A GitHub URL split across a PDF line wrap (hyphen+spaces+newline) is recovered whole."""
    import scan_data_availability as sda
    text = (
        "Data availability The code is publicly available at: "
        "https://github.com/benny817817/early-life-sugar - \n"
        "rationing-respiratory-health . Appendix A."
    )
    links = sda._extract_accession_links(text)
    assert links == ["https://github.com/benny817817/early-life-sugar-rationing-respiratory-health"]
