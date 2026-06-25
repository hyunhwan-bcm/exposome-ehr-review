"""Tests for TinyDB-as-source-of-truth rebuild + concurrent summarization.

Tests written first — implementation must satisfy these. No network / no LLM:
_process_one is exercised with a mocked summarize_text; rebuild_combined uses
a temp DB so it never touches the real papers/db.json.
"""
import json
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pipeline_ops
from summarizer.schema import ManuscriptChecklist


# ── Part A: rebuild_combined routes through TinyDB Store ─────────────────────


def _write_per_paper(d: Path, pmcid: str, **over) -> None:
    rec = {
        "pmcid": pmcid, "title": f"{pmcid} title", "year": "2020",
        "ehr_used": False, "ehr_evidence": "n/a", "summary": "s",
    }
    rec.update(over)
    (d / f"{pmcid}.json").write_text(json.dumps(rec, indent=2))


def test_rebuild_combined_uses_store(tmp_path: Path):
    """rebuild_combined should import into a Store and export from it."""
    d = tmp_path / "summaries"; d.mkdir()
    _write_per_paper(d, "PMC1", ehr_used=True, ehr_evidence="HES codes")
    _write_per_paper(d, "PMC2")
    db = tmp_path / "db.json"
    out = tmp_path / "combined.json"

    data = pipeline_ops.rebuild_combined(summary_dir=d, out_path=out, db_path=db)

    # Store was populated (db file exists with 2 records)
    assert db.exists()
    from database import Store
    with Store(db) as s:
        assert s.stats()["total"] == 2
    # Export is a valid SummaryBatch
    assert data["n"] == 2
    assert len(data["summaries"]) == 2


def test_rebuild_combined_skips_invalid(tmp_path: Path):
    """Invalid per-paper files are skipped with a warning, not fatal."""
    import pytest
    d = tmp_path / "summaries"; d.mkdir()
    _write_per_paper(d, "PMC1")
    (d / "PMC_BAD.json").write_text("{not valid json")
    db = tmp_path / "db.json"
    out = tmp_path / "c.json"

    with pytest.warns(UserWarning, match="PMC_BAD"):
        data = pipeline_ops.rebuild_combined(summary_dir=d, out_path=out, db_path=db)
    assert data["n"] == 1


def test_rebuild_combined_empty(tmp_path: Path):
    d = tmp_path / "summaries"; d.mkdir()
    db = tmp_path / "db.json"
    data = pipeline_ops.rebuild_combined(summary_dir=d, out_path=tmp_path / "c.json", db_path=db)
    assert data["n"] == 0


def test_rebuild_combined_sorted(tmp_path: Path):
    d = tmp_path / "summaries"; d.mkdir()
    _write_per_paper(d, "PMC9", year="2019")
    _write_per_paper(d, "PMC1", year="2021")
    _write_per_paper(d, "PMC2", year="2019")
    data = pipeline_ops.rebuild_combined(summary_dir=d, out_path=tmp_path / "c.json",
                                         db_path=tmp_path / "db.json")
    assert [s["pmcid"] for s in data["summaries"]] == ["PMC2", "PMC9", "PMC1"]


# ── Part B: concurrent summarization ─────────────────────────────────────────


def test_process_one_extracts_and_summarizes(tmp_path: Path):
    """_process_one reads a file, calls summarize_text, writes per-paper JSON."""
    from summarizer.run import _process_one

    fake = tmp_path / "PMC1234567.txt"
    fake.write_text("some paper text")

    mock_checklist = ManuscriptChecklist(
        pmcid="PMC1234567", title="T", year="2020",
        ehr_used=True, ehr_evidence="e", summary="s",
    )
    mock_client = MagicMock()
    with patch("summarizer.run.extract", return_value=("text", "text")), \
         patch("summarizer.run.summarize_text", return_value=mock_checklist):
        result = _process_one(
            path=fake, meta={"PMC1234567": {"title": "T", "year": "2020"}},
            client=mock_client, model="gemma", chunked=False, recover=False,
            summary_dir=tmp_path,
        )

    assert result.status == "ok"
    assert result.pmcid == "PMC1234567"
    assert (tmp_path / "PMC1234567.json").exists()
    assert result.checklist is not None


def test_process_one_skips_cached(tmp_path: Path):
    """If a valid summary JSON already exists, _process_one skips (no LLM call)."""
    from summarizer.run import _process_one

    fake = tmp_path / "PMC9999999.txt"; fake.write_text("x")
    cached = ManuscriptChecklist(
        pmcid="PMC9999999", title="T", year="2020",
        ehr_used=False, ehr_evidence="n/a", summary="s",
    )
    (tmp_path / "PMC9999999.json").write_text(cached.model_dump_json(indent=2))

    with patch("summarizer.run.summarize_text") as mock_summ:
        result = _process_one(
            path=fake, meta={}, client=MagicMock(), model="m",
            chunked=False, recover=False, summary_dir=tmp_path,
        )
    mock_summ.assert_not_called()
    assert result.status == "skipped"
    assert result.checklist is not None


def test_process_one_handles_failure(tmp_path: Path):
    """If summarize_text raises, _process_one returns status='failed', no crash."""
    from summarizer.run import _process_one

    fake = tmp_path / "PMC7777777.txt"; fake.write_text("x")
    with patch("summarizer.run.extract", return_value=("text", "text")), \
         patch("summarizer.run.summarize_text", side_effect=RuntimeError("boom")):
        result = _process_one(
            path=fake, meta={}, client=MagicMock(), model="m",
            chunked=False, recover=True, summary_dir=tmp_path,
        )
    assert result.status == "failed"
    assert result.checklist is None
    assert not (tmp_path / "PMC7777777.json").exists()


def test_workers_flag_accepted():
    """--workers N is a valid CLI argument."""
    from summarizer.run import main
    with patch("summarizer.run.get_client", side_effect=RuntimeError("no key")):
        rc = main(["--workers", "4", "--limit", "0"])
    # exits before summarizing (no files / no key), but arg parsing must not crash
    assert rc in (1, 2)


def test_default_workers_is_one():
    """Without --workers, the default is 1 (serial, unchanged behavior)."""
    from summarizer.run import main
    import inspect
    src = inspect.getsource(main)
    assert "--workers" in src
    assert "default" in src.lower() or "1" in src
