"""Tests for the Dagster asset orchestration layer (feature/dagster-orchestration).

Tests written first — implementation must satisfy these. No network / no LLM: the
subprocess wrappers are exercised with mocked subprocess; the combined-file rebuild
is pure Python.
"""
import json
import sys
from pathlib import Path
from unittest.mock import patch

# allow `pytest` from repo root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pipeline
import pipeline_ops
from summarizer.schema import SummaryBatch


REPO_ROOT = Path(__file__).resolve().parent.parent


# ── asset wiring / lineage ───────────────────────────────────────────────────


def test_definitions_loads_and_has_four_assets():
    from dagster import Definitions
    assert isinstance(pipeline.defs, Definitions)
    graph = pipeline.defs.resolve_asset_graph()
    keys = {".".join(k.path) for k in graph.get_all_asset_keys()}
    assert keys == {"download_log", "per_paper_summaries",
                    "manuscript_summaries", "results"}


def test_asset_lineage_is_wired():
    deps = {name: {".".join(k.path) for k in fn.dependency_keys}
            for name, fn in [
                ("download_log", pipeline.download_log),
                ("per_paper_summaries", pipeline.per_paper_summaries),
                ("manuscript_summaries", pipeline.manuscript_summaries),
                ("results", pipeline.results),
            ]}
    assert deps["download_log"] == set()
    assert deps["per_paper_summaries"] == {"download_log"}
    assert deps["manuscript_summaries"] == {"per_paper_summaries"}
    assert deps["results"] == {"manuscript_summaries"}


# ── subprocess wrappers (mocked — no network / no LLM) ───────────────────────


def test_fetch_papers_invokes_fetcher_module():
    with patch("pipeline_ops.subprocess.run") as run:
        run.return_value.returncode = 0
        pipeline_ops.fetch_papers()
    assert run.called
    cmd = run.call_args[0][0]
    assert cmd[:3] == [sys.executable, "-m", "fetch_pmc_papers"]
    assert run.call_args[1]["cwd"] == str(REPO_ROOT)


def test_summarize_papers_invokes_summarizer_run():
    with patch("pipeline_ops.subprocess.run") as run:
        run.return_value.returncode = 0
        pipeline_ops.summarize_papers()
    cmd = run.call_args[0][0]
    assert cmd[:3] == [sys.executable, "-m", "summarizer.run"]
    assert "--recover" in cmd  # resume-only mode


def test_build_results_invokes_build_results_script():
    with patch("pipeline_ops.subprocess.run") as run:
        run.return_value.returncode = 0
        pipeline_ops.build_results()
    cmd = run.call_args[0][0]
    assert cmd[:2] == [sys.executable, str(REPO_ROOT / "build_results.py")]


def test_ops_propagate_nonzero_exit_as_runtimeerror():
    with patch("pipeline_ops.subprocess.run") as run:
        run.return_value.returncode = 1
        import pytest
        with pytest.raises(RuntimeError, match="fetch_pmc_papers"):
            pipeline_ops.fetch_papers()


# ── rebuild_combined (pure python; mirrors TinyDB export without the dep) ──────


def _write_per_paper(d: Path, pmcid: str, **over) -> None:
    rec = {
        "pmcid": pmcid, "title": f"{pmcid} title", "year": "2020",
        "ehr_used": False, "ehr_evidence": "n/a", "summary": "s",
    }
    rec.update(over)
    (d / f"{pmcid}.json").write_text(json.dumps(rec, indent=2))


def test_rebuild_combined_roundtrips_through_summarybatch(tmp_path: Path):
    d = tmp_path / "summaries"; d.mkdir()
    _write_per_paper(d, "PMC1", ehr_used=True, ehr_evidence="HES codes")
    _write_per_paper(d, "PMC2")
    (d / "README.txt").write_text("ignore")  # non-json skipped

    out = tmp_path / "combined.json"
    data = pipeline_ops.rebuild_combined(summary_dir=d, out_path=out, db_path=tmp_path/"db.json")
    assert data["n"] == 2
    assert len(data["summaries"]) == 2
    batch = SummaryBatch.model_validate(data)  # re-validates every record
    assert batch.n == 2
    assert {s.pmcid for s in batch.summaries} == {"PMC1", "PMC2"}
    assert out.read_text()  # written to disk


def test_rebuild_combined_sorted_by_year_then_pmcid(tmp_path: Path):
    d = tmp_path / "summaries"; d.mkdir()
    _write_per_paper(d, "PMC9", year="2019")
    _write_per_paper(d, "PMC1", year="2021")
    _write_per_paper(d, "PMC2", year="2019")
    data = pipeline_ops.rebuild_combined(summary_dir=d, out_path=tmp_path/"c.json", db_path=tmp_path/"db.json")
    assert [s["pmcid"] for s in data["summaries"]] == ["PMC2", "PMC9", "PMC1"]


def test_rebuild_combined_empty_is_valid(tmp_path: Path):
    d = tmp_path / "summaries"; d.mkdir()
    data = pipeline_ops.rebuild_combined(summary_dir=d, out_path=tmp_path/"c.json", db_path=tmp_path/"db.json")
    assert data["n"] == 0
    SummaryBatch.model_validate(data)  # empty batch is valid


def test_rebuild_combined_skips_invalid_per_paper_file(tmp_path: Path):
    d = tmp_path / "summaries"; d.mkdir()
    _write_per_paper(d, "PMC1")
    (d / "PMC_BAD.json").write_text("{not valid json")  # corrupt → skipped, not fatal
    import pytest
    with pytest.warns(UserWarning, match="PMC_BAD"):
        data = pipeline_ops.rebuild_combined(summary_dir=d, out_path=tmp_path/"c.json", db_path=tmp_path/"db.json")
    assert data["n"] == 1  # only the good one kept