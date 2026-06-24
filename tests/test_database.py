"""Tests for the TinyDB-backed manuscript summary store (feature/tinydb-manuscript-store).

Tests written first — implementation must satisfy these. No live API calls.
"""
import json
import sys
from pathlib import Path

# allow `pytest` from repo root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest

from summarizer.schema import ManuscriptChecklist, SummaryBatch
from database import Store


# ── fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture
def store(tmp_path: Path) -> Store:
    return Store(tmp_path / "db.json")


def _valid_record(pmcid: str = "PMC1234567", **overrides) -> dict:
    base = {
        "pmcid": pmcid,
        "title": "Childhood T1DM EWAS across England",
        "year": "2020",
        "ehr_used": True,
        "ehr_evidence": "Hospital Episode Statistics ICD codes",
        "summary": "EWAS of childhood T1DM incidence across England.",
        "key_findings": ["15 environmental factors associated with T1DM."],
        "captured_features": ["HES ICD codes"],
        "pathologies_diseases": ["type 1 diabetes"],
        "study_design": "ecological EWAS",
        "data_source_type": "EHR",
        "population": "children 0-9 yrs across England",
        "exposure_domain": "air pollution",
        "limitations": ["ecological design"],
        "confidence": "medium",
    }
    base.update(overrides)
    return base


# ── add / get ────────────────────────────────────────────────────────────────


def test_add_validates_and_stores(store: Store):
    rec = store.add(_valid_record())
    assert rec["pmcid"] == "PMC1234567"  # normalized to full form
    got = store.get("1234567")  # lookup tolerates bare digits
    assert got is not None
    assert got["ehr_used"] is True
    assert got["pathologies_diseases"] == ["type 1 diabetes"]


def test_add_normalizes_pmcid(store: Store):
    rec = store.add(_valid_record(pmcid="7654321"))
    assert rec["pmcid"] == "PMC7654321"
    assert store.get("PMC7654321") is not None


def test_add_rejects_invalid_record(store: Store):
    bad = _valid_record()
    bad["ehr_used"] = 5  # not coercible to bool (schema coerces strings, not ints)
    with pytest.raises(ValueError):
        store.add(bad)


def test_add_duplicate_rejects(store: Store):
    store.add(_valid_record())
    with pytest.raises(ValueError, match="exists"):
        store.add(_valid_record())


def test_get_missing_returns_none(store: Store):
    assert store.get("PMC0000000") is None


# ── list / stats ─────────────────────────────────────────────────────────────


def test_list_sorted_by_year_then_pmcid(store: Store):
    store.add(_valid_record("PMC0000001", year="2019"))
    store.add(_valid_record("PMC0000003", year="2019"))
    store.add(_valid_record("PMC0000002", year="2021"))
    rows = store.list_all()
    pmcids = [r["pmcid"] for r in rows]
    assert pmcids == ["PMC0000001", "PMC0000003", "PMC0000002"]


def test_stats_counts(store: Store):
    store.add(_valid_record("PMC1", ehr_used=True))
    store.add(_valid_record("PMC2", ehr_used=False))
    s = store.stats()
    assert s["total"] == 2
    assert s["ehr_used_true"] == 1
    assert s["ehr_used_false"] == 1


# ── update ───────────────────────────────────────────────────────────────────


def test_update_merges_and_revalidates(store: Store):
    store.add(_valid_record())
    store.update("PMC1234567", {"ehr_used": False, "ehr_evidence": "n/a", "confidence": "low"})
    got = store.get("PMC1234567")
    assert got["ehr_used"] is False
    assert got["ehr_evidence"] == "n/a"
    assert got["confidence"] == "low"
    # untouched fields preserved
    assert got["title"] == "Childhood T1DM EWAS across England"


def test_update_rejects_invalid(store: Store):
    store.add(_valid_record())
    with pytest.raises(ValueError):
        store.update("PMC1234567", {"pmcid": "PMC999"})  # can't rekey
    with pytest.raises(ValueError):
        store.update("PMC1234567", {"ehr_used": 5})  # not coercible to bool


def test_update_missing_raises(store: Store):
    with pytest.raises(KeyError):
        store.update("PMC0000000", {"summary": "x"})


# ── delete ───────────────────────────────────────────────────────────────────


def test_delete_removes(store: Store):
    store.add(_valid_record())
    assert store.delete("PMC1234567") is True
    assert store.get("PMC1234567") is None
    assert store.delete("PMC1234567") is False  # already gone


# ── find ──────────────────────────────────────────────────────────────────────


def test_find_by_ehr_and_disease(store: Store):
    store.add(_valid_record("PMC1", ehr_used=True, pathologies_diseases=["asthma"]))
    store.add(_valid_record("PMC2", ehr_used=False, pathologies_diseases=["type 1 diabetes"]))
    assert [r["pmcid"] for r in store.find(ehr_used=True)] == ["PMC1"]
    assert [r["pmcid"] for r in store.find(disease="asthma")] == ["PMC1"]
    assert [r["pmcid"] for r in store.find(ehr_used=False)] == ["PMC2"]
    assert len(store.find(exposure="air pollution")) == 2  # both share exposure_domain


def test_find_by_query_matches_title_or_summary(store: Store):
    store.add(_valid_record("PMC1", title="Mercury exposure in kids", summary="lead study"))
    # 'mercury' only in a title
    hits = store.find(query="mercury")
    assert [r["pmcid"] for r in hits] == ["PMC1"]
    # 'lead' only in summary
    hits = store.find(query="lead")
    assert [r["pmcid"] for r in hits] == ["PMC1"]


# ── export ────────────────────────────────────────────────────────────────────


def test_export_combined_roundtrips_through_summarybatch(store: Store, tmp_path: Path):
    store.add(_valid_record("PMC1"))
    store.add(_valid_record("PMC2"))
    out = tmp_path / "manuscript_summaries.json"
    store.export_combined(out)
    data = json.loads(out.read_text())
    # build_results.load() expects a top-level "summaries" list
    assert set(data.keys()) >= {"n", "model", "summaries"}
    assert data["n"] == 2
    assert len(data["summaries"]) == 2
    # the export must re-validate as a real SummaryBatch
    batch = SummaryBatch.model_validate(data)
    assert batch.n == 2
    assert all(isinstance(s, ManuscriptChecklist) for s in batch.summaries)


# ── import ───────────────────────────────────────────────────────────────────


def test_import_from_records_upserts(store: Store):
    inserted, updated = store.import_records([_valid_record("PMC1"), _valid_record("PMC2")])
    assert inserted == 2 and updated == 0
    assert store.stats()["total"] == 2
    # re-import -> idempotent: only changed rows count as updates, none inserted
    inserted, updated = store.import_records(
        [_valid_record("PMC1", ehr_used=False, ehr_evidence="n/a"), _valid_record("PMC2")]
    )
    assert inserted == 0  # no new rows
    assert updated == 1    # only PMC1 changed; unchanged PMC2 is not an update
    assert store.stats()["total"] == 2  # no growth
    assert store.get("PMC1")["ehr_used"] is False


def test_import_from_per_paper_dir(store: Store, tmp_path: Path):
    d = tmp_path / "summaries"
    d.mkdir()
    (d / "PMC1111111.json").write_text(json.dumps(_valid_record("PMC1111111"), indent=2))
    (d / "PMC2222222.json").write_text(json.dumps(_valid_record("PMC2222222"), indent=2))
    # a stray non-json file should be skipped without error
    (d / "README.txt").write_text("ignore me")
    inserted, updated = store.import_records_from_dir(d)
    assert inserted == 2
    assert store.stats()["total"] == 2


# ── sanity: pip dep present ──────────────────────────────────────────────────


def test_tinydb_importable():
    import tinydb  # noqa: F401