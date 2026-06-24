"""TinyDB-backed store for manuscript summary JSONs.

Per-paper summaries (``summarizer.schema.ManuscriptChecklist``) live as
documents in a single TinyDB JSON file. Every ``add`` / ``update`` validates
through the existing Pydantic schema, so the store is never left in a state
that ``build_results.py`` can't read. ``export_combined`` regenerates the
combined ``SummaryBatch`` JSON that the rest of the pipeline consumes.

The DB file defaults to ``papers/db.json`` (gitignored, like the rest of
``papers/``); the tracked review artifact is the ``results/`` copy written by
``db.py export``.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from tinydb import TinyDB
from tinydb.table import Document

from summarizer.schema import ManuscriptChecklist, SummaryBatch

DEFAULT_DB_PATH = Path("papers/db.json")
DEFAULT_COMBINED = Path("papers/manuscript_summaries.json")
DEFAULT_RESULTS = Path("results/manuscript_summaries.json")
DEFAULT_SUMMARY_DIR = Path("papers/summaries")

# The combined file historically carried the LLM model used for a summarizer
# run. The DB aggregates records that may span runs, so the export derives the
# batch ``model`` from the records themselves (most common non-empty value),
# falling back to an empty string. ``build_results.py`` does not read it.
EXPORT_MODEL = ""


def _normalize_pmcid(pmcid: str) -> str:
    """Same normalization as ``ManuscriptChecklist._normalize_pmcid``."""
    v = (pmcid or "").strip().upper()
    if not v.startswith("PMC"):
        v = f"PMC{v}"
    return v


def _most_common_model(rows: list[Mapping[str, Any]]) -> str:
    """Most common non-empty per-record ``model``; else ``EXPORT_MODEL``."""
    from collections import Counter

    counts = Counter(str(r.get("model", "")) for r in rows if str(r.get("model", "")).strip())
    return counts.most_common(1)[0][0] if counts else EXPORT_MODEL


def _validated_dump(data: Mapping[str, Any]) -> dict:
    """Validate ``data`` as a ManuscriptChecklist and return the normalized dict."""
    return ManuscriptChecklist(**dict(data)).model_dump()


class Store:
    """A validated CRUD layer over a TinyDB JSON document store."""

    def __init__(self, path: str | Path = DEFAULT_DB_PATH) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # indent keeps the JSON human-readable; each row is one Document.
        self._db = TinyDB(str(self.path), indent=2, sort_keys=True)

    # ── lifecycle ───────────────────────────────────────────────────────────
    def close(self) -> None:
        self._db.close()

    def __enter__(self) -> "Store":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # ── query helpers ───────────────────────────────────────────────────────
    def _find_doc(self, pmcid: str) -> Document | None:
        target = _normalize_pmcid(pmcid)
        docs = self._db.search(lambda d: d.get("pmcid") == target)
        return docs[0] if docs else None

    # ── create ──────────────────────────────────────────────────────────────
    def add(self, data: Mapping[str, Any]) -> dict:
        """Validate and insert a new record. Raises if the pmcid already exists."""
        validated = _validated_dump(data)
        if self._find_doc(validated["pmcid"]) is not None:
            raise ValueError(f"pmcid {validated['pmcid']} already exists; use update() instead")
        self._db.insert(validated)
        return validated

    # ── read ───────────────────────────────────────────────────────────────
    def get(self, pmcid: str) -> dict | None:
        doc = self._find_doc(pmcid)
        return dict(doc) if doc is not None else None

    def list_all(self) -> list[dict]:
        rows = [dict(d) for d in self._db.all()]
        rows.sort(key=lambda r: (str(r.get("year", "")), str(r.get("pmcid", ""))))
        return rows

    def stats(self) -> dict:
        rows = self.list_all()
        return {
            "total": len(rows),
            "ehr_used_true": sum(1 for r in rows if r.get("ehr_used") is True),
            "ehr_used_false": sum(1 for r in rows if r.get("ehr_used") is False),
        }

    # ── update ─────────────────────────────────────────────────────────────
    def update(self, pmcid: str, partial: Mapping[str, Any]) -> dict:
        """Merge ``partial`` onto the existing record, re-validate, and replace.

        The ``pmcid`` cannot be changed (it is the key); pass it and it must match.
        """
        existing = self.get(pmcid)
        if existing is None:
            raise KeyError(f"pmcid {_normalize_pmcid(pmcid)} not found")
        if "pmcid" in partial and _normalize_pmcid(str(partial["pmcid"])) != existing["pmcid"]:
            raise ValueError("cannot change pmcid via update(); delete + add instead")
        merged = {**existing, **dict(partial)}
        validated = _validated_dump(merged)  # raises on bad data
        self._db.remove(lambda d: d.get("pmcid") == validated["pmcid"])
        self._db.insert(validated)
        return validated

    # ── delete ──────────────────────────────────────────────────────────────
    def delete(self, pmcid: str) -> bool:
        target = _normalize_pmcid(pmcid)
        before = len(self._db.all())
        self._db.remove(lambda d: d.get("pmcid") == target)
        return len(self._db.all()) < before

    # ── find ───────────────────────────────────────────────────────────────
    def find(
        self,
        *,
        ehr_used: bool | None = None,
        disease: str | None = None,
        exposure: str | None = None,
        query: str | None = None,
        design: str | None = None,
        data_source: str | None = None,
        confidence: str | None = None,
    ) -> list[dict]:
        """Filter records in memory (the catalog is small). All filters are AND."""
        def _has(needle: str, hay: str | Iterable[str]) -> bool:
            n = needle.lower()
            if isinstance(hay, str):
                return n in hay.lower()
            return any(n in str(x).lower() for x in hay)

        out: list[dict] = []
        for r in self.list_all():
            if ehr_used is not None and bool(r.get("ehr_used")) is not ehr_used:
                continue
            if disease is not None and not _has(disease, r.get("pathologies_diseases", [])):
                continue
            if exposure is not None and not _has(exposure, r.get("exposure_domain", "")):
                continue
            if design is not None and not _has(design, r.get("study_design", "")):
                continue
            if data_source is not None and not _has(data_source, r.get("data_source_type", "")):
                continue
            if confidence is not None and str(r.get("confidence", "")).lower() != confidence.lower():
                continue
            if query is not None and not _has(
                query, f"{r.get('title', '')} {r.get('summary', '')}"
            ):
                continue
            out.append(r)
        return out

    # ── export ──────────────────────────────────────────────────────────────
    def export_combined(self, path: str | Path = DEFAULT_COMBINED) -> dict:
        """Write a ``SummaryBatch`` JSON; re-using the schema re-validates everything."""
        summaries = self.list_all()
        batch = SummaryBatch(n=len(summaries), model=_most_common_model(summaries), summaries=summaries)
        payload = batch.model_dump_json(indent=2)
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(payload)
        return json.loads(payload)

    # ── import ──────────────────────────────────────────────────────────────
    def import_records(self, records: Iterable[Mapping[str, Any]]) -> tuple[int, int]:
        """Upsert records by pmcid. Returns (inserted, updated)."""
        inserted = updated = 0
        for rec in records:
            validated = _validated_dump(rec)
            existing = self._find_doc(validated["pmcid"])
            if existing is None:
                self._db.insert(validated)
                inserted += 1
            elif dict(existing) != validated:
                self._db.remove(lambda d: d.get("pmcid") == validated["pmcid"])
                self._db.insert(validated)
                updated += 1
        return inserted, updated

    def import_records_from_dir(self, dir_path: str | Path = DEFAULT_SUMMARY_DIR) -> tuple[int, int]:
        """Import every ``<pmcid>.json`` under ``dir_path`` (idempotent upsert)."""
        d = Path(dir_path)
        records: list[dict] = []
        for f in sorted(d.glob("*.json")):
            records.append(json.loads(f.read_text()))
        return self.import_records(records)