#!/usr/bin/env python3
"""CLI for the TinyDB-backed manuscript summary store.

Examples
  python db.py add PMC1234567 --title "..." --year 2020 \
      --ehr-used --ehr-evidence "HES ICD codes" --summary "..." \
      --disease asthma --exposure-domain "air pollution" --confidence medium
  python db.py add PMC1234567 --file path/to/checklist.json
  python db.py get PMC1234567
  python db.py update PMC1234567 --no-ehr-used --ehr-evidence "n/a"
  python db.py list
  python db.py find --ehr-used --disease asthma
  python db.py import --from-dir papers/summaries
  python db.py export                       # -> papers/manuscript_summaries.json + results/...
  python db.py stats
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from database import (
    DEFAULT_COMBINED,
    DEFAULT_DB_PATH,
    DEFAULT_RESULTS,
    DEFAULT_SUMMARY_DIR,
    Store,
)


# ── argument helpers ────────────────────────────────────────────────────────

SCALAR_FIELDS = {
    "title": "--title",
    "year": "--year",
    "ehr_evidence": "--ehr-evidence",
    "summary": "--summary",
    "study_design": "--study-design",
    "data_source_type": "--data-source",
    "population": "--population",
    "exposure_domain": "--exposure-domain",
    "confidence": "--confidence",
    "source_format": "--source-format",
    "model": "--model",
}
LIST_FIELDS = {
    "key_findings": "--key-finding",
    "captured_features": "--captured-feature",
    "pathologies_diseases": "--disease",
    "limitations": "--limitation",
}


def _add_record_args(p: argparse.ArgumentParser) -> None:
    """Field flags shared by `add` and `update`."""
    for field, flag in SCALAR_FIELDS.items():
        p.add_argument(flag, dest=field, default=None, help=f"{field}")
    p.add_argument(
        "--ehr-used", dest="ehr_used", action="store_true", default=None,
        help="set ehr_used true",
    )
    p.add_argument(
        "--no-ehr-used", dest="ehr_used", action="store_false",
        help="set ehr_used false",
    )
    for field, flag in LIST_FIELDS.items():
        p.add_argument(flag, dest=field, action="append", default=None, help=f"{field} (repeatable)")
    p.add_argument("--json", dest="json_blob", default=None, help="inline JSON record / patch")
    p.add_argument("--file", dest="json_file", default=None, help="path to a JSON record / patch")


def _record_from_args(args: argparse.Namespace) -> dict:
    """Build a record dict from --file, --json, and field flags (flags win)."""
    rec: dict = {}
    if args.json_file:
        rec.update(json.loads(Path(args.json_file).read_text()))
    if args.json_blob:
        rec.update(json.loads(args.json_blob))
    for field in SCALAR_FIELDS:
        v = getattr(args, field, None)
        if v is not None:
            rec[field] = v
    if getattr(args, "ehr_used", None) is not None:
        rec["ehr_used"] = args.ehr_used
    for field in LIST_FIELDS:
        v = getattr(args, field, None)
        if v is not None:
            rec[field] = v
    return rec


# ── output helpers ──────────────────────────────────────────────────────────

def _row(r: dict) -> str:
    ehr = "EHR" if r.get("ehr_used") else "  -"
    dis = "; ".join(r.get("pathologies_diseases") or []) or "—"
    return f"{r.get('pmcid','?'):<13} {str(r.get('year','?'))[:4]:<4} {ehr}  {r.get('title','')}"


# ── commands ────────────────────────────────────────────────────────────────

def cmd_add(args: argparse.Namespace) -> int:
    rec = _record_from_args(args)
    rec["pmcid"] = args.pmcid
    with Store(args.db) as store:
        out = store.add(rec)
    print(f"✓ added {out['pmcid']}")
    return 0


def cmd_get(args: argparse.Namespace) -> int:
    with Store(args.db) as store:
        rec = store.get(args.pmcid)
    if rec is None:
        print(f"not found: {_norm(args.pmcid)}", file=sys.stderr)
        return 1
    print(json.dumps(rec, indent=2, sort_keys=True))
    return 0


def cmd_update(args: argparse.Namespace) -> int:
    patch = _record_from_args(args)
    patch.pop("pmcid", None)  # pmcid is the key, set via the positional
    with Store(args.db) as store:
        out = store.update(args.pmcid, patch)
    print(f"✓ updated {out['pmcid']}")
    return 0


def cmd_delete(args: argparse.Namespace) -> int:
    if not args.yes:
        print(f"refusing to delete {_norm(args.pmcid)} without --yes", file=sys.stderr)
        return 2
    with Store(args.db) as store:
        ok = store.delete(args.pmcid)
    print("✓ deleted" if ok else "— not found")
    return 0 if ok else 1


def cmd_list(args: argparse.Namespace) -> int:
    with Store(args.db) as store:
        rows = store.list_all()
    print(f"# {len(rows)} records\n")
    for r in rows:
        print(_row(r))
    return 0


def cmd_find(args: argparse.Namespace) -> int:
    filters = {}
    if args.ehr_used is not None:
        filters["ehr_used"] = args.ehr_used
    for key, flag in [
        ("disease", args.disease), ("exposure", args.exposure),
        ("query", args.query), ("design", args.design),
        ("data_source", args.data_source), ("confidence", args.confidence),
    ]:
        if flag is not None:
            filters[key] = flag
    with Store(args.db) as store:
        rows = store.find(**filters)
    print(f"# {len(rows)} match(es)\n")
    for r in rows:
        print(_row(r))
    return 0


def cmd_import(args: argparse.Namespace) -> int:
    with Store(args.db) as store:
        if args.from_dir:
            ins, upd = store.import_records_from_dir(args.from_dir)
        elif args.from_combined:
            data = json.loads(Path(args.from_combined).read_text())
            summaries = data.get("summaries", data if isinstance(data, list) else [])
            ins, upd = store.import_records(summaries)
        else:
            print("import needs --from-dir or --from-combined", file=sys.stderr)
            return 2
    print(f"✓ imported: {ins} new, {upd} updated")
    return 0


def cmd_export(args: argparse.Namespace) -> int:
    with Store(args.db) as store:
        data = store.export_combined(args.combined)
        if args.results:
            Path(args.results).parent.mkdir(parents=True, exist_ok=True)
            Path(args.results).write_text(json.dumps(data, indent=2))
    print(f"✓ exported {data['n']} records")
    print(f"  combined: {args.combined}")
    if args.results:
        print(f"  results:  {args.results}")
    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    with Store(args.db) as store:
        s = store.stats()
    for k, v in s.items():
        print(f"{k:<16} {v}")
    return 0


def _norm(pmcid: str) -> str:
    v = (pmcid or "").strip().upper()
    return v if v.startswith("PMC") else f"PMC{v}"


# ── argparse ────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="TinyDB manuscript summary store (CRUD CLI).")
    p.add_argument("--db", default=str(DEFAULT_DB_PATH), help=f"DB path (default {DEFAULT_DB_PATH})")
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("add", help="add a new record")
    a.add_argument("pmcid"); _add_record_args(a); a.set_defaults(func=cmd_add)

    g = sub.add_parser("get", help="print a record as JSON")
    g.add_argument("pmcid"); g.set_defaults(func=cmd_get)

    u = sub.add_parser("update", help="merge fields onto a record (re-validates)")
    u.add_argument("pmcid"); _add_record_args(u); u.set_defaults(func=cmd_update)

    d = sub.add_parser("delete", help="delete a record")
    d.add_argument("pmcid"); d.add_argument("--yes", action="store_true"); d.set_defaults(func=cmd_delete)

    l = sub.add_parser("list", help="list all records")
    l.set_defaults(func=cmd_list)

    f = sub.add_parser("find", help="filter records")
    f.add_argument("--ehr-used", dest="ehr_used", action="store_true", default=None)
    f.add_argument("--no-ehr-used", dest="ehr_used", action="store_false")
    f.add_argument("--disease", default=None)
    f.add_argument("--exposure", default=None)
    f.add_argument("--query", default=None, help="substring of title or summary")
    f.add_argument("--design", default=None)
    f.add_argument("--source", dest="data_source", default=None)
    f.add_argument("--confidence", default=None)
    f.set_defaults(func=cmd_find)

    i = sub.add_parser("import", help="ingest existing JSONs into the DB (upsert)")
    i.add_argument("--from-dir", default=None, help=f"dir of <pmcid>.json (default {DEFAULT_SUMMARY_DIR})")
    i.add_argument("--from-combined", default=None, help="path to a combined SummaryBatch JSON")
    i.set_defaults(func=cmd_import)

    e = sub.add_parser("export", help="export combined SummaryBatch JSON")
    e.add_argument("--combined", default=str(DEFAULT_COMBINED), help=f"(default {DEFAULT_COMBINED})")
    e.add_argument("--results", default=str(DEFAULT_RESULTS), help=f"tracked copy (default {DEFAULT_RESULTS})")
    e.add_argument("--no-results", action="store_true", help="skip the results/ copy")
    e.set_defaults(func=cmd_export)

    s = sub.add_parser("stats", help="record counts")
    s.set_defaults(func=cmd_stats)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if getattr(args, "no_results", False):
        args.results = None
    try:
        return args.func(args)
    except ValueError as e:  # validation / business-rule errors
        print(f"error: {e}", file=sys.stderr)
        return 2
    except KeyError as e:  # missing record
        print(f"error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())