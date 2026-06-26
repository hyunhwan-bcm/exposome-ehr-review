#!/usr/bin/env python3
"""
Summarize collected manuscripts into Pydantic-validated JSON checklists via
Gemma 4 12B.

Usage:
    python -m summarizer.run                 # all downloaded papers
    python -m summarizer.run --pmcid PMC7145790   # single paper
    python -m summarizer.run --limit 5       # pilot on first N
    python -m summarizer.run --force          # re-summarize even if cached

Env:
    GEMMA_BASE_URL  default http://bioinfolder.com:8000/v1
    GEMMA_API_KEY    required
    GEMMA_MODEL      default unsloth/gemma-4-12B-it-qat-GGUF
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

# Load .env if present (optional dependency — skipped if unavailable).
try:
    from dotenv import load_dotenv
    load_dotenv()
except Exception:
    pass

from .extract import extract, pmcid_from_filename
from .llm_client import get_client, summarize_text, summarize_chunked
from .schema import ManuscriptChecklist, SummaryBatch


@dataclass
class PaperResult:
    """Result of processing a single paper."""
    pmcid: str
    status: str  # "ok" | "skipped" | "failed"
    checklist: ManuscriptChecklist | None = None

PAPERS_DIR = Path("papers")
SUMMARY_DIR = PAPERS_DIR / "summaries"
COMBINED_PATH = PAPERS_DIR / "manuscript_summaries.json"
LOG_PATH = Path("papers/download_log.json")


def load_metadata() -> dict[str, dict]:
    """pmcid -> {title, year, journal, authors} from the download log."""
    if not LOG_PATH.exists():
        return {}
    log = json.loads(LOG_PATH.read_text())
    return {p["pmcid"]: p for p in log.get("papers", [])}


def _fetch_missing_metadata(pmids: list[str]) -> dict[str, dict]:
    """Fetch PMCID title/year from PMC ESummary when the local download log lacks it."""
    out: dict[str, dict] = {}
    for i in range(0, len(pmids), 20):
        ids = [p.removeprefix("PMC") for p in pmids[i:i + 20]]
        qs = urllib.parse.urlencode({"db": "pmc", "id": ",".join(ids), "retmode": "json"})
        try:
            with urllib.request.urlopen(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{qs}", timeout=30) as r:
                data = json.loads(r.read())
        except Exception:
            continue
        for uid in data.get("result", {}).get("uids", []):
            item = data["result"].get(uid, {})
            out[f"PMC{uid}"] = {
                "pmcid": f"PMC{uid}",
                "title": item.get("title", ""),
                "journal": item.get("source", ""),
                "year": item.get("pubdate", "")[:4],
                "authors": ", ".join(a.get("name", "") for a in item.get("authors", [])[:2]),
            }
        time.sleep(0.4)
    return out


def _filename_like_title(title: str, pmcid: str) -> bool:
    return bool(re.match(r"^\d{4}_" + re.escape(pmcid) + r"_", title or ""))


def discover_files() -> list[Path]:
    """All downloaded full-text files (PDF + XML), pmcid-sorted."""
    files = sorted(
        list(PAPERS_DIR.glob("*.pdf")) + list(PAPERS_DIR.glob("*.xml")),
        key=lambda p: pmcid_from_filename(p),
    )
    return files


def find_failed() -> list[Path]:
    """Files whose PMCID has no summary JSON yet (the failed/missing set)."""
    done = {p.stem for p in SUMMARY_DIR.glob("*.json")}
    return [f for f in discover_files() if pmcid_from_filename(f) not in done]


def _process_one(
    path: Path,
    meta: dict,
    client,
    model: str,
    chunked: bool,
    recover: bool,
    summary_dir: Path,
) -> PaperResult:
    """Process a single paper: extract -> summarize -> write per-paper JSON.

    Thread-safe: each paper writes to its own file; the OpenAI client is
    safe for concurrent use. Returns a PaperResult (never raises).
    """
    pmcid = pmcid_from_filename(path)
    out_path = summary_dir / f"{pmcid}.json"
    m = meta.get(pmcid, {})
    title = m.get("title", path.stem)
    year = m.get("year", "")

    # Skip if a valid cached summary exists. If older cached JSON used filename
    # fallback identity fields, repair those fields from reproducible PMC metadata.
    if out_path.exists():
        try:
            checklist = ManuscriptChecklist.model_validate_json(out_path.read_text())
            updates = {}
            if m.get("title") and _filename_like_title(checklist.title, pmcid):
                updates["title"] = m["title"]
            if m.get("year") and not checklist.year:
                updates["year"] = m["year"]
            if updates:
                checklist = checklist.model_copy(update=updates)
                out_path.write_text(checklist.model_dump_json(indent=2))
            return PaperResult(pmcid=pmcid, status="skipped", checklist=checklist)
        except Exception:
            pass  # cached file invalid → re-summarize

    try:
        text, src_fmt = extract(path)
        if chunked:
            checklist = summarize_chunked(
                text=text, pmcid=pmcid, title=title, year=year,
                source_format=src_fmt, client=client, model=model,
            )
        else:
            try:
                checklist = summarize_text(
                    text=text, pmcid=pmcid, title=title, year=year,
                    source_format=src_fmt, client=client, model=model,
                )
            except Exception as single_err:
                if not recover:
                    raise
                checklist = summarize_chunked(
                    text=text, pmcid=pmcid, title=title, year=year,
                    source_format=src_fmt, client=client, model=model,
                )
        out_path.write_text(checklist.model_dump_json(indent=2))
        return PaperResult(pmcid=pmcid, status="ok", checklist=checklist)
    except Exception:
        return PaperResult(pmcid=pmcid, status="failed", checklist=None)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pmcid", help="Summarize a single paper by PMCID.")
    ap.add_argument("--limit", type=int, help="Only process the first N papers.")
    ap.add_argument("--force", action="store_true",
                    help="Re-summarize even if a cached JSON exists.")
    ap.add_argument("--recover", action="store_true",
                    help="Only process papers without a summary yet; fall back to "
                         "chunked (dissect-into-two) extraction when single-shot fails.")
    ap.add_argument("--chunked", action="store_true",
                    help="Always use chunked extraction (dissect into chunks + merge).")
    ap.add_argument("--workers", type=int, default=1,
                    help="Number of concurrent LLM workers (default 1 = serial). "
                         "The OpenAI client is thread-safe; useful when the "
                         "bottleneck is network-bound LLM calls.")
    args = ap.parse_args(argv)

    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    meta = load_metadata()

    if args.recover:
        files = find_failed()
    elif args.pmcid:
        pmcid = args.pmcid.upper()
        if not pmcid.startswith("PMC"):
            pmcid = f"PMC{pmcid}"
        files = [f for f in discover_files() if pmcid_from_filename(f) == pmcid]
        if not files:
            print(f"No downloaded file found for {pmcid}", file=sys.stderr)
            return 1
    else:
        files = discover_files()
    if args.limit:
        files = files[: args.limit]

    if not files:
        print("No papers found under papers/. Run `make download` first.", file=sys.stderr)
        return 1

    pmcids = [pmcid_from_filename(f) for f in files]
    missing_meta = [p for p in pmcids if not meta.get(p, {}).get("title") or not meta.get(p, {}).get("year")]
    if missing_meta:
        fetched = _fetch_missing_metadata(missing_meta)
        meta.update(fetched)
        if fetched:
            print(f"Fetched missing PMC metadata: {len(fetched)}")

    try:
        client, model = get_client()
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        return 2

    print(f"Model: {model}")
    print(f"Papers to summarize: {len(files)}")
    print(f"Workers: {args.workers}")
    print("=" * 60)

    ok = failed = skipped = 0
    checklists: list[ManuscriptChecklist] = []

    def _handle(result: PaperResult) -> None:
        nonlocal ok, failed, skipped
        if result.status == "skipped":
            skipped += 1
            if result.checklist:
                checklists.append(result.checklist)
            print(f"  {result.pmcid}  [skip — cached]")
        elif result.status == "ok":
            ok += 1
            checklists.append(result.checklist)
            ehr = "EHR" if result.checklist.ehr_used else "no-EHR"
            print(f"  {result.pmcid}  ✓ {ehr} | "
                  f"{len(result.checklist.pathologies_diseases)} disease(s) | "
                  f"{len(result.checklist.captured_features)} feature(s)")
        else:
            failed += 1
            print(f"  {result.pmcid}  ✗ failed")

    if args.workers <= 1:
        for idx, path in enumerate(files, 1):
            result = _process_one(
                path=path, meta=meta, client=client, model=model,
                chunked=args.chunked, recover=args.recover,
                summary_dir=SUMMARY_DIR,
            )
            _handle(result)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(
                    _process_one, path=p, meta=meta, client=client, model=model,
                    chunked=args.chunked, recover=args.recover,
                    summary_dir=SUMMARY_DIR,
                ): p for p in files
            }
            for fut in as_completed(futures):
                _handle(fut.result())

    # ── combined file ───────────────────────────────────────────────────
    checklists.sort(key=lambda c: (c.year, c.pmcid))
    batch = SummaryBatch(n=len(checklists), model=model, summaries=checklists)
    COMBINED_PATH.write_text(batch.model_dump_json(indent=2))
    n_ehr = sum(1 for c in checklists if c.ehr_used)

    print("=" * 60)
    print(f"  Summarized : {ok}")
    print(f"  Skipped    : {skipped}  (cached)")
    print(f"  Failed     : {failed}")
    print(f"  Total in combined file: {len(checklists)}")
    print(f"  EHR-based  : {n_ehr} / {len(checklists)}")
    print(f"  Per-paper  : {SUMMARY_DIR}/")
    print(f"  Combined   : {COMBINED_PATH}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
