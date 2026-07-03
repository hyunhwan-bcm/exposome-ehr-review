"""Thin operations layer for the Dagster asset pipeline.

Each op wraps an EXISTING script so the pipeline is orchestrated (lineage,
caching/resume, UI) without rewriting the fetcher, summarizer, or results
builder. ``rebuild_combined`` is the one pure-Python step: it derives the
combined ``SummaryBatch`` from the per-paper summary files (no LLM, no TinyDB
dependency), mirroring how the TinyDB store exports.

Everything runs with ``cwd`` at the repo root so the scripts' relative
``papers/`` paths resolve.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import warnings
from pathlib import Path

from summarizer.schema import ManuscriptChecklist, SummaryBatch
from database import Store, DEFAULT_DB_PATH

REPO_ROOT = Path(__file__).resolve().parent
PAPERS_DIR = REPO_ROOT / "papers"
SUMMARY_DIR = PAPERS_DIR / "summaries"
COMBINED_PATH = PAPERS_DIR / "manuscript_summaries.json"
DOWNLOAD_LOG = PAPERS_DIR / "download_log.json"

# Default summarizer mode: resume-only (skip cached papers, chunked recovery on
# failure) + concurrent workers (the bottleneck is network-bound LLM calls).
# Override worker counts via SUMMARIZE_WORKERS / SCAN_WORKERS (e.g. from CI
# workflow_dispatch inputs) without touching the Dagster asset graph.
SUMMARIZE_ARGS = ["--recover", "--workers", os.environ.get("SUMMARIZE_WORKERS", "4")]

# Focused data-availability pass. This updates the already-created per-paper
# JSONs with data_availability/accession fields before the combined artifact is
# rebuilt. Keeping it as a separate Dagster asset makes the lineage explicit.
DATA_AVAILABILITY_ARGS = ["--workers", os.environ.get("SCAN_WORKERS", "4")]


def _run(cmd: list[str], *, label: str) -> int:
    """Run ``cmd`` in the repo root; raise RuntimeError on non-zero exit."""
    proc = subprocess.run(cmd, cwd=str(REPO_ROOT), text=True, capture_output=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-6:]
        raise RuntimeError(f"{label} failed (exit {proc.returncode})\n" + "\n".join(tail))
    return proc.returncode


def fetch_papers() -> int:
    """Run the PMC fetch/download stage -> papers/ + download_log.json."""
    return _run([sys.executable, "-m", "fetch_pmc_papers"], label="fetch_pmc_papers")


def _limit_args() -> list[str]:
    """``--limit N`` for a pilot run, sourced from the PIPELINE_LIMIT env var."""
    limit = os.environ.get("PIPELINE_LIMIT", "").strip()
    return ["--limit", limit] if limit else []


def summarize_papers() -> int:
    """Run the summarizer over missing papers (resume-only, chunked recovery)."""
    return _run([sys.executable, "-m", "summarizer.run"] + SUMMARIZE_ARGS + _limit_args(),
                label="summarizer.run")


def scan_data_availability() -> int:
    """Run the focused data-availability scan over downloaded papers."""
    return _run([sys.executable, str(REPO_ROOT / "scan_data_availability.py")] + DATA_AVAILABILITY_ARGS + _limit_args(),
                label="scan_data_availability")


def build_results() -> int:
    """Run build_results.py -> results/ (SUMMARY.md, checklist.md, combined copy)."""
    return _run([sys.executable, str(REPO_ROOT / "build_results.py")],
                label="build_results")


def rebuild_combined(
    *,
    summary_dir: Path = SUMMARY_DIR,
    out_path: Path = COMBINED_PATH,
    db_path: Path = DEFAULT_DB_PATH,
) -> dict:
    """Rebuild the combined ``SummaryBatch`` via the TinyDB Store.

    Imports every ``<pmcid>.json`` under ``summary_dir`` into the Store (which
    validates each through ``ManuscriptChecklist``; invalid files are skipped
    with a warning), then exports the combined file from the Store. The Store
    is the source of truth — sorting and batch-model derivation happen there.
    """
    store = Store(db_path)
    try:
        for f in sorted(summary_dir.glob("*.json")):
            try:
                rec = json.loads(f.read_text())
                store.import_records([rec])
            except Exception as e:
                warnings.warn(f"skipping {f.name}: {type(e).__name__}", UserWarning,
                              stacklevel=2)
        return store.export_combined(out_path)
    finally:
        store.close()