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
import subprocess
import sys
import warnings
from collections import Counter
from pathlib import Path

from summarizer.schema import ManuscriptChecklist, SummaryBatch

REPO_ROOT = Path(__file__).resolve().parent
PAPERS_DIR = REPO_ROOT / "papers"
SUMMARY_DIR = PAPERS_DIR / "summaries"
COMBINED_PATH = PAPERS_DIR / "manuscript_summaries.json"
DOWNLOAD_LOG = PAPERS_DIR / "download_log.json"

# Default summarizer mode: resume-only (skip cached papers, chunked recovery on
# failure) — the safe, idempotent choice for incremental materialization.
SUMMARIZE_ARGS = ["--recover"]


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


def summarize_papers() -> int:
    """Run the summarizer over missing papers (resume-only, chunked recovery)."""
    return _run([sys.executable, "-m", "summarizer.run"] + SUMMARIZE_ARGS,
                label="summarizer.run")


def build_results() -> int:
    """Run build_results.py -> results/ (SUMMARY.md, checklist.md, combined copy)."""
    return _run([sys.executable, str(REPO_ROOT / "build_results.py")],
                label="build_results")


def _most_common_model(rows: list[ManuscriptChecklist]) -> str:
    counts = Counter(r.model for r in rows if str(r.model).strip())
    return counts.most_common(1)[0][0] if counts else ""


def rebuild_combined(
    *,
    summary_dir: Path = SUMMARY_DIR,
    out_path: Path = COMBINED_PATH,
) -> dict:
    """Rebuild the combined ``SummaryBatch`` from per-paper summary files.

    Reads every ``<pmcid>.json`` under ``summary_dir``, validates each as a
    ``ManuscriptChecklist`` (invalid files are skipped with a warning), sorts by
    ``(year, pmcid)``, derives the batch ``model`` from the records, and writes
    the combined JSON to ``out_path``. Pure-Python — no LLM, no TinyDB.
    """
    checklists: list[ManuscriptChecklist] = []
    for f in sorted(summary_dir.glob("*.json")):
        try:
            checklists.append(
                ManuscriptChecklist.model_validate_json(f.read_text())
            )
        except (json.JSONDecodeError, ValueError) as e:
            warnings.warn(f"skipping {f.name}: {type(e).__name__}", UserWarning,
                          stacklevel=2)

    checklists.sort(key=lambda c: (str(c.year), c.pmcid))
    batch = SummaryBatch(n=len(checklists), model=_most_common_model(checklists),
                         summaries=checklists)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(batch.model_dump_json(indent=2))
    return json.loads(out_path.read_text())