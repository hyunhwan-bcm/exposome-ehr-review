"""Dagster asset pipeline for the pediatric EWAS / EHR literature review.

A four-asset lineage that orchestrates the existing scripts via thin wrappers
(see ``pipeline_ops``):

    download_log  ->  per_paper_summaries  ->  manuscript_summaries  ->  results

Run it:

    make dagster                       # `dagster dev -f pipeline.py` (UI + lineage)
    # or materialize the whole graph headlessly:
    make materialize

Each asset materializes by running the corresponding stage; the summarizer runs
in ``--recover`` (resume-only, chunked recovery) mode so re-materializing only
fills in the gaps. No rewrite of the fetcher, summarizer, or build_results.
"""
from __future__ import annotations

import json
from pathlib import Path

from dagster import asset, Definitions

import pipeline_ops

PAPERS_DIR = pipeline_ops.PAPERS_DIR
SUMMARY_DIR = pipeline_ops.SUMMARY_DIR
DOWNLOAD_LOG = pipeline_ops.DOWNLOAD_LOG


@asset
def download_log() -> dict:
    """Fetch PMC papers + the download log (fetch_pmc_papers.py)."""
    pipeline_ops.fetch_papers()
    if not DOWNLOAD_LOG.exists():
        return {"n": 0}
    log = json.loads(DOWNLOAD_LOG.read_text())
    return {"n": len(log.get("papers", []))}


@asset
def per_paper_summaries(download_log: dict) -> dict:
    """Summarize every missing paper (summarizer.run --recover) -> papers/summaries/."""
    pipeline_ops.summarize_papers()
    n = len(list(SUMMARY_DIR.glob("*.json")))
    return {"n": n, "downloaded": download_log.get("n", 0)}


@asset
def manuscript_summaries(per_paper_summaries: dict) -> dict:
    """Rebuild the combined SummaryBatch from per-paper files (pure python, no LLM)."""
    data = pipeline_ops.rebuild_combined()
    return {"n": data["n"], "model": data.get("model", ""),
            "per_paper": per_paper_summaries.get("n", 0)}


@asset
def results(manuscript_summaries: dict) -> dict:
    """Export readable results/ (SUMMARY.md, checklist.md, combined copy)."""
    pipeline_ops.build_results()
    return {"n": manuscript_summaries.get("n", 0)}


defs = Definitions(assets=[download_log, per_paper_summaries,
                           manuscript_summaries, results])