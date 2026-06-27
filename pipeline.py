"""Dagster asset pipeline for the pediatric EWAS / EHR literature review.

A five-asset lineage that orchestrates the existing scripts via thin wrappers
(see ``pipeline_ops``):

    download_log  ->  per_paper_summaries  ->  data_availability_scan  ->
    manuscript_summaries  ->  results

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

from dagster import asset, Definitions, MetadataValue

import pipeline_ops

PAPERS_DIR = pipeline_ops.PAPERS_DIR
SUMMARY_DIR = pipeline_ops.SUMMARY_DIR
DOWNLOAD_LOG = pipeline_ops.DOWNLOAD_LOG


def summarize_data_availability_counts(summary_dir: Path = SUMMARY_DIR) -> dict:
    """Count data-availability categories and JSON parse failures in summaries."""
    counts: dict[str, int] = {}
    parse_failures = 0
    for f in summary_dir.glob("*.json"):
        try:
            rec = json.loads(f.read_text())
        except Exception:
            parse_failures += 1
            continue
        cat = rec.get("data_availability", "not-stated") or "not-stated"
        counts[cat] = counts.get(cat, 0) + 1
    return {"n": sum(counts.values()), "counts": counts, "parse_failures": parse_failures}


@asset
def download_log(context) -> dict:
    """Fetch PMC papers + the download log (fetch_pmc_papers.py)."""
    pipeline_ops.fetch_papers()
    if not DOWNLOAD_LOG.exists():
        data = {"n": 0}
    else:
        log = json.loads(DOWNLOAD_LOG.read_text())
        data = {"n": len(log.get("papers", []))}
    context.add_output_metadata({
        "papers": data["n"],
        "path": MetadataValue.path(DOWNLOAD_LOG),
    })
    return data


@asset
def per_paper_summaries(context, download_log: dict) -> dict:
    """Summarize every missing paper (summarizer.run --recover) -> papers/summaries/."""
    pipeline_ops.summarize_papers()
    n = len(list(SUMMARY_DIR.glob("*.json")))
    data = {"n": n, "downloaded": download_log.get("n", 0)}
    context.add_output_metadata({
        "summaries": n,
        "downloaded_candidates": data["downloaded"],
        "path": MetadataValue.path(SUMMARY_DIR),
    })
    return data


@asset
def data_availability_scan(context, per_paper_summaries: dict) -> dict:
    """Update per-paper summaries with data-availability fields."""
    pipeline_ops.scan_data_availability()
    data = summarize_data_availability_counts()
    data["per_paper"] = per_paper_summaries.get("n", 0)
    parse_failures = data["parse_failures"]
    counts = data["counts"]
    context.add_output_metadata({
        "summaries_scanned": data["n"],
        "per_paper_summaries": data["per_paper"],
        "parse_failures": parse_failures,
        "data_availability_counts": MetadataValue.json(counts),
    })
    return data


@asset
def manuscript_summaries(context, data_availability_scan: dict) -> dict:
    """Rebuild the combined SummaryBatch from enriched per-paper files."""
    data = pipeline_ops.rebuild_combined()
    out = {"n": data["n"], "model": data.get("model", ""),
           "per_paper": data_availability_scan.get("per_paper", data_availability_scan.get("n", 0))}
    context.add_output_metadata({
        "summaries": out["n"],
        "model": out["model"],
        "path": MetadataValue.path(pipeline_ops.COMBINED_PATH),
    })
    return out


@asset
def results(context, manuscript_summaries: dict) -> dict:
    """Export readable results/ (SUMMARY.md, checklist.md, combined copy)."""
    pipeline_ops.build_results()
    data = {"n": manuscript_summaries.get("n", 0)}
    context.add_output_metadata({
        "summaries": data["n"],
        "summary_md": MetadataValue.path(pipeline_ops.REPO_ROOT / "results" / "SUMMARY.md"),
        "checklist_md": MetadataValue.path(pipeline_ops.REPO_ROOT / "results" / "checklist.md"),
    })
    return data


defs = Definitions(assets=[download_log, per_paper_summaries,
                           data_availability_scan, manuscript_summaries,
                           results])
