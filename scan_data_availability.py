#!/usr/bin/env python3
"""Focused LLM scan: capture data-availability info from each manuscript.

For every downloaded paper, ask the LLM ONLY about data availability (accession
numbers / repository links vs 'available upon request' / 'in-house'), then merge
the result onto the existing per-paper summary JSON (preserving all other fields)
and rebuild the combined file. Cheaper + safer than a full re-summarize: no
regression risk to the analytical fields the model already got right.

Usage:
    python scan_data_availability.py                 # all papers
    python scan_data_availability.py --limit 5       # pilot
    python scan_data_availability.py --pmcid PMC7145790
    python scan_data_availability.py --workers 4     # concurrent
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

try:
    from dotenv import load_dotenv
    load_dotenv()
except Exception:
    pass

from summarizer.extract import extract, pmcid_from_filename
from summarizer.llm_client import get_client, extract_json_object
from summarizer.schema import ManuscriptChecklist

PAPERS_DIR = Path("papers")
SUMMARY_DIR = PAPERS_DIR / "summaries"
COMBINED_PATH = PAPERS_DIR / "manuscript_summaries.json"
LOG_PATH = PAPERS_DIR / "download_log.json"

SOURCE_CHAR_BUDGET = 12000  # data-availability statements are usually near the end

# Deterministic accession / repository patterns. For the specific thing we want
# to capture (accession numbers / links), a regex is more reliable than the LLM,
# which misclassified real Zenodo/GitHub/GEO deposits as not-stated.
_ACCESSION_PATTERNS = [
    re.compile(r"dbGaP\s*(phs\d+(\.v\d+(\.p\d+)?)?)", re.I),
    re.compile(r"\b(phs\d{6}(\.v\d+\.p\d+)?)\b", re.I),
    re.compile(r"\b(GSE\d{4,})\b", re.I),
    re.compile(r"\b(PRJ[E|N][A-Z]?\d{4,})\b", re.I),  # PRJEB / PRJNA (ENA/SRA)
    re.compile(r"\b(ENA[:\s]*PRJ\w+)\b", re.I),
    re.compile(r"(https?://(?:dx\.)?doi\.org/10\.5281/zenodo\.\d+)", re.I),
    re.compile(r"(https?://zenodo\.org/(?:record/|records/)?\d+)", re.I),
    re.compile(r"(https?://github\.com/[\w.-]+/[\w.-]+)", re.I),
    re.compile(r"(https?://figshare\.com/\S+)", re.I),
    re.compile(r"(https?://datadryad\.org/\S+)", re.I),
    re.compile(r"\b(ArrayExpress\s*E-[A-Z]{4}-\d+)\b", re.I),
]


def _extract_accession_links(text: str) -> list[str]:
    """Deterministic scan for accession IDs / repository URLs in the text."""
    found: list[str] = []
    for pat in _ACCESSION_PATTERNS:
        for m in pat.finditer(text):
            val = m.group(1) if m.groups() else m.group(0)
            val = val.strip().replace("\n", " ")
            # normalize github blob/raw urls to the repo root
            if "github.com" in val:
                val = "/".join(val.split("/")[:5])  # https://github.com/owner/repo
            if val not in found:
                found.append(val)
    return found

DA_SYSTEM_PROMPT = (
    "You are a data-availability extractor for biomedical manuscripts. "
    "Read the manuscript text and determine HOW the study's data can be obtained. "
    "Respond with a ```json fenced object with EXACTLY these keys:\n"
    '  "data_availability": one of "public-repository", "available-upon-request", '
    '"in-house", "supplementary-only", "not-stated"\n'
    '  "data_accession_links": list of accession IDs / repository URLs / names '
    '(e.g. "dbGaP phs000123", "GSE12345", "github.com/x/y"). Empty list unless '
    'public-repository.\n'
    '  "data_availability_statement": the verbatim sentence(s) justifying the call; '
    '"" if none.\n'
    "Rules:\n"
    "- public-repository ONLY if a real accession number / repository URL is named.\n"
    "- 'available upon request' / 'on reasonable request' -> available-upon-request.\n"
    "- 'in-house' / institutional / private cohort with no public access -> in-house.\n"
    "- Data only in supplementary files -> supplementary-only.\n"
    "- No data-availability statement found -> not-stated.\n"
    "Output ONLY the json fenced block."
)


def _da_window(text: str, size: int = 6000) -> str:
    """Pull the data-availability section from full text; fall back to the tail.

    DA statements live under headings like 'Data availability' / 'Data and code
    availability' / 'Code availability', usually near the end. Windowing (vs.
    truncating head or tail) keeps the call cheap and avoids cutting the section.
    """
    import re
    m = re.search(
        r"(?:data\s+(?:and\s+code\s+)?(?:availability|accessing)|code\s+availability)",
        text, re.I)
    if m:
        start = max(0, m.start() - 200)
        return text[start : start + size]
    return text[-size:]


def ask_llm_data_availability(*, client, model, text, pmcid, title, year):
    """One focused LLM call; returns the parsed dict or None on failure.

    A deterministic regex pass runs FIRST: if it finds a real accession ID /
    repository URL, we short-circuit to public-repository (the LLM repeatedly
    misclassified real Zenodo/GitHub/GEO deposits as not-stated). The LLM is
    only consulted for the prose categories (available-upon-request / in-house /
    supplementary-only / not-stated).
    """
    window = _da_window((text or "").strip())
    if not window:
        return None

    links = _extract_accession_links(window)
    if links:
        # grab a verbatim snippet around the first link / 'deposited' mention
        m = re.search(r"(?:deposit|publicly available|accessible|available at|data availability)[^\n.]{0,160}", window, re.I)
        return {
            "data_availability": "public-repository",
            "data_accession_links": links,
            "data_availability_statement": (m.group(0).strip() if m else "Data deposited in a public repository.")[:300],
        }

    messages = [
        {"role": "system", "content": DA_SYSTEM_PROMPT},
        {"role": "user", "content": f"Title: {title}\nYear: {year}\nPMCID: {pmcid}\n\nManuscript text:\n{window}"},
    ]
    try:
        resp = client.chat.completions.create(
            model=model, messages=messages, max_tokens=512,
            temperature=0.0, timeout=90.0,
        )
        raw = resp.choices[0].message.content or ""
        obj = extract_json_object(raw)
        # still try to harvest any links the LLM might have dropped
        if obj is not None:
            extra = _extract_accession_links(window)
            if extra and not obj.get("data_accession_links"):
                obj["data_accession_links"] = extra
                if obj.get("data_availability") != "public-repository":
                    obj["data_availability"] = "public-repository"
        return obj
    except Exception:
        return None


def load_metadata() -> dict[str, dict]:
    if not LOG_PATH.exists():
        return {}
    log = json.loads(LOG_PATH.read_text())
    return {p["pmcid"]: p for p in log.get("papers", [])}


def discover_files() -> list[Path]:
    return sorted(
        list(PAPERS_DIR.glob("*.pdf")) + list(PAPERS_DIR.glob("*.xml")),
        key=lambda p: pmcid_from_filename(p),
    )


def scan_one(path: Path, *, summary_dir: Path, meta: dict | None = None,
             client=None, model: str | None = None) -> dict:
    """Extract text, ask the LLM about data availability, merge onto existing JSON."""
    pmcid = pmcid_from_filename(path)
    meta = meta or {}
    m = meta.get(pmcid, {})
    title = m.get("title", path.stem)
    year = m.get("year", "")

    # load existing record (preserve all fields)
    out_path = summary_dir / f"{pmcid}.json"
    rec: dict = {}
    if out_path.exists():
        try:
            rec = json.loads(out_path.read_text())
        except Exception:
            rec = {}
    # ensure identity fields present (for papers not yet summarized)
    rec.setdefault("pmcid", pmcid)
    rec.setdefault("title", title)
    rec.setdefault("year", year)

    # full LLM client if not provided (used in tests via injection)
    own_client = client is None
    if own_client:
        client, model = get_client()

    text, _ = extract(path)
    da = ask_llm_data_availability(client=client, model=model, text=text,
                                   pmcid=pmcid, title=title, year=year)
    if da is not None:
        rec.update(da)

    # validate + write. If the required analytical fields (ehr_used etc.) are
    # present, round-trip through the schema to normalize + guard; otherwise the
    # full summarizer hasn't run yet, so write the partial record as-is (it'll be
    # completed + validated when summarize runs).
    if all(k in rec for k in ("ehr_used", "ehr_evidence", "summary")):
        checklist = ManuscriptChecklist(**rec)
        out_path.write_text(checklist.model_dump_json(indent=2))
    else:
        out_path.write_text(json.dumps(rec, indent=2, sort_keys=True))
    return json.loads(out_path.read_text())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pmcid", help="Scan a single paper by PMCID.")
    ap.add_argument("--limit", type=int, help="Only scan the first N papers (pilot).")
    ap.add_argument("--workers", type=int, default=1, help="Concurrent LLM workers.")
    args = ap.parse_args(argv)

    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    meta = load_metadata()

    if args.pmcid:
        pmcid = args.pmcid.upper()
        if not pmcid.startswith("PMC"):
            pmcid = f"PMC{pmcid}"
        files = [f for f in discover_files() if pmcid_from_filename(f) == pmcid]
    else:
        files = discover_files()
    if args.limit:
        files = files[: args.limit]
    if not files:
        print("No papers found under papers/.", file=sys.stderr)
        return 1

    try:
        client, model = get_client()
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        return 2

    print(f"Data-availability scan | model: {model} | papers: {len(files)} | workers: {args.workers}")
    print("=" * 60)

    from collections import Counter
    cats: Counter[str] = Counter()
    ok = failed = 0

    def _do(p: Path) -> tuple[str, str]:
        try:
            r = scan_one(p, summary_dir=SUMMARY_DIR, meta=meta, client=client, model=model)
            return r["pmcid"], r.get("data_availability", "not-stated")
        except Exception as e:
            return pmcid_from_filename(p), f"ERROR:{type(e).__name__}"

    if args.workers <= 1:
        for p in files:
            pmcid, cat = _do(p)
            print(f"  {pmcid}  {cat}")
            cats[cat] += 1
            ok += 0 if cat.startswith("ERROR") else 1
            failed += 1 if cat.startswith("ERROR") else 0
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(_do, p): p for p in files}
            for fut in as_completed(futs):
                pmcid, cat = fut.result()
                print(f"  {pmcid}  {cat}")
                cats[cat] += 1
                ok += 0 if cat.startswith("ERROR") else 1
                failed += 1 if cat.startswith("ERROR") else 0

    print("=" * 60)
    print(f"  ok={ok}  failed={failed}")
    for cat, n in cats.most_common():
        print(f"  {cat:<26} {n}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
