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
from typing import Literal
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

try:
    from dotenv import load_dotenv
    load_dotenv()
except Exception:
    pass

from summarizer.extract import extract, pmcid_from_filename
from summarizer.llm_client import get_client, extract_json_object, MAX_OUTPUT_TOKENS
from summarizer.schema import ManuscriptChecklist
from pydantic import BaseModel, Field, field_validator, model_validator

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
    re.compile(r"(https?://github\.com/[\w.-]+/[\w.-]+(?:-[\w.-]+)*)", re.I),
    re.compile(r"(https?://figshare\.com/\S+)", re.I),
    re.compile(r"(https?://datadryad\.org/\S+)", re.I),
    re.compile(r"\b(ArrayExpress\s*E-[A-Z]{4}-\d+)\b", re.I),
]


def _normalize_url_text(s: str) -> str:
    """Repair URLs broken across PDF line wraps before matching.

    pypdf often splits a URL like 'early-life-sugar-rationing-respiratory-health'
    across lines as 'early-life-sugar -\nrationing-respiratory-health' (a hyphen
    that belongs to the URL, surrounded by spaces and a newline). Collapse those
    line-wrap artifacts so the URL matches as one token.
    """
    # 'word - \n word' inside a URL -> 'word-word' (hyphen was a URL char; drop only spaces+newline)
    s = re.sub(r"([\w])\s*-\s*\n\s*([\w])", r"\1-\2", s)
    # 'word-\n word' (hyphen immediately before newline) -> 'wordword' (hyphen was a wrap artifact)
    s = re.sub(r"-\s*\n\s*", "", s)
    return s


def _extract_accession_links(text: str) -> list[str]:
    """Deterministic scan for accession IDs / repository URLs in the text."""
    # Repair line-wrapped URLs across the WHOLE text first (cheap, one pass), so
    # patterns match full URLs that pypdf split across lines.
    text = _normalize_url_text(text)
    found: list[str] = []
    for pat in _ACCESSION_PATTERNS:
        for m in pat.finditer(text):
            val = m.group(1) if m.groups() else m.group(0)
            val = val.strip()
            # drop a trailing blob/raw/tree/commit path segment -> repo root, but keep
            # the full repo name (incl. hyphenated suffixes like 'early-life-sugar-...').
            if "github.com" in val:
                parts = val.split("/")
                # https://github.com/owner/repo  -> keep first 5 segments (scheme+host+owner+repo)
                # but strip trailing blob/raw/tree/etc. paths
                if len(parts) > 5 and parts[5] in {"blob", "raw", "tree", "commit", "releases", "archive"}:
                    parts = parts[:5]
                val = "/".join(parts)
            if val not in found:
                found.append(val)
    return found


_DA_VALUES = {
    "public-repository",
    "available-upon-request",
    "in-house",
    "supplementary-only",
    "not-stated",
}


def _normalize_da_value(value: str) -> str | None:
    val = value.strip().lower().replace("_", "-").replace(" ", "-")
    return val if val in _DA_VALUES else None


def _extract_markdown_da_response(raw: str) -> dict | None:
    """Recover data-availability fields from Gemma's markdown prose output.

    The prompt asks for JSON, but the deployed model sometimes emits bullet
    prose like ``* `data_availability`: "public-repository"``. Treat this as a
    salvage path only; proper JSON remains the preferred contract.
    """
    if not raw:
        return None

    category = None
    for m in re.finditer(r"`?data_availability`?\s*:\s*[\"'`* ]*([A-Za-z_-]+(?:\s+[A-Za-z_-]+)*)", raw):
        candidate = _normalize_da_value(m.group(1).strip(' "\'`*.'))
        if candidate is not None:
            category = candidate

    links = _extract_accession_links(raw)
    if links and category in {None, "not-stated"}:
        category = "public-repository"

    if category is None:
        return None

    if category != "public-repository":
        links = []

    quoted = [
        s.strip()
        for s in re.findall(r'"([^"]*)"', raw)
        if 20 <= len(s) <= 500
        and re.search(r"data|deposit|repository|available|accession|zenodo|github|supplement", s, re.I)
    ]
    statement = " ".join(quoted[:3])

    return {
        "data_availability": category,
        "data_accession_links": links,
        "data_availability_statement": statement,
    }


class DataAvailabilityResult(BaseModel):
    """Pydantic-validated data-availability extraction from one LLM response."""

    data_availability: Literal[
        "public-repository", "available-upon-request", "in-house",
        "supplementary-only", "not-stated"
    ] = Field(default="not-stated",
              description="How the study's data can be obtained.")
    data_accession_links: list[str] = Field(
        default_factory=list,
        description="Accession IDs / repository URLs / names.",
    )
    data_availability_statement: str = Field(
        default="", description="Verbatim sentence(s) justifying the call.")
    data_mentioned: bool = Field(default=False)
    data_available: bool = Field(default=False)
    data_open_access: bool = Field(default=False)

    @model_validator(mode="after")
    def _derive_data_access_flags(self):
        da = self.data_availability
        has_statement = bool(self.data_availability_statement.strip())
        has_links = bool(self.data_accession_links)
        self.data_mentioned = da != "not-stated" or has_statement or has_links
        self.data_available = da in {"public-repository", "supplementary-only", "available-upon-request"} or has_links
        self.data_open_access = da in {"public-repository", "supplementary-only"} or has_links
        return self

    @field_validator("data_availability", mode="before")
    @classmethod
    def _norm_cat(cls, v):
        if isinstance(v, str):
            v2 = v.strip().lower().replace("_", "-").replace(" ", "-")
            valid = {"public-repository", "available-upon-request", "in-house",
                     "supplementary-only", "not-stated"}
            return v2 if v2 in valid else "not-stated"
        return v

    @field_validator("data_accession_links", mode="before")
    @classmethod
    def _coerce_links(cls, v):
        if isinstance(v, str):
            return [v] if v.strip() else []
        if isinstance(v, list):
            return [str(x).strip() for x in v if str(x).strip()]
        return v


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
    '  "data_mentioned": true if any data availability/sharing/access statement exists.\n'
    '  "data_available": true if repository, supplement, or request path exists.\n'
    '  "data_open_access": true only if accessible without permission/application/authors.\n'
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

    The LLM owns the job: it classifies data_availability and emits
    data_accession_links, validated through DataAvailabilityResult (pydantic).
    The deterministic regex is only a SAFETY NET — it runs AFTER the LLM to
    supplement any real accession/URL the model dropped, and to upgrade a
    not-stated/available-upon-request call to public-repository when a real
    accession pattern is plainly in the text. It never short-circuits the LLM.
    """
    window = _da_window((text or "").strip())
    if not window:
        return None

    messages = [
        {"role": "system", "content": DA_SYSTEM_PROMPT},
        {"role": "user", "content": f"Title: {title}\nYear: {year}\nPMCID: {pmcid}\n\nManuscript text:\n{window}"},
    ]
    result = None
    try:
        resp = client.chat.completions.create(
            model=model, messages=messages, max_tokens=MAX_OUTPUT_TOKENS,
            temperature=0.0, timeout=180.0,
        )
        raw = resp.choices[0].message.content or ""
        obj = extract_json_object(raw)
        if obj is None:
            obj = _extract_markdown_da_response(raw)
        if obj is not None:
            # validate through pydantic — enforces the enum + coerces links,
            # and rejects garbage the model might emit.
            try:
                da = DataAvailabilityResult.model_validate(obj)
                result = da.model_dump()
            except Exception:
                result = None
    except Exception:
        result = None

    # SAFETY NET: harvest any real accession/URL the LLM dropped. If the model
    # missed a link that's plainly in the text, add it; and if it said
    # not-stated/available-upon-request while a real accession is present,
    # upgrade to public-repository. (The LLM — now with 32k tokens — usually
    # gets this right; this only catches its occasional misses.)
    extra = _extract_accession_links(window)
    if extra:
        if result is None:
            m = re.search(r"(?:deposit|publicly available|accessible|available at|data availability)[^\n.]{0,160}", window, re.I)
            result = {
                "data_availability": "public-repository",
                "data_accession_links": extra,
                "data_availability_statement": (m.group(0).strip() if m else "Data deposited in a public repository.")[:300],
            }
        else:
            existing = set(result.get("data_accession_links") or [])
            result["data_accession_links"] = list(existing | set(extra))
            if result.get("data_availability") != "public-repository":
                result["data_availability"] = "public-repository"
    return result


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

    total = len(files)
    completed = 0
    started = 0
    lock = __import__("threading").Lock()

    def _print_start(pmcid: str) -> None:
        nonlocal started
        with lock:
            started += 1
            print(f"  start {started:>3}/{total} {pmcid} ...", flush=True)

    def _print_progress(pmcid: str, cat: str) -> None:
        nonlocal completed
        with lock:
            completed += 1
            print(f"  done  [{completed:>3}/{total}] {pmcid:<13} {cat}", flush=True)
            cats[cat] += 1

    def _do(p: Path) -> tuple[str, str]:
        _print_start(pmcid_from_filename(p))
        try:
            r = scan_one(p, summary_dir=SUMMARY_DIR, meta=meta, client=client, model=model)
            return r["pmcid"], r.get("data_availability", "not-stated")
        except Exception as e:
            return pmcid_from_filename(p), f"ERROR:{type(e).__name__}"

    if args.workers <= 1:
        for p in files:
            pmcid, cat = _do(p)
            _print_progress(pmcid, cat)
            ok += 0 if cat.startswith("ERROR") else 1
            failed += 1 if cat.startswith("ERROR") else 0
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(_do, p): p for p in files}
            for fut in as_completed(futs):
                pmcid, cat = fut.result()
                _print_progress(pmcid, cat)
                ok += 0 if cat.startswith("ERROR") else 1
                failed += 1 if cat.startswith("ERROR") else 0

    print("=" * 60)
    print(f"  ok={ok}  failed={failed}")
    for cat, n in cats.most_common():
        print(f"  {cat:<26} {n}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
