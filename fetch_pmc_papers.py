#!/usr/bin/env python3
"""
Fetch and download open-access PMC papers related to EWAS using EHR data.
Uses NCBI E-utilities + PMC OA API for reliable PDF resolution.
FTP URLs are converted to HTTPS; tar.gz packages are extracted to retrieve the PDF.
"""

import io
import json
import os
import re
import tarfile
import time
import requests
import xml.etree.ElementTree as ET
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
OUTPUT_DIR = Path("papers")
NCBI_BASE  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PMC_OA_API = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"
DELAY      = 0.4   # seconds between API calls (NCBI limit ≈ 3/sec without key)

# Filters applied to every query:
#   - open access[filter]      : only OA full-text articles
#   - NOT Review[Publication Type] : exclude review articles
#   - NOT systematic[Title]    : exclude systematic reviews
#   - NOT meta-analysis[Title] : exclude meta-analyses
#   - NOT bibliometric[Title]  : exclude bibliometric surveys
#   - NOT protocol[Title]      : exclude study protocols
# "environment-wide" is spelled out to avoid confusion with
# "epigenome-wide" (also abbreviated EWAS in the methylation literature).
_FILTERS = (
    'open access[filter] '
    'NOT Review[Publication Type] '
    'NOT systematic[Title] '
    'NOT meta-analysis[Title] '
    'NOT bibliometric[Title] '
    'NOT protocol[Title]'
)

SEARCH_QUERIES = [
    f'"environment-wide association" electronic health record {_FILTERS}',
    f'"environment-wide association" EHR cohort {_FILTERS}',
    f'"environment-wide association" biobank phenome {_FILTERS}',
    f'"exposome-wide association" electronic health record {_FILTERS}',
    f'"exposome-wide association" EHR {_FILTERS}',
    f'"environment-wide association" clinical data epidemiology {_FILTERS}',
    f'"environment-wide association study" {_FILTERS}',
]

MAX_PER_QUERY = 30


# ── Helpers ───────────────────────────────────────────────────────────────────

def ftp_to_https(url: str) -> str:
    """NCBI FTP and HTTPS share the same path — swap the scheme."""
    return url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov", 1)


def sanitize(text: str, maxlen: int = 80) -> str:
    text = re.sub(r'[^\w\s-]', '', text)
    text = re.sub(r'\s+', '_', text.strip())
    return text[:maxlen]


def search_pmc(query: str) -> list:
    try:
        r = requests.get(f"{NCBI_BASE}/esearch.fcgi", timeout=15, params={
            "db": "pmc", "term": query, "retmax": MAX_PER_QUERY,
            "retmode": "json", "sort": "relevance",
        })
        r.raise_for_status()
        data = r.json().get("esearchresult", {})
        if "count" not in data:
            print(f"  API error: {r.text[:200]}")
            return []
        print(f"  hits: {data['count']:>5}  |  retrieved: {len(data['idlist'])}")
        return data["idlist"]
    except Exception as e:
        print(f"  Search failed: {e}")
        return []


def fetch_summaries(ids: list) -> dict:
    if not ids:
        return {}
    r = requests.get(f"{NCBI_BASE}/esummary.fcgi", timeout=15,
                     params={"db": "pmc", "id": ",".join(ids), "retmode": "json"})
    r.raise_for_status()
    return r.json().get("result", {})


def get_oa_links(pmcid: str) -> tuple:
    """
    Returns (pdf_url, tgz_url) from the PMC OA API.
    Both are converted from ftp:// → https://.
    """
    r = requests.get(PMC_OA_API, params={"id": pmcid}, timeout=15)
    if r.status_code != 200:
        return None, None
    try:
        root = ET.fromstring(r.text)
    except ET.ParseError:
        return None, None

    pdf_url = tgz_url = None
    for link in root.iter("link"):
        fmt  = link.attrib.get("format", "")
        href = ftp_to_https(link.attrib.get("href", ""))
        if fmt == "pdf":
            pdf_url = href
        elif fmt == "tgz":
            tgz_url = href
    return pdf_url, tgz_url


def download_bytes(url: str) -> bytes | None:
    """Download raw bytes from a URL; return None on failure."""
    try:
        r = requests.get(url, timeout=60,
                         headers={"User-Agent": "Mozilla/5.0 (academic research)"})
        if r.status_code == 200 and len(r.content) > 5_000:
            return r.content
        print(f"    Bad response: {r.status_code}, {len(r.content)} bytes")
    except Exception as e:
        print(f"    Download error: {e}")
    return None


SUPP_PATTERN = re.compile(
    r'(suppl?e?m?e?n?t?|supp?\d|_s\d+[\._]|[-_]s\d+\.pdf$'
    r'|app\d+|appendix|fig(ure)?\d*|table\d*)',
    re.IGNORECASE
)

def pdf_from_tgz(data: bytes) -> bytes | None:
    """
    Extract the main article PDF from a tar.gz archive.
    Strategy:
      1. Exclude files whose basenames match supplementary/appendix patterns.
      2. From the remaining candidates, pick the largest (= main article).
      3. If no main candidates exist (archive is supplementary-only), return None.
    """
    try:
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
            pdf_members = [m for m in tf.getmembers() if m.name.endswith(".pdf")]
            if not pdf_members:
                return None

            main_candidates = [m for m in pdf_members
                               if not SUPP_PATTERN.search(os.path.basename(m.name))]

            if not main_candidates:
                names = [os.path.basename(m.name) for m in pdf_members]
                print(f"    Skipped: tarball contains only supplementary/appendix PDFs")
                print(f"    Files: {names}")
                return None

            best = max(main_candidates, key=lambda m: m.size)
            print(f"    Extracted main article: {os.path.basename(best.name)}")
            return tf.extractfile(best).read()
    except Exception as e:
        print(f"    tar.gz extraction error: {e}")
    return None


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    log_path = OUTPUT_DIR / "download_log.json"
    log = json.loads(log_path.read_text()) if log_path.exists() else {}
    already_done = set(log.get("downloaded", []))

    # ── 1. Search ──────────────────────────────────────────────────────────
    print("=" * 65)
    print("  SEARCHING PubMed Central (open-access filter)")
    print("=" * 65)
    all_ids: set = set()
    for q in SEARCH_QUERIES:
        label = q.replace(" open access[filter]", "")
        print(f"\n  Query: {label!r}")
        all_ids.update(search_pmc(q))
        time.sleep(DELAY)
    print(f"\nUnique PMC IDs: {len(all_ids)}")

    # ── 2. Metadata ────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("  FETCHING METADATA")
    print("=" * 65)
    id_list = list(all_ids)
    summaries: dict = {}
    for i in range(0, len(id_list), 20):
        res = fetch_summaries(id_list[i:i + 20])
        for uid in res.get("uids", []):
            summaries[uid] = res[uid]
        time.sleep(DELAY)

    # ── 3. Filter out reviews / non-primary articles ───────────────────────
    REVIEW_TITLE_KW = re.compile(
        r'\b(review|systematic review|meta.analysis|bibliometric|'
        r'narrative review|scoping review|overview|commentary|'
        r'perspective|editorial|protocol|roadmap)\b',
        re.IGNORECASE
    )
    EPIGENOME_KW = re.compile(
        r'\b(epigenome.wide|methylation|DNA methyl|CpG|histone)\b',
        re.IGNORECASE
    )

    # ── 4. Print table ─────────────────────────────────────────────────────
    print(f"\n{'#':<4} {'PMCID':<14} {'Yr':<5} {'Flag':<9} {'Title':<52} Journal")
    print("-" * 115)
    papers = []       # primary research papers
    skipped_reviews = []
    for i, (uid, item) in enumerate(summaries.items(), 1):
        title   = item.get("title",   "N/A")
        journal = item.get("source",  "N/A")
        year    = item.get("pubdate", "")[:4]
        authors = ", ".join(a.get("name","") for a in item.get("authors", [])[:2])
        pub_type = " ".join(item.get("pubtype", []))

        is_review   = bool(REVIEW_TITLE_KW.search(title)) or "Review" in pub_type
        is_epigenome = bool(EPIGENOME_KW.search(title)) and \
                       "environment-wide" not in title.lower() and \
                       "exposome" not in title.lower()

        if is_review:
            flag = "[REVIEW]"
            skipped_reviews.append({"pmcid": f"PMC{uid}", "title": title, "reason": "review"})
        elif is_epigenome:
            flag = "[EPIOME]"
            skipped_reviews.append({"pmcid": f"PMC{uid}", "title": title, "reason": "epigenome-wide (not environment-wide)"})
        else:
            flag = "[OK]"
            papers.append((uid, title, journal, year, authors))

        print(f"{i:<4} PMC{uid:<11} {year:<5} {flag:<9} {title[:51]:<52} {journal[:25]}")

    print(f"\n  Primary research papers : {len(papers)}")
    print(f"  Excluded (reviews)      : {sum(1 for s in skipped_reviews if s['reason']=='review')}")
    print(f"  Excluded (epigenome)    : {sum(1 for s in skipped_reviews if 'epigenome' in s['reason'])}")
    log["excluded"] = skipped_reviews

    # ── 5. Download ────────────────────────────────────────────────────────
    print(f"\n{'=' * 65}")
    print(f"  DOWNLOADING PDFs  →  ./{OUTPUT_DIR}/")
    print("=" * 65 + "\n")

    downloaded = skipped = no_oa = failed = 0

    for idx, (uid, title, journal, year, authors) in enumerate(papers, 1):
        pmcid    = f"PMC{uid}"
        out_stem = sanitize(f"{year}_{pmcid}_{title}")
        pdf_path = OUTPUT_DIR / f"{out_stem}.pdf"
        pfx      = f"  [{idx:>2}/{len(papers)}] {pmcid} ({year})"

        if uid in already_done or pdf_path.exists():
            print(f"{pfx}  [SKIP]")
            skipped += 1
            continue

        print(f"{pfx}  {title[:55]}")

        pdf_url, tgz_url = get_oa_links(pmcid)
        time.sleep(DELAY)

        if not pdf_url and not tgz_url:
            print(f"           → Not in OA subset")
            log.setdefault("no_oa", []).append({"pmcid": pmcid, "title": title})
            no_oa += 1
            continue

        pdf_data: bytes | None = None

        if pdf_url:
            print(f"           → PDF: {pdf_url[:78]}")
            pdf_data = download_bytes(pdf_url)

        if pdf_data is None and tgz_url:
            print(f"           → TAR: {tgz_url[:78]}")
            raw = download_bytes(tgz_url)
            if raw:
                pdf_data = pdf_from_tgz(raw)

        if pdf_data:
            pdf_path.write_bytes(pdf_data)
            size_kb = pdf_path.stat().st_size // 1024
            print(f"           ✓ {pdf_path.name}  ({size_kb} KB)")
            log.setdefault("downloaded", []).append(uid)
            already_done.add(uid)
            downloaded += 1
        else:
            print(f"           ✗ Could not retrieve PDF")
            log.setdefault("failed", []).append({"pmcid": pmcid, "title": title})
            failed += 1

        time.sleep(DELAY)

    # ── 6. Save log ────────────────────────────────────────────────────────
    log["papers"] = [
        {"pmcid": f"PMC{uid}", "title": t, "journal": j, "year": y, "authors": a}
        for uid, t, j, y, a in papers
    ]
    log_path.write_text(json.dumps(log, indent=2))

    print(f"\n{'=' * 50}")
    print(f"  Downloaded : {downloaded}")
    print(f"  Skipped    : {skipped}  (already on disk)")
    print(f"  Not OA     : {no_oa}  (subscription-only)")
    print(f"  Failed     : {failed}")
    print(f"  Log        : {log_path}")
    print(f"  Output     : {OUTPUT_DIR.resolve()}")
    print("=" * 50)


if __name__ == "__main__":
    main()
