#!/usr/bin/env python3
"""
Fetch and download open-access PMC papers on pediatric/childhood environmental
exposure (exposome / EWAS) studies using EHR / administrative health data.

Pipeline:
  1. Search PubMed Central with pediatric-constrained EWAS/exposome + EHR queries.
  2. Fetch metadata, filter out reviews / conference abstracts / epigenome-only.
  3. Resolve a full-text PDF: NCBI PMC OA (PDF or tar.gz) -> Europe PMC fallback.

Europe PMC is used as a fallback for papers that are open access but have no
resolvable PDF link on the NCBI OA service (e.g. project-design / late-deposited
articles such as the EXPOsOMICS project paper, PMC6192011).
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
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
DELAY      = 0.4   # seconds between API calls (NCBI limit ≈ 3/sec without key)

# ── Search strategy ──────────────────────────────────────────────────────────
# Goal: primary research studies that examine ENVIRONMENTAL FACTOR EXPOSURES
# using EHR/administrative health data as the primary data source, plus
# pediatric vaccine/immunization as an exposure (Tier 5).
#
# All terms use [Title/Abstract] so the concepts must appear in the abstract,
# not just in references or acknowledgements.
#
# Excluded via publication type + title keywords:
#   Review, systematic review, meta-analysis, bibliometric, protocol, narrative review
#
_FILTERS = (
    'open access[filter] '
    'NOT Review[Publication Type] '
    'NOT Congress[Publication Type] '
    'NOT "Published Erratum"[Publication Type] '
    'NOT systematic[Title/Abstract] '
    'NOT meta-analysis[Title/Abstract] '
    'NOT bibliometric[Title/Abstract] '
    'NOT protocol[Title/Abstract] '
    'NOT narrative review[Title/Abstract] '
    # Pediatric scope guard: exclude clearly adult-only outcomes that leak in
    # because "child"/"children" appears in their abstracts incidentally.
    'NOT ("dementia"[Title/Abstract] OR "Alzheimer"[Title/Abstract] '
    'OR "midlife"[Title/Abstract] OR "mid-life"[Title/Abstract] '
    'OR "older adults"[Title/Abstract] OR "menopause"[Title/Abstract] '
    'OR "postmenopausal"[Title/Abstract] OR "late-onset"[Title/Abstract])'
)

# EHR synonyms — with [Title/Abstract] to require it in core text
_EHR_TA = (
    '"electronic health record"[Title/Abstract] OR '
    '"electronic medical record"[Title/Abstract] OR '
    '"EHR"[Title/Abstract] OR '
    '"EMR"[Title/Abstract] OR '
    '"claims data"[Title/Abstract] OR '
    '"administrative health data"[Title/Abstract] OR '
    '"health records"[Title/Abstract]'
)

# Pediatric / childhood population constraint — required in title/abstract.
# Every query below ANDs this in, so all retained hits are pediatric by construction.
_PED_TERMS = (
    '"pediatric"[Title/Abstract] OR "paediatric"[Title/Abstract] OR '
    '"child"[Title/Abstract] OR "children"[Title/Abstract] OR '
    '"childhood"[Title/Abstract] OR "infant"[Title/Abstract] OR '
    '"newborn"[Title/Abstract] OR "neonatal"[Title/Abstract] OR '
    '"adolescent"[Title/Abstract] OR "adolescence"[Title/Abstract] OR '
    '"youth"[Title/Abstract] OR "early life"[Title/Abstract] OR '
    '"pediatrics"[Title/Abstract]'
)
_PEDIATRIC_TA = f'({_PED_TERMS})'

SEARCH_QUERIES = [
    # ── Tier 1: explicit EWAS / exposome-wide in pediatric populations ─────
    f'("environment-wide association"[Title/Abstract] OR "exposome-wide association"[Title/Abstract]) {_PEDIATRIC_TA} {_FILTERS}',
    f'exposome[Title/Abstract] {_PEDIATRIC_TA} (association[Title/Abstract] OR risk[Title/Abstract]) {_FILTERS}',

    # ── Tier 2: environmental exposure + EHR + pediatric (all in abstract) ──
    f'({ _EHR_TA }) "environmental exposure"[Title/Abstract] {_PEDIATRIC_TA} (association[Title/Abstract] OR health[Title/Abstract]) {_FILTERS}',
    f'({ _EHR_TA }) "air pollution"[Title/Abstract] {_PEDIATRIC_TA} (cohort[Title/Abstract] OR association[Title/Abstract]) {_FILTERS}',
    f'({ _EHR_TA }) "chemical exposure"[Title/Abstract] {_PEDIATRIC_TA} health[Title/Abstract] {_FILTERS}',
    f'({ _EHR_TA }) "neighborhood environment"[Title/Abstract] {_PEDIATRIC_TA} {_FILTERS}',
    f'({ _EHR_TA }) "built environment"[Title/Abstract] {_PEDIATRIC_TA} health[Title/Abstract] {_FILTERS}',
    f'({ _EHR_TA }) "social determinants of health"[Title/Abstract] {_PEDIATRIC_TA} {_FILTERS}',
    f'({ _EHR_TA }) "prenatal"[Title/Abstract] ({_PED_TERMS}) exposure[Title/Abstract] {_FILTERS}',

    # ── Tier 3: geospatial / contextual exposure linkage to pediatric EHR ──
    f'({ _EHR_TA }) ("deprivation index"[Title/Abstract] OR "area deprivation"[Title/Abstract]) {_PEDIATRIC_TA} {_FILTERS}',
    f'({ _EHR_TA }) ("geospatial"[Title/Abstract] OR "geocod"[Title/Abstract]) {_PEDIATRIC_TA} exposure[Title/Abstract] {_FILTERS}',

    # ── Tier 4: birth-cohort / linked-data pediatric exposome (EHR term optional).
    # Pediatric exposome research is heavily birth-cohort based; the data source
    # is often linked administrative/registry data that does not say "EHR" in the
    # abstract, so requiring it (as Tier 2 does) misses most of the literature.
    f'(exposome[Title/Abstract] OR "environment-wide association"[Title/Abstract]) ({_PED_TERMS}) ("birth cohort"[Title/Abstract] OR cohort[Title/Abstract] OR "linked data"[Title/Abstract] OR "administrative data"[Title/Abstract]) {_FILTERS}',
    f'("prenatal exposure"[Title/Abstract] OR "prenatal exposome"[Title/Abstract] OR "early-life exposure"[Title/Abstract] OR "early life exposome"[Title/Abstract]) ({_PED_TERMS}) (association[Title/Abstract] OR risk[Title/Abstract]) {_FILTERS}',
    f'("air pollution"[Title/Abstract] OR "particulate matter"[Title/Abstract] OR PM2.5[Title/Abstract]) ({_PED_TERMS}) ("birth cohort"[Title/Abstract] OR cohort[Title/Abstract] OR "linked data"[Title/Abstract]) (asthma[Title/Abstract] OR respiratory[Title/Abstract] OR birth[Title/Abstract]) {_FILTERS}',
    f'("blood lead"[Title/Abstract] OR "chemical exposure"[Title/Abstract] OR "endocrine disruptor"[Title/Abstract]) ({_PED_TERMS}) (cohort[Title/Abstract] OR "linked data"[Title/Abstract] OR "administrative data"[Title/Abstract]) {_FILTERS}',
    f'({ _EHR_TA }) ("birth cohort"[Title/Abstract] OR "linked data"[Title/Abstract]) ({_PED_TERMS}) exposure[Title/Abstract] {_FILTERS}',

    # ── Tier 5: pediatric vaccine / immunization as the exposure ──
    # Vaccination is treated as an exposure (vaccine type / schedule / timing)
    # predicting a pediatric health/biomarker outcome (safety, febrile seizure,
    # asthma, BMI, infection, fever, neurodevelopment, autoimmune, SIDS).
    # No EHR term required: vaccine studies are often registry/claims/cohort
    # based and don't name 'EHR' in the abstract (same rationale as Tier 4).
    f'(vaccine[Title/Abstract] OR vaccination[Title/Abstract] OR immunization[Title/Abstract] OR immunisation[Title/Abstract]) ({_PED_TERMS}) ("adverse event"[Title/Abstract] OR "febrile seizure"[Title/Abstract] OR "vaccine safety"[Title/Abstract] OR "vaccine-associated"[Title/Abstract]) {_FILTERS}',
    f'(MMR[Title/Abstract] OR DTaP[Title/Abstract] OR "BCG vaccine"[Title/Abstract] OR "rotavirus vaccine"[Title/Abstract] OR "HPV vaccine"[Title/Abstract] OR "human papillomavirus vaccine"[Title/Abstract] OR "human papillomavirus vaccination"[Title/Abstract] OR "influenza vaccine"[Title/Abstract]) ({_PED_TERMS}) (febrile seizure[Title/Abstract] OR fever[Title/Abstract] OR asthma[Title/Abstract] OR infection[Title/Abstract] OR neurodevelopment[Title/Abstract] OR autoimmune[Title/Abstract]) {_FILTERS}',
    f'("vaccine safety"[Title/Abstract] OR "vaccine schedule"[Title/Abstract] OR "vaccination schedule"[Title/Abstract]) ({_PED_TERMS}) (cohort[Title/Abstract] OR "linked data"[Title/Abstract] OR "administrative data"[Title/Abstract] OR claims[Title/Abstract]) {_FILTERS}',
]

MAX_PER_QUERY = 200  # high enough to capture every hit from the broadest query


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
        count = int(data['count'])
        ids = data["idlist"]
        cap = f"  *** CAPPED: raise MAX_PER_QUERY ({count} > {MAX_PER_QUERY})" if count > len(ids) else ""
        print(f"  hits: {count:>5}  |  retrieved: {len(ids)}{cap}")
        return ids
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


def europepmc_pdf_url(pmcid: str) -> str | None:
    """
    Fallback PDF resolver via Europe PMC.
    Some open-access articles (e.g. late-deposited or project-design papers)
    expose no direct PDF link on the NCBI OA service; Europe PMC indexes the
    same PMC corpus and often carries a resolvable PDF URL. Returns None if no
    PDF-style full-text URL is found.
    """
    try:
        r = requests.get(EPMC_SEARCH, timeout=20, params={
            "query": f"PMCID:{pmcid}",
            "resultType": "core",
            "format": "json",
        })
        r.raise_for_status()
        results = r.json().get("resultList", {}).get("result", [])
        if not results:
            return None
        for entry in results[0].get("fullTextUrlList", {}).get("fullTextUrl", []):
            if entry.get("documentStyle", "").lower() == "pdf":
                return entry.get("url")
    except Exception as e:
        print(f"    Europe PMC lookup error: {e}")
    return None


EPMC_FULLTEXT_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"

def europepmc_fulltext_xml(pmcid: str) -> bytes | None:
    """
    Last-resort full-text retriever: download the JATS XML full text from
    Europe PMC. Used when no PDF is obtainable anywhere (NCBI OA tgz is a
    dead link AND Europe PMC has no resolvable PDF). Saves the complete
    article text + references, which is sufficient for review extraction.
    """
    try:
        r = requests.get(EPMC_FULLTEXT_XML.format(pmcid=pmcid), timeout=30,
                         headers={"User-Agent": "Mozilla/5.0 (academic research)"})
        if r.status_code == 200 and len(r.content) > 5_000:
            return r.content
        print(f"    Europe PMC XML: HTTP {r.status_code}, {len(r.content)} bytes")
    except Exception as e:
        print(f"    Europe PMC XML error: {e}")
    return None


def validate_fulltext(path: Path, is_xml: bool) -> bool:
    """
    Return True only if the saved file is genuine full text.

    - PDF: must start with the %PDF- magic and be non-trivially large.
    - XML: must be a JATS *article* (not an abstract-only record) with a real
      <body>. Conference-supplement records (e.g. Alzheimer's & Dementia,
      J Endocr Soc abstracts) come back as article-type="abstract" with no
      body — those are NOT full papers and are rejected here.
    """
    try:
        data = path.read_bytes()
    except OSError:
        return False
    if not is_xml:
        return data[:5] == b"%PDF-" and len(data) > 20_000
    if b'article-type="abstract"' in data[:3000]:
        return False
    try:
        root = ET.fromstring(data)
    except ET.ParseError:
        return False
    if root.attrib.get("article-type", "") == "abstract":
        return False
    body = root.find(".//body")
    if body is None:
        return False
    text = "".join(body.itertext()).strip()
    return len(text) > 2000


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

def _load_log(log_path: Path) -> dict:
    """Load the shared download log with an advisory file lock (fcntl) so parallel
    per-query jobs don't clobber each other. Returns the parsed dict (or {})."""
    import fcntl
    log_path.parent.mkdir(parents=True, exist_ok=True)
    # ponytail: fcntl.LOCK_SH on the log file serializes reads across per-query jobs.
    with open(log_path, "a+") as fh:
        fcntl.flock(fh, fcntl.LOCK_SH)
        fh.seek(0)
        try:
            return json.loads(fh.read() or "{}")
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def _save_log(log_path: Path, log: dict) -> None:
    """Write the shared log under an exclusive lock (merges per-query results)."""
    import fcntl
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a+") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        fh.seek(0); existing = fh.read().strip()
        base = json.loads(existing) if existing else {}
        # union each list/set field so concurrent per-query runs merge
        for k in ("downloaded","xml_only","excluded","failed","abstract_only","papers"):
            if k in log:
                merged = []
                seen = set()
                for item in list(base.get(k, [])) + list(log[k] if isinstance(log[k], list) else []):
                    key = json.dumps(item, sort_keys=True)
                    if key not in seen:
                        seen.add(key); merged.append(item)
                base[k] = merged
        base.update({k: v for k, v in log.items() if k not in ("downloaded","xml_only","excluded","failed","abstract_only","papers")})
        fh.seek(0); fh.truncate(); fh.write(json.dumps(base, indent=2))
        fcntl.flock(fh, fcntl.LOCK_UN)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--query-index", type=int, default=None,
                    help="Run ONLY SEARCH_QUERIES[N] (0-based). Lets each query run as its "
                         "own isolated job; the shared download log converges via file locks.")
    args = ap.parse_args(argv)

    OUTPUT_DIR.mkdir(exist_ok=True)
    log_path = OUTPUT_DIR / "download_log.json"
    log = _load_log(log_path)
    already_done = set(log.get("downloaded", []))

    # ── 1. Search ──────────────────────────────────────────────────────────
    print("=" * 65)
    print("  SEARCHING PubMed Central (open-access filter)")
    print("=" * 65)
    all_ids: set = set()
    queries = SEARCH_QUERIES if args.query_index is None else [SEARCH_QUERIES[args.query_index]]
    for q in queries:
        label = q.split(" open access[filter]")[0].strip()
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

    # ── 3. Filter out reviews / conference abstracts / off-topic articles ────
    REVIEW_TITLE_KW = re.compile(
        r'\b(review|systematic review|meta.analysis|bibliometric|'
        r'narrative review|scoping review|overview|commentary|'
        r'perspective|editorial|protocol|roadmap|journal club|'
        r'educational|tutorial|primer|poster)\b',
        re.IGNORECASE
    )
    EPIGENOME_KW = re.compile(
        r'\b(epigenome.wide|methylation|DNA methyl|CpG|histone)\b',
        re.IGNORECASE
    )
    # Conference abstract title patterns: "P-1234.", "SAT-123", "123 Title", poster codes
    CONF_ABSTRACT_KW = re.compile(
        r'^(P-\d+\.|SAT-\d+|SUN-\d+|MON-\d+|TUE-\d+|THU-\d+|\d{4,5}\s)',
        re.IGNORECASE
    )
    # Journals that are exclusively conference supplement/abstract publications
    CONF_JOURNALS = {
        'innov aging',        # Gerontological Society abstract supplement
        'j clin transl sci',  # ACTS conference abstracts
        'open forum infect dis',  # IDSA conference supplement abstracts
        'j endocr soc',       # ENDO abstract supplement
    }

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

        is_review    = bool(REVIEW_TITLE_KW.search(title)) or "Review" in pub_type
        is_epigenome = bool(EPIGENOME_KW.search(title)) and \
                       "environment-wide" not in title.lower() and \
                       "exposome" not in title.lower()
        is_conf_abs  = bool(CONF_ABSTRACT_KW.search(title)) or \
                       journal.lower() in CONF_JOURNALS

        if is_review:
            flag = "[REVIEW]"
            skipped_reviews.append({"pmcid": f"PMC{uid}", "title": title, "reason": "review"})
        elif is_epigenome:
            flag = "[EPIOME]"
            skipped_reviews.append({"pmcid": f"PMC{uid}", "title": title, "reason": "epigenome-wide (not environment-wide)"})
        elif is_conf_abs:
            flag = "[CONF]"
            skipped_reviews.append({"pmcid": f"PMC{uid}", "title": title, "reason": "conference abstract"})
        else:
            flag = "[OK]"
            papers.append((uid, title, journal, year, authors))

        print(f"{i:<4} PMC{uid:<11} {year:<5} {flag:<9} {title[:51]:<52} {journal[:25]}")

    print(f"\n  Primary research papers : {len(papers)}")
    print(f"  Excluded (reviews)      : {sum(1 for s in skipped_reviews if s['reason']=='review')}")
    print(f"  Excluded (epigenome)    : {sum(1 for s in skipped_reviews if 'epigenome' in s['reason'])}")
    print(f"  Excluded (conf abstract): {sum(1 for s in skipped_reviews if s['reason']=='conference abstract')}")
    log["excluded"] = skipped_reviews

    # ── 5. Download ────────────────────────────────────────────────────────
    print(f"\n{'=' * 65}")
    print(f"  DOWNLOADING PDFs  →  ./{OUTPUT_DIR}/")
    print("=" * 65 + "\n")

    downloaded = skipped = failed = abstract_only = 0

    def record_success(uid_, pmcid_, out_stem_, data_: bytes, fmt: str, title_: str) -> bool:
        """Write + validate a downloaded file. Returns True if kept as full text."""
        nonlocal downloaded, abstract_only
        is_xml = fmt != "pdf"
        ext = "xml" if is_xml else "pdf"
        path = OUTPUT_DIR / f"{out_stem_}.{ext}"
        path.write_bytes(data_)
        size_kb = path.stat().st_size // 1024
        if not validate_fulltext(path, is_xml=is_xml):
            path.unlink(missing_ok=True)
            print(f"           ✗ {ext.upper()} not full text (abstract-only record) — discarded")
            log.setdefault("abstract_only", []).append({"pmcid": pmcid_, "title": title_})
            abstract_only += 1
            return False
        tag = "" if not is_xml else " (xml)"
        print(f"           ✓{tag} {path.name}  ({size_kb} KB)")
        log.setdefault("downloaded", []).append(uid_)
        if is_xml:
            log.setdefault("xml_only", []).append(pmcid_)
        already_done.add(uid_)
        downloaded += 1
        return True

    for idx, (uid, title, journal, year, authors) in enumerate(papers, 1):
        pmcid    = f"PMC{uid}"
        out_stem = sanitize(f"{year}_{pmcid}_{title}")
        pdf_path = OUTPUT_DIR / f"{out_stem}.pdf"
        xml_path = OUTPUT_DIR / f"{out_stem}.xml"
        pfx      = f"  [{idx:>2}/{len(papers)}] {pmcid} ({year})"

        if uid in already_done or pdf_path.exists() or xml_path.exists():
            print(f"{pfx}  [SKIP]")
            skipped += 1
            continue

        print(f"{pfx}  {title[:55]}")

        pdf_url, tgz_url = get_oa_links(pmcid)
        time.sleep(DELAY)

        data: bytes | None = None
        fmt = "pdf"

        # 1. NCBI OA direct PDF.
        if pdf_url:
            print(f"           → PDF: {pdf_url[:78]}")
            data = download_bytes(pdf_url)
        # 2. NCBI OA tar.gz → extract main article PDF.
        if data is None and tgz_url:
            print(f"           → TAR: {tgz_url[:78]}")
            raw = download_bytes(tgz_url)
            if raw:
                data = pdf_from_tgz(raw)
        # 3. Europe PMC PDF (often resolvable where NCBI OA is a dead link).
        if data is None:
            ep_url = europepmc_pdf_url(pmcid)
            if ep_url:
                print(f"           → EuropePMC PDF: {ep_url[:78]}")
                data = download_bytes(ep_url)
                time.sleep(DELAY)
        # 4. Europe PMC JATS XML full text — last resort to capture full paper.
        if data is None:
            xml_data = europepmc_fulltext_xml(pmcid)
            if xml_data:
                print(f"           → EuropePMC XML (full text)")
                data = xml_data
                fmt = "xml"

        if data is not None:
            record_success(uid, pmcid, out_stem, data, fmt, title)
        else:
            print(f"           ✗ Could not retrieve full text")
            log.setdefault("failed", []).append({"pmcid": pmcid, "title": title})
            failed += 1

        time.sleep(DELAY)

    # ── 6. Save log ────────────────────────────────────────────────────────
    log["papers"] = [
        {"pmcid": f"PMC{uid}", "title": t, "journal": j, "year": y, "authors": a}
        for uid, t, j, y, a in papers
    ]
    _save_log(log_path, log)

    print(f"\n{'=' * 50}")
    print(f"  Downloaded    : {downloaded}")
    print(f"  Skipped       : {skipped}  (already on disk)")
    print(f"  Abstract-only : {abstract_only}  (not full text — discarded)")
    print(f"  Failed        : {failed}")
    print(f"  Log           : {log_path}")
    print(f"  Output        : {OUTPUT_DIR.resolve()}")
    print("=" * 50)


if __name__ == "__main__":
    main()
