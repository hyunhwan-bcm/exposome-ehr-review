"""Full-text extraction from downloaded PDFs and JATS-XML records."""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path


def extract_pdf(path: Path, max_pages: int = 60) -> str:
    """Extract text from a PDF via pypdf, capped to the first ``max_pages``."""
    from pypdf import PdfReader

    reader = PdfReader(str(path))
    pages = []
    for i, page in enumerate(reader.pages[:max_pages]):
        try:
            pages.append(page.extract_text() or "")
        except Exception:
            continue
    return "\n".join(pages).strip()


def extract_jats_xml(path: Path) -> str:
    """Extract readable text from a JATS-XML full-text record.

    Pulls <article-title>, <abstract>, the <body> paragraphs, and any
    data-availability section (which in JATS lives in <back>, not <body> —
    without this the data-availability scan would never see it).
    Conference / abstract-only records (article-type="abstract") have no
    <body> and return just title + abstract.
    """
    try:
        root = ET.parse(str(path)).getroot()
    except ET.ParseError:
        return Path(path).stem

    def text_of(tag: str) -> str:
        chunks = []
        for el in root.iter(tag):
            chunks.append("".join(el.itertext()))
        return " ".join(chunks).strip()

    title = text_of("article-title")
    abstract = text_of("abstract")
    body = text_of("body")

    # data-availability: <data-availability> element or a <sec> titled for it.
    da_chunks = []
    for el in root.iter("data-availability"):
        da_chunks.append("".join(el.itertext()))
    for sec in root.iter("sec"):
        sec_type = sec.attrib.get("sec-type", "")
        title_el = sec.find("title")
        heading = "".join(title_el.itertext()) if title_el is not None else ""
        if "data-availability" in sec_type.lower() or "data availability" in heading.lower():
            da_chunks.append("".join(sec.itertext()))
    data_availability = " ".join(da_chunks).strip()

    parts = [p for p in (title, abstract, body, data_availability) if p]
    return "\n\n".join(parts).strip()


def extract(path: Path) -> tuple[str, str]:
    """Return (text, source_format) for a downloaded paper file."""
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        return extract_pdf(path), "pdf"
    if suffix == ".xml":
        return extract_jats_xml(path), "xml"
    return path.read_text(errors="ignore"), "text"


# ── filename → pmcid ──────────────────────────────────────────────────────────
_PMCID_RE = re.compile(r"(PMC\d+)")


def pmcid_from_filename(path: Path) -> str:
    m = _PMCID_RE.search(path.stem)
    return m.group(1) if m else path.stem
