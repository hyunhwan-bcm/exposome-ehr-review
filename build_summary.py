#!/usr/bin/env python3
"""
Regenerate paper_summary.md as a pediatric-focused inventory from download_log.json.

Only papers that were actually retrieved as full text (PDF or JATS XML) are
listed; abstract-only / failed records are excluded. Papers are bucketed by
exposure domain and by health outcome using title-keyword heuristics so the
collection can be scanned quickly for a systematic-review worktable.
"""
import json
import re
from pathlib import Path

LOG = Path("papers/download_log.json")
OUT = Path("paper_summary.md")

# ── keyword → bucket ──────────────────────────────────────────────────────────
EXPOSURE_DOMAINS = {
    "Air pollution":   [r"air pollution", r"particulate", r"PM2\.5", r"PM10",
                        r"black carbon", r"traffic", r"nitrogen dioxide", r"NO2",
                        r"wildfire", r"ultra low emission", r"ULEZ", r"ozone",
                        r"endotoxin", r"smok(?!(e|ing)\s+(cessation|quit))"],
    "Metals / lead":   [r"\blead\b", r"blood lead", r"mercury", r"arsenic",
                        r"cadmium", r"manganese", r"metals?"],
    "Chemicals / EDC": [r"perfluoro", r"PFAS", r"PFOA", r"PFOS", r"alkylphenol",
                        r"bisphenol", r"phthalate", r"persistent organic pollut",
                        r"\bPOP[s]?\b", r"pesticide", r"insecticide",
                        r"organochlor", r"DDE", r"DDT", r"endocrine disrupt",
                        r"pyrethroid", r"paraben", r"triclosan"],
    "Prenatal / perinatal": [r"prenatal", r"gestational", r"maternal",
                        r"in utero", r"fetal", r"perinatal", r"pregnancy"],
    "SDOH / neighborhood": [r"social determinants", r"deprivation", r"neighborhood",
                        r"built environment", r"socioeconomic", r"foster care",
                        r"adverse childhood experience", r"\bACE[s]?\b",
                        r"inequality", r"disparit"],
    "Tobacco / e-cig": [r"electronic cigarette", r"e-cig", r"tobacco",
                        r"secondhand smoke", r"nicotine", r"maternal smoking"],
    "Other / mixed":   [],
}

OUTCOMES = {
    "Respiratory":     [r"asthma", r"respiratory", r"lung function", r"wheez",
                        r"FEV", r"airway", r"pulmonary"],
    "Neurodevelopment": [r"neurodevelop", r"cognit", r"intelligence", r"IQ",
                        r"intellectual", r"autism", r"behavior", r"behaviour",
                        r"neurobehavior", r"mental", r"psychomotor"],
    "Growth / birth":   [r"birth weight", r"fetal growth", r"preterm", r"gestational age",
                        r"small for gestational", r"low birth", r"BMI", r"obesity",
                        r"adipos", r"growth", r"birth outcome"],
    "Metabolic / endocrine": [r"diabet", r"glucose", r"insulin", r"cardiometabolic",
                        r"metabolic", r"thyroid", r"puberty", r"lipid"],
    "Renal / other":    [r"renal", r"glomerular", r"kidney"],
    "Behavioral / social": [r"suicid", r"adverse childhood", r"foster care",
                        r"wear-and-tear", r"stress"],
    "Methods / characterization": [r"metabolome", r"exposome", r"reference",
                        r"characterization", r"quality control", r"biomarker"],
    "Other / mixed":    [],
}


def first_bucket(title: str, table: dict) -> str:
    t = title.lower()
    for label, patterns in table.items():
        if any(re.search(p, t) for p in patterns):
            return label
    return "Other / mixed"


def main():
    log = json.loads(LOG.read_text())
    downloaded_uids = {uid.replace("PMC", "") for uid in log.get("downloaded", [])}
    # downloaded[] stores bare PMC numeric ids; papers[] uses "PMC<n>"
    downloaded_pmcids = {f"PMC{u}" for u in downloaded_uids}

    papers = [p for p in log.get("papers", [])
              if p["pmcid"] in downloaded_pmcids]
    xml_only = set(log.get("xml_only", []))
    papers.sort(key=lambda p: (p.get("year", ""), p["pmcid"]))

    # ── bucketing ─────────────────────────────────────────────────────────
    by_exposure = {}
    by_outcome = {}
    for p in papers:
        ex = first_bucket(p["title"], EXPOSURE_DOMAINS)
        oc = first_bucket(p["title"], OUTCOMES)
        p["exposure"] = ex
        p["outcome"] = oc
        by_exposure.setdefault(ex, []).append(p)
        by_outcome.setdefault(oc, []).append(p)

    n_pdf = sum(1 for p in papers if p["pmcid"] not in xml_only)
    n_xml = sum(1 for p in papers if p["pmcid"] in xml_only)
    years = sorted({p.get("year", "?") for p in papers})

    # ── markdown ───────────────────────────────────────────────────────────
    L = []
    L.append("# Pediatric Exposome / Environmental-Exposure Studies with EHR & Linked Health Data")
    L.append("")
    L.append(f"> **{len(papers)} open-access full-text papers** retrieved from PubMed Central.")
    L.append("> Scope: **childhood / pediatric** environmental-exposure (exposome / EWAS) studies")
    L.append("> using EHR, administrative/claims, or linked cohort data. Pivot from the prior")
    L.append("> adult-focused set; broadened to birth-cohort / linked-data literature where the")
    L.append("> abstract does not name \"EHR\" explicitly.")
    L.append("")
    L.append(f"> {n_pdf} PDFs · {n_xml} JATS-XML full texts · span {years[0]}–{years[-1]}.  ")
    L.append(f"> {len(log.get('abstract_only', []))} abstract-only conference records discarded; "
             f"{len(log.get('excluded', []))} reviews/conference abstracts excluded.")
    L.append("")
    L.append("---")
    L.append("")

    # main table
    L.append("## Inventory")
    L.append("")
    L.append("| # | PMCID | Year | Journal | Title | Exposure domain | Outcome | Format |")
    L.append("|---|-------|------|---------|-------|------------------|---------|--------|")
    for i, p in enumerate(papers, 1):
        fmt = "XML" if p["pmcid"] in xml_only else "PDF"
        title = p["title"].replace("|", "\\|")
        journal = (p.get("journal") or "").replace("|", "\\|")
        L.append(f"| {i} | {p['pmcid']} | {p.get('year','')} | {journal} | "
                 f"{title} | {p['exposure']} | {p['outcome']} | {fmt} |")
    L.append("")
    L.append("---")
    L.append("")

    # by exposure
    L.append("## Quick Reference by Exposure Domain")
    L.append("")
    L.append("| Domain | Count | PMCIDs |")
    L.append("|--------|-------|--------|")
    for label in EXPOSURE_DOMAINS:
        grp = by_exposure.get(label, [])
        if not grp:
            continue
        ids = ", ".join(p["pmcid"] for p in grp)
        L.append(f"| {label} | {len(grp)} | {ids} |")
    L.append("")
    L.append("---")
    L.append("")

    # by outcome
    L.append("## Quick Reference by Health Outcome")
    L.append("")
    L.append("| Outcome | Count | PMCIDs |")
    L.append("|---------|-------|--------|")
    for label in OUTCOMES:
        grp = by_outcome.get(label, [])
        if not grp:
            continue
        ids = ", ".join(p["pmcid"] for p in grp)
        L.append(f"| {label} | {len(grp)} | {ids} |")
    L.append("")
    L.append("---")
    L.append("")
    L.append("*Regenerated by `build_summary.py` from `papers/download_log.json`. "
    "Run `make summary` to reprint, `make fresh` to re-collect.*")

    OUT.write_text("\n".join(L) + "\n")
    print(f"✓ wrote {OUT} — {len(papers)} papers, "
          f"{len(by_exposure)} exposure domains, {len(by_outcome)} outcome buckets")


if __name__ == "__main__":
    main()
