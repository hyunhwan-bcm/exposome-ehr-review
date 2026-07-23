# Data collection

How the corpus is searched, filtered, downloaded, and validated. Implemented in
[`fetch_pmc_papers.py`](../fetch_pmc_papers.py).

← back to the [README](../README.md)

## Process (seven stages)

1. **Search PubMed Central** (NCBI ESearch) — Tier 1–5 queries, `retmax = 200`,
   union of unique PMC IDs across tiers (a CAPPED warning is logged if any
   query's hit count exceeds `retmax`).
2. **Fetch metadata** (NCBI ESummary).
3. **Filter candidates** by title / pubtype — drop reviews / posters /
   educational, epigenome-only (non-environment), and conference abstracts
   (`P-1234` / `SAT-xxx` / conference-only journals). Adult-outcome negation
   (dementia, Alzheimer, midlife, menopause, older adults) is applied at the
   query level.
4. **Resolve full text** via cascading fallbacks — NCBI OA direct PDF → NCBI OA
   tar.gz (main PDF extracted, supplementary/appendix filtered out) → Europe PMC
   PDF → Europe PMC JATS XML.
5. **Validate full text** — PDFs by `%PDF-` magic + size > 20 KB; XML by not
   being `article-type="abstract"` and having a `<body>` with > 2 000 chars.
   Abstract-only records are discarded.
6. **Write** `papers/*.pdf | *.xml` and append to `download_log.json`
   (`downloaded` / `xml_only` / `failed` / `excluded` / `abstract_only`).
7. **Regenerate the inventory** — `build_summary.py` → `paper_summary.md`,
   grouped by exposure domain & health outcome.

## Search strategy

Every query ANDs in a pediatric population constraint
(`pediatric` / `paediatric` / `child` / `children` / `childhood` / `infant` /
`newborn` / `neonatal` / `adolescent` / `youth` / `early life` / `pediatrics`,
all `[Title/Abstract]`) so every retained hit is pediatric by construction.
Adult-only outcomes are negated at the query level to stop adult studies
leaking in via incidental "child" mentions.

| Tier | Rationale |
|------|-----------|
| 1 | Explicit `environment-wide` / `exposome-wide` association in pediatric populations |
| 2 | Environmental exposure × EHR/claims/admin data × pediatric (all in abstract) |
| 3 | Geospatial / area-deprivation exposure linked to pediatric EHR |
| 4 | Birth-cohort / linked-data pediatric exposome — broadened because most pediatric exposome research does not name "EHR" in the abstract |
| 5 | Pediatric **vaccine / immunization as the exposure** → health outcome (safety, febrile seizure, fever, asthma, infection, neurodevelopment, autoimmune); named vaccines MMR/DTaP/BCG/rotavirus/HPV/flu; exposure_domain tag `vaccine / immunization`. No EHR term required — vaccine studies are often registry/claims/cohort based |

## Full-text resolution & validation

Because many open-access papers have dead NCBI OA links, retrieval uses a
cascade of fallbacks: NCBI direct PDF → NCBI tar.gz extraction → Europe PMC
PDF → Europe PMC JATS XML. Every downloaded file is then **validated** as
genuine full text — PDFs by the `%PDF-` magic + minimum size, XML by checking
it is not an `article-type="abstract"` record and that it has a real `<body>`
(plus `<data-availability>` / `sec-type="data-availability"` from `<back>`) —
so abstract-only conference / supplement records that slipped past the title
filter are discarded rather than counted as papers.
