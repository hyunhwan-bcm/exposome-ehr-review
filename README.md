# Pediatric Exposome / EWAS Literature Collection

A reproducible pipeline that searches PubMed Central for **pediatric / childhood
environmental-exposure (exposome / EWAS) studies** using EHR, administrative /
claims, or linked cohort data, then downloads the open-access full text and
produces a structured inventory for systematic-review work.

> **Current collection: 116 full-text papers** (99 PDFs + 17 JATS-XML),
> 1990–2026. See [`paper_summary.md`](./paper_summary.md) for the inventory.

## Quick start

```bash
make setup      # create .venv and install requests
make download   # run the pediatric PMC fetcher (incremental — skips what's on disk)
make summary    # regenerate paper_summary.md + print a compact inventory
make fresh      # clean + re-download from scratch
make clean      # remove downloaded PDFs/XML + download log
```

`papers/` (PDFs, XML, `download_log.json`) is gitignored — only the scripts,
Makefile, and `paper_summary.md` are tracked.

## Data-collection process

The pipeline is implemented in [`fetch_pmc_papers.py`](./fetch_pmc_papers.py)
and runs in seven stages: a pediatric-constrained PMC search → metadata fetch
→ candidate filtering → a cascading full-text resolver → full-text
validation → on-disk output + audit log → summary regeneration.

```mermaid
flowchart TD
    Start([Pediatric exposome / EWAS literature collection]) --> S1

    subgraph S1["1. Search PubMed Central — NCBI ESearch"]
        direction TB
        T1["Tier 1 — EWAS / exposome-wide × pediatric"]
        T2["Tier 2 — environmental exposure × EHR × pediatric"]
        T3["Tier 3 — geospatial / area deprivation × pediatric EHR"]
        T4["Tier 4 — birth cohort / linked data × pediatric"]
        T1 & T2 & T3 & T4 --> ES["esearch.fcgi<br/>retmax = 200"]
        ES --> Cap{"hits &gt; retmax?"}
        Cap -- yes --> Warn["⚠ CAPPED warning in log"]
        Cap -- no --> IDs
        Warn --> IDs["unique PMC IDs<br/>(union across tiers)"]
    end

    IDs --> S2["2. Fetch metadata — NCBI ESummary"]

    S2 --> S3
    subgraph S3["3. Filter candidates (title / pubtype)"]
        direction TB
        F1{"Review / journal club /<br/>poster / educational?"}
        F2{"Epigenome-only<br/>(not environment)?"}
        F3{"Conference abstract?<br/>P-1234 / SAT-xxx / conf journal"}
        F1 -- yes --> XR[Excluded — review / poster]
        F1 -- no --> F2
        F2 -- yes --> XE[Excluded — epigenome]
        F2 -- no --> F3
        F3 -- yes --> XC[Excluded — conference]
        F3 -- no --> Cand[Candidate paper]
    end

    NoteA["Adult-outcome negation<br/>(dementia, Alzheimer, midlife,<br/>menopause, older adults)<br/>applied at query level"]
    NoteA -.-> S3

    Cand --> S4
    subgraph S4["4. Resolve full text — cascading fallbacks"]
        direction TB
        R1["a. NCBI OA direct PDF"]
        R2["b. NCBI OA tar.gz → extract main PDF<br/>(supplementary / appendix filtered out)"]
        R3["c. Europe PMC PDF"]
        R4["d. Europe PMC JATS XML full text"]
        R1 -- no PDF --> R2
        R2 -- no PDF --> R3
        R3 -- no PDF --> R4
    end

    S4 --> Got{Got bytes?}
    Got -- no --> Fail[Logged as failed]
    Got -- yes --> S5

    subgraph S5["5. Validate full text"]
        direction TB
        VF{"Format?"}
        VP{"PDF: %PDF- magic<br/>and size &gt; 20 KB?"}
        VX{"XML: not article-type=abstract,<br/>has &lt;body&gt; with &gt; 2 000 chars?"}
        VF -- PDF --> VP
        VF -- XML --> VX
        VP -- no --> AO[Discard — abstract-only]
        VX -- no --> AO
        VP -- yes --> Keep
        VX -- yes --> Keep
    end

    Keep --> S6["6. Write papers/*.pdf | *.xml<br/>+ append to download_log.json<br/>(downloaded / xml_only / failed /<br/>excluded / abstract_only)"]
    S6 --> S7["7. build_summary.py<br/>→ paper_summary.md<br/>grouped by exposure domain &<br/>health outcome"]
    S7 --> Done([Done])
```

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

## Full-text resolution & validation

Because many open-access papers have dead NCBI OA links, retrieval uses a
cascade of fallbacks: NCBI direct PDF → NCBI tar.gz extraction → Europe PMC
PDF → Europe PMC JATS XML. Every downloaded file is then **validated** as
genuine full text — PDFs by the `%PDF-` magic + minimum size, XML by checking
it is not an `article-type="abstract"` record and that it has a real `<body>`
— so abstract-only conference / supplement records that slipped past the title
filter are discarded rather than counted as papers.

## Output

| Path | Contents |
|------|----------|
| `papers/*.pdf` | Downloaded full-text PDFs |
| `papers/*.xml` | JATS-XML full text (where no PDF was resolvable) |
| `papers/download_log.json` | Audit log: `downloaded`, `excluded`, `failed`, `abstract_only`, `xml_only`, `papers` |
| `paper_summary.md` | Human-readable inventory grouped by exposure domain and health outcome |

## Files

| File | Purpose |
|------|---------|
| `fetch_pmc_papers.py` | Search + filter + download pipeline |
| `build_summary.py` | Regenerates `paper_summary.md` from the download log |
| `Makefile` | `setup` / `download` / `clean` / `fresh` / `summary` targets |
| `paper_summary.md` | Generated inventory (do not hand-edit) |
