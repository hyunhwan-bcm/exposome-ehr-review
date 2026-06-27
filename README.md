# Pediatric Exposome / EWAS Literature Collection

A reproducible pipeline that searches PubMed Central for **pediatric / childhood
environmental-exposure (exposome / EWAS) studies** and **pediatric
vaccine/immunization-exposure studies**, using EHR, administrative / claims, or
linked cohort data, downloads the open-access full text, summarizes each
manuscript with **Gemma 4 12B**, and captures **data-availability** (accession
numbers / repository links) for systematic-review work.

> **Current collection: 122 full-text papers**, 1990–2026, **22 EHR-based**.
> See [`paper_summary.md`](./paper_summary.md) for the inventory.

## Quick start

```bash
make setup      # create .venv + install deps (requests, openai, pypdf, pydantic, tinydb, dagster)
make download   # pediatric PMC fetcher (incremental — skips what's on disk)
make summarize  # Gemma 4 12B -> per-paper + combined JSON (concurrent; --workers via SUMMARIZE_ARGS)
make results    # export readable results/ (SUMMARY.md, checklist.md, combined JSON)
make test       # unit tests (no live API calls)
```

`papers/` (PDFs, XML, `download_log.json`, TinyDB `db.json`) is gitignored —
only scripts, the Makefile, `results/`, and `paper_summary.md` are tracked.

## Pipeline overview

```mermaid
flowchart LR
    D[make download<br/>fetch_pmc_papers.py] --> S[make summarize<br/>summarizer/run.py]
    S -->|per-paper JSON| DB[(TinyDB store<br/>papers/db.json)]
    S --> A[scan_data_availability.py<br/>Gemma 4 12B + regex]
    A -->|data_availability, accession links| DB
    DB --> C[make db-export<br/>combined JSON]
    C --> R[make results<br/>build_results.py]
    R --> OUT[results/SUMMARY.md<br/>results/checklist.md<br/>results/manuscript_summaries.json]
```

The same graph is orchestrated as **Dagster assets** (`pipeline.py`) — run the
UI with `make dagster` or materialize headlessly with `make materialize`:

```mermaid
flowchart LR
    download_log --> per_paper_summaries --> data_availability_scan --> manuscript_summaries --> results
```

## Manuscript summarization (Gemma 4 12B → Pydantic JSON)

Each manuscript is summarized into a structured **checklist** by **Gemma 4 12B**
via an external OpenAI-compatible endpoint, validated with a Pydantic schema.
`--workers N` runs papers concurrently (the OpenAI client is thread-safe;
bottleneck is network-bound LLM calls).

```bash
cp .env.example .env        # fill in GEMMA_API_KEY (never committed)
make summarize              # all manuscripts (resume-only, chunked recovery)
make summarize-paper PMC=PMC7145790   # single paper
python -m summarizer.run --recover --workers 4   # 4 concurrent
```

**Output:** `papers/summaries/<pmcid>.json` (one per paper) +
`papers/manuscript_summaries.json` (combined).

**Checklist schema** (`summarizer/schema.py`, `ManuscriptChecklist`):

| Field | Type | Description |
|-------|------|-------------|
| `pmcid` / `title` / `year` | str | identity (from download log) |
| `ehr_used` | bool | does the study use EHR/EMR/claims/admin data? |
| `ehr_evidence` | str | sentence(s) justifying `ehr_used` |
| `summary` | str | 2-4 sentence summary |
| `key_findings` | list[str] | main results |
| `captured_features` | list[str] | EHR features/variables captured |
| `pathologies_diseases` | list[str] | disease(s)/outcome(s) |
| `study_design` / `data_source_type` / `population` / `exposure_domain` | str | review fields |
| `limitations` | list[str] | stated limitations |
| `confidence` | high\|medium\|low\|unclear | fit for a pediatric EHR/exposome review |
| **`data_availability`** | enum | public-repository\|available-upon-request\|in-house\|supplementary-only\|not-stated |
| **`data_accession_links`** | list[str] | accession IDs / repository URLs / names |
| **`data_availability_statement`** | str | verbatim sentence(s) justifying the call |
| `source_format` / `model` | str | provenance |

The model wraps JSON in chain-of-thought, so the client uses a strict
fixed-key prompt, robust fenced/balanced-JSON extraction with light repair for
truncated responses, lenient field validators, and a retry-with-nudge loop.

## Data-availability scan

[`scan_data_availability.py`](./scan_data_availability.py) asks Gemma 4 12B
(32k-token budget, pydantic-validated) **only** about how a study's data can be
obtained, in a focused window around the "Data availability" section. A
deterministic regex safety-net supplements any accession/URL (dbGaP / GSE /
PRJEB / Zenodo / GitHub / figshare / Dryad) the model drops, and repairs URLs
broken across PDF line wraps.

```bash
python scan_data_availability.py                 # all papers
python scan_data_availability.py --limit 5       # pilot
python scan_data_availability.py --workers 4     # concurrent
```

Latest scan (122 papers): 96 not-stated, 10 supplementary-only,
7 available-upon-request, **7 public-repository** (with accession/links), 2 in-house.

## TinyDB store

[`database.py`](./database.py) + [`db.py`](./db.py): a TinyDB-backed document
store that validates every record through `ManuscriptChecklist` (the catalog
can never drift into a state `build_results.py` can't read). The DB is the
source of truth for the combined file.

```bash
make db-import    # papers/summaries/*.json -> papers/db.json
make db-export    # store -> combined + results/manuscript_summaries.json
make db-stats     # record counts
python db.py find --ehr-used --disease asthma
python db.py update PMC1234567 --no-ehr-used --ehr-evidence "n/a"
```

## Dagster orchestration

[`pipeline.py`](./pipeline.py) wraps the existing scripts as assets with
lineage: `download_log` → `per_paper_summaries` → `data_availability_scan` →
`manuscript_summaries` → `results`. The summarizer runs in `--recover` +
`--workers 4` mode, and the focused data-availability scan is a first-class
asset so `make materialize` reproduces the same enriched outputs as the manual
workflow.

```bash
make dagster       # UI + lineage browser (dagster dev -m pipeline)
make materialize   # materialize the whole graph headlessly
```

## Data-collection process

Implemented in [`fetch_pmc_papers.py`](./fetch_pmc_papers.py) in seven stages:
pediatric-constrained PMC search → metadata fetch → candidate filtering →
cascading full-text resolver → full-text validation → on-disk output + audit log
→ summary regeneration.

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

## Output

| Path | Contents |
|------|----------|
| `papers/*.pdf` | Downloaded full-text PDFs (gitignored) |
| `papers/*.xml` | JATS-XML full text where no PDF was resolvable (gitignored) |
| `papers/download_log.json` | Audit log: `downloaded`, `excluded`, `failed`, `abstract_only`, `xml_only`, `papers` (gitignored) |
| `papers/summaries/<pmcid>.json` | One checklist per paper (gitignored) |
| `papers/db.json` | TinyDB document store (gitignored, source of truth) |
| `papers/manuscript_summaries.json` | Combined `SummaryBatch` JSON (gitignored, regenerated by the store) |
| `results/manuscript_summaries.json` | **Tracked** combined copy |
| `results/SUMMARY.md` / `results/checklist.md` | **Tracked** readable inventory (word-boundary cell clipping) |
| `paper_summary.md` | Human-readable inventory grouped by exposure domain and health outcome |

## Files

| File | Purpose |
|------|---------|
| `fetch_pmc_papers.py` | Search + filter + download pipeline |
| `build_summary.py` | Regenerates `paper_summary.md` from the download log |
| `build_results.py` | Exports readable `results/` from the combined JSON |
| `summarizer/` | Manuscript summarization (Pydantic schema + Gemma client + extractor + runner) |
| `scan_data_availability.py` | Focused LLM + regex scan for data-availability / accession links |
| `database.py` / `db.py` | TinyDB store + CRUD CLI |
| `pipeline.py` / `pipeline_ops.py` | Dagster asset orchestration (fetch → summarize → combined → results) |
| `Makefile` | `setup` / `download` / `summarize` / `results` / `db-*` / `dagster` / `materialize` / `test` targets |
| `paper_summary.md` | Generated inventory (do not hand-edit) |