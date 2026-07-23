# Manuscript summarization & data-availability

How each downloaded paper is turned into a structured, validated checklist and
how its data-availability is captured.

← back to the [README](../README.md)

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

### Checklist schema (`summarizer/schema.py`, `ManuscriptChecklist`)

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

[`scan_data_availability.py`](../scan_data_availability.py) asks Gemma 4 12B
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

Latest scan (166 papers): 127 not-stated, 16 supplementary-only,
11 available-upon-request, **9 public-repository** (with accession/links), 3 in-house.
