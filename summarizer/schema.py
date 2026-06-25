"""Pydantic checklist schema for manuscript summarization.

Each collected paper is summarized by Gemma 4 12B into a structured record
conforming to :class:`ManuscriptChecklist`. The schema is the single source of
truth for the JSON written under ``papers/summaries/`` and the combined
``papers/manuscript_summaries.json``.
"""
from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, field_validator


class ManuscriptChecklist(BaseModel):
    """Structured checklist extracted from a single manuscript."""

    # ── identity (filled from the download log, not the LLM) ──────────────
    pmcid: str = Field(..., description="PubMed Central ID, e.g. PMC7145790")
    title: str = Field(..., description="Article title")
    year: str = Field(..., description="Publication year")

    # ── existence-of-EHR checklist (the core question) ────────────────────
    ehr_used: bool = Field(
        ...,
        description="True if the study actually uses EHR / EMR / claims / "
                    "administrative health records as a data source.",
    )
    ehr_evidence: str = Field(
        ...,
        description="The sentence(s) from the manuscript that justify the "
                    "ehr_used decision; 'n/a' if not EHR-based.",
    )

    # ── narrative summary + findings ───────────────────────────────────────
    summary: str = Field(
        ..., description="2-4 sentence plain-language summary of the study."
    )
    key_findings: list[str] = Field(
        default_factory=list,
        description="Bullet list of the main results / conclusions.",
    )

    # ── EHR features / variables captured ─────────────────────────────────
    captured_features: list[str] = Field(
        default_factory=list,
        description="EHR features / variables the study captured (e.g. "
                    "'BMI percentile', 'ICD-9 asthma codes', 'lab panels'). "
                    "Empty list if the study is not EHR-based.",
    )

    # ── pathology / disease ───────────────────────────────────────────────
    pathologies_diseases: list[str] = Field(
        default_factory=list,
        description="Disease(s) / pathology / health outcome studied.",
    )

    # ── extended review fields ─────────────────────────────────────────────
    study_design: str = Field(
        default="",
        description="e.g. retrospective cohort, cross-sectional, EWAS, "
                    "case-control, birth cohort, ecological.",
    )
    data_source_type: str = Field(
        default="",
        description="Primary data source, e.g. 'EHR', 'claims', 'birth cohort', "
                    "'NHANES', 'linked administrative data', 'biospecimen cohort'.",
    )
    population: str = Field(
        default="",
        description="Study population, e.g. 'children 0-9 yrs across England'.",
    )
    exposure_domain: str = Field(
        default="",
        description="e.g. 'air pollution', 'metals/lead', 'chemicals/EDC', "
                    "'prenatal', 'SDOH/neighborhood'.",
    )
    limitations: list[str] = Field(
        default_factory=list,
        description="Stated or apparent study limitations.",
    )

    # ── data availability ─────────────────────────────────────────────────
    data_availability: Literal[
        "public-repository", "available-upon-request", "in-house",
        "supplementary-only", "not-stated"
    ] = Field(
        default="not-stated",
        description="How the study's data can be obtained. "
                    "public-repository = accession number / repository link "
                    "(dbGaP, GEO, ArrayExpress, Zenodo, GitHub, etc.); "
                    "available-upon-request = authors offer data on request; "
                    "in-house = institutional/private dataset, no public access; "
                    "supplementary-only = data only in supplementary files; "
                    "not-stated = no data-availability statement found.",
    )
    data_accession_links: list[str] = Field(
        default_factory=list,
        description="Accession IDs, repository URLs, or repository names "
                    "(e.g. 'dbGaP phs000123', 'GSE12345', 'github.com/x/y'). "
                    "Empty unless data_availability is public-repository.",
    )
    data_availability_statement: str = Field(
        default="",
        description="Verbatim sentence(s) from the manuscript that justify the "
                    "data_availability call; '' if not stated.",
    )
    confidence: Literal["high", "medium", "low", "unclear"] = Field(
        default="unclear",
        description="Confidence that the study fits a pediatric EHR/exposome "
                    "review, based on how clearly EHR + pediatric + exposure "
                    "are established.",
    )

    # ── provenance ─────────────────────────────────────────────────────────
    source_format: Literal["pdf", "xml", "text"] = Field(
        default="text",
        description="Format of the source full text that was summarized.",
    )
    model: str = Field(default="", description="LLM model id used for extraction.")

    @field_validator("data_availability", mode="before")
    @classmethod
    def _normalize_data_availability(cls, v):
        if isinstance(v, str):
            v2 = v.strip().lower().replace("_", "-").replace(" ", "-")
            valid = {"public-repository", "available-upon-request", "in-house",
                     "supplementary-only", "not-stated"}
            if v2 in valid:
                return v2
            return "not-stated"
        return v

    @field_validator("data_accession_links", mode="before")
    @classmethod
    def _coerce_accession_links(cls, v):
        if isinstance(v, str):
            return [v] if v.strip() else []
        if isinstance(v, list):
            return [str(x).strip() for x in v if str(x).strip()]
        return v

    @field_validator("pmcid")
    @classmethod
    def _normalize_pmcid(cls, v: str) -> str:
        v = v.strip().upper()
        if not v.startswith("PMC"):
            v = f"PMC{v}"
        return v

    @field_validator("ehr_used", mode="before")
    @classmethod
    def _coerce_ehr_used(cls, v):
        """Models often emit 'Yes'/'No'/'Y'/'N'/'true'/'false' for booleans."""
        if isinstance(v, bool):
            return v
        if isinstance(v, str):
            v2 = v.strip().lower()
            if v2 in {"true", "yes", "y", "1", "t", "used"}:
                return True
            if v2 in {"false", "no", "n", "0", "f", "not used"}:
                return False
        return v  # let Pydantic raise on anything else

    @field_validator("confidence", mode="before")
    @classmethod
    def _normalize_confidence(cls, v):
        if isinstance(v, str):
            v2 = v.strip().lower()
            if v2 in {"high", "medium", "low", "unclear"}:
                return v2
            return "unclear"
        return v

    @field_validator("key_findings", "captured_features",
                     "pathologies_diseases", "limitations", mode="before")
    @classmethod
    def _coerce_str_list(cls, v):
        """Accept a single string as a one-element list; drop empties."""
        if isinstance(v, str):
            return [v] if v.strip() else []
        if isinstance(v, list):
            return [str(x).strip() for x in v if str(x).strip()]
        return v


class SummaryBatch(BaseModel):
    """Combined file: a list of per-paper checklists + run metadata."""

    n: int
    model: str
    summaries: list[ManuscriptChecklist]


# ── JSON schema for the LLM prompt ──────────────────────────────────────────
# A trimmed schema handed to the model so it knows the expected shape. Identity
# / provenance fields are filled by us, so the model is only asked for the
# analytical fields.
LLM_FIELDS_SCHEMA = {
    "type": "object",
    "properties": {
        "ehr_used": {"type": "boolean"},
        "ehr_evidence": {"type": "string"},
        "summary": {"type": "string"},
        "key_findings": {"type": "array", "items": {"type": "string"}},
        "captured_features": {"type": "array", "items": {"type": "string"}},
        "pathologies_diseases": {"type": "array", "items": {"type": "string"}},
        "study_design": {"type": "string"},
        "data_source_type": {"type": "string"},
        "population": {"type": "string"},
        "exposure_domain": {"type": "string"},
        "limitations": {"type": "array", "items": {"type": "string"}},
        "confidence": {"type": "string", "enum": ["high", "medium", "low", "unclear"]},
    },
    "required": [
        "ehr_used", "ehr_evidence", "summary", "key_findings",
        "captured_features", "pathologies_diseases", "study_design",
        "data_source_type", "population", "exposure_domain", "limitations",
        "confidence",
    ],
    "additionalProperties": False,
}
