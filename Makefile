# Makefile for the pediatric EWAS / EHR literature collection + summarizer.
#
#   make help        show available targets
#   make setup       create venv + install all deps
#   make download    run the pediatric-focused PMC fetcher
#   make clean       remove downloaded PDFs + download log
#   make fresh       clean + download (start over)
#   make summary     print a compact inventory + regenerate paper_summary.md
#   make summarize   summarize all manuscripts via Gemma 4 12B -> JSON
#   make results     export readable results/ (SUMMARY.md, checklist.md, combined JSON)
#   make summarize-paper PMC=PMC7145790   summarize one paper
#   make test        run the summarizer unit tests
#
# Override the interpreter via env vars if needed:
#   make download PYTHON=python3.11

PYTHON        ?= python3
SCRIPT        := fetch_pmc_papers.py
PAPERS_DIR    := papers
DOWNLOAD_LOG  := $(PAPERS_DIR)/download_log.json
COMBINED_JSON := $(PAPERS_DIR)/manuscript_summaries.json
VENV          := .venv
PIP           := $(VENV)/bin/pip
VENV_PYTHON   := $(VENV)/bin/python

.PHONY: help setup download clean fresh summary summarize summarize-paper test results summary db-import db-export db-stats dagster materialize

help: ## Show available targets
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-16s\033[0m %s\n", $$1, $$2}'

setup: $(VENV) ## Create virtualenv and install all deps
	@$(PIP) install --quiet --upgrade pip
	@$(PIP) install --quiet requests openai pypdf pydantic python-dotenv pytest tinydb dagster
	@echo "✓ venv ready at $(VENV)"

$(VENV):
	@$(PYTHON) -m venv $(VENV)
	@echo "✓ created venv at $(VENV)"

download: $(DOWNLOAD_LOG) ## Run the pediatric-focused PMC fetcher

$(DOWNLOAD_LOG): $(SCRIPT) $(VENV)
	@$(VENV_PYTHON) $(SCRIPT)

clean: ## Remove downloaded PDFs/XML + download log + summaries
	@rm -f $(PAPERS_DIR)/*.pdf $(PAPERS_DIR)/*.xml $(DOWNLOAD_LOG)
	@rm -rf $(PAPERS_DIR)/summaries $(COMBINED_JSON)
	@echo "✓ cleaned $(PAPERS_DIR) (PDFs/XML + log + summaries removed)"

fresh: clean download ## Clean then re-download from scratch

summarize: $(COMBINED_JSON) ## Summarize all manuscripts via Gemma 4 12B -> JSON

$(COMBINED_JSON): summarizer $(VENV)
	@if [ ! -f .env ] && [ -z "$$GEMMA_API_KEY" ]; then \
		echo "⚠ No GEMMA_API_KEY: copy .env.example to .env and fill in the key."; \
		exit 2; \
	fi
	@$(VENV_PYTHON) -m summarizer.run

summarize-paper: ## Summarize a single paper: make summarize-paper PMC=PMC7145790
	@test -n "$(PMC)" || (echo "Usage: make summarize-paper PMC=PMC7145790"; exit 2)
	@$(VENV_PYTHON) -m summarizer.run --pmcid $(PMC)

test: $(VENV) ## Run summarizer unit tests (no live API calls)
	@$(VENV_PYTHON) -m pytest tests/ -q

results: $(COMBINED_JSON) ## Export readable results/ (SUMMARY.md, checklist.md, combined JSON)
	@$(VENV_PYTHON) build_results.py
	@echo "✓ wrote results/ (SUMMARY.md, checklist.md, manuscript_summaries.json)"

db-import: $(VENV) ## Import existing per-paper JSONs (papers/summaries) into the TinyDB store
	@$(VENV_PYTHON) db.py import --from-dir papers/summaries

db-export: $(VENV) ## Export the TinyDB store -> combined + results/manuscript_summaries.json
	@$(VENV_PYTHON) db.py export
	@echo "✓ exported from TinyDB store (papers/db.json)"

db-stats: $(VENV) ## Print TinyDB store record counts
	@$(VENV_PYTHON) db.py stats

dagster: $(VENV) ## Launch the Dagster asset UI + lineage browser
	@$(VENV)/bin/dagster dev -m pipeline

materialize: $(VENV) ## Materialize the whole asset graph headlessly (fetch -> summarize -> results)
	@$(VENV)/bin/dagster asset materialize -m pipeline --select "*"

web: $(VENV) ## Serve the browsable review web app on http://localhost:8010
	@$(VENV_PYTHON) webapp.py

summary: $(DOWNLOAD_LOG) ## Print a compact inventory + regenerate paper_summary.md
	@$(VENV_PYTHON) build_summary.py
	@echo ""
	@$(VENV_PYTHON) -c "\
import json, pathlib; \
log = json.loads(pathlib.Path('$(DOWNLOAD_LOG)').read_text()); \
papers = log.get('papers', []); \
print(f'{\"#\":<4} {\"PMCID\":<14} {\"Yr\":<5} {\"Journal\":<22} Title'); \
print('-' * 100); \
[print(f'{i:<4} {p[\"pmcid\"]:<14} {p[\"year\"]:<5} {p[\"journal\"][:21]:<22} {p[\"title\"][:60]}') \
 for i, p in enumerate(papers, 1)]; \
print(f'\nTotal candidates: {len(papers)} | Downloaded: {len(log.get(\"downloaded\",[]))} | ' \
      f'Abstract-only: {len(log.get(\"abstract_only\",[]))} | Excluded: {len(log.get(\"excluded\",[]))}')"
