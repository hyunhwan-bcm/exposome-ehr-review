# Makefile for the pediatric EWAS / EHR literature collection.
#
#   make help        show available targets
#   make setup       create venv + install requests
#   make download    run the PMC fetcher (pediatric-focused)
#   make clean       remove downloaded PDFs + download log
#   make fresh       clean + download (start over)
#   make summary     print a compact inventory from download_log.json
#
# Override the interpreter / script via env vars if needed:
#   make download PYTHON=python3.11

PYTHON       ?= python3
SCRIPT       := fetch_pmc_papers.py
PAPERS_DIR   := papers
DOWNLOAD_LOG := $(PAPERS_DIR)/download_log.json
VENV         := .venv
PIP          := $(VENV)/bin/pip
VENV_PYTHON  := $(VENV)/bin/python

.PHONY: help setup download clean fresh summary check-venv

help: ## Show available targets
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

setup: $(VENV) ## Create virtualenv and install requests
	@$(PIP) install --quiet --upgrade pip
	@$(PIP) install --quiet requests
	@echo "✓ venv ready at $(VENV)"

$(VENV):
	@$(PYTHON) -m venv $(VENV)
	@echo "✓ created venv at $(VENV)"

download: $(DOWNLOAD_LOG) ## Run the pediatric-focused PMC fetcher

$(DOWNLOAD_LOG): $(SCRIPT) $(VENV)
	@$(VENV_PYTHON) $(SCRIPT)

clean: ## Remove downloaded PDFs/XML + download log
	@rm -f $(PAPERS_DIR)/*.pdf $(PAPERS_DIR)/*.xml $(DOWNLOAD_LOG)
	@echo "✓ cleaned $(PAPERS_DIR) (PDFs/XML + log removed)"

fresh: clean download ## Clean then re-download from scratch

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
print(f'\\nTotal candidates: {len(papers)} | Downloaded: {len(log.get(\"downloaded\",[]))} | ' \
      f'Abstract-only: {len(log.get(\"abstract_only\",[]))} | Excluded: {len(log.get(\"excluded\",[]))}')"
