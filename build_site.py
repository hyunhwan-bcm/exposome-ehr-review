#!/usr/bin/env python3
"""Generate the static GitHub Pages site (``docs/index.html``) from the combined
manuscript summaries.

The page is a self-contained, Tailwind-styled inventory of every summarized
paper — searchable, filterable (EHR / public-repository), sortable, with an
expandable detail row per paper — plus headline stats and the pipeline figure.
Data is embedded inline so the page works when GitHub Pages serves only the
``docs/`` folder (no runtime fetch / CORS).

Usage:
    python build_site.py            # results/manuscript_summaries.json -> docs/index.html
"""
from __future__ import annotations

import collections
import html
import json
from pathlib import Path

REPO = "hyunhwan-bcm/exposome-ehr-review"
COMBINED = Path("results/manuscript_summaries.json")
OUT = Path("docs/index.html")

# Fields carried into the page (keep the payload lean).
KEEP = [
    "pmcid", "title", "year", "ehr_used", "summary", "key_findings",
    "pathologies_diseases", "data_source_type", "exposure_domain",
    "data_availability", "data_accession_links", "confidence",
]


def load() -> dict:
    return json.loads(COMBINED.read_text())


def build_records(summaries: list[dict]) -> list[dict]:
    recs = [{k: r.get(k) for k in KEEP} for r in summaries]
    recs.sort(key=lambda r: (r.get("year") or "", r.get("pmcid") or ""))
    return recs


def stats(summaries: list[dict]) -> dict:
    n = len(summaries)
    ehr = sum(1 for r in summaries if r.get("ehr_used"))
    avail = collections.Counter(r.get("data_availability") or "not-stated" for r in summaries)
    years = [r.get("year") for r in summaries if r.get("year")]
    return {
        "n": n, "ehr": ehr, "public": avail.get("public-repository", 0),
        "year_min": min(years) if years else "", "year_max": max(years) if years else "",
    }


# Placeholders use %%TOKEN%% so the literal {} in the embedded CSS/JS are safe.
PAGE = """<!doctype html>
<html lang="en" class="scroll-smooth">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pediatric Exposome / EWAS Literature Collection</title>
<meta name="description" content="Searchable inventory of %%N%% pediatric exposome / EWAS manuscripts summarized with Gemma 4 12B, %%EHR%% EHR-based, spanning %%YMIN%%-%%YMAX%%.">
<link rel="stylesheet" href="tailwind.css">
</head>
<body class="bg-slate-50 text-slate-800 dark:bg-slate-950 dark:text-slate-200 antialiased">

<header class="border-b border-slate-200 dark:border-slate-800 bg-white/80 dark:bg-slate-900/60 backdrop-blur sticky top-0 z-10">
  <div class="max-w-6xl mx-auto px-4 py-4 flex items-center justify-between gap-4">
    <div>
      <h1 class="text-lg sm:text-xl font-bold tracking-tight">Pediatric Exposome / EWAS Literature Collection</h1>
      <p class="text-xs sm:text-sm text-slate-500 dark:text-slate-400">Summarized with Gemma 4&nbsp;12B · %%YMIN%%–%%YMAX%%</p>
    </div>
    <a href="https://github.com/%%REPO%%" class="shrink-0 text-sm font-medium text-blue-600 dark:text-blue-400 hover:underline">GitHub&nbsp;↗</a>
  </div>
</header>

<main class="max-w-6xl mx-auto px-4 py-8 space-y-10">

  <section class="grid grid-cols-2 sm:grid-cols-4 gap-4">
    <div class="rounded-xl bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-5">
      <div class="text-3xl font-bold text-blue-600 dark:text-blue-400">%%N%%</div>
      <div class="text-sm text-slate-500 dark:text-slate-400 mt-1">manuscripts</div>
    </div>
    <div class="rounded-xl bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-5">
      <div class="text-3xl font-bold text-emerald-600 dark:text-emerald-400">%%EHR%%</div>
      <div class="text-sm text-slate-500 dark:text-slate-400 mt-1">EHR-based</div>
    </div>
    <div class="rounded-xl bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-5">
      <div class="text-3xl font-bold text-violet-600 dark:text-violet-400">%%PUBLIC%%</div>
      <div class="text-sm text-slate-500 dark:text-slate-400 mt-1">public repository</div>
    </div>
    <div class="rounded-xl bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-5">
      <div class="text-3xl font-bold text-amber-600 dark:text-amber-400">%%YMIN%%–%%YMAX%%</div>
      <div class="text-sm text-slate-500 dark:text-slate-400 mt-1">year span</div>
    </div>
  </section>

  <section class="rounded-xl bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-6">
    <h2 class="text-base font-semibold mb-4">Pipeline</h2>
    <figure class="flex flex-col items-center">
      <img src="pipeline.png" alt="End-to-end literature pipeline diagram" class="max-w-full sm:max-w-lg rounded-lg border border-slate-200 dark:border-slate-800">
      <figcaption class="text-sm text-slate-500 dark:text-slate-400 mt-4 max-w-2xl text-center">
        <strong>Figure 1.</strong> PubMed Central → <code>make download</code> → <code>make summarize</code>
        (⇄ Gemma&nbsp;4&nbsp;12B) → <code>scan_data_availability.py</code> → TinyDB → <code>make results</code>.
      </figcaption>
    </figure>
  </section>

  <section>
    <div class="flex flex-wrap items-center gap-3 mb-4">
      <h2 class="text-base font-semibold mr-auto">Inventory</h2>
      <input id="q" type="search" placeholder="Search title, PMCID, disease, exposure…"
        class="w-full sm:w-80 rounded-lg border border-slate-300 dark:border-slate-700 bg-white dark:bg-slate-900 px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500">
      <div class="flex gap-1 text-sm" id="filters">
        <button data-f="all"    class="filter px-3 py-2 rounded-lg bg-blue-600 text-white">All</button>
        <button data-f="ehr"    class="filter px-3 py-2 rounded-lg bg-slate-200 dark:bg-slate-800">EHR</button>
        <button data-f="public" class="filter px-3 py-2 rounded-lg bg-slate-200 dark:bg-slate-800">Public repo</button>
      </div>
    </div>

    <div class="overflow-x-auto rounded-xl border border-slate-200 dark:border-slate-800">
      <table class="min-w-full text-sm">
        <thead class="bg-slate-100 dark:bg-slate-900 text-slate-600 dark:text-slate-400">
          <tr>
            <th class="text-left font-medium px-4 py-3 cursor-pointer select-none" data-sort="year">Year ▾</th>
            <th class="text-left font-medium px-4 py-3">PMCID</th>
            <th class="text-left font-medium px-4 py-3">Title</th>
            <th class="text-left font-medium px-4 py-3">Exposure</th>
            <th class="text-left font-medium px-4 py-3">EHR</th>
            <th class="text-left font-medium px-4 py-3">Data availability</th>
          </tr>
        </thead>
        <tbody id="rows" class="divide-y divide-slate-200 dark:divide-slate-800"></tbody>
      </table>
    </div>
    <p id="count" class="text-xs text-slate-500 dark:text-slate-400 mt-3"></p>
  </section>

  <section class="rounded-xl bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 p-6">
    <h2 class="text-base font-semibold mb-3">Raw results</h2>
    <ul class="text-sm space-y-1 text-blue-600 dark:text-blue-400">
      <li><a class="hover:underline" href="https://github.com/%%REPO%%/blob/main/results/SUMMARY.md">results/SUMMARY.md</a> <span class="text-slate-500 dark:text-slate-400">— headline stats &amp; top diseases/features</span></li>
      <li><a class="hover:underline" href="https://github.com/%%REPO%%/blob/main/results/checklist.md">results/checklist.md</a> <span class="text-slate-500 dark:text-slate-400">— full per-paper checklist</span></li>
      <li><a class="hover:underline" href="https://github.com/%%REPO%%/blob/main/results/manuscript_summaries.json">results/manuscript_summaries.json</a> <span class="text-slate-500 dark:text-slate-400">— combined machine-readable JSON</span></li>
    </ul>
  </section>

</main>

<footer class="border-t border-slate-200 dark:border-slate-800 py-6 text-center text-xs text-slate-500 dark:text-slate-400">
  Generated by <code>build_site.py</code> from the Gemma&nbsp;4&nbsp;12B Pydantic checklists.
</footer>

<script id="data" type="application/json">%%DATA%%</script>
<script>
const DATA = JSON.parse(document.getElementById('data').textContent);
const rowsEl = document.getElementById('rows');
const countEl = document.getElementById('count');
const qEl = document.getElementById('q');
let filter = 'all', sortDir = -1, sortKey = 'year', query = '';

const esc = s => (s ?? '').toString().replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const availStyle = {
  'public-repository':      'bg-violet-100 text-violet-800 dark:bg-violet-900/40 dark:text-violet-300',
  'available-upon-request': 'bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300',
  'supplementary-only':     'bg-sky-100 text-sky-800 dark:bg-sky-900/40 dark:text-sky-300',
  'in-house':               'bg-rose-100 text-rose-800 dark:bg-rose-900/40 dark:text-rose-300',
  'not-stated':             'bg-slate-100 text-slate-600 dark:bg-slate-800 dark:text-slate-400',
};

function matches(r) {
  if (filter === 'ehr' && !r.ehr_used) return false;
  if (filter === 'public' && r.data_availability !== 'public-repository') return false;
  if (!query) return true;
  const hay = [r.pmcid, r.title, r.exposure_domain, (r.pathologies_diseases||[]).join(' ')].join(' ').toLowerCase();
  return hay.includes(query);
}

function render() {
  let list = DATA.filter(matches);
  list.sort((a, b) => ((a[sortKey]||'') > (b[sortKey]||'') ? 1 : -1) * sortDir);
  rowsEl.innerHTML = list.map((r, i) => {
    const av = r.data_availability || 'not-stated';
    const findings = (r.key_findings||[]).map(f => `<li>${esc(f)}</li>`).join('');
    const links = (r.data_accession_links||[]).map(l => `<code class="text-xs">${esc(l)}</code>`).join(' · ');
    return `
    <tr class="hover:bg-slate-50 dark:hover:bg-slate-900/50 cursor-pointer align-top" data-i="${i}">
      <td class="px-4 py-3 whitespace-nowrap tabular-nums">${esc(r.year)}</td>
      <td class="px-4 py-3 whitespace-nowrap"><a href="https://www.ncbi.nlm.nih.gov/pmc/articles/${esc(r.pmcid)}/" target="_blank" rel="noopener" class="text-blue-600 dark:text-blue-400 hover:underline" onclick="event.stopPropagation()">${esc(r.pmcid)}</a></td>
      <td class="px-4 py-3 max-w-md">${esc(r.title)}</td>
      <td class="px-4 py-3 text-slate-500 dark:text-slate-400 whitespace-nowrap">${esc(r.exposure_domain||'—')}</td>
      <td class="px-4 py-3">${r.ehr_used ? '<span class="inline-block px-2 py-0.5 rounded-full text-xs bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300">EHR</span>' : '<span class="text-slate-400">—</span>'}</td>
      <td class="px-4 py-3"><span class="inline-block px-2 py-0.5 rounded-full text-xs ${availStyle[av]}">${esc(av)}</span></td>
    </tr>
    <tr class="detail hidden bg-slate-50 dark:bg-slate-900/40" data-detail="${i}">
      <td colspan="6" class="px-4 py-4">
        <p class="text-sm mb-2">${esc(r.summary)}</p>
        ${findings ? `<ul class="list-disc list-inside text-sm text-slate-600 dark:text-slate-400 space-y-0.5">${findings}</ul>` : ''}
        <div class="mt-2 text-xs text-slate-500 dark:text-slate-400 space-y-0.5">
          ${(r.pathologies_diseases||[]).length ? `<div><span class="font-medium">Diseases:</span> ${esc((r.pathologies_diseases||[]).join(', '))}</div>` : ''}
          ${r.data_source_type ? `<div><span class="font-medium">Data source:</span> ${esc(r.data_source_type)}</div>` : ''}
          ${r.confidence ? `<div><span class="font-medium">Confidence:</span> ${esc(r.confidence)}</div>` : ''}
          ${links ? `<div><span class="font-medium">Accession:</span> ${links}</div>` : ''}
        </div>
      </td>
    </tr>`;
  }).join('');
  countEl.textContent = `${list.length} of ${DATA.length} papers`;
}

rowsEl.addEventListener('click', e => {
  const tr = e.target.closest('tr[data-i]');
  if (!tr) return;
  const d = rowsEl.querySelector(`tr[data-detail="${tr.dataset.i}"]`);
  if (d) d.classList.toggle('hidden');
});
qEl.addEventListener('input', () => { query = qEl.value.trim().toLowerCase(); render(); });
document.querySelectorAll('.filter').forEach(b => b.addEventListener('click', () => {
  filter = b.dataset.f;
  document.querySelectorAll('.filter').forEach(x => { x.className = 'filter px-3 py-2 rounded-lg bg-slate-200 dark:bg-slate-800'; });
  b.className = 'filter px-3 py-2 rounded-lg bg-blue-600 text-white';
  render();
}));
document.querySelector('[data-sort="year"]').addEventListener('click', () => { sortDir *= -1; render(); });

render();
</script>
</body>
</html>
"""


def main() -> None:
    data = load()
    summaries = data["summaries"] if isinstance(data, dict) else data
    st = stats(summaries)
    recs = build_records(summaries)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    page = (PAGE
            .replace("%%REPO%%", REPO)
            .replace("%%N%%", str(st["n"]))
            .replace("%%EHR%%", str(st["ehr"]))
            .replace("%%PUBLIC%%", str(st["public"]))
            .replace("%%YMIN%%", html.escape(st["year_min"]))
            .replace("%%YMAX%%", html.escape(st["year_max"]))
            .replace("%%DATA%%", json.dumps(recs, ensure_ascii=False)))
    OUT.write_text(page)
    print(f"Wrote {OUT} ({OUT.stat().st_size // 1024} KB) — {st['n']} papers, "
          f"{st['ehr']} EHR, {st['public']} public-repo.")


if __name__ == "__main__":
    main()
