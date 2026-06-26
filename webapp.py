#!/usr/bin/env python3
"""Minimal stdlib-only web app to browse the exposome review collection.

No dependencies. Serves the tracked combined JSON + the TinyDB store live.
Reading data on every request is fine at this scale (hundreds of records).

  python webapp.py                 # http://localhost:8010
  PORT=9000 python webapp.py       # custom port
  OPEN=0 python webapp.py          # don't auto-open a browser

Routes:
  /                - browsable HTML table of all papers
  /api/summaries   - raw combined JSON
  /api/stats       - counts (by ehr_used, data_availability, exposure, confidence)
  /api/summary/PMC1234567 - one paper's record (from the TinyDB store)
  /search?q=asthma          - filter the table by title/summary/pmcid/disease
"""
from __future__ import annotations

import json
import os
import sys
import html
import webbrowser
from http.server import HTTPServer, BaseHTTPRequestHandler
from pathlib import Path
from urllib.parse import urlparse, parse_qs

REPO_ROOT = Path(__file__).resolve().parent
COMBINED = REPO_ROOT / "results" / "manuscript_summaries.json"
DB_PATH = REPO_ROOT / "papers" / "db.json"

PORT = int(os.environ.get("PORT", "8010"))
AUTO_OPEN = os.environ.get("OPEN", "1") != "0"


def load_combined() -> dict:
    """Read the tracked combined JSON (the review artifact)."""
    if not COMBINED.exists():
        return {"n": 0, "model": "", "summaries": []}
    return json.loads(COMBINED.read_text())


def store_record(pmcid: str) -> dict | None:
    """One record from the live TinyDB store (source of truth)."""
    if not DB_PATH.exists():
        return None
    # ponytail: lazy import — only needed when asked, keeps startup cheap
    sys.path.insert(0, str(REPO_ROOT))
    from database import Store  # type: ignore
    with Store(DB_PATH) as s:
        rec = s.get(pmcid)
    return rec


def _is_vaccine(r: dict) -> bool:
    haystack = " ".join(str(r.get(k, "")) for k in ("title", "summary", "exposure_domain")).lower()
    return "vacc" in haystack or "immuniz" in haystack or "immunis" in haystack


def stats(combined: dict) -> dict:
    rows = combined.get("summaries", [])
    from collections import Counter
    return {
        "n": len(rows),
        "ehr_used_true": sum(1 for r in rows if r.get("ehr_used") is True),
        "ehr_used_false": sum(1 for r in rows if r.get("ehr_used") is False),
        "vaccine_like": sum(1 for r in rows if _is_vaccine(r)),
        "by_data_availability": dict(Counter(r.get("data_availability", "not-stated") for r in rows)),
        "by_exposure_domain": dict(Counter(r.get("exposure_domain", "unclear") for r in rows)),
        "by_confidence": dict(Counter(r.get("confidence", "unclear") for r in rows)),
    }


ESCAPE = html.escape


def _chip(cat: str) -> str:
    color = {
        "public-repository": "#1a7f37",
        "available-upon-request": "#bf8700",
        "in-house": "#8250df",
        "supplementary-only": "#0969da",
        "not-stated": "#636c76",
    }.get(cat, "#636c76")
    return f'<span style="color:#fff;background:ឧ{color};padding:2px 6px;border-radius:3px;font-size:11px">{ESCAPE(cat)}</span>'.replace("ឧ", "")


def render_table(combined: dict, q: str = "") -> str:
    """Browsable HTML table, filtered by q (title/summary/pmcid/disease)."""
    rows = combined.get("summaries", [])
    rows.sort(key=lambda r: (str(r.get("year", "")), str(r.get("pmcid", ""))))
    if q:
        ql = q.lower()
        rows = [r for r in rows if ql in (
            f"{r.get('title','')} {r.get('summary','')} "
            f"{r.get('pmcid','')} "
            f"{' '.join(r.get('pathologies_diseases',[]))} "
            f"{r.get('data_source_type','')} "
            f"{r.get('exposure_domain','')} "
            f"{' '.join(r.get('data_accession_links',[]))}"
        ).lower()]

    def links_cell(r: dict) -> str:
        links = r.get("data_accession_links", [])
        if not links:
            return "—"
        return "<br>".join(f'<a href="{ESCAPE(l)}" target="_blank" style="font-size:11px">{ESCAPE(l)}</a>' for l in links)

    def diseases(r: dict) -> str:
        return ESCAPE("; ".join(r.get("pathologies_diseases", []))) or "—"

    body_rows = "\n".join(
        f"<tr>"
        f'<td style="font-family:monospace"><a href="/api/summary/{r.get("pmcid","")}">{ESCAPE(str(r.get("pmcid","")))}</a></td>'
        f"<td>{ESCAPE(str(r.get('year','')))}</td>"
        f'<td>{r.get("ehr_used") and "✓ EHR" or ""}</td>'
        f"<td>{ESCAPE(str(r.get('title','')))[:90]}</td>"
        f"<td>{_chip(r.get('data_availability','not-stated'))}</td>"
        f"<td>{links_cell(r)}</td>"
        f"<td>{diseases(r)}</td>"
        f"<td>{ESCAPE(str(r.get('exposure_domain','')))[:24]}</td>"
        f"<td>{ESCAPE(str(r.get('confidence','')))}</td>"
        f"</tr>"
        for r in rows
    )

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Exposome / EWAS Literature Collection</title>
<style>
 body {{ font-family: -apple-system, system-ui, sans-serif; margin: 0; color: #1f2328; background: #fff; }}
 header {{ background: #1f2328; color: #fff; padding: 18px 28px; }}
 header h1 {{ margin: 0 0 4px; font-size: 19px; }}
 header p {{ margin: 0; opacity: .8; font-size: 13px; }}
 nav {{ padding: 10px 28px; border-bottom: 1px solid #d0d7de; background: #f6f8fa; }}
 nav a {{ margin-right: 16px; color: #0969da; text-decoration: none; font-size: 13px; }}
 nav a:hover {{ text-decoration: underline; }}
 main {{ padding: 16px 28px; }}
 .stats {{ display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 16px; }}
 .stat {{ background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px; padding: 8px 12px; font-size: 13px; }}
 .stat b {{ font-size: 18px; display: block; }}
 form {{ margin: 12px 0 16px; }}
 input {{ padding: 6px 10px; width: 380px; border: 1px solid #d0d7de; border-radius: 6px; font-size: 13px; }}
 button {{ padding: 6px 14px; border: 1px solid #d0d7de; background: #fff; border-radius: 6px; cursor: pointer; font-size: 13px; }}
 table {{ border-collapse: collapse; width: 100%; font-size: 12.5px; }}
 th, td {{ border: 1px solid #d0d7de; padding: 6px 8px; text-align: left; vertical-align: top; }}
 th {{ background: #f6f8fa; position: sticky; top: 0; }}
 tr:hover {{ background: #f6f8fa; }}
 code {{ font-family: monospace; }}
 pre {{ background: #f6f8fa; padding: 14px; border-radius: 6px; overflow: auto; font-size: 12px; }}
</style></head><body>
<header>
  <h1>Pediatric Exposome / EWAS Literature Collection</h1>
  <p>{len(rows)} papers shown · {sum(1 for r in rows if r.get('ehr_used'))} EHR-based · {sum(1 for r in rows if _is_vaccine(r))} vaccine/immunization · combined JSON at <code>/api/summaries</code></p>
</header>
<nav>
  <a href="/">Table</a>
  <a href="/api/summaries">Combined JSON</a>
  <a href="/api/stats">Stats</a>
  <a href="/search?q=vaccine">Vaccines</a>
</nav>
<main>
  <form method="get" action="/search">
    <input name="q" placeholder="filter: asthma, lead, PMC1234567, github, zenodo…" value="{ESCAPE(q)}">
    <button>Filter</button>
  </form>
  <table>
    <thead><tr>
      <th>PMCID</th><th>Year</th><th>EHR</th><th>Title</th>
      <th>Data avail.</th><th>Accession / links</th><th>Diseases</th><th>Exposure</th><th>Conf.</th>
    </tr></thead>
    <tbody>{body_rows or '<tr><td colspan="9">no papers</td></tr>'}</tbody>
  </table>
</main>
</body></html>"""


class Handler(BaseHTTPRequestHandler):
    def _json(self, obj, code=200):
        body = json.dumps(obj, indent=2).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _html(self, s, code=200):
        body = s.encode()
        self.send_response(code)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        u = urlparse(self.path)
        qs = parse_qs(u.query)
        if u.path == "/" or u.path == "/search":
            q = qs.get("q", [""])[0]
            self._html(render_table(load_combined(), q))
        elif u.path == "/api/summaries":
            self._json(load_combined())
        elif u.path == "/api/stats":
            self._json(stats(load_combined()))
        elif u.path.startswith("/api/summary/"):
            pmcid = u.path.split("/api/summary/", 1)[1]
            rec = store_record(pmcid)
            if rec is None:
                self._json({"error": "not found", "pmcid": pmcid}, 404)
            else:
                self._html(f"<pre>{ESCAPE(json.dumps(rec, indent=2, sort_keys=True))}</pre>")
        else:
            self._json({"error": "not found"}, 404)

    def log_message(self, fmt, *a):
        sys.stderr.write(f"[{self.log_date_time_string()}] {fmt % a}\n")


def main() -> int:
    srv = HTTPServer(("0.0.0.0", PORT), Handler)
    print(f"exposome review web app -> http://localhost:{PORT}")
    print(f"  data: {COMBINED}")
    print("  ^C to stop")
    if AUTO_OPEN:
        webbrowser.open(f"http://localhost:{PORT}")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())