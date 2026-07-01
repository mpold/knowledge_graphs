#!/usr/bin/env python3
"""Download PDFs for the high-impact PMC articles whose XML has **no <body>**.

This script is a replica of ``bacs/doi_openalex.py`` -- same download strategy and
helpers -- repointed at a new input and output:

  * INPUT  : the ``PMC*.xml`` files in ``high_impact_xmls/`` that lack a
             ``<body>`` element (front-matter / abstract-only records, per
             xml_structure.py). Those carry no full-text JATS, so we fetch the PDF
             instead, to be converted to XML downstream by GROBID.
  * OUTPUT : ``ncbi_pdfs_grobid/PMC*.pdf`` + a download-summary TSV.

Strategy (try in order, stop on first success) -- identical to doi_openalex.py:
1. PMC PoW   -- scrape the PMC article page for the PDF filename, solve the
   Proof-of-Work challenge, and download with the solved cookie via curl
2. OpenAlex  -- look up OA PDF URLs from the OpenAlex API, download via curl
3. CrossRef  -- follow full-text PDF links from CrossRef metadata via curl

All HTTP downloads use subprocess+curl to avoid TLS fingerprint blocking. PMCIDs
are resolved to PMIDs + DOIs via the NCBI ID Converter (seeded from the PMID/DOI
already present in each XML's front matter). Finished PDFs on disk are skipped, so
re-running resumes.

AUTO-RECOVERY (resume after interruption / time-out / crash)
------------------------------------------------------------
* Resume: finished ``<PMCID>.pdf`` files are skipped, so a re-run continues where
  it stopped.
* Append-only journals updated *live* (flushed + fsync'd after every article):
    - ``_downloaded.tsv`` -- PMCIDs that succeeded (+ which source worked),
    - ``_failed.tsv``     -- PMCIDs that failed every source.
  Both downloaded and failed PMCIDs are skipped on the next run (re-try the
  failures with ``RETRY_FAILED=1``), so no work is repeated after a restart.
* Atomic writes: every PDF is downloaded to a ``.part`` file, validated
  (``%PDF-`` header, >1 KB), then ``os.replace``'d into place -- a killed curl
  never leaves a truncated PDF that looks complete. Stray ``.part`` files are
  swept on start-up.
* Wall-clock cap: ``TIME_BUDGET=<seconds>`` makes a run self-exit cleanly before a
  CI/command time-out; re-run to finish (0 = unlimited, the default).
* Ctrl-C / unexpected exit is caught: the summary TSV is always (re)written in a
  ``finally`` block, reconstructed from the journals + the PDFs on disk, so it is
  correct and complete no matter how many times the run was interrupted/resumed.

NOT executed here -- run it yourself.

CAVEAT: many of these records are publisher-restricted and not openly available,
so a PDF may be unobtainable through all three sources; those rows are marked
``failed`` in the summary TSV.
"""

import sys
import os
import re
import csv
import json
import glob
import hashlib
import subprocess
import time

import requests

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

# Paths are resolved relative to this file so the script works from any cwd.
BASE = os.path.dirname(os.path.abspath(__file__))
XML_DIR = os.path.join(BASE, "high_impact_xmls")          # source XMLs
OUTPUT_DIR = os.path.join(BASE, "ncbi_pdfs_grobid")       # PDF export target
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "ncbi_pdf_download_summary.tsv")
DONE_LOG = os.path.join(OUTPUT_DIR, "_downloaded.tsv")     # append-only success journal
FAILED_LOG = os.path.join(OUTPUT_DIR, "_failed.tsv")       # append-only failure journal
SUMMARY_DIR = os.path.join(BASE, "summaries")
SUMMARY_HTML = os.path.join(SUMMARY_DIR, "ncbi_pdf.html")
BATCH_SIZE = 200
EMAIL = "your-email@example.com"
DELAY = 0.35

# Auto-recovery knobs (env-overridable).
TIME_BUDGET = float(os.environ.get("TIME_BUDGET", "0"))    # seconds; 0 = unlimited
RETRY_FAILED = os.environ.get("RETRY_FAILED", "").strip() in ("1", "true", "yes")

CURL_UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
           "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36")

RE_BODY = re.compile(r"<body\b")
LOG_HEADER = ["PMCID", "PMID", "DOI", "Status", "Source"]


# ---------------------------------------------------------------------------
# Auto-recovery journals (append-only, flushed live so progress survives a kill)
# ---------------------------------------------------------------------------

def append_log(path, row):
    """Append one row to an append-only TSV journal, flushed + fsync'd."""
    new = not os.path.exists(path)
    with open(path, "a", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        if new:
            w.writerow(LOG_HEADER)
        w.writerow(row)
        f.flush()
        os.fsync(f.fileno())


def load_log(path):
    """Return {pmcid: row} from a journal written by append_log (skip header)."""
    out = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            for row in csv.reader(f, delimiter="\t"):
                if row and row[0] not in ("", "PMCID", "Total_attempted", "Total"):
                    out[row[0]] = row
    return out


def sweep_partials(output_dir):
    """Delete stray *.part files left by a killed curl/download."""
    for p in glob.glob(os.path.join(output_dir, "*.part")):
        try:
            os.remove(p)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _xml_article_id(text, idtype):
    """Pull an <article-id pub-id-type="idtype"> value from XML front matter."""
    m = re.search(r'<article-id\b[^>]*pub-id-type="%s"[^>]*>([^<]+)</article-id>' % idtype,
                  text)
    return m.group(1).strip() if m else ""


def load_no_body_records(xml_dir):
    """Select the PMC*.xml files lacking a <body> and read their ids.

    Replaces doi_openalex.py's ``load_pmcid_records`` (which parsed a markdown
    table). Returns ``(records, xml_meta)`` where ``records`` is a list of
    ``(pmid, pmcid)`` tuples (PMID parsed from the XML when present, else "") and
    ``xml_meta`` maps ``pmcid -> {"pmid":..., "doi":...}`` from the front matter.
    """
    records = []
    xml_meta = {}
    seen = set()
    for path in sorted(glob.glob(os.path.join(xml_dir, "PMC*.xml"))):
        try:
            with open(path, encoding="utf-8", errors="ignore") as f:
                text = f.read()
        except OSError:
            continue
        if RE_BODY.search(text):
            continue                                       # has full text -> not a target
        pmcid = os.path.splitext(os.path.basename(path))[0]   # 'PMC1014204'
        if pmcid in seen:
            continue
        seen.add(pmcid)
        pmid = _xml_article_id(text, "pmid")
        doi = _xml_article_id(text, "doi")
        records.append((pmid, pmcid))
        xml_meta[pmcid] = {"pmid": pmid, "doi": doi}
    return records, xml_meta


def get_already_downloaded(output_dir):
    done = set()
    if os.path.isdir(output_dir):
        for fn in os.listdir(output_dir):
            m = re.match(r"(PMC\d+)\.pdf$", fn)
            if m:
                done.add(m.group(1))
    return done


def is_valid_pdf_file(filepath):
    try:
        if not os.path.exists(filepath):
            return False
        if os.path.getsize(filepath) < 1000:
            os.remove(filepath)
            return False
        with open(filepath, "rb") as f:
            header = f.read(5)
        if header != b"%PDF-":
            os.remove(filepath)
            return False
        return True
    except Exception:
        return False


def curl_download(url, filepath, extra_headers=None, cookies=None):
    """Download a URL to filepath using curl. Returns True if valid PDF.

    Auto-recovery: curl writes to a ``.part`` temp; only a validated PDF is
    atomically moved into place, so an interrupted download never leaves a
    truncated file that a later resume would mistake for "done".
    """
    tmp = filepath + ".part"
    cmd = ["curl", "-sL", "--max-time", "120", "-H", f"User-Agent: {CURL_UA}"]
    if extra_headers:
        for k, v in extra_headers.items():
            cmd.extend(["-H", f"{k}: {v}"])
    if cookies:
        cookie_str = "; ".join(f"{k}={v}" for k, v in cookies.items())
        cmd.extend(["-b", cookie_str])
    cmd.extend(["--output", tmp, url])
    try:
        subprocess.run(cmd, capture_output=True, timeout=130)
        if is_valid_pdf_file(tmp):          # validates header/size, removes if bad
            os.replace(tmp, filepath)
            return True
        if os.path.exists(tmp):
            os.remove(tmp)
        return False
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        return False


def curl_get_text(url, extra_headers=None):
    """GET a URL with curl and return response body as string."""
    cmd = ["curl", "-sL", "--max-time", "30", "-H", f"User-Agent: {CURL_UA}"]
    if extra_headers:
        for k, v in extra_headers.items():
            cmd.extend(["-H", f"{k}: {v}"])
    cmd.append(url)
    try:
        result = subprocess.run(cmd, capture_output=True, timeout=35)
        return result.stdout.decode("utf-8", errors="replace")
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Step 0: Batch resolve PMCIDs → PMIDs + DOIs via NCBI ID Converter
# ---------------------------------------------------------------------------

def resolve_ids(pmcids):
    """Convert PMCIDs to PMIDs and DOIs using the NCBI ID Converter API."""
    pmcid_to_pmid = {}
    pmcid_to_doi = {}
    for i in range(0, len(pmcids), BATCH_SIZE):
        batch = pmcids[i:i + BATCH_SIZE]
        ids_str = ",".join(batch)
        url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
        params = {"ids": ids_str, "format": "json", "tool": "ncbi_pdf", "email": EMAIL}
        try:
            resp = requests.get(url, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            for record in data.get("records", []):
                pmcid = record.get("pmcid", "")
                pmid = str(record.get("pmid", ""))
                doi = record.get("doi", "")
                if pmcid and pmid:
                    pmcid_to_pmid[pmcid] = pmid
                if pmcid and doi:
                    pmcid_to_doi[pmcid] = doi
        except Exception as e:
            print(f"  Error in ID conversion batch {i + 1}: {e}")
        pmid_found = sum(1 for p in batch if p in pmcid_to_pmid)
        doi_found = sum(1 for p in batch if p in pmcid_to_doi)
        print(f"  Batch {i + 1}-{min(i + BATCH_SIZE, len(pmcids))}: "
              f"{pmid_found} PMIDs, {doi_found} DOIs")
        time.sleep(0.5)
    return pmcid_to_pmid, pmcid_to_doi


# ---------------------------------------------------------------------------
# Source 1: PMC with Proof-of-Work challenge solver
# ---------------------------------------------------------------------------

def solve_pow(challenge, difficulty):
    """Find nonce where sha256(challenge + nonce) starts with '0' * difficulty."""
    prefix = "0" * difficulty
    nonce = 0
    while True:
        h = hashlib.sha256((challenge + str(nonce)).encode()).hexdigest()
        if h.startswith(prefix):
            return nonce
        nonce += 1


def get_pmc_pdf_filename(pmcid):
    """Scrape the PMC article page to find the actual PDF filename."""
    html = curl_get_text(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/")
    if not html:
        return None
    matches = re.findall(r'href="([^"]*\.pdf)"', html)
    for m in matches:
        if m.startswith("pdf/") or "/pdf/" in m:
            if m.startswith("pdf/"):
                return m
            parts = m.split("/pdf/")
            if len(parts) == 2:
                return "pdf/" + parts[1]
    if matches:
        m = matches[0]
        if m.startswith("/"):
            return None
        return m
    return None


def try_pmc_pow(pmcid, filepath):
    """Download PDF from PMC by solving the Proof-of-Work challenge."""
    pdf_filename = get_pmc_pdf_filename(pmcid)
    if not pdf_filename:
        pdf_filename = "pdf/"

    pdf_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/{pdf_filename}"

    html = curl_get_text(pdf_url)
    if not html:
        return False, "pmc_pow_no_response"

    if html.startswith("%PDF-"):
        tmp = filepath + ".part"
        with open(tmp, "wb") as f:
            f.write(html.encode("latin-1"))
        if is_valid_pdf_file(tmp):          # removes tmp if not a valid PDF
            os.replace(tmp, filepath)
            return True, "pmc_direct"
        if os.path.exists(tmp):
            os.remove(tmp)

    challenge_m = re.search(r'POW_CHALLENGE\s*=\s*"([^"]+)"', html)
    difficulty_m = re.search(r'POW_DIFFICULTY\s*=\s*"(\d+)"', html)
    cookie_name_m = re.search(r'POW_COOKIE_NAME\s*=\s*"([^"]+)"', html)

    if not challenge_m:
        if "403" in html[:200]:
            return False, "pmc_403"
        return False, "pmc_no_challenge"

    challenge = challenge_m.group(1)
    difficulty = int(difficulty_m.group(1))
    cookie_name = cookie_name_m.group(1)

    nonce = solve_pow(challenge, difficulty)
    cookie_value = f"{challenge},{nonce}"

    if curl_download(pdf_url, filepath, cookies={cookie_name: cookie_value}):
        return True, "pmc_pow"

    if pdf_filename != "pdf/":
        pdf_url2 = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/pdf/"
        if curl_download(pdf_url2, filepath, cookies={cookie_name: cookie_value}):
            return True, "pmc_pow_generic"

    return False, "pmc_pow_download_failed"


# ---------------------------------------------------------------------------
# Source 2: OpenAlex → curl download
# ---------------------------------------------------------------------------

def try_openalex(doi, pmid, pmcid, filepath):
    """Look up OA PDF URL from OpenAlex API and download via curl."""
    identifiers = []
    if doi:
        identifiers.append(f"https://doi.org/{doi}")
    if pmid:
        identifiers.append(f"pmid:{pmid}")
    if pmcid:
        identifiers.append(f"pmcid:{pmcid}")

    for ident in identifiers:
        url = f"https://api.openalex.org/works/{ident}?mailto={EMAIL}"
        try:
            text = curl_get_text(url)
            if not text or text.startswith("<!"):
                continue
            data = json.loads(text)
        except Exception:
            continue

        candidates = []
        oa = data.get("open_access", {})
        oa_url = oa.get("oa_url")
        if oa_url:
            candidates.append(oa_url)

        best_loc = data.get("best_oa_location") or {}
        for key in ("pdf_url", "landing_page_url"):
            u = best_loc.get(key)
            if u and u not in candidates:
                candidates.append(u)

        for loc in data.get("locations", []):
            if not loc.get("is_oa"):
                continue
            for key in ("pdf_url", "landing_page_url"):
                u = loc.get(key)
                if u and u not in candidates:
                    candidates.append(u)

        for cand in candidates:
            if curl_download(cand, filepath):
                return True, "openalex"
            time.sleep(0.2)

        time.sleep(DELAY)

    return False, "openalex_no_pdf"


# ---------------------------------------------------------------------------
# Source 3: CrossRef full-text links → curl download
# ---------------------------------------------------------------------------

def try_crossref(doi, filepath):
    if not doi:
        return False, "crossref_no_doi"

    url = f"https://api.crossref.org/works/{doi}"
    try:
        text = curl_get_text(url, extra_headers={"User-Agent": f"mailto:{EMAIL}"})
        if not text:
            return False, "crossref_no_response"
        data = json.loads(text).get("message", {})
    except Exception:
        return False, "crossref_parse_error"

    candidates = []
    for link in data.get("link", []):
        ct = link.get("content-type", "")
        link_url = link.get("URL", "")
        if "pdf" in ct.lower() and link_url:
            candidates.append(link_url)

    resource = data.get("resource", {}).get("primary", {}).get("URL")
    if resource:
        candidates.append(resource)

    for cand in candidates:
        if curl_download(cand, filepath):
            return True, "crossref"
        time.sleep(0.3)

    return False, "crossref_no_working_pdf"


# ---------------------------------------------------------------------------
# Cascade orchestrator
# ---------------------------------------------------------------------------

def download_one(pmcid, pmid, doi, filepath):
    ok, src = try_pmc_pow(pmcid, filepath)
    if ok:
        return True, src
    time.sleep(DELAY)

    ok, src = try_openalex(doi, pmid, pmcid, filepath)
    if ok:
        return True, src
    time.sleep(DELAY)

    if doi:
        ok, src = try_crossref(doi, filepath)
        if ok:
            return True, src

    return False, "all_failed"


# ---------------------------------------------------------------------------
# Summary (reconstructed from the journals + the PDFs on disk -> always correct)
# ---------------------------------------------------------------------------

def write_summary(xml_meta):
    """Rebuild SUMMARY_FILE from the append-only journals and the files on disk.

    Idempotent and resume-safe: callable at any point (incl. after an
    interruption) and always reflects the true current state.
    """
    done = load_log(DONE_LOG)
    failed = load_log(FAILED_LOG)
    on_disk = get_already_downloaded(OUTPUT_DIR)

    rows = []
    for pmcid in sorted(on_disk):
        if pmcid in done:
            rows.append(done[pmcid][:5])
        else:                                # PDF present but journal entry lost
            m = xml_meta.get(pmcid, {})
            rows.append([pmcid, m.get("pmid", ""), m.get("doi", ""), "downloaded", "prior"])
    for pmcid in sorted(failed):
        if pmcid not in on_disk:             # a later success supersedes an old failure
            rows.append(failed[pmcid][:5])

    success = sum(1 for r in rows if r[3] == "downloaded")
    fail = sum(1 for r in rows if r[3] == "failed")
    with open(SUMMARY_FILE, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(LOG_HEADER)
        for row in sorted(rows, key=lambda x: x[3]):
            writer.writerow(row)
        writer.writerow(["Total", "", len(rows),
                         f"downloaded={success}", f"failed={fail}"])
    print(f"Summary saved to {SUMMARY_FILE}  (downloaded={success}, failed={fail})")


# ---------------------------------------------------------------------------
# HTML summary of the PDF export (same look as summaries/ncbi_xml.html)
# ---------------------------------------------------------------------------
CSS = """
 body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        margin: 2rem auto; max-width: 1100px; color: #222; line-height: 1.45; padding: 0 1rem; }
 h1 { margin-bottom: .25rem; } h2 { margin-top: 2.25rem; border-bottom: 1px solid #ddd; padding-bottom: .25rem; }
 .meta { color: #555; margin-bottom: 1rem; }
 .stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: .75rem; margin: 1rem 0 1.5rem; }
 .stat { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: .7rem 1rem; }
 .stat .v { font-size: 1.35rem; font-weight: 600; } .stat .k { color: #555; font-size: .82rem; }
 table { border-collapse: collapse; width: 100%; margin: .5rem 0 1rem; font-size: .9rem; }
 th, td { border: 1px solid #e1e4e8; padding: .3rem .55rem; text-align: left; vertical-align: middle; }
 th { background: #f6f8fa; }
 td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
 code { background: #f6f8fa; padding: 1px 4px; border-radius: 3px; font-size: .88em; }
 pre { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: .9rem 1rem; overflow-x: auto; font-size: .85rem; line-height: 1.4; white-space: pre-wrap; }
 .bar { display:inline-block; height:.72em; background:#3b7dd8; border-radius:2px; }
 .bar.g { background:#2da44e; } .bar.o { background:#bf8700; }
 .dim { color: #888; font-size: .85em; }
 .key  { background: #ddf4ff; border-left: 4px solid #0969da; padding: .6rem .9rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 ol.strategy > li { margin: .45rem 0; }
"""

STRATEGY_HTML = """
<ol class="strategy">
  <li><strong>Target set.</strong> The articles whose <code>high_impact_xmls/</code> XML has no
      <code>&lt;body&gt;</code> (front-matter / abstract-only records flagged by
      <code>xml_structure.py</code>) &mdash; their full text is unavailable as JATS, so the PDF is
      fetched for downstream GROBID conversion.</li>
  <li><strong>Cascade (replicates <code>bacs/doi_openalex.py</code>).</strong> Per article, in order,
      stop on first success: (1) <strong>PMC</strong> &mdash; scrape the article page for the PDF name,
      solve the proof-of-work challenge, download with the solved cookie via <code>curl</code>;
      (2) <strong>OpenAlex</strong> open-access PDF URLs; (3) <strong>CrossRef</strong> full-text PDF
      links. PMCIDs are resolved to PMID+DOI via the NCBI ID Converter, seeded from each XML's front
      matter.</li>
  <li><strong>Auto-recovery.</strong> Finished PDFs are skipped; successes/failures are journalled live
      (<code>_downloaded.tsv</code> / <code>_failed.tsv</code>) and skipped on re-run; PDFs are written
      atomically via a <code>.part</code> temp; <code>TIME_BUDGET</code> allows a clean self-exit; and
      this summary plus the TSV are always rewritten in a <code>finally</code> block, reconstructed from
      the journals + the files on disk.</li>
  <li><strong>This page</strong> reports the export outcome: downloaded vs failed vs pending, the source
      that worked, and the bytes on disk.</li>
</ol>
"""


def _fmt(x):
    return format(x, ",")


def write_html_summary(xml_meta, total_targets):
    """Write summaries/pdf_exports_summary.html from the journals + the PDFs on disk.

    Recovery-safe and idempotent: reflects the true current state whenever called.
    """
    done = load_log(DONE_LOG)
    failed = load_log(FAILED_LOG)
    on_disk = get_already_downloaded(OUTPUT_DIR)

    downloaded = len(on_disk)
    failed_n = sum(1 for p in failed if p not in on_disk)
    attempted = downloaded + failed_n
    pending = max(0, total_targets - attempted)
    rate = (downloaded / attempted * 100) if attempted else 0

    # which source produced each downloaded PDF
    src_counts = {}
    for pmcid in on_disk:
        row = done.get(pmcid)
        src = row[4] if (row and len(row) > 4 and row[4]) else "prior/unknown"
        src_counts[src] = src_counts.get(src, 0) + 1

    # bytes on disk
    total_bytes = 0
    for pmcid in on_disk:
        try:
            total_bytes += os.path.getsize(os.path.join(OUTPUT_DIR, "%s.pdf" % pmcid))
        except OSError:
            pass
    mb = total_bytes / (1024 * 1024)
    mean_kb = (total_bytes / 1024 / downloaded) if downloaded else 0

    import html as _html
    H = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<title>PDF export &mdash; summary</title>",
         "<style>" + CSS + "</style></head><body>"]
    H.append("<h1>PDF export &mdash; summary</h1>")
    H.append("<p class='meta'>Generated by <code>ncbi_pdf.py</code> (replica of "
             "<code>bacs/doi_openalex.py</code>) &middot; source: the <strong>%s</strong> "
             "<code>high_impact_xmls/</code> records with no <code>&lt;body&gt;</code> &middot; "
             "output: <code>ncbi_pdfs_grobid/PMC*.pdf</code></p>" % _fmt(total_targets))

    cards = [
        (_fmt(total_targets), "no-&lt;body&gt; targets"),
        (_fmt(downloaded), "PDFs downloaded"),
        (_fmt(failed_n), "failed (all sources)"),
        (_fmt(pending), "pending (not yet tried)"),
        ("%.1f%%" % rate, "success rate (of attempted)"),
        ("%.1f MB" % mb, "downloaded on disk"),
    ]
    H.append("<div class='stat-grid'>")
    for v, k in cards:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")

    H.append("<div class='key'>Of <strong>%s</strong> articles without a <code>&lt;body&gt;</code>, "
             "<strong>%s</strong> PDFs were downloaded and <strong>%s</strong> failed every source "
             "(<strong>%s</strong> not yet attempted). Mean PDF size <strong>%.0f KB</strong>. "
             "Row-level status: <code>ncbi_pdfs_grobid/ncbi_pdf_download_summary.tsv</code>.</div>"
             % (_fmt(total_targets), _fmt(downloaded), _fmt(failed_n), _fmt(pending), mean_kb))

    H.append("<h2>1. Strategy</h2>" + STRATEGY_HTML)

    H.append("<h2>2. Export outcome</h2>")
    H.append("<table><thead><tr><th>outcome</th><th class='num'>articles</th><th>dist.</th>"
             "<th class='num'>%</th></tr></thead><tbody>")
    omax = max(downloaded, failed_n, pending) or 1
    for name, c, cls in (("downloaded", downloaded, "bar g"),
                         ("failed (all sources)", failed_n, "bar o"),
                         ("pending (not yet tried)", pending, "bar")):
        w = max(2, round(c / omax * 260))
        pc = (c / total_targets * 100) if total_targets else 0
        H.append("<tr><td>%s</td><td class='num'>%s</td>"
                 "<td><span class='%s' style='width:%dpx'></span></td>"
                 "<td class='num dim'>%.1f%%</td></tr>" % (name, _fmt(c), cls, w, pc))
    H.append("</tbody></table>")

    H.append("<h2>3. Successful source</h2>")
    if src_counts:
        H.append("<p class='dim'>Which step in the cascade produced each downloaded PDF.</p>")
        H.append("<table><thead><tr><th>source</th><th class='num'>PDFs</th><th>dist.</th>"
                 "<th class='num'>%</th></tr></thead><tbody>")
        smax = max(src_counts.values()) or 1
        for src, c in sorted(src_counts.items(), key=lambda kv: kv[1], reverse=True):
            w = max(2, round(c / smax * 260))
            pc = (c / downloaded * 100) if downloaded else 0
            H.append("<tr><td>%s</td><td class='num'>%s</td>"
                     "<td><span class='bar g' style='width:%dpx'></span></td>"
                     "<td class='num dim'>%.1f%%</td></tr>"
                     % (_html.escape(src), _fmt(c), w, pc))
        H.append("</tbody></table>")
    else:
        H.append("<p class='dim'>No PDFs downloaded yet &mdash; run <code>ncbi_pdf.py</code>.</p>")

    H.append("""<h2>4. Caveats</h2><ul>
  <li><strong>Counts reflect the journals + files on disk at generation time.</strong> The page is
      rewritten on every run (incl. after an interruption or <code>TIME_BUDGET</code> exit), so it
      stays current as the export resumes.</li>
  <li><strong>Failures are expected.</strong> Many no-<code>&lt;body&gt;</code> records are
      publisher-restricted and not openly available, so no PDF is obtainable through PMC, OpenAlex or
      CrossRef; those are logged and skipped on re-run.</li>
</ul>""")
    H.append("</body></html>")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(SUMMARY_HTML, "w", encoding="utf-8") as f:
        f.write("\n".join(H))
    print(f"HTML summary saved to {SUMMARY_HTML}")


def write_all_summaries(xml_meta, total_targets):
    """Write both the TSV and the HTML summary (called live + in the finally block)."""
    write_summary(xml_meta)
    write_html_summary(xml_meta, total_targets)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _empty_content_dir(path):
    """Empty an output directory of generated content before a fresh run so each
    run starts clean. Python scripts are never deleted -- only generated content."""
    if os.path.isdir(path):
        for root, dirs, files in os.walk(path, topdown=False):
            for name in files:
                if not name.endswith(".py"):
                    os.remove(os.path.join(root, name))
            for name in dirs:
                try:
                    os.rmdir(os.path.join(root, name))   # remove only if now empty
                except OSError:
                    pass
    os.makedirs(path, exist_ok=True)


def main():
    os.makedirs(XML_DIR, exist_ok=True)      # create input dir if missing
    records, xml_meta = load_no_body_records(XML_DIR)
    total_targets = len(records)
    print(f"Loaded {total_targets} no-<body> PMC*.xml records from {XML_DIR}")
    _empty_content_dir(OUTPUT_DIR)           # fresh run: clear stale PDFs before writing
    sweep_partials(OUTPUT_DIR)               # auto-recovery: clear half-written downloads

    # Resume: skip PDFs already on disk and PMCIDs already journalled as failed.
    already = get_already_downloaded(OUTPUT_DIR)
    failed_prior = set() if RETRY_FAILED else set(load_log(FAILED_LOG))
    remaining = [(pmid, pmcid) for pmid, pmcid in records
                 if pmcid not in already and pmcid not in failed_prior]
    print(f"Already have {len(already)} PDFs; {len(failed_prior)} previously failed "
          f"(skipped{' - RETRY_FAILED set, not skipped' if RETRY_FAILED else ''}); "
          f"{len(remaining)} remaining to try")

    if not remaining:
        print("Nothing to download!")
        write_all_summaries(xml_meta, total_targets)
        return

    pmcids_remaining = [pmcid for _, pmcid in remaining]

    print(f"\nStep 1: Resolving {len(pmcids_remaining)} PMCIDs -> PMIDs + DOIs...")
    pmcid_to_pmid, pmcid_to_doi = resolve_ids(pmcids_remaining)
    # Seed any gaps with the PMID/DOI parsed directly from the XML front matter.
    for pmcid, m in xml_meta.items():
        if m.get("pmid") and pmcid not in pmcid_to_pmid:
            pmcid_to_pmid[pmcid] = m["pmid"]
        if m.get("doi") and pmcid not in pmcid_to_doi:
            pmcid_to_doi[pmcid] = m["doi"]
    print(f"Found {len(pmcid_to_pmid)} PMIDs, {len(pmcid_to_doi)} DOIs")

    print(f"\nStep 2: Downloading PDFs (PMC PoW -> OpenAlex -> CrossRef)...")
    success = failed = 0
    total = len(remaining)
    deadline = (time.time() + TIME_BUDGET) if TIME_BUDGET > 0 else 0
    interrupted = False

    try:
        for idx, (xml_pmid, pmcid) in enumerate(remaining, 1):
            if deadline and time.time() >= deadline:
                print(f"\n[recovery] TIME_BUDGET reached - stopping cleanly at "
                      f"{idx - 1}/{total}; re-run to resume")
                interrupted = True
                break

            # Prefer the PMID parsed from the XML; fall back to the ID-converter response.
            pmid = xml_pmid or pmcid_to_pmid.get(pmcid, "")
            doi = pmcid_to_doi.get(pmcid, "")
            filepath = os.path.join(OUTPUT_DIR, f"{pmcid}.pdf")

            id_str = " ".join(s for s in (
                pmcid,
                f"PMID:{pmid}" if pmid else "",
                f"doi:{doi}" if doi else "no DOI",
            ) if s)
            print(f"  [{idx}/{total}] {id_str}...", end=" ", flush=True)

            ok, source = download_one(pmcid, pmid, doi, filepath)

            if ok:
                size_kb = os.path.getsize(filepath) // 1024
                print(f"OK ({size_kb} KB) [{source}]")
                append_log(DONE_LOG, [pmcid, pmid, doi, "downloaded", source])
                success += 1
            else:
                print(f"FAIL [{source}]")
                append_log(FAILED_LOG, [pmcid, pmid, doi, "failed", source])
                failed += 1

            time.sleep(DELAY)
    except KeyboardInterrupt:
        interrupted = True
        print("\n[recovery] interrupted by user - progress journalled; re-run to resume")
    finally:
        # Auto-recovery: always (re)write the TSV + HTML summaries from journals + disk.
        write_all_summaries(xml_meta, total_targets)

    print(f"\n{'=' * 50}")
    print(f"This run -- downloaded: {success}  failed: {failed}"
          f"{'  (INCOMPLETE - re-run to resume)' if interrupted else ''}")


if __name__ == "__main__":
    main()
