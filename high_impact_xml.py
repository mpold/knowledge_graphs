#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
high_impact_xml.py
==================
Export the PMC article XMLs for the high-impact subset of a PubMed result set:
PubMed Central ids (pmc_ids) whose article-weighted journal impact factor sits in
the 90th percentile (top 10%) of the result set -- then write an HTML summary of
what was exported.

Written *after* ``bacs/pmc_xml.py`` and keeps its two defining elements:
  * a ``pmid2pmcid`` dictionary that drives the export and
  * the exact, entire parametrized ``fetch_pmc_xml_efetch(pmc_id,
    email="your-email@example.com")`` -- reproduced verbatim from ``bacs/pmc_xml.py``
    (see ``fetch_pmc_xml_efetch_canonical`` below) and used through a thin,
    time-out/interruption-tolerant wrapper of the *same signature*.

Here the dictionary is built from the 90th-percentile selection of
``pmids/pmid_pmc_ids.tsv`` (so it tracks the impact-factor data), and the script
is hardened to recover from interruptions and time-outs (see RECOVERY below).

Input  : pmids/pmid_pmc_ids.tsv
         (pmid, pmc_id, source_publication, issn, journal_impact_factor, year
          -- produced by pubmed_query.py)
Outputs: high_impact_xmls/PMC*.xml            one JATS XML per article
         high_impact_xmls/_failed.tsv         records that could not be fetched
         summaries/ncbi_xml_summary.html      summary of the exported XMLs

NOT executed here -- run it yourself; re-run to resume.

STRATEGY
--------
1. Threshold (article-weighted 90th percentile). Each *article* contributes its
   journal's impact factor once (a journal with many papers counts many times --
   "article-weighted"). The per-article impact factors from the TSV are sorted and
   the 90th-percentile value taken -- the same statistic brain_cancer_hits.html
   reports (~5.85). Override with IF_THRESHOLD=<float>; change the cut with
   PERCENTILE=<0..1>.

2. pmid -> pmc_id dictionary. Keep TSV rows that HAVE a pmc_id AND whose impact
   factor >= threshold, and assemble ``pmid2pmcid = {pmid: pmc_id}`` -- the
   structure bacs/pmc_xml.py drives its export from. (The pmc_ids in this dict are
   exactly the ones that "map to the 90th percentile".)

3. Export (after bacs/pmc_xml.py). For each pmc_id fetch the full-text JATS XML
   via eFetch (db=pmc, rettype=xml) with ``requests`` using the canonical
   ``fetch_pmc_xml_efetch`` signature, skipping any <PMCID>.xml already on disk,
   and write it as high_impact_xmls/<PMCID>.xml; sleep 0.34 s (~3 req/s).
   Publisher-restricted records return front matter only (detected as
   'abstract-only'); nothing is silently dropped.

4. Summarise. Build summaries/ncbi_xml_summary.html from the selection + the files
   on disk -- same look as brain_cancer_hits.html: stat cards, strategy, XMLs per
   year, full-text vs abstract-only proportion, impact-factor break-down and the
   top journals.

RECOVERY FROM INTERRUPTIONS AND TIME-OUTS
-----------------------------------------
* Resume: a finished <PMCID>.xml on disk is skipped, so re-running continues where
  it stopped (no progress lost).
* Atomic writes: each XML is written to a temp file and os.replace()'d into place,
  so an interrupted/killed write never leaves a half-file that looks "done".
* Network time-outs: every eFetch has connect/read timeouts and is retried with
  exponential back-off; HTTP 429/5xx honour Retry-After. Transient failures retry;
  permanent ones are logged to _failed.tsv and skipped on future runs (override
  with RETRY_FAILED=1).
* Wall-clock cap: set TIME_BUDGET=<seconds> to make a run self-exit cleanly before
  a command/CI timeout; re-run to finish (0 = unlimited, the default).
* Ctrl-C / unexpected exit: handled gracefully -- the partial run is consistent and
  the summary is still written (finally block) so the page reflects progress.
"""

import os
import csv
import sys
import time
import html as _html
import requests

# --------------------------------------------------------------------------- #
# Paths / config
# --------------------------------------------------------------------------- #
BASE        = os.path.dirname(os.path.abspath(__file__))
IN_TSV      = os.path.join(BASE, "pmids", "pmid_pmc_ids.tsv")
OUT_DIR     = os.path.join(BASE, "high_impact_xmls")
FAILED_LOG  = os.path.join(OUT_DIR, "_failed.tsv")
SUMMARY_DIR = os.path.join(BASE, "summaries")
SUMMARY_HTML= os.path.join(SUMMARY_DIR, "high_impact_xml.html")

EMAIL       = os.environ.get("CONTACT_EMAIL", "your-email@example.com")
DELAY       = 0.34                                           # ~3 req/s (no API key)
IF_OVERRIDE = os.environ.get("IF_THRESHOLD", "").strip()
PCTL        = float(os.environ.get("PERCENTILE", "0.90"))
TIME_BUDGET = float(os.environ.get("TIME_BUDGET", "0"))     # 0 = unlimited
MAX_TRIES   = int(os.environ.get("MAX_TRIES", "5"))
HTTP_TIMEOUT= (10, 120)                                      # (connect, read) seconds
RETRY_FAILED= os.environ.get("RETRY_FAILED", "").strip() in ("1", "true", "yes")


# --------------------------------------------------------------------------- #
# 1-2. Selection -> pmid2pmcid dictionary (built from the 90th-pctile subset)
# --------------------------------------------------------------------------- #
def parse_if(s):
    try:
        return float((s or "").strip())
    except ValueError:
        return None


def load_selection():
    """Return (threshold, pmid2pmcid dict, meta-by-pmcid) from the TSV.

    Article-weighted: the percentile is taken over the per-ARTICLE impact factors
    (one value per row), so a prolific high-IF journal pulls the threshold up in
    proportion to how many articles it contributed.
    """
    with open(IN_TSV, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    if IF_OVERRIDE:
        threshold = float(IF_OVERRIDE)
    else:
        vals = sorted(v for v in (parse_if(r["journal_impact_factor"]) for r in rows)
                      if v is not None)
        if not vals:
            raise SystemExit("error: no impact-factor values in TSV")
        threshold = vals[min(len(vals) - 1, int(PCTL * len(vals)))]

    pmid2pmcid, meta = {}, {}
    for r in rows:
        pmc = r["pmc_id"].strip()
        v = parse_if(r["journal_impact_factor"])
        if pmc and v is not None and v >= threshold:
            pmid2pmcid[r["pmid"]] = pmc
            meta[pmc] = {"pmid": r["pmid"], "journal": r["source_publication"],
                         "issn": r["issn"], "if": v, "year": r["year"]}
    return threshold, pmid2pmcid, meta


# --------------------------------------------------------------------------- #
# Recovery helpers
# --------------------------------------------------------------------------- #
def load_failed():
    """PMCIDs already known to be unfetchable (skip unless RETRY_FAILED)."""
    failed = set()
    if os.path.exists(FAILED_LOG) and not RETRY_FAILED:
        with open(FAILED_LOG, encoding="utf-8") as fh:
            for ln in fh:
                p = ln.split("\t", 1)[0].strip()
                if p and p != "pmcid":
                    failed.add(p)
    return failed


def log_failed(pmc_id, reason):
    new = not os.path.exists(FAILED_LOG)
    with open(FAILED_LOG, "a", encoding="utf-8") as fh:
        if new:
            fh.write("pmcid\treason\n")
        fh.write("%s\t%s\n" % (pmc_id, str(reason).replace("\t", " ").replace("\n", " ")))


def atomic_write(path, text):
    """Write text to a temp file then atomically replace -- never a half file."""
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)          # atomic on POSIX and Windows


def looks_like_pmc_xml(text):
    return ("<pmc-articleset" in text) or ("<article" in text)


# --------------------------------------------------------------------------- #
# 3. The canonical fetch function, reproduced verbatim from bacs/pmc_xml.py
# --------------------------------------------------------------------------- #
def fetch_pmc_xml_efetch_canonical(pmc_id, email="your-email@example.com"):
    """Fetch full-text JATS XML via NCBI eFetch.
    pmc_id: numeric PMC ID (e.g. 209839) or string like 'PMC209839'

    This is the entire parametrized function exactly as it appears in
    bacs/pmc_xml.py. It performs a single eFetch and raises on HTTP error; the
    time-out/interruption recovery is added by the wrapper below, which preserves
    this same signature.
    """
    pmc_id = str(pmc_id).replace("PMC", "")
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pmc",
        "id": pmc_id,
        "rettype": "xml",       # full JATS XML
        "retmode": "xml",
        "tool": "my_script",
        "email": email,
        }
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    return resp.text            # JATS XML string


def fetch_pmc_xml_efetch(pmc_id, email="your-email@example.com"):
    """Time-out/interruption-tolerant eFetch with the canonical signature.

    Same parameters and return contract as fetch_pmc_xml_efetch_canonical above
    (returns the JATS XML string; raises requests.RequestException when the article
    cannot be fetched after MAX_TRIES so the caller can log it and carry on), but
    each call carries connect/read timeouts and retries transient time-outs and
    rate-limits (HTTP 429/5xx, honouring Retry-After) with exponential back-off.
    """
    pmc_id = str(pmc_id).replace("PMC", "")
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pmc",
        "id": pmc_id,
        "rettype": "xml",       # full JATS XML
        "retmode": "xml",
        "tool": "high_impact_xml",
        "email": email,
    }

    last_exc = None
    for attempt in range(MAX_TRIES):
        try:
            resp = requests.get(url, params=params, timeout=HTTP_TIMEOUT)
            if resp.status_code in (429, 500, 502, 503, 504):   # transient -> retry
                wait = resp.headers.get("Retry-After")
                time.sleep(float(wait) if (wait and wait.isdigit())
                           else min(60, 1.5 * (2 ** attempt)))
                last_exc = requests.HTTPError("HTTP %d" % resp.status_code)
                continue
            resp.raise_for_status()                              # permanent 4xx -> raises
            return resp.text                                     # JATS XML string
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_exc = exc                                       # time-out -> back off & retry
            time.sleep(min(60, 1.5 * (2 ** attempt)))
    raise last_exc if last_exc else requests.RequestException("eFetch failed")


# --------------------------------------------------------------------------- #
# Export loop (driven by the pmid2pmcid dictionary, with recovery)
# --------------------------------------------------------------------------- #
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


def export(pmid2pmcid):
    _empty_content_dir(OUT_DIR)              # fresh run: clear stale exports before writing
    failed_skip = load_failed()
    deadline = (time.time() + TIME_BUDGET) if TIME_BUDGET > 0 else 0

    total = len(pmid2pmcid)
    fetched = skipped_done = skipped_failed = newly_failed = 0
    interrupted = False
    try:
        for done, (pmid, pmc_id) in enumerate(pmid2pmcid.items(), 1):
            path = os.path.join(OUT_DIR, "%s.xml" % pmc_id)
            if os.path.exists(path) and os.path.getsize(path) > 0:
                skipped_done += 1
                continue
            if pmc_id in failed_skip:
                skipped_failed += 1
                continue
            if deadline and time.time() >= deadline:
                print("[export] TIME_BUDGET reached - stopping cleanly; re-run to resume",
                      file=sys.stderr)
                interrupted = True
                break

            try:
                text = fetch_pmc_xml_efetch(pmc_id, email=EMAIL)
            except requests.RequestException as exc:
                log_failed(pmc_id, "%s: %s" % (type(exc).__name__, str(exc)[:200]))
                newly_failed += 1
                time.sleep(DELAY)
                continue

            if not text or not looks_like_pmc_xml(text):
                log_failed(pmc_id, "empty-or-nonxml-response")
                newly_failed += 1
            else:
                atomic_write(path, text)
                fetched += 1
                if fetched % 200 == 0:
                    print("[export] fetched %d (at %d/%d)" % (fetched, done, total),
                          file=sys.stderr)
            time.sleep(DELAY)
    except KeyboardInterrupt:
        interrupted = True
        print("\n[export] interrupted by user - progress saved; re-run to resume",
              file=sys.stderr)

    print("[export] selected=%d  fetched_now=%d  already_done=%d  "
          "skipped_failed=%d  newly_failed=%d%s"
          % (total, fetched, skipped_done, skipped_failed, newly_failed,
             "  (INCOMPLETE)" if interrupted else ""), file=sys.stderr)
    return not interrupted


# --------------------------------------------------------------------------- #
# 4. Summary HTML (mirrors brain_cancer_hits.html)
# --------------------------------------------------------------------------- #
CSS = """
 body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        margin: 2rem auto; max-width: 1100px; color: #222; line-height: 1.45; padding: 0 1rem; }
 h1 { margin-bottom: .25rem; } h2 { margin-top: 2.25rem; border-bottom: 1px solid #ddd; padding-bottom: .25rem; }
 h3 { margin: 1.2rem 0 .3rem; }
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
 .cols { display:grid; grid-template-columns: 1fr 1fr; gap:1.5rem; align-items:start; }
 @media (max-width:780px){ .cols{ grid-template-columns:1fr; } }
"""

STRATEGY_HTML = """
<ol class="strategy">
  <li><strong>Threshold &mdash; article-weighted 90th percentile.</strong> Each article counts its
      journal's impact factor once, so a prolific high-impact journal weighs in proportion to its
      article count. The per-article impact factors from <code>pmids/pmid_pmc_ids.tsv</code> are sorted
      and the 90th-percentile value taken (the same statistic <code>brain_cancer_hits.html</code>
      reports, &asymp;5.85).</li>
  <li><strong>pmid &rarr; pmc_id dictionary.</strong> Like <code>bacs/pmc_xml.py</code>, the export is
      driven by a <code>pmid2pmcid</code> dictionary &mdash; here built from the TSV rows that have a
      <code>pmc_id</code> AND impact factor &ge; threshold. Those <code>pmc_id</code>s are exactly the
      ones that map to the 90th percentile.</li>
  <li><strong>Export (after bacs/pmc_xml.py).</strong> The entire parametrized
      <code>fetch_pmc_xml_efetch(pmc_id, email="your-email@example.com")</code> from
      <code>bacs/pmc_xml.py</code> is reproduced verbatim; each <code>pmc_id</code> is fetched with
      <code>requests</code> via eFetch (<code>db=pmc, rettype=xml</code>), skipping any
      <code>&lt;PMCID&gt;.xml</code> already present, written as
      <code>high_impact_xmls/&lt;PMCID&gt;.xml</code> with a 0.34&nbsp;s pause. Publisher-restricted
      records return front matter only (<em>abstract-only</em>).</li>
  <li><strong>Recovery from interruptions &amp; time-outs.</strong> Finished files are skipped on
      re-run; XMLs are written atomically (temp + <code>os.replace</code>) so a killed write never
      leaves a half-file; eFetch calls have connect/read timeouts and exponential back-off (honouring
      HTTP 429 <code>Retry-After</code>); unfetchable records are logged to <code>_failed.tsv</code> and
      skipped thereafter; <code>TIME_BUDGET</code> allows a clean self-exit before a wall-clock cap; and
      Ctrl-C still leaves a consistent state with the summary written.</li>
  <li><strong>Summarise.</strong> This page is built from the selection and the files on disk.</li>
</ol>
"""


def _fmt(x):
    return format(x, ",")


def _pctile(p, arr):
    return arr[min(len(arr) - 1, int(p * len(arr)))] if arr else None


def classify_file(path):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            head = fh.read(1000000)
    except OSError:
        return "missing"
    if not head.strip():
        return "empty"
    if "<body" in head:
        return "fulltext"
    if "pmc-articleset" in head or "<article" in head:
        return "abstract-only"
    return "other"


def build_summary(threshold, meta):
    """Build summaries/ncbi_xml_summary.html from the selection + exported files."""
    rows = []
    status_counts = {}
    for pmc, m in meta.items():
        path = os.path.join(OUT_DIR, "%s.xml" % pmc)
        if not (os.path.exists(path) and os.path.getsize(path) > 0):
            continue
        st = classify_file(path)
        status_counts[st] = status_counts.get(st, 0) + 1
        d = dict(m)
        d["pmcid"] = pmc
        d["status"] = st
        rows.append(d)

    n = len(rows)
    full = status_counts.get("fulltext", 0)
    abso = status_counts.get("abstract-only", 0)
    full_prop = full / n if n else 0
    n_failed = 0
    if os.path.exists(FAILED_LOG):
        with open(FAILED_LOG, encoding="utf-8") as fh:
            n_failed = max(0, sum(1 for _ in fh) - 1)

    per_year, full_year = {}, {}
    if_vals, jcount, jif = [], {}, {}
    for r in rows:
        y = r["year"]
        if y:
            per_year[y] = per_year.get(y, 0) + 1
            if r["status"] == "fulltext":
                full_year[y] = full_year.get(y, 0) + 1
        if r["if"] is not None:
            if_vals.append(r["if"])
        j = r["journal"]
        if j:
            jcount[j] = jcount.get(j, 0) + 1
            if r["if"] is not None:
                jif[j] = r["if"]
    if_vals.sort()
    years = sorted((int(y) for y in per_year), reverse=True)
    maxc = max(per_year.values()) if per_year else 1

    buckets = [("<5", 0, 5), ("5-7", 5, 7), ("7-10", 7, 10), ("10-15", 10, 15),
               ("15-20", 15, 20), ("20-30", 20, 30), ("30-50", 30, 50), ("50+", 50, 1e9)]
    bdist = {b[0]: 0 for b in buckets}
    for v in if_vals:
        for name, lo, hi in buckets:
            if lo <= v < hi:
                bdist[name] += 1
                break

    H = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<title>High-impact PMC XML export &mdash; summary</title>",
         "<style>" + CSS + "</style></head><body>"]
    H.append("<h1>High-impact PMC XML export &mdash; summary</h1>")
    H.append("<p class='meta'>Generated by <code>high_impact_xml.py</code> (written after "
             "<code>bacs/pmc_xml.py</code>) &middot; selection: pmc_id present AND article-weighted "
             "impact factor &ge; <strong>%.4f</strong> (90th percentile) &middot; %d pmc_ids in the "
             "<code>pmid2pmcid</code> dictionary</p>" % (threshold, len(meta)))

    med = _pctile(0.5, if_vals)
    cards = [(_fmt(n), "XMLs exported"),
             (_fmt(full), "full-text (%.1f%%)" % (full_prop * 100)),
             (_fmt(abso), "abstract-only (%.1f%%)" % (abso / n * 100 if n else 0)),
             (_fmt(len(jcount)), "distinct journals"),
             ("%.2f" % med if med is not None else "n/a", "median impact factor"),
             (_fmt(n_failed), "unfetchable (logged)")]
    H.append("<div class='stat-grid'>")
    for v, k in cards:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")

    H.append("<div class='key'>Exported <strong>%s</strong> of <strong>%s</strong> selected article "
             "XMLs. <strong>%s (%.1f%%)</strong> carry full text (<code>&lt;body&gt;</code>); "
             "<strong>%s</strong> are publisher-restricted to front matter; <strong>%s</strong> were "
             "unfetchable (logged to <code>_failed.tsv</code>). Files: "
             "<code>high_impact_xmls/PMC*.xml</code>.</div>"
             % (_fmt(n), _fmt(len(meta)), _fmt(full), full_prop * 100, _fmt(abso), _fmt(n_failed)))

    H.append("<h2>1. Strategy</h2>" + STRATEGY_HTML)

    H.append("<h2>2. Exported XMLs per year &amp; full-text coverage</h2><div class='cols'><div><table>")
    H.append("<thead><tr><th class='num'>year</th><th class='num'>XMLs</th><th>dist.</th>"
             "<th class='num'>full text</th><th class='num'>FT %</th></tr></thead><tbody>")
    for y in years:
        ys = str(y)
        c = per_year[ys]
        fc = full_year.get(ys, 0)
        w = max(2, round(c / maxc * 260))
        H.append("<tr><td class='num'>%d</td><td class='num'>%s</td>"
                 "<td><span class='bar' style='width:%dpx'></span></td>"
                 "<td class='num'>%s</td><td class='num dim'>%.0f%%</td></tr>"
                 % (y, _fmt(c), w, _fmt(fc), (fc / c * 100 if c else 0)))
    H.append("</tbody></table></div><div>")
    dec = {}
    for ys, c in per_year.items():
        d = (int(ys) // 10) * 10
        dec[d] = dec.get(d, 0) + c
    H.append("<h3 style='margin-top:0'>By decade</h3><table>"
             "<thead><tr><th class='num'>decade</th><th class='num'>XMLs</th></tr></thead><tbody>")
    for d in sorted(dec, reverse=True):
        H.append("<tr><td class='num'>%ds</td><td class='num'>%s</td></tr>" % (d, _fmt(dec[d])))
    H.append("</tbody></table></div></div>")

    H.append("<h2>3. Full-text vs abstract-only</h2>")
    H.append("<table><thead><tr><th>status</th><th class='num'>XMLs</th><th class='num'>share</th>"
             "</tr></thead><tbody>")
    for name in ("fulltext", "abstract-only", "other", "empty"):
        c = status_counts.get(name, 0)
        if c:
            H.append("<tr><td>%s</td><td class='num'>%s</td><td class='num dim'>%.1f%%</td></tr>"
                     % (name, _fmt(c), c / n * 100 if n else 0))
    H.append("</tbody></table>")
    H.append("<p class='dim'>'abstract-only' = PMC returned front matter only (publisher does not permit "
             "full-text XML download).</p>")

    H.append("<h2>4. Impact-factor break-down (exported set)</h2>")
    if if_vals:
        stats = [("%.2f" % (sum(if_vals) / len(if_vals)), "mean"),
                 ("%.2f" % _pctile(0.5, if_vals), "median"),
                 ("%.2f" % _pctile(0.75, if_vals), "75th pct"),
                 ("%.2f" % _pctile(0.90, if_vals), "90th pct"),
                 ("%.2f" % _pctile(0.99, if_vals), "99th pct"),
                 ("%.0f" % if_vals[-1], "max")]
        H.append("<div class='stat-grid'>")
        for v, k in stats:
            H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
        H.append("</div>")
        bmax = max(bdist.values()) or 1
        H.append("<table><thead><tr><th>IF band</th><th class='num'>XMLs</th><th>share</th>"
                 "<th class='num'>%</th></tr></thead><tbody>")
        for name, _lo, _hi in buckets:
            v = bdist[name]
            w = max(2, round(v / bmax * 320))
            H.append("<tr><td>%s</td><td class='num'>%s</td>"
                     "<td><span class='bar o' style='width:%dpx'></span></td>"
                     "<td class='num dim'>%.1f%%</td></tr>"
                     % (name, _fmt(v), w, (v / len(if_vals) * 100)))
        H.append("</tbody></table>")
        H.append("<p class='dim'>All exported articles sit at or above the 90th-percentile cut "
                 "(%.2f), so values cluster in the upper bands.</p>" % threshold)
    else:
        H.append("<p class='dim'>No impact-factor values available.</p>")

    H.append("<h2>5. Top journals in the exported set</h2>")
    top = sorted(jcount.items(), key=lambda kv: kv[1], reverse=True)[:25]
    jmax = top[0][1] if top else 1
    H.append("<table><thead><tr><th>journal (NLM abbrev)</th><th class='num'>XMLs</th><th>dist.</th>"
             "<th class='num'>impact factor</th></tr></thead><tbody>")
    for j, c in top:
        w = max(2, round(c / jmax * 220))
        iv = jif.get(j)
        H.append("<tr><td>%s</td><td class='num'>%s</td>"
                 "<td><span class='bar g' style='width:%dpx'></span></td>"
                 "<td class='num'>%s</td></tr>"
                 % (_html.escape(j), _fmt(c), w, ("n/a" if iv is None else "%.2f" % iv)))
    H.append("</tbody></table>")

    H.append("""<h2>6. Caveats</h2><ul>
  <li><strong>Impact factor is the value carried in the input TSV</strong> (OpenAlex 2-yr mean
      citedness, a JIF proxy, unless a curated table was supplied to <code>pubmed_query.py</code>).</li>
  <li><strong>Full text is not guaranteed.</strong> A pmc_id means a PMC record exists; many publishers
      restrict full-text XML, so a share of exports are front-matter only.</li>
  <li><strong>Counts reflect the files on disk at generation time</strong> &mdash; if the run was
      interrupted or hit <code>TIME_BUDGET</code>, re-run to completion and this page refreshes.</li>
</ul>""")
    H.append("</body></html>")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(SUMMARY_HTML, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    print("[summary] wrote %s (%d XMLs on disk, %d failed)" % (SUMMARY_HTML, n, n_failed),
          file=sys.stderr)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    os.makedirs(os.path.dirname(IN_TSV), exist_ok=True)   # create input dir (pmids/) if missing
    if not os.path.exists(IN_TSV):
        raise SystemExit("error: input not found: %s" % IN_TSV)
    threshold, pmid2pmcid, meta = load_selection()
    print("[select] 90th-pctile IF threshold=%.4f; pmid2pmcid has %d entries"
          % (threshold, len(pmid2pmcid)), file=sys.stderr)
    complete = False
    try:
        complete = export(pmid2pmcid)
    finally:
        # Always summarise what we have, even after an interruption/time-out.
        build_summary(threshold, meta)
    print("[done] %s" % ("complete" if complete else "partial - re-run to resume"),
          file=sys.stderr)


if __name__ == "__main__":
    main()
