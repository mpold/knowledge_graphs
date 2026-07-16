#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pubmed_query.py
===============
Retrieve ALL PMIDs for an arbitrary PubMed query (read from STDIN), annotate each
with its PubMed Central id, source journal, journal ISSN and a journal impact
factor, write everything to a TSV, and emit an HTML summary.

USAGE
-----
    python pubmed_query.py
        -> prompts:  enter command line argument:
           type the query and press Enter (reads one line from STDIN).
    echo "brain cancer" | python pubmed_query.py   # piped STDIN also works

    The query is read strictly from STDIN; there is no command-line argument
    fallback, so the query recorded in the report is always exactly what was
    entered on STDIN.

OUTPUTS  (paths are resolved relative to this script's directory)
-------
    pmids/pmid_pmc_ids.tsv          pmid, pmc_id, source_publication, issn,
                                     journal_impact_factor, year
    summaries/query_hits.html  strategy + per-year counts,
                                     PMC-id proportion, impact-factor breakdown, ...

Caches under pmids/ (_pmids.txt, _annot.tsv, _journal_if.tsv) make long runs
resumable: re-running continues where it stopped.

STRATEGY (also embedded into the HTML report)
---------------------------------------------
1. Beat the 10,000-record export cap. PubMed's esearch returns at most ~10,000
   UIDs per query and cannot page beyond that index even with the history server.
   We partition the query by publication year (datetype=pdat) and recursively
   subdivide any window still above the cap (year -> month -> day) until every
   leaf window is retrievable in one page; the de-duplicated union is complete.

2. Annotate PMID -> PMCID + journal + ISSN + year in one pass. esummary (POST,
   batches of 200) returns a docsum whose articleids carry the PMC id, whose
   `source`/`issn`/`essn` give the journal and ISSNs, and whose `sortpubdate`
   gives the publication year. No PMC id present => not in PubMed Central.

3. Journal impact factor. The official Journal Impact Factor (Clarivate JCR) is
   proprietary with no free API, so by default we use OpenAlex
   `summary_stats.2yr_mean_citedness` -- mean citations in the last two years to
   a journal's papers, the same construction as JIF -- as an open, reproducible
   proxy, looked up per ISSN (batched 50/request). If a curated table is supplied
   (env JIF_TABLE=path to a "issn<TAB>jif" file), its values take precedence, so
   true JCR figures can be plugged in without code changes.

4. Polite & resumable. Requests are throttled (~3/s NCBI, ~7/s OpenAlex) with
   retries, back-off and tolerant JSON parsing; every phase checkpoints to a TSV
   after each batch. Set TIME_BUDGET=<seconds> to make each run self-exit before
   a wall-clock cap (0 = run to completion, the default).
"""

import os
import sys
import json
import time
import http.client
import urllib.request
import urllib.parse
import urllib.error

# --------------------------------------------------------------------------- #
# Configuration / paths (project-relative)
# --------------------------------------------------------------------------- #
BASE        = os.path.dirname(os.path.abspath(__file__))
PMID_DIR    = os.path.join(BASE, "pmids")
SUMMARY_DIR = os.path.join(BASE, "summaries")
OUT_TSV     = os.path.join(PMID_DIR, "pmid_pmc_ids.tsv")
SUMMARY_HTML= os.path.join(SUMMARY_DIR, "pubmed_query.html")

CACHE_PMIDS = os.path.join(PMID_DIR, "_pmids.txt")
CACHE_ANNOT = os.path.join(PMID_DIR, "_annot.tsv")        # pmid pmcid year source issn essn
CACHE_JIF   = os.path.join(PMID_DIR, "_journal_if.tsv")   # issn if h_index name issn_l
CACHE_META  = os.path.join(PMID_DIR, "_query.json")       # signature of the query these caches belong to

ESEARCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
OA_SRC   = "https://api.openalex.org/sources"

API_KEY    = os.environ.get("NCBI_API_KEY", "").strip()
EMAIL      = os.environ.get("CONTACT_EMAIL", "your-email@example.com")
JIF_TABLE  = os.environ.get("JIF_TABLE", "").strip()      # optional curated issn->jif
NCBI_DELAY = 0.12 if API_KEY else 0.34
OA_DELAY   = 0.15
DATETYPE   = "pdat"
PAGE       = 9999            # < 10k esearch ceiling
ESUM_BATCH = 200            # smaller POST -> less chance of a truncated chunked reply
OA_BATCH   = 50
YEAR_HI    = int(os.environ.get("YEAR_HI", "2026"))
YEAR_LO    = int(os.environ.get("YEAR_LO", "1900"))
TIME_BUDGET= float(os.environ.get("TIME_BUDGET", "0"))    # 0 = unlimited


# --------------------------------------------------------------------------- #
# HTTP helper (retries, back-off, tolerant JSON)
# --------------------------------------------------------------------------- #
def http_json(url, data=None, tries=6):
    """GET (or POST if data) and parse JSON; tolerant of stray control chars.

    NCBI/OpenAlex occasionally drop a chunked response mid-stream
    (http.client.IncompleteRead): the partial body is truncated JSON and is
    unusable, so we discard it and retry with exponential back-off rather than
    surfacing the transport error to the caller.
    """
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=data,
                                         headers={"Accept-Encoding": "identity"})
            with urllib.request.urlopen(req, timeout=120) as r:
                raw = r.read().decode("utf-8", "replace")
            return json.loads(raw, strict=False)
        except (urllib.error.URLError, urllib.error.HTTPError,
                http.client.IncompleteRead, http.client.HTTPException,
                ConnectionError, TimeoutError, OSError,
                json.JSONDecodeError) as exc:
            if attempt == tries - 1:
                raise
            wait = min(30.0, 1.5 * (2 ** attempt))     # 1.5, 3, 6, 12, 24, 30...
            print("[http] retry %d/%d (%s): %s"
                  % (attempt + 1, tries, type(exc).__name__, url[:80]),
                  file=sys.stderr)
            time.sleep(wait)


def clean(s):
    return (s or "").replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


# --------------------------------------------------------------------------- #
# Phase 1 -- retrieve ALL PMIDs (date partitioning + recursive subdivision)
# --------------------------------------------------------------------------- #
def esearch(term, d1, d2, retstart=0):
    """One esearch over [d1, d2]; returns (total_count, idlist <= PAGE).

    NCBI intermittently returns an error payload instead of a result (rate-limit
    {"error": ...} with no esearchresult, or {"esearchresult": {"ERROR": ...}}),
    which has no "count" field. Treat that as transient: back off and retry, and
    only raise once it persists, rather than crashing on KeyError.
    """
    params = {"db": "pubmed", "term": term, "datetype": DATETYPE,
              "mindate": d1, "maxdate": d2, "retmode": "json",
              "retmax": PAGE, "retstart": retstart}
    if API_KEY:
        params["api_key"] = API_KEY
    url = ESEARCH + "?" + urllib.parse.urlencode(params)
    tries = 8
    for attempt in range(tries):
        try:
            doc = http_json(url)
        except Exception as exc:                       # network/JSON gave up -> retry here too
            doc = {"error": "http_json failed: %s" % exc}
        res = doc.get("esearchresult", {}) or {}
        time.sleep(NCBI_DELAY)
        # NCBI returns the count as a string; accept it whenever present and numeric.
        cnt = res.get("count")
        if cnt is not None and str(cnt).isdigit():
            return int(cnt), res.get("idlist", [])
        err = (res.get("ERROR") or res.get("errorlist") or
               doc.get("error") or doc)
        if attempt == tries - 1:
            raise RuntimeError("esearch returned no count for [%s..%s] after %d tries: %s"
                               % (d1, d2, tries, err))
        # Exponential back-off (capped) to ride out rate-limit / transient errors.
        wait = min(30.0, 2.0 * (2 ** attempt))
        print("[esearch] retry %d/%d for [%s..%s] in %.0fs: %s"
              % (attempt + 1, tries, d1, d2, wait, err), file=sys.stderr)
        time.sleep(wait)


def _days_in_month(year, month):
    if month == 12:
        nxt = (year + 1, 1)
    else:
        nxt = (year, month + 1)
    import datetime
    return (datetime.date(nxt[0], nxt[1], 1) - datetime.date(year, month, 1)).days


def fetch_day(term, dd):
    count, ids = esearch(term, dd, dd)
    while len(ids) < count and len(ids) < PAGE:        # single day > cap: implausible
        _, more = esearch(term, dd, dd, retstart=len(ids))
        if not more:
            break
        ids += more
    return ids


def fetch_month(term, year, month):
    last = _days_in_month(year, month)
    count, ids = esearch(term, "%d/%02d/01" % (year, month),
                         "%d/%02d/%02d" % (year, month, last))
    if count <= PAGE:
        return ids
    out = []
    for day in range(1, last + 1):
        out += fetch_day(term, "%d/%02d/%02d" % (year, month, day))
    return out


def fetch_year(term, year):
    """Return (true_year_count, [all pmids for that year]) subdividing as needed."""
    count, ids = esearch(term, "%d/01/01" % year, "%d/12/31" % year)
    if count <= PAGE:
        return count, ids
    out = []
    for month in range(1, 13):
        out += fetch_month(term, year, month)
    return count, out


def cache_signature(term):
    """Identity of the query a cache belongs to: the term plus the year window.

    The pmids/ caches are keyed on PMIDs, not on the query that produced them, so
    reusing them across a *different* query silently returns the wrong corpus. The
    signature makes that reuse conditional: same query + same YEAR_LO/HI => resume;
    anything different => the caches are stale and must be rebuilt.
    """
    return {"query": term, "year_lo": YEAR_LO, "year_hi": YEAR_HI}


def load_cache_signature():
    """Return the signature stored for the current caches, or None if absent/unreadable."""
    try:
        with open(CACHE_META, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def invalidate_caches(reason):
    """Drop the PMID/annotation/JIF caches so a changed query starts clean."""
    print("[pmids] cache invalidated (%s) -- refetching" % reason, file=sys.stderr)
    for path in (CACHE_PMIDS, CACHE_ANNOT, CACHE_JIF, CACHE_META):
        try:
            os.remove(path)
        except OSError:
            pass


def get_all_pmids(term):
    """Complete, de-duplicated PMID set for the query (cached to CACHE_PMIDS).

    The cache is reused only when its stored signature matches the current query
    and year window (see cache_signature); on a mismatch every pmids/ cache is
    wiped first so a new query can never inherit a previous run's PMIDs.
    """
    sig = cache_signature(term)
    if os.path.exists(CACHE_PMIDS):
        stored = load_cache_signature()
        if stored != sig:
            was = (stored or {}).get("query") if stored else None
            invalidate_caches("query changed from %r to %r" % (was, term))
        else:
            with open(CACHE_PMIDS, encoding="utf-8") as fh:
                cached = [ln.strip() for ln in fh if ln.strip()]
            if cached:
                print("[pmids] using cache: %d" % len(cached), file=sys.stderr)
                return cached
    total, _ = esearch(term, "%d/01/01" % YEAR_LO, "%d/12/31" % YEAR_HI)  # warm + sanity
    if total == 0:
        # PubMed matched nothing over the whole window. Return now instead of
        # sweeping all YEAR_HI..YEAR_LO years -- each per-year esearch would also
        # return 0, print nothing (the progress line is gated on a non-zero count),
        # and the run would look frozen for a minute-plus before yielding an empty
        # corpus. A frequent cause is a NOT phrase that PubMed tokenizes as a
        # substring of the main phrase -- e.g. "small cell lung cancer" is found
        # inside "non-small cell lung cancer" (the hyphen splits "non-small"), so
        # `... NOT "small cell lung cancer"[Title/Abstract]` excludes every record.
        print("[pmids] PubMed returned 0 records for this query over %d-%d -- nothing "
              "to fetch.\n        Check the query: a NOT phrase that is a token-substring "
              "of the main\n        phrase (e.g. \"small cell lung cancer\" within "
              "\"non-small cell lung cancer\")\n        excludes every record."
              % (YEAR_LO, YEAR_HI), file=sys.stderr)
        return []
    all_ids = set()
    for year in range(YEAR_HI, YEAR_LO - 1, -1):
        count, ids = fetch_year(term, year)
        if count:
            all_ids.update(ids)
            print("[pmids] %d: %6d  (unique so far %d)" % (year, count, len(all_ids)),
                  file=sys.stderr)
    ordered = sorted(all_ids, key=int)
    with open(CACHE_PMIDS, "w", encoding="utf-8") as fh:
        fh.write("\n".join(ordered))
    with open(CACHE_META, "w", encoding="utf-8") as fh:
        json.dump(sig, fh)
    print("[pmids] total unique %d" % len(ordered), file=sys.stderr)
    return ordered


# --------------------------------------------------------------------------- #
# Phase 2 -- annotate PMID -> pmcid, journal, issn, essn, year (esummary)
# --------------------------------------------------------------------------- #
def annotate(pmids, deadline):
    done = set()
    if os.path.exists(CACHE_ANNOT):
        with open(CACHE_ANNOT, encoding="utf-8") as fh:
            for ln in fh:
                p = ln.split("\t", 1)[0]
                if p:
                    done.add(p)
    todo = [p for p in pmids if p not in done]
    print("[annot] %d done, %d remaining" % (len(done), len(todo)), file=sys.stderr)
    out = open(CACHE_ANNOT, "a", encoding="utf-8")
    i = 0
    while i < len(todo):
        if deadline and time.time() >= deadline:
            out.close()
            return False
        batch = todo[i:i + ESUM_BATCH]
        i += ESUM_BATCH
        params = {"db": "pubmed", "retmode": "json", "id": ",".join(batch)}
        if API_KEY:
            params["api_key"] = API_KEY
        res = http_json(ESUMMARY, data=urllib.parse.urlencode(params).encode()).get("result", {})
        for uid in res.get("uids", []):
            a = res.get(uid, {})
            pmcid = ""
            for x in a.get("articleids", []):
                v = str(x.get("value", ""))
                if x.get("idtype") in ("pmc", "pmcid") and v.startswith("PMC"):
                    pmcid = v
                    break
            spd = a.get("sortpubdate") or a.get("pubdate") or a.get("epubdate") or ""
            year = spd[:4] if spd[:4].isdigit() else ""
            out.write("\t".join([uid, pmcid, year, clean(a.get("source")),
                                  clean(a.get("issn")), clean(a.get("essn"))]) + "\n")
        out.flush()
        print("[annot] +%d (~%d/%d this run)" % (len(batch), i, len(todo)), file=sys.stderr)
        time.sleep(NCBI_DELAY)
    out.close()
    return True


def load_annotations():
    rows = []
    with open(CACHE_ANNOT, encoding="utf-8") as fh:
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if len(c) >= 6:
                rows.append(dict(pmid=c[0], pmcid=c[1], year=c[2],
                                 source=c[3], issn=c[4], essn=c[5]))
    return rows


# --------------------------------------------------------------------------- #
# Phase 3 -- journal impact factor per ISSN (curated table or OpenAlex proxy)
# --------------------------------------------------------------------------- #
def load_curated_jif():
    table = {}
    if JIF_TABLE and os.path.exists(JIF_TABLE):
        with open(JIF_TABLE, encoding="utf-8") as fh:
            for ln in fh:
                c = ln.rstrip("\n").split("\t")
                if len(c) >= 2:
                    try:
                        table[c[0].strip()] = float(c[1])
                    except ValueError:
                        pass
        print("[jif] curated table: %d issns" % len(table), file=sys.stderr)
    return table


def impact_factors(issns, deadline):
    """OpenAlex 2yr_mean_citedness per ISSN; cached to CACHE_JIF (proxy for JIF)."""
    done = set()
    if os.path.exists(CACHE_JIF):
        with open(CACHE_JIF, encoding="utf-8") as fh:
            for ln in fh:
                p = ln.split("\t", 1)[0]
                if p:
                    done.add(p)
    todo = sorted(issns - done)
    print("[jif] %d issns done, %d remaining (of %d)" % (len(done), len(todo), len(issns)),
          file=sys.stderr)
    out = open(CACHE_JIF, "a", encoding="utf-8")
    i = 0
    while i < len(todo):
        if deadline and time.time() >= deadline:
            out.close()
            return False
        batch = todo[i:i + OA_BATCH]
        i += OA_BATCH
        url = (OA_SRC + "?filter=issn:" + "|".join(batch) +
               "&per_page=%d&mailto=%s" % (OA_BATCH, EMAIL))
        d = http_json(url)
        found = {}
        for s in d.get("results", []):
            ss = s.get("summary_stats", {}) or {}
            row = (ss.get("2yr_mean_citedness"), ss.get("h_index"),
                   clean(s.get("display_name")), s.get("issn_l") or "")
            for isn in (s.get("issn") or []):
                found[isn] = row
            if s.get("issn_l"):
                found[s["issn_l"]] = row
        for isn in batch:
            r = found.get(isn)
            if r:
                if2, h, name, isnl = r
                out.write("\t".join([isn, "" if if2 is None else "%.4f" % if2,
                                     "" if h is None else str(h), name, isnl]) + "\n")
            else:
                out.write("\t".join([isn, "", "", "", ""]) + "\n")
        out.flush()
        print("[jif] +%d (~%d/%d this run)" % (len(batch), i, len(todo)), file=sys.stderr)
        time.sleep(OA_DELAY)
    out.close()
    return True


def load_jif_map():
    jif = dict(load_curated_jif())          # curated first (takes precedence)
    if os.path.exists(CACHE_JIF):
        with open(CACHE_JIF, encoding="utf-8") as fh:
            for ln in fh:
                c = ln.rstrip("\n").split("\t")
                if len(c) >= 2 and c[1] and c[0] not in jif:
                    jif[c[0]] = float(c[1])
    return jif


# --------------------------------------------------------------------------- #
# Output -- TSV
# --------------------------------------------------------------------------- #
def write_tsv(rows, jif):
    with open(OUT_TSV, "w", encoding="utf-8") as fh:
        fh.write("pmid\tpmc_id\tsource_publication\tissn\tjournal_impact_factor\tyear\n")
        for r in rows:
            issn = r["essn"] or r["issn"]                       # one canonical ISSN
            v = jif.get(r["essn"]) if r["essn"] in jif else jif.get(r["issn"])
            fh.write("\t".join([
                r["pmid"], r["pmcid"], r["source"], issn,
                "" if v is None else "%.4f" % v, r["year"]]) + "\n")
    print("[tsv] wrote %d rows -> %s" % (len(rows), OUT_TSV), file=sys.stderr)


# --------------------------------------------------------------------------- #
# Output -- HTML summary
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
 .warn { background: #fff8c5; border-left: 4px solid #d4a72c; padding: .6rem .9rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 ol.strategy > li { margin: .45rem 0; }
 .cols { display:grid; grid-template-columns: 1fr 1fr; gap:1.5rem; align-items:start; }
 @media (max-width:780px){ .cols{ grid-template-columns:1fr; } }
"""

STRATEGY_HTML = """
<ol class="strategy">
  <li><strong>Retrieve all PMIDs (beat the 10k cap).</strong> PubMed's esearch returns at most ~10,000
      UIDs and cannot page beyond that index even via the history server. The query is partitioned by
      <em>publication year</em> (<code>datetype=pdat</code>) and any window still over 9,999 is
      recursively subdivided (year &rarr; month &rarr; day) until each leaf is fully retrievable; the
      de-duplicated union is the complete set.</li>
  <li><strong>Annotate PMID &rarr; PMCID + journal + ISSN + year in one pass.</strong> <code>esummary</code>
      (POST, 200 ids/call) returns a docsum whose <code>articleids</code> carry the PMC id, whose
      <code>source</code>/<code>issn</code>/<code>essn</code> give the journal and ISSNs, and whose
      <code>sortpubdate</code> gives the year. No PMC id &rArr; not in PubMed Central.</li>
  <li><strong>Journal impact factor.</strong> The official Journal Impact Factor (Clarivate JCR) is
      proprietary with no free API, so the default metric is <strong>OpenAlex
      <code>2yr_mean_citedness</code></strong> &mdash; mean citations over the last two years to a
      journal's papers, the same construction as JIF &mdash; an open, reproducible proxy looked up per
      ISSN (batched 50/request). A curated <code>issn&rarr;jif</code> table (env <code>JIF_TABLE</code>)
      overrides the proxy when supplied, so true JCR figures can be plugged in.</li>
  <li><strong>Polite &amp; resumable.</strong> Requests are throttled (~3/s NCBI, ~7/s OpenAlex) with
      retries and tolerant JSON parsing; each phase checkpoints after every batch so runs resume
      cleanly.</li>
</ol>
"""


def _fmt(x):
    return format(x, ",")


def build_html(query, rows, jif):
    import html as _html

    n = len(rows)
    per_year, pmc_year = {}, {}
    pmc_total = 0
    if_vals = []
    jcount, jif_by_journal = {}, {}
    for r in rows:
        y = r["year"]
        if y:
            per_year[y] = per_year.get(y, 0) + 1
        if r["pmcid"]:
            pmc_total += 1
            if y:
                pmc_year[y] = pmc_year.get(y, 0) + 1
        v = jif.get(r["essn"]) if r["essn"] in jif else jif.get(r["issn"])
        if r["source"]:
            jcount[r["source"]] = jcount.get(r["source"], 0) + 1
            if v is not None:
                jif_by_journal[r["source"]] = v
        if v is not None:
            if_vals.append(v)

    years = sorted((int(y) for y in per_year), reverse=True)
    maxc = max(per_year.values()) if per_year else 1
    if_vals.sort()

    def pctile(p):
        if not if_vals:
            return None
        return if_vals[min(len(if_vals) - 1, int(p * len(if_vals)))]

    buckets = [("0-1", 0, 1), ("1-2", 1, 2), ("2-3", 2, 3), ("3-5", 3, 5),
               ("5-10", 5, 10), ("10-20", 10, 20), ("20-50", 20, 50), ("50+", 50, 1e9)]
    bdist = {b[0]: 0 for b in buckets}
    for v in if_vals:
        for name, lo, hi in buckets:
            if lo <= v < hi:
                bdist[name] += 1
                break

    nif = len(if_vals)
    pmc_prop = pmc_total / n if n else 0
    ifcov = nif / n if n else 0
    qd = _html.escape(query)

    H = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<title>PubMed query results &mdash; PMIDs, PMC coverage &amp; impact factor</title>",
         "<style>" + CSS + "</style></head><body>"]
    H.append("<h1>PubMed query results &mdash; PMIDs, PMC coverage &amp; impact factor</h1>")
    H.append("<p class='meta'>Query (read strictly from &lt;STDIN&gt;): <code>%s</code> "
             "&middot; Source: NCBI E-utilities (esearch/esummary) + OpenAlex "
             "&middot; impact-factor metric: %s</p>"
             % (qd, "curated JIF table + OpenAlex 2-yr mean citedness" if JIF_TABLE
                else "OpenAlex 2-yr mean citedness (JIF proxy)"))

    # stat cards
    cards = [(_fmt(n), "unique PMIDs"),
             (_fmt(pmc_total), "with PMC id (%.1f%%)" % (pmc_prop * 100)),
             (_fmt(nif), "with impact factor (%.1f%%)" % (ifcov * 100)),
             (_fmt(len(jcount)), "distinct journals"),
             ("%.2f" % pctile(0.5) if if_vals else "n/a", "median impact factor"),
             (str(len(years)), "years covered")]
    H.append("<div class='stat-grid'>")
    for v, k in cards:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")

    H.append("<div class='key'>Retrieved <strong>%s</strong> unique PMIDs for "
             "<code>%s</code>, annotated each with its PubMed Central id where one exists "
             "(<strong>%s</strong>, <strong>%.1f%%</strong>), the source journal + ISSN, and a "
             "journal impact factor (<strong>%.1f%%</strong> of articles). Row-level data: "
             "<code>pmids/pmid_pmc_ids.tsv</code>.</div>"
             % (_fmt(n), qd, _fmt(pmc_total), pmc_prop * 100, ifcov * 100))

    H.append("<h2>1. Strategy</h2>" + STRATEGY_HTML)

    # per-year + pmc
    H.append("<h2>2. PMIDs per year &amp; PMC coverage</h2><div class='cols'><div><table>")
    H.append("<thead><tr><th class='num'>year</th><th class='num'>PMIDs</th><th>dist.</th>"
             "<th class='num'>in PMC</th><th class='num'>PMC %</th></tr></thead><tbody>")
    for y in years:
        ys = str(y)
        c = per_year[ys]
        pc = pmc_year.get(ys, 0)
        w = max(2, round(c / maxc * 260))
        H.append("<tr><td class='num'>%d</td><td class='num'>%s</td>"
                 "<td><span class='bar' style='width:%dpx'></span></td>"
                 "<td class='num'>%s</td><td class='num dim'>%.0f%%</td></tr>"
                 % (y, _fmt(c), w, _fmt(pc), (pc / c * 100 if c else 0)))
    H.append("</tbody></table></div><div>")
    # decade rollup
    dec = {}
    for ys, c in per_year.items():
        d = (int(ys) // 10) * 10
        dec[d] = dec.get(d, 0) + c
    H.append("<h3 style='margin-top:0'>By decade</h3><table>"
             "<thead><tr><th class='num'>decade</th><th class='num'>PMIDs</th></tr></thead><tbody>")
    for d in sorted(dec, reverse=True):
        H.append("<tr><td class='num'>%ds</td><td class='num'>%s</td></tr>" % (d, _fmt(dec[d])))
    H.append("</tbody></table><p class='dim'>Year = publication year of the record "
             "(esummary <code>sortpubdate</code>).</p></div></div>")

    # pmc proportion
    H.append("<h2>3. PubMed Central id proportion</h2>")
    H.append("<p><strong>%s of %s articles (%.1f%%)</strong> carry a PMC id. Coverage is strongly "
             "time-dependent &mdash; near zero before the PMC era, rising in recent years &mdash; "
             "reflecting open-access mandates rather than indexing of older work.</p>"
             % (_fmt(pmc_total), _fmt(n), pmc_prop * 100))

    # impact factor breakdown
    H.append("<h2>4. Impact-factor breakdown</h2>")
    if if_vals:
        stats = [("%.2f" % (sum(if_vals) / nif), "mean"),
                 ("%.2f" % pctile(0.5), "median"),
                 ("%.2f" % pctile(0.25), "25th pct"),
                 ("%.2f" % pctile(0.75), "75th pct"),
                 ("%.2f" % pctile(0.90), "90th pct"),
                 ("%.0f" % if_vals[-1], "max")]
        H.append("<div class='stat-grid'>")
        for v, k in stats:
            H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
        H.append("</div>")
        bmax = max(bdist.values()) or 1
        H.append("<table><thead><tr><th>IF band</th><th class='num'>articles</th><th>share</th>"
                 "<th class='num'>%</th></tr></thead><tbody>")
        for name, _lo, _hi in buckets:
            v = bdist[name]
            w = max(2, round(v / bmax * 320))
            H.append("<tr><td>%s</td><td class='num'>%s</td>"
                     "<td><span class='bar o' style='width:%dpx'></span></td>"
                     "<td class='num dim'>%.1f%%</td></tr>"
                     % (name, _fmt(v), w, (v / nif * 100 if nif else 0)))
        H.append("</tbody></table>")
        H.append("<p class='dim'>Article-weighted over the %s articles whose journal matched a "
                 "value (%.1f%% of all).</p>" % (_fmt(nif), ifcov * 100))
    else:
        H.append("<p class='dim'>No impact-factor values were resolved.</p>")

    # top journals
    H.append("<h2>5. Top journals by article count</h2>")
    top = sorted(jcount.items(), key=lambda kv: kv[1], reverse=True)[:25]
    jmax = top[0][1] if top else 1
    H.append("<table><thead><tr><th>journal (NLM abbrev)</th><th class='num'>articles</th><th>dist.</th>"
             "<th class='num'>impact factor</th></tr></thead><tbody>")
    for j, c in top:
        w = max(2, round(c / jmax * 220))
        iv = jif_by_journal.get(j)
        H.append("<tr><td>%s</td><td class='num'>%s</td>"
                 "<td><span class='bar g' style='width:%dpx'></span></td>"
                 "<td class='num'>%s</td></tr>"
                 % (_html.escape(j), _fmt(c), w, ("n/a" if iv is None else "%.2f" % iv)))
    H.append("</tbody></table>")

    # caveats
    H.append("""<h2>6. Caveats</h2><ul>
  <li><strong>Impact factor is a proxy unless a curated table is supplied.</strong> OpenAlex
      <code>2yr_mean_citedness</code> is built like JIF but on OpenAlex's citation graph, so absolute
      values differ from Clarivate JCR; set <code>JIF_TABLE</code> to use official figures.</li>
  <li><strong>ISSN&rarr;journal matching is imperfect.</strong> Print vs electronic ISSN can resolve to
      a secondary/legacy OpenAlex source, occasionally giving an implausibly low value for a major
      journal. The <em>distribution</em> is more robust than any single journal's number.</li>
  <li><strong>PMC id = availability flag</strong>, not a guarantee of open full text (some records are
      embargoed or author-manuscript only).</li>
  <li><strong>Live data.</strong> PubMed, PMC and OpenAlex update continuously, so re-running shifts
      recent years slightly.</li>
</ul>""")
    H.append("<p class='dim'>Reproduce: <code>echo \"%s\" | python pubmed_query.py</code></p>" % qd)
    H.append("</body></html>")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(SUMMARY_HTML, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    print("[html] wrote %s" % SUMMARY_HTML, file=sys.stderr)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def read_query():
    """Prompt for the query with the exact phrase, then read one line from STDIN.

    The query comes strictly from STDIN -- there is deliberately no argv fallback,
    so the reported query is always exactly what was piped/typed in. Using input()
    (one line, returns on Enter) instead of sys.stdin.read() (which blocks until
    EOF) avoids the 'stuck waiting for input' behaviour.
    """
    try:
        q = input("enter command line argument: ").strip()
    except EOFError:
        q = ""                       # empty/closed STDIN
    return q


def main():
    query = read_query()
    if not query:
        sys.exit("error: no query provided on STDIN")
    os.makedirs(PMID_DIR, exist_ok=True)
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    deadline = (time.time() + TIME_BUDGET) if TIME_BUDGET > 0 else 0

    pmids = get_all_pmids(query)
    if not pmids:
        # No PMIDs -> no corpus. Stop with a non-zero exit so the orchestrator
        # aborts at step 1 with a clear reason, rather than writing an empty TSV
        # that makes step 2 (high_impact_xml.py) fail with "no impact-factor values".
        sys.exit("error: query matched 0 PubMed records -- nothing to build a corpus "
                 "from; refine the query and re-run")

    if not annotate(pmids, deadline):
        print("ANNOTATE_INCOMPLETE - re-run to resume", file=sys.stderr)
        return
    rows = load_annotations()

    issns = set()
    for r in rows:
        if r["issn"]:
            issns.add(r["issn"])
        if r["essn"]:
            issns.add(r["essn"])
    if not impact_factors(issns, deadline):
        print("IMPACTFACTOR_INCOMPLETE - re-run to resume", file=sys.stderr)
        return

    jif = load_jif_map()
    write_tsv(rows, jif)
    build_html(query, rows, jif)
    print("ALL DONE: %d pmids, %s -> %s & %s"
          % (len(rows), query, OUT_TSV, SUMMARY_HTML), file=sys.stderr)


if __name__ == "__main__":
    main()
