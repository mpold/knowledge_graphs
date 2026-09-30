#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
from_archive.py
===============
Serve the part of a PubMed result set that is **already on disk** from a local
XML archive, instead of re-downloading it from NCBI.

This is stage **1b** of the Step 1 publications pipeline (see
``step_1_publications.html``): it sits between ``pubmed_query.py`` (stage 1,
which resolves the query to ``pmids/pmid_pmc_ids.tsv``) and
``high_impact_xml.py`` (stage 2, which eFetches the article XML). It answers one
question -- *which publications returned by the upstream query do we already
have?* -- and acts on the answer in three ways:

  1. it **copies** each archive hit into ``gpu_bundle/experimental_ner/`` (the
     corpus the GPU bundle consumes) and into ``archive_xmls/`` (a durable
     staging copy of the same set, see RE-SEEDING below);
  2. it writes ``pmids/from_archive_pmcids.txt``, the list of PMC ids that
     ``high_impact_xml.py`` should therefore **not** download; and
  3. it writes ``summaries/from_archive.html`` describing the match.

Input  : pmids/pmid_pmc_ids.tsv                    (from pubmed_query.py)
         <archive>/PMC*.xml, PMC*.grobid.tei.xml   (default: nearest ancestor xmls/)
Outputs: archive_xmls/PMC*.{xml,grobid.tei.xml}    durable copy of the hit set
         gpu_bundle/experimental_ner/PMC*.{...}    seeded into the NER corpus
         pmids/from_archive_pmcids.txt             download skip-list for stage 2
         summaries/from_archive.html               summary of the match

NOT executed here -- run it yourself; re-run freely (it is idempotent).

STRATEGY
--------
1. Identity = PMC id. The archive and every pipeline stage name files by PMC id
   -- JATS as ``PMC<digits>.xml``, GROBID TEI as ``PMC<digits>.grobid.tei.xml``
   -- so the filename stem is the shared key, exactly as in
   ``named_entity_xml.py``. The archive is indexed by that stem in one scan.

2. Match against the upstream query. Every row of ``pmids/pmid_pmc_ids.tsv``
   that carries a ``pmc_id`` is a candidate; the intersection with the archive
   index is the hit set. The impact-factor percentile is deliberately *not*
   applied here -- the archive is free, so anything the query returned and we
   already own is worth keeping, and stage 2's percentile still governs what
   gets *downloaded*. (Set ``MIN_IF`` to filter the archive hits by impact
   factor anyway.)

3. One representative per publication, best text wins. A paper can be in the
   archive as JATS, as GROBID TEI, or both. GROBID TEI only ever exists because
   the JATS export had no ``<body>``, so the preference is: **JATS with a body**
   > **GROBID TEI** > **abstract-only JATS**. That yields the richest available
   full text per paper and never emits two files for one PMC id.

4. Copy, never move. The archive is a shared repository that other projects read
   -- it is opened read-only and nothing in it is renamed, moved or deleted.
   Copies are ``shutil.copy2`` via a temp file + ``os.replace`` so an interrupted
   run never leaves a half-file that looks complete.

5. Tell stage 2 what to skip. ``pmids/from_archive_pmcids.txt`` lists the hit
   PMC ids; ``high_impact_xml.py`` drops them from its ``pmid2pmcid`` export
   dictionary, so the eFetch loop only fetches what the archive did not supply.
   Delete that file (or set ``USE_ARCHIVE_SKIP=0`` on stage 2) to download
   everything regardless.

6. Summarise. ``summaries/from_archive.html`` reports coverage (query hits vs
   archive hits vs still-to-download), which scheme served each hit, full-text
   vs abstract-only, the year and journal profile of the served set, and the
   article-type breakdown of what was seeded into ``experimental_ner/``.

RE-SEEDING -- WHY THIS SCRIPT RUNS TWICE
----------------------------------------
Stage 6 (``named_entity_xml.py``) **clean-rebuilds** ``gpu_bundle/experimental_ner/``
(``shutil.rmtree`` + recreate) from the freshly downloaded corpus, which would
discard the archive contribution. So ``from_archive.py`` is run at two points:

    stage 1b  (after pubmed_query.py)   -- match, stage, seed, write the skip-list
    stage 6b  (after named_entity_xml.py) -- re-seed experimental_ner/

The second run is cheap: the hit set is already in ``archive_xmls/``, so it only
re-copies from there. Unlike the other bulk stages this script **never empties**
``gpu_bundle/experimental_ner/`` -- it merges into it, because that directory
holds the union of the downloaded corpus and the archive contribution.
``archive_xmls/`` *is* pruned to the current hit set, so changing the query does
not leave stale papers behind.

ENVIRONMENT
-----------
ARCHIVE_DIR   archive to read (default: the nearest ``xmls/`` in an ancestor of this
              script's directory, e.g. ``../xmls``)
MIN_IF        keep only archive hits with journal impact factor >= this value
ORIGINAL_ONLY 1 = seed experimental_ner/ with original-results article types only
              (same RESULTS_CATEGORY mapping as named_entity_xml.py); default 0
              copies every hit, since experimental_ner/ is the corpus the user
              asked the archive to be copied into
NO_SKIP_LIST  1 = do not write pmids/from_archive_pmcids.txt (stage 2 then
              re-downloads the archive hits as well)
DRY_RUN       1 = plan + write the summary without copying anything
"""

import os
import csv
import sys
import shutil
import html as _html
from collections import Counter

# --------------------------------------------------------------------------- #
# Paths / config -- all script-relative, so this runs from any working directory
# --------------------------------------------------------------------------- #
BASE         = os.path.dirname(os.path.abspath(__file__))
IN_TSV       = os.path.join(BASE, "pmids", "pmid_pmc_ids.tsv")


def _find_default_archive(start):
    """Nearest ``xmls/`` directory in an ancestor of ``start`` (``start`` excluded).

    The shared archive sits beside the project checkout (``PythonProject/xmls``
    for ``PythonProject/anaerobic``) or, in the older layout, beside
    ``relationship_graphs/`` two levels up -- so walk upwards rather than hard-code
    the depth. Falls back to ``../xmls`` (reported as missing) when none exists.
    """
    here = os.path.abspath(start)
    parent = os.path.dirname(here)
    while parent != here:
        candidate = os.path.join(parent, "xmls")
        if os.path.isdir(candidate):
            return candidate
        here, parent = parent, os.path.dirname(parent)
    return os.path.join(os.path.dirname(os.path.abspath(start)), "xmls")


ARCHIVE_DIR  = os.path.abspath(os.environ.get("ARCHIVE_DIR", "").strip()
                               or _find_default_archive(BASE))
STAGE_DIR    = os.path.join(BASE, "archive_xmls")                        # durable copy
EXPERIMENTAL_DIR = os.path.join(BASE, "gpu_bundle", "experimental_ner")  # seeded corpus
SKIP_LIST    = os.path.join(BASE, "pmids", "from_archive_pmcids.txt")
SUMMARY_DIR  = os.path.join(BASE, "summaries")
SUMMARY_HTML = os.path.join(SUMMARY_DIR, "from_archive.html")

GROBID_SUFFIX = ".grobid.tei.xml"
JATS_SUFFIX   = ".xml"

MIN_IF        = os.environ.get("MIN_IF", "").strip()
ORIGINAL_ONLY = os.environ.get("ORIGINAL_ONLY", "0") not in ("0", "", "false", "False")
NO_SKIP_LIST  = os.environ.get("NO_SKIP_LIST", "0") not in ("0", "", "false", "False")
DRY_RUN       = os.environ.get("DRY_RUN", "0") not in ("0", "", "false", "False")

HEAD_BYTES    = 1000000        # how much of a file to read when classifying it

PROMPT = (
    "design 'from_archive.py' that determines the publications from the upstream "
    "query step that already exist in the '*\\PythonProject\\xmls', and the copy "
    "them to 'gpu_bundle/experimental_ner'; then proceed with xml downloads and "
    "rest of the pipeline described in 'step_1_publications.html'"
)

# Article-type -> results bucket, mirroring named_entity_xml.py so a hit seeded
# from the archive is classified by exactly the same rule as a downloaded one.
RESULTS_ORIGINAL    = "original-results"
RESULTS_SECONDARY   = "secondary-synthesis"
RESULTS_NONRESEARCH = "non-research"

RESULTS_CATEGORY = {
    "research-article":    RESULTS_ORIGINAL,
    "brief-report":        RESULTS_ORIGINAL,
    "case-report":         RESULTS_ORIGINAL,
    "rapid-communication": RESULTS_ORIGINAL,
    "report":              RESULTS_ORIGINAL,
    "review-article":      RESULTS_SECONDARY,
    "systematic-review":   RESULTS_SECONDARY,
    "meeting-report":      RESULTS_SECONDARY,
}


def results_category(article_type):
    return RESULTS_CATEGORY.get(article_type, RESULTS_NONRESEARCH)


# --------------------------------------------------------------------------- #
# 1. Index the archive by PMC id
# --------------------------------------------------------------------------- #
def index_archive(directory):
    """{pmcid: {"jats": path|None, "grobid": path|None}} for a flat XML archive.

    One ``os.scandir`` pass over the directory -- no file is opened here, so the
    index is cheap even for a six-figure archive. ``*.grobid.tei.xml`` also ends
    in ``.xml``, so TEI is tested first.
    """
    index = {}
    if not os.path.isdir(directory):
        return index
    with os.scandir(directory) as it:
        for entry in it:
            if not entry.is_file():
                continue
            name = entry.name
            if name.endswith(GROBID_SUFFIX):
                pmc, scheme = name[:-len(GROBID_SUFFIX)], "grobid"
            elif name.endswith(JATS_SUFFIX):
                pmc, scheme = name[:-len(JATS_SUFFIX)], "jats"
            else:
                continue
            index.setdefault(pmc, {"jats": None, "grobid": None})[scheme] = entry.path
    return index


# --------------------------------------------------------------------------- #
# 2. The upstream query's publications
# --------------------------------------------------------------------------- #
def parse_if(s):
    try:
        return float((s or "").strip())
    except (TypeError, ValueError):
        return None


def load_query_rows():
    """Rows of pmids/pmid_pmc_ids.tsv that carry a pmc_id, keyed by PMC id.

    Rows without a pmc_id have no PMC record at all -- they can be neither
    archived nor downloaded, so they are counted for the summary and dropped.
    """
    with open(IN_TSV, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    with_pmc, floor = {}, parse_if(MIN_IF) if MIN_IF else None
    for r in rows:
        pmc = (r.get("pmc_id") or "").strip()
        if not pmc:
            continue
        v = parse_if(r.get("journal_impact_factor"))
        if floor is not None and (v is None or v < floor):
            continue
        with_pmc[pmc] = {"pmid": (r.get("pmid") or "").strip(),
                         "journal": (r.get("source_publication") or "").strip(),
                         "issn": (r.get("issn") or "").strip(),
                         "if": v,
                         "year": (r.get("year") or "").strip()}
    return len(rows), with_pmc


# --------------------------------------------------------------------------- #
# 3. Pick one representative file per hit
# --------------------------------------------------------------------------- #
_ARTICLE_TYPE_MARK = 'article-type="'


def read_head(path):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            return fh.read(HEAD_BYTES)
    except OSError:
        return ""


def has_body(head):
    """True when a JATS document carries a full-text <body> (not front matter)."""
    return "<body" in head


def jats_article_type(head):
    """article-type of the opening <article> tag, or 'unknown'."""
    i = head.find("<article")
    if i < 0:
        return "unknown"
    end = head.find(">", i)
    tag = head[i:end if end > 0 else i + 2000]
    j = tag.find(_ARTICLE_TYPE_MARK)
    if j < 0:
        return "unknown"
    j += len(_ARTICLE_TYPE_MARK)
    k = tag.find('"', j)
    return tag[j:k] if k > 0 else "unknown"


def choose_representative(entry):
    """Best single archive file for one publication.

    GROBID TEI exists in the archive only because that paper's JATS export had
    no ``<body>``, so full-text JATS is the richer source when it exists and the
    TEI is the fallback; an abstract-only JATS is the last resort. Returns
    (path, scheme, status, article_type).
    """
    jats, grobid = entry.get("jats"), entry.get("grobid")
    jats_head = read_head(jats) if jats else ""
    atype = jats_article_type(jats_head) if jats_head else "unknown (no JATS twin)"

    if jats and has_body(jats_head):
        return jats, "jats", "fulltext", atype
    if grobid:
        return grobid, "grobid", "fulltext", atype
    if jats:
        return jats, "jats", "abstract-only", atype
    return None, None, "missing", atype


# --------------------------------------------------------------------------- #
# 4. Copy (never move) -- atomic, idempotent
# --------------------------------------------------------------------------- #
def copy_into(src, dst_dir):
    """Copy src into dst_dir. Returns 'copied', 'present', 'planned' or 'missing'.

    Written to ``<name>.part`` then ``os.replace``d, so an interrupted copy never
    leaves a truncated file that a later run would mistake for done. A
    destination that already matches the source size is left alone.
    """
    if not src or not os.path.exists(src):
        return "missing"
    dst = os.path.join(dst_dir, os.path.basename(src))
    if os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(src):
        return "present"
    if DRY_RUN:
        return "planned"
    tmp = dst + ".part"
    shutil.copy2(src, tmp)
    os.replace(tmp, dst)
    return "copied"


def prune_stage_dir(keep_names):
    """Drop files from archive_xmls/ that are no longer in the hit set.

    Only ``archive_xmls/`` is pruned -- it is this script's own staging area and
    must track the current query. ``gpu_bundle/experimental_ner/`` is never
    emptied here: it holds the union of the downloaded corpus and the archive
    contribution.
    """
    removed = 0
    if not os.path.isdir(STAGE_DIR):
        return removed
    with os.scandir(STAGE_DIR) as it:
        stale = [e.path for e in it
                 if e.is_file() and e.name not in keep_names and not e.name.endswith(".py")]
    for path in stale:
        if not DRY_RUN:
            os.remove(path)
        removed += 1
    return removed


# --------------------------------------------------------------------------- #
# 5. Skip-list for high_impact_xml.py
# --------------------------------------------------------------------------- #
def write_skip_list(pmcids):
    """PMC ids stage 2 must not re-download, one per line (# comments allowed)."""
    os.makedirs(os.path.dirname(SKIP_LIST), exist_ok=True)
    tmp = SKIP_LIST + ".part"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write("# PMC ids served from the local archive by from_archive.py -- "
                 "high_impact_xml.py skips these.\n")
        fh.write("# Archive: %s\n" % ARCHIVE_DIR)
        for pmc in sorted(pmcids):
            fh.write(pmc + "\n")
    os.replace(tmp, SKIP_LIST)


# --------------------------------------------------------------------------- #
# 6. Summary HTML (house style of the other stage summaries)
# --------------------------------------------------------------------------- #
CSS = """
 body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        margin: 2rem auto; max-width: 1100px; color: #222; line-height: 1.45; padding: 0 1rem; }
 h1 { margin-bottom: .25rem; } h2 { margin-top: 2.25rem; border-bottom: 1px solid #ddd; padding-bottom: .25rem; }
 h3 { margin: 1.2rem 0 .3rem; }
 .meta { color: #555; margin-bottom: 1rem; }
 .stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: .75rem; margin: 1rem 0 1.5rem; }
 .stat { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: .75rem 1rem; }
 .stat .v { font-size: 1.4rem; font-weight: 600; } .stat .k { color: #555; font-size: .85rem; }
 table { border-collapse: collapse; width: 100%; margin: .5rem 0 1rem; font-size: .9rem; }
 th, td { border: 1px solid #e1e4e8; padding: .35rem .55rem; text-align: left; vertical-align: middle; }
 th { background: #f6f8fa; }
 td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
 code { background: #f6f8fa; padding: 1px 4px; border-radius: 3px; font-size: .88em; }
 pre { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: .9rem 1rem; overflow-x: auto; font-size: .85rem; line-height: 1.4; white-space: pre-wrap; }
 .bar { display:inline-block; height:.72em; background:#3b7dd8; border-radius:2px; vertical-align: middle; }
 .bar.g { background:#2da44e; } .bar.o { background:#bf8700; }
 .dim { color: #888; font-size: .85em; }
 .key  { background: #ddf4ff; border-left: 4px solid #0969da; padding: .6rem .9rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 .note { background: #fff8c5; border-left: 4px solid #d4a72c; padding: .5rem .75rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 ol.strategy > li { margin: .45rem 0; }
"""

STRATEGY_HTML = """
<ol class="strategy">
  <li><strong>Identity = PMC id.</strong> The archive and every pipeline stage name files by PMC id
      &mdash; JATS as <code>PMC&lt;digits&gt;.xml</code>, GROBID TEI as
      <code>PMC&lt;digits&gt;.grobid.tei.xml</code> &mdash; so the filename stem is the shared key
      (the same identity <code>named_entity_xml.py</code> uses). The archive is indexed by that stem
      in a single <code>os.scandir</code> pass, opening no files.</li>
  <li><strong>Match against the upstream query.</strong> Every row of
      <code>pmids/pmid_pmc_ids.tsv</code> carrying a <code>pmc_id</code> is a candidate; the
      intersection with the archive index is the hit set. The impact percentile is deliberately
      <em>not</em> applied &mdash; archive files are free, and stage 2's percentile still governs what
      is <em>downloaded</em> (set <code>MIN_IF</code> to filter the archive hits too).</li>
  <li><strong>One representative per publication, best text wins.</strong> A paper may be in the
      archive as JATS, as GROBID TEI, or both. TEI only ever exists because that paper's JATS export
      had no <code>&lt;body&gt;</code>, so the preference is <strong>full-text JATS &rarr; GROBID TEI
      &rarr; abstract-only JATS</strong>: the richest text available, and never two files for one
      PMC id.</li>
  <li><strong>Copy, never move.</strong> The archive is a shared repository read by other projects; it
      is opened read-only and nothing in it is renamed, moved or deleted. Each copy goes to a
      <code>.part</code> temp file and is <code>os.replace</code>d into position, so an interrupted run
      leaves no truncated file that a later run would treat as done.</li>
  <li><strong>Tell stage 2 what to skip.</strong> <code>pmids/from_archive_pmcids.txt</code> lists the
      hit PMC ids; <code>high_impact_xml.py</code> drops them from its <code>pmid2pmcid</code> export
      dictionary, so the eFetch loop only downloads what the archive did not supply. Delete that file
      (or run stage 2 with <code>USE_ARCHIVE_SKIP=0</code>) to download everything regardless.</li>
  <li><strong>Summarise.</strong> This page is built from the match and the files on disk.</li>
</ol>
"""


def _fmt(x):
    return format(x, ",")


def build_summary(stats):
    """Write summaries/from_archive.html from the match + the files on disk."""
    n_rows      = stats["n_rows"]
    n_with_pmc  = stats["n_with_pmc"]
    hits        = stats["hits"]              # {pmc: {...}}
    n_hits      = len(hits)
    n_to_fetch  = n_with_pmc - n_hits
    scheme      = stats["scheme_counts"]
    status      = stats["status_counts"]
    types       = stats["type_counts"]
    copy_stage  = stats["copy_stage"]
    copy_exp    = stats["copy_exp"]
    pruned      = stats["pruned"]
    n_seeded    = stats["n_seeded"]
    n_exp_now   = stats["n_exp_now"]

    cover = (100.0 * n_hits / n_with_pmc) if n_with_pmc else 0.0
    n_original = sum(c for t, c in types.items() if results_category(t) == RESULTS_ORIGINAL)

    per_year, jcount = Counter(), Counter()
    for m in hits.values():
        if m["year"]:
            per_year[m["year"]] += 1
        if m["journal"]:
            jcount[m["journal"]] += 1

    H = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<title>from_archive &mdash; publications served from the local XML archive</title>",
         "<style>" + CSS + "</style></head><body>"]
    H.append("<h1>Publications already in the local archive</h1>")
    H.append("<p class='meta'>Generated by <code>from_archive.py</code> (stage&nbsp;1b of the Step&nbsp;1 "
             "publications pipeline) &middot; archive: <code>%s</code> &middot; identity: "
             "<strong>PMC id</strong> &middot; the archive is read <strong>read-only</strong> and "
             "copied, never moved%s</p>"
             % (_html.escape(ARCHIVE_DIR), " &middot; <strong>DRY_RUN</strong>" if DRY_RUN else ""))

    cards = [(_fmt(n_rows), "publications from the query"),
             (_fmt(n_with_pmc), "with a PMC id (fetchable)"),
             (_fmt(n_hits), "already in the archive (%.1f%%)" % cover),
             (_fmt(n_to_fetch), "left for stage 2 to download"),
             (_fmt(stats["n_archive"]), "publications in the archive"),
             (_fmt(n_exp_now), "files now in experimental_ner/")]
    H.append("<div class='stat-grid'>")
    for v, k in cards:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")

    H.append("<div class='key'>Of the <strong>%s</strong> query publications that have a PMC id, "
             "<strong>%s (%.1f%%)</strong> were already present in <code>%s</code> and were copied "
             "into <code>gpu_bundle/experimental_ner/</code> (and staged in <code>archive_xmls/</code>). "
             "The remaining <strong>%s</strong> are left to <code>high_impact_xml.py</code>, which "
             "reads <code>pmids/from_archive_pmcids.txt</code> and skips the archive-served ids &mdash; "
             "so no publication is fetched twice.</div>"
             % (_fmt(n_with_pmc), _fmt(n_hits), cover, _html.escape(ARCHIVE_DIR), _fmt(n_to_fetch)))

    H.append("<h2>1. Prompt</h2><pre>%s</pre>" % _html.escape(PROMPT))
    H.append("<h2>2. Strategy</h2>" + STRATEGY_HTML)

    # 3. Coverage
    H.append("<h2>3. Coverage of the upstream query</h2><table>")
    H.append("<thead><tr><th>set</th><th class='num'>publications</th><th class='num'>share</th>"
             "<th>disposition</th></tr></thead><tbody>")
    base = n_rows or 1
    for label, count, note in [
            ("query result (<code>pmid_pmc_ids.tsv</code>)", n_rows, "everything the query returned"),
            ("&mdash; without a PMC id", n_rows - n_with_pmc,
             "no PMC record exists; neither archivable nor fetchable"),
            ("&mdash; with a PMC id", n_with_pmc, "candidates for the archive match"),
            ("&nbsp;&nbsp;&mdash; served from the archive", n_hits,
             "copied to <code>experimental_ner/</code> &amp; <code>archive_xmls/</code>"),
            ("&nbsp;&nbsp;&mdash; not in the archive", n_to_fetch,
             "downloaded by <code>high_impact_xml.py</code> (subject to its percentile)")]:
        H.append("<tr><td>%s</td><td class='num'>%s</td><td class='num dim'>%.1f%%</td>"
                 "<td class='dim'>%s</td></tr>" % (label, _fmt(count), 100.0 * count / base, note))
    H.append("</tbody></table>")

    # 4. Which scheme served each hit
    H.append("<h2>4. Which archive file served each publication</h2>")
    H.append("<p>GROBID TEI exists in the archive only for papers whose JATS export carried no "
             "<code>&lt;body&gt;</code>, so full-text JATS is preferred, TEI is the fallback that "
             "recovers the body, and an abstract-only JATS is the last resort.</p>")
    H.append("<table><thead><tr><th>representative</th><th>scheme</th><th class='num'>publications</th>"
             "<th class='num'>share</th><th>text recovered</th></tr></thead><tbody>")
    rep_rows = [("full-text JATS", "jats", "fulltext", "<code>&lt;body&gt;</code> present in the NCBI export"),
                ("GROBID TEI", "grobid", "fulltext", "body recovered from the PDF by GROBID"),
                ("abstract-only JATS", "jats", "abstract-only", "front matter only &mdash; publisher-restricted")]
    for label, sch, st, note in rep_rows:
        c = stats["rep_counts"].get((sch, st), 0)
        H.append("<tr><td><strong>%s</strong></td><td><code>%s</code></td><td class='num'>%s</td>"
                 "<td class='num dim'>%.1f%%</td><td class='dim'>%s</td></tr>"
                 % (label, sch, _fmt(c), (100.0 * c / n_hits if n_hits else 0), note))
    H.append("<tfoot><tr><th>Total</th><th></th><th class='num'>%s</th><th class='num'>100%%</th>"
             "<th></th></tr></tfoot></table>" % _fmt(n_hits))
    H.append("<p class='dim'>By scheme: %s &middot; by text status: %s</p>"
             % (", ".join("%s=%s" % (k, _fmt(v)) for k, v in sorted(scheme.items())) or "&mdash;",
                ", ".join("%s=%s" % (k, _fmt(v)) for k, v in sorted(status.items())) or "&mdash;"))

    # 5. Article types seeded
    H.append("<h2>5. Article types seeded into <code>experimental_ner/</code></h2>")
    H.append("<div class='note'>Types are read from each hit's JATS <code>article-type</code>; a "
             "GROBID-served paper with no JATS twin in the archive is <code>unknown (no JATS twin)</code>. "
             "The mapping to results buckets is the same <code>RESULTS_CATEGORY</code> table "
             "<code>named_entity_xml.py</code> applies to the downloaded corpus. This stage seeds "
             "<strong>%s</strong> by default; set <code>ORIGINAL_ONLY=1</code> to seed only the "
             "original-results bucket.</div>"
             % ("only the original-results bucket" if ORIGINAL_ONLY else "every hit, of any type"))
    tmax = max(types.values()) if types else 1
    H.append("<table><thead><tr><th>article-type</th><th>results bucket</th><th class='num'>files</th>"
             "<th class='num'>share</th><th>dist.</th></tr></thead><tbody>")
    for t, c in sorted(types.items(), key=lambda kv: (-kv[1], kv[0])):
        b = results_category(t)
        cls = "bar g" if b == RESULTS_ORIGINAL else "bar"
        H.append("<tr><td><code>%s</code></td><td class='dim'>%s</td><td class='num'>%s</td>"
                 "<td class='num dim'>%.1f%%</td>"
                 "<td><span class='%s' style='width:%dpx'></span></td></tr>"
                 % (_html.escape(t), b, _fmt(c), (100.0 * c / n_hits if n_hits else 0),
                    cls, max(1, int(round(280 * c / tmax)))))
    H.append("</tbody></table>")
    H.append("<p class='dim'><strong>%s of %s</strong> archive hits (%.1f%%) are original-results "
             "article types.</p>"
             % (_fmt(n_original), _fmt(n_hits), (100.0 * n_original / n_hits if n_hits else 0)))

    # 6. Years & journals
    H.append("<h2>6. Year &amp; journal profile of the archive-served set</h2><table>")
    H.append("<thead><tr><th class='num'>year</th><th class='num'>publications</th><th>dist.</th>"
             "</tr></thead><tbody>")
    ymax = max(per_year.values()) if per_year else 1
    for y in sorted(per_year, key=lambda s: (-int(s) if s.isdigit() else 0, s))[:25]:
        c = per_year[y]
        H.append("<tr><td class='num'>%s</td><td class='num'>%s</td>"
                 "<td><span class='bar' style='width:%dpx'></span></td></tr>"
                 % (_html.escape(y), _fmt(c), max(1, int(round(260 * c / ymax)))))
    H.append("</tbody></table>")
    H.append("<h3>Top journals</h3><table><thead><tr><th>journal (NLM abbrev)</th>"
             "<th class='num'>publications</th><th>dist.</th></tr></thead><tbody>")
    top = jcount.most_common(20)
    jmax = top[0][1] if top else 1
    for j, c in top:
        H.append("<tr><td>%s</td><td class='num'>%s</td>"
                 "<td><span class='bar g' style='width:%dpx'></span></td></tr>"
                 % (_html.escape(j), _fmt(c), max(1, int(round(220 * c / jmax)))))
    H.append("</tbody></table>")

    # 7. What this run did
    H.append("<h2>7. What this run did</h2><table>")
    H.append("<thead><tr><th>action</th><th>destination</th><th class='num'>copied</th>"
             "<th class='num'>already present</th><th class='num'>missing</th></tr></thead><tbody>")
    for label, dest, c in [("stage the hit set", "<code>archive_xmls/</code>", copy_stage),
                           ("seed the NER corpus", "<code>gpu_bundle/experimental_ner/</code>", copy_exp)]:
        H.append("<tr><td>%s</td><td>%s</td><td class='num'>%s</td><td class='num'>%s</td>"
                 "<td class='num'>%s</td></tr>"
                 % (label, dest, _fmt(c.get("copied", 0) + c.get("planned", 0)),
                    _fmt(c.get("present", 0)), _fmt(c.get("missing", 0))))
    H.append("</tbody></table>")
    H.append("<p class='dim'>%s stale file(s) pruned from <code>archive_xmls/</code> (no longer in the "
             "hit set). <strong>%s</strong> files were seeded into <code>experimental_ner/</code> this "
             "run; the directory now holds <strong>%s</strong> files in total &mdash; it is "
             "<em>merged into, never emptied</em>, because it carries the union of the downloaded corpus "
             "and the archive contribution.</p>" % (_fmt(pruned), _fmt(n_seeded), _fmt(n_exp_now)))

    # 8. Where this sits in the pipeline
    H.append("<h2>8. Where this sits in the pipeline</h2>")
    H.append("<p>Stage&nbsp;1b, between <code>pubmed_query.py</code> and <code>high_impact_xml.py</code> "
             "&mdash; and again as stage&nbsp;6b, after <code>named_entity_xml.py</code>:</p>")
    H.append("<pre>1.  pubmed_query.py          query (STDIN)        -&gt; pmids/pmid_pmc_ids.tsv\n"
             "1b. from_archive.py          pmid_pmc_ids.tsv + archive\n"
             "                                                  -&gt; archive_xmls/, experimental_ner/,\n"
             "                                                     pmids/from_archive_pmcids.txt\n"
             "2.  high_impact_xml.py       tsv + percentile     -&gt; high_impact_xmls/   (skips 1b's ids)\n"
             "3.  xml_structure.py         high_impact_xmls/    -&gt; summaries/xml_structure.html\n"
             "4.  ncbi_pdf.py              no-&lt;body&gt; records     -&gt; ncbi_pdfs_grobid/\n"
             "5.  grobid_xml.py            PDFs (Docker+GROBID) -&gt; grobid_xmls/\n"
             "6.  named_entity_xml.py      grobid + high_impact -&gt; experimental_ner/  (CLEAN REBUILD)\n"
             "6b. from_archive.py          archive_xmls/        -&gt; experimental_ner/  (re-seed)\n"
             "7.  pre_ner_xml_structure.py experimental_ner/    -&gt; summaries/pre_ner_xml_structure.html</pre>")
    H.append("<div class='note'><strong>Why 6b exists.</strong> <code>named_entity_xml.py</code> "
             "clean-rebuilds <code>gpu_bundle/experimental_ner/</code> from the freshly downloaded "
             "corpus, which would discard the archive contribution. Re-running "
             "<code>from_archive.py</code> afterwards restores it from <code>archive_xmls/</code> &mdash; "
             "a cheap, idempotent copy. Run it a third time and nothing changes.</div>")

    # 9. Caveats
    H.append("<h2>9. Notes &amp; caveats</h2><ul>")
    H.append("<li><strong>The archive is never modified.</strong> Files are copied out with "
             "<code>shutil.copy2</code>; nothing in <code>%s</code> is written, moved or deleted.</li>"
             % _html.escape(ARCHIVE_DIR))
    H.append("<li><strong>Matching is by PMC id, not by content.</strong> An archived file is trusted to "
             "be that publication's XML. A truncated or stale archive entry would be copied as-is; the "
             "full-text/abstract-only column above is the sanity check.</li>")
    H.append("<li><strong>No impact-factor filter by default.</strong> Archive hits cost nothing to "
             "reuse, so all of them are served; stage 2's percentile still decides what is downloaded. "
             "Set <code>MIN_IF=&lt;float&gt;</code> to apply a floor here as well.</li>")
    H.append("<li><strong>Article type is advisory.</strong> A GROBID-served paper with no JATS twin in "
             "the archive has no type field at all and is reported as "
             "<code>unknown (no JATS twin)</code> &mdash; under <code>ORIGINAL_ONLY=1</code> such papers "
             "would be dropped, which is why the default seeds every hit.</li>")
    H.append("<li><strong>Idempotent.</strong> Re-running copies only what is missing and prunes "
             "<code>archive_xmls/</code> to the current hit set; <code>DRY_RUN=1</code> plans and writes "
             "this page without touching any file.</li>")
    H.append("</ul></body></html>")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(SUMMARY_HTML, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    print("[summary] wrote %s" % SUMMARY_HTML, file=sys.stderr)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    os.makedirs(os.path.dirname(IN_TSV), exist_ok=True)
    if not os.path.exists(IN_TSV):
        raise SystemExit("error: input not found: %s\n"
                         "       run pubmed_query.py (stage 1) first" % IN_TSV)
    if not os.path.isdir(ARCHIVE_DIR):
        raise SystemExit("error: archive not found: %s\n"
                         "       set ARCHIVE_DIR=<path> to point at it" % ARCHIVE_DIR)

    n_rows, query_pmcs = load_query_rows()
    archive = index_archive(ARCHIVE_DIR)
    print("[scan] archive %s -> %d publications" % (ARCHIVE_DIR, len(archive)), file=sys.stderr)
    print("[query] %d rows, %d with a pmc_id%s"
          % (n_rows, len(query_pmcs), (" (MIN_IF=%s)" % MIN_IF) if MIN_IF else ""), file=sys.stderr)

    # 3. Resolve one representative file per hit.
    hits, rep_counts = {}, Counter()
    scheme_counts, status_counts, type_counts = Counter(), Counter(), Counter()
    for pmc, meta in query_pmcs.items():
        entry = archive.get(pmc)
        if not entry:
            continue
        path, scheme, st, atype = choose_representative(entry)
        if not path:
            continue
        if ORIGINAL_ONLY and results_category(atype) != RESULTS_ORIGINAL:
            continue
        m = dict(meta)
        m.update(path=path, scheme=scheme, status=st, atype=atype)
        hits[pmc] = m
        rep_counts[(scheme, st)] += 1
        scheme_counts[scheme] += 1
        status_counts[st] += 1
        type_counts[atype] += 1
    print("[match] %d of %d fetchable publications are already in the archive"
          % (len(hits), len(query_pmcs)), file=sys.stderr)

    # 4. Copy: archive -> archive_xmls/ (staging) and -> experimental_ner/ (corpus).
    if not DRY_RUN:
        os.makedirs(STAGE_DIR, exist_ok=True)
        os.makedirs(EXPERIMENTAL_DIR, exist_ok=True)
    keep_names = {os.path.basename(m["path"]) for m in hits.values()}
    pruned = prune_stage_dir(keep_names)

    copy_stage, copy_exp = Counter(), Counter()
    try:
        for i, m in enumerate(hits.values(), 1):
            copy_stage[copy_into(m["path"], STAGE_DIR)] += 1
            # Seed from the staged copy when it exists, so stage 6b works even if
            # the archive is offline; fall back to the archive on the first run.
            staged = os.path.join(STAGE_DIR, os.path.basename(m["path"]))
            copy_exp[copy_into(staged if os.path.exists(staged) else m["path"],
                               EXPERIMENTAL_DIR)] += 1
            if i % 500 == 0:
                print("[copy] %d/%d" % (i, len(hits)), file=sys.stderr)
    except KeyboardInterrupt:
        print("\n[copy] interrupted -- copies made so far are complete; re-run to finish",
              file=sys.stderr)

    # 5. Skip-list for stage 2.
    if NO_SKIP_LIST:
        print("[skip-list] NO_SKIP_LIST=1 -- not written; stage 2 will download these too",
              file=sys.stderr)
    elif not DRY_RUN:
        write_skip_list(hits)
        print("[skip-list] wrote %s (%d ids)" % (SKIP_LIST, len(hits)), file=sys.stderr)

    n_exp_now = len([e for e in os.scandir(EXPERIMENTAL_DIR)
                     if e.is_file()]) if os.path.isdir(EXPERIMENTAL_DIR) else 0
    print("[copy] archive_xmls: %s | experimental_ner: %s | pruned %d stale"
          % (dict(copy_stage), dict(copy_exp), pruned), file=sys.stderr)

    build_summary({"n_rows": n_rows, "n_with_pmc": len(query_pmcs), "n_archive": len(archive),
                   "hits": hits, "rep_counts": rep_counts, "scheme_counts": scheme_counts,
                   "status_counts": status_counts, "type_counts": type_counts,
                   "copy_stage": copy_stage, "copy_exp": copy_exp, "pruned": pruned,
                   "n_seeded": copy_exp.get("copied", 0) + copy_exp.get("planned", 0),
                   "n_exp_now": n_exp_now})
    print("[done] %d publications served from the archive; %d left for high_impact_xml.py"
          % (len(hits), len(query_pmcs) - len(hits)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
