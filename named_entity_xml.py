#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
named_entity_xml.py
===================
Assemble the de-duplicated corpus for named-entity recognition: take the two
XML sources -- the GROBID full-text TEI in ``grobid_xmls/`` and the high-impact
JATS XML in ``high_impact_xmls/`` -- identify the *unique publications* across
them and **move** one representative file per publication into
``named_entity_xmls/`` so that **every** XML in ``grobid_xmls/`` is included.
Then determine each moved file's article type and write an HTML summary.

Written *after* ``summaries/bacs/xmls_for_ner.html`` (the predecessor assembly of
``ner_xmls/``) and keeps its two defining ideas:
  * publication identity = **PMC ID** (the filename stem, shared by both schemes), and
  * on overlap the **GROBID full text wins** and the JATS twin is excluded, so
    each paper appears exactly once.

Input  : grobid_xmls/PMC*.grobid.tei.xml                :GROBID TEI (full text from PDFs)
         high_impact_xmls/PMC*.xml                      :high-impact JATS XML (eFetch)
Outputs: named_entity_xmls/PMC*.{xml,grobid.tei.xml}    :one file per publication
         gpu_bundle/experimental_ner/PMC*.{xml,grobid.tei.xml}  :original-results subset (clean-rebuilt in the GPU bundle)
         summaries/named_entity_xmls.html               :summary of both moves + types

NOT executed here -- run it yourself. Set DRY_RUN=1 to plan + (re)write the
summary without moving anything; the page is also rewritten on every run from the
state on disk, so re-running after a partial move simply refreshes it.

STRATEGY
--------
1. Equivalence key = PMC ID. Both directories name files by PMC ID -- JATS as
   ``PMC<digits>.xml``, GROBID as ``PMC<digits>.grobid.tei.xml``. GROBID TEI
   carries no PMID/PMCID *inside* the document (only DOI/ISSN/ORCID), so the
   filename stem is the reliable, shared identity for a publication. The two
   naming schemes are disjoint, so a JATS file and a GROBID file for the same
   paper never collide in the destination.

2. Unique publications, GROBID-preferred. The unique set is the union of the two
   sources by PMC ID. To guarantee "all of grobid_xmls is included", GROBID is the
   winner on any overlap: move *every* GROBID file, and move a JATS file *only if
   its PMC ID is absent from the GROBID set*. The overlapping JATS records (the
   2,609 papers also in grobid_xmls -- publisher-restricted, often abstract-only)
   are therefore excluded in favour of the richer GROBID full text. Net: each
   publication appears once, and the GROBID corpus is fully represented.

3. Move & verify. ``shutil.move`` each selected file into ``named_entity_xmls/``;
   a file already present in the destination is treated as done (resume-safe), and
   moves are verified by source/destination counts. Distinct filename schemes mean
   no collisions.

4. Determine article type of the files now in named_entity_xmls. JATS files
   (``*.xml``) carry an authoritative ``<article article-type="...">`` -- read it
   directly. GROBID TEI (``*.grobid.tei.xml``) has *no* document-level type field
   (every ``type="..."`` in it is structural: figure / bibr / main / published),
   so each GROBID file is typed by joining on its PMC ID to the matching JATS
   record (the overlapping JATS twins stay in ``high_impact_xmls/`` after the move,
   exactly because GROBID won the overlap) -- the same method as
   ``summaries/bacs/ncbi_grobid_article_types.html``. A GROBID file with no JATS
   twin is recorded as ``unknown (no JATS twin)``.

5. Classify which types carry original results. Each article type is mapped to
   one of three buckets -- original-results (primary research: own methods, data &
   findings, incl. case reports and clinical trials, the latter tagged
   research-article by JATS), secondary-synthesis (reviews/meta-analyses/meeting
   reports of prior results) and non-research (editorial, correspondence, notices,
   fragments). JATS has no dedicated clinical-trial type, so clinical trials sit
   inside research-article.

6. Isolate the original-results subset. Every file in the original-results bucket
   is moved a SECOND time, out of named_entity_xmls/ into experimental_ner/ --
   collecting the primary-research papers (research articles, brief reports, case
   reports, clinical trials, ...) for experimental NER and leaving the secondary /
   non-research files in named_entity_xmls/. GROBID files are categorised by
   PMC-ID join, JATS files by their own article-type; both schemes keep their
   filename, so no paper is duplicated. Resume-safe like the first move.

7. Summarise. ``summaries/named_entity_xmls.html`` is built from disk: what moved
   by source, the directory state before/after, the article-type distribution of
   the assembled corpus with an explicit GROBID vs non-GROBID breakdown per type,
   the original-results classification (also broken down GROBID vs non-GROBID), and
   the second move to experimental_ner/. All three prompts and this strategy are
   embedded in the page. The corpus universe spans named_entity_xmls/ AND
   experimental_ner/, so counts stay coherent after the second move. The summary is
   written in a ``finally`` block so an interrupted run still leaves an accurate page.
"""

import os
import re
import sys
import glob
import shutil
import html as _html
from collections import Counter

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
HERE          = os.path.dirname(os.path.abspath(__file__))
GROBID_DIR    = os.path.join(HERE, "grobid_xmls")        # *.grobid.tei.xml (TEI)
HIGH_IMPACT_DIR = os.path.join(HERE, "high_impact_xmls") # *.xml (JATS)
OUT_DIR       = os.path.join(HERE, "named_entity_xmls")  # destination (assembled corpus)
EXPERIMENTAL_DIR = os.path.join(HERE, "gpu_bundle", "experimental_ner")  # original-results subset (written into the GPU bundle)
SUMMARY_DIR   = os.path.join(HERE, "summaries")
SUMMARY_HTML  = os.path.join(SUMMARY_DIR, "named_entity_xml.html")

GROBID_SUFFIX = ".grobid.tei.xml"
JATS_SUFFIX   = ".xml"

DRY_RUN = os.environ.get("DRY_RUN", "0") not in ("0", "", "false", "False")

PROMPT = (
    "design a strategy and write a script that identifies unique publications in "
    "'grobid_xmls' and 'high_impact_xmls' and moves the corresponding xmls to "
    "'named_entity_xmls' so that all xmls in 'grobid_xmls' are included; summarize "
    "the moved xmls in 'summaries/named_entity_xmls.html'; determine the article "
    "type of xmls in 'named_entity_xmls'; summarize the article types in "
    "'summaries/named_entity_xmls.html'; do not execute the script; include this "
    "prompt and the strategy in 'summaries/named_entity_xmls.html'"
)

PROMPT2 = (
    "update named_entity_xml.py so that the breakdown of the grobid and non-grobid "
    "xmls is included in 'named_entity_xmls.html'; also, determine the article types "
    "that contain the original results based on methods, case reports, clinical "
    "trials etc.; include this prompt in the 'named_entity_xmls.html'; do not "
    "execute 'named_entity_xml.py'"
)

PROMPT3 = (
    "add the following capacity to 'named_entity_xmls.html' - move the xmls that "
    "fall under the 'original-results'-category in '6. Article types that contain "
    "original results' in '/summaries/named_entity_xmls.html' from "
    "'named_entity_xmls' to 'experimental_ner'; summarize this move in "
    "'/summaries/named_entity_xmls.html'"
)

# Opening <article ...> tag of a JATS document carries the authoritative type.
_ARTICLE_TYPE_RE = re.compile(r'<article\b[^>]*?\barticle-type="([^"]+)"')

# --------------------------------------------------------------------------- #
# Original-results classification of JATS article types.
#
# Which article types report ORIGINAL primary results -- i.e. a study with its
# own methods and findings (original experiments, case reports, clinical trials,
# etc.) -- versus those that synthesise prior work or carry no results at all.
# Note: JATS has no dedicated "clinical-trial" article-type; clinical trials are
# tagged "research-article", so they fall inside that bucket.
# --------------------------------------------------------------------------- #
RESULTS_ORIGINAL    = "original-results"      # primary research: own methods + findings
RESULTS_SECONDARY   = "secondary-synthesis"   # reviews / meta-analyses of prior results
RESULTS_NONRESEARCH = "non-research"          # editorial / correspondence / notices / fragments

RESULTS_CATEGORY = {
    # --- original primary research (own methods, data & findings) -------------
    "research-article":    RESULTS_ORIGINAL,   # incl. clinical trials (JATS tags them so)
    "brief-report":        RESULTS_ORIGINAL,   # short original study
    "case-report":         RESULTS_ORIGINAL,   # original clinical observation
    "rapid-communication": RESULTS_ORIGINAL,   # fast-tracked original study
    "report":              RESULTS_ORIGINAL,   # original data report
    # --- secondary synthesis of prior results --------------------------------
    "review-article":      RESULTS_SECONDARY,
    "systematic-review":   RESULTS_SECONDARY,
    "meeting-report":      RESULTS_SECONDARY,
    # --- everything else falls through to non-research (see results_category) -
}

# Human-readable description + ordering of the three buckets.
RESULTS_BUCKETS = [
    (RESULTS_ORIGINAL,
     "Original results &mdash; primary research presenting its own methods, data and "
     "findings (original experiments, case reports and clinical trials &mdash; the latter "
     "tagged <code>research-article</code> by JATS)."),
    (RESULTS_SECONDARY,
     "Secondary synthesis &mdash; reviews, systematic reviews/meta-analyses and meeting "
     "reports that analyse or summarise prior results rather than generate new primary data."),
    (RESULTS_NONRESEARCH,
     "Non-research &mdash; editorials, correspondence (letters/replies/commentary), news, "
     "corrections/retractions and other notices or fragments that carry no study results."),
]


def results_category(article_type):
    """Map a JATS article-type to one of the three results buckets."""
    return RESULTS_CATEGORY.get(article_type, RESULTS_NONRESEARCH)


# --------------------------------------------------------------------------- #
# Identity & typing helpers
# --------------------------------------------------------------------------- #
def pmcid_of(name):
    """PMC-ID stem of a file in either scheme (PMC10006201)."""
    base = os.path.basename(name)
    if base.endswith(GROBID_SUFFIX):
        return base[:-len(GROBID_SUFFIX)]
    if base.endswith(JATS_SUFFIX):
        return base[:-len(JATS_SUFFIX)]
    return base


def list_pmcids(directory, suffix):
    out = {}
    for path in glob.glob(os.path.join(directory, "*" + suffix)):
        out[pmcid_of(path)] = path
    return out


def grobid_files(directory):
    """GROBID TEI files (PMC*.grobid.tei.xml) in a directory."""
    return glob.glob(os.path.join(directory, "*" + GROBID_SUFFIX))


def jats_files(directory):
    """True JATS files (PMC*.xml) in a directory, excluding GROBID TEI.

    Note: ``*.grobid.tei.xml`` also ends in ``.xml``, so a plain ``*.xml`` glob
    would wrongly include GROBID files -- filter them out here.
    """
    return [p for p in glob.glob(os.path.join(directory, "*" + JATS_SUFFIX))
            if not p.endswith(GROBID_SUFFIX)]


def jats_article_type(path):
    """Read article-type from a JATS file's opening <article> tag, or None."""
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            head = fh.read(20000)
    except OSError:
        return None
    m = _ARTICLE_TYPE_RE.search(head)
    return m.group(1) if m else "unknown"


def build_jats_type_map():
    """PMC ID -> JATS article-type, scanning BOTH high_impact_xmls/ and
    named_entity_xmls/.

    Scanning both means a publication's type resolves whether or not its JATS file
    has already been moved: the overlapping JATS twins stay in high_impact_xmls/
    (so every GROBID file can be typed by PMC-ID join), while the JATS files that
    were moved are read from named_entity_xmls/. Robust to pre-, mid- and
    post-move states.
    """
    type_by_pmc = {}
    for directory in (HIGH_IMPACT_DIR, OUT_DIR, EXPERIMENTAL_DIR):
        for path in jats_files(directory):
            type_by_pmc.setdefault(pmcid_of(path), jats_article_type(path))
    return type_by_pmc


# --------------------------------------------------------------------------- #
# 1-3. Plan & move
# --------------------------------------------------------------------------- #
def plan_move():
    """Decide what to move. Returns (grobid_paths, jats_paths, excluded_pmcids).

    * grobid_paths : GROBID files still in grobid_xmls/ that need moving.
    * jats_paths   : JATS files whose PMC ID is NOT in the GROBID winner set.
    * excluded     : PMC IDs present in both sources (JATS twin dropped).

    The GROBID winner set is everything in grobid_xmls/ PLUS any GROBID file
    already in named_entity_xmls/ from a prior run -- so once GROBID has moved out,
    its excluded JATS twins are still recognised and never moved in as duplicates
    (the two filename schemes would not collide on their own).
    """
    grobid = list_pmcids(GROBID_DIR, GROBID_SUFFIX)
    jats   = {pmcid_of(p): p for p in jats_files(HIGH_IMPACT_DIR)}
    grobid_ids = set(grobid) | {pmcid_of(p) for p in grobid_files(OUT_DIR)}
    excluded = sorted(grobid_ids & set(jats))
    jats_to_move = {p: path for p, path in jats.items() if p not in grobid_ids}
    return grobid, jats_to_move, excluded


def safe_move(src, dst_dir):
    """Move src into dst_dir; skip (return 'present') if already there."""
    dst = os.path.join(dst_dir, os.path.basename(src))
    if os.path.exists(dst):
        # Already moved on a previous run -- resume-safe no-op.
        if os.path.exists(src) and os.path.abspath(src) != os.path.abspath(dst):
            # Stray duplicate left in the source; leave both, report it.
            return "duplicate-in-source"
        return "present"
    if not os.path.exists(src):
        return "missing"
    if DRY_RUN:
        return "planned"
    shutil.move(src, dst)
    return "moved"


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


def do_move(grobid, jats_to_move):
    if not DRY_RUN and (grobid or jats_to_move):
        # Fresh run: clear named_entity_xmls/ before re-assembling it. Guarded by
        # having files to move so a no-op re-run never wipes the corpus.
        _empty_content_dir(OUT_DIR)
    elif not DRY_RUN:
        os.makedirs(OUT_DIR, exist_ok=True)
    counts = Counter()
    for path in list(grobid.values()) + list(jats_to_move.values()):
        counts[safe_move(path, OUT_DIR)] += 1
    return counts


# --------------------------------------------------------------------------- #
# 1b. Second move: original-results files -> experimental_ner/
# --------------------------------------------------------------------------- #
def file_category(path, jats_type_map):
    """results_category for a corpus file, by its scheme.

    JATS (*.xml) are typed from their own opening <article> tag; GROBID TEI
    (*.grobid.tei.xml) carry no type, so they are typed by PMC-ID join to the
    JATS map -- the same rule used for the article-type summary.
    """
    base = os.path.basename(path)
    if base.endswith(GROBID_SUFFIX):
        atype = jats_type_map.get(pmcid_of(base), "unknown (no JATS twin)")
    else:
        atype = jats_article_type(path) or "unknown"
    return results_category(atype)


def plan_experimental_move(jats_type_map):
    """Original-results files currently in named_entity_xmls/ to move out.

    Returns {pmcid: path} for every file in OUT_DIR whose article type falls in
    the original-results bucket (section 5) -- both GROBID and JATS schemes.
    """
    to_move = {}
    for path in grobid_files(OUT_DIR) + jats_files(OUT_DIR):
        if file_category(path, jats_type_map) == RESULTS_ORIGINAL:
            to_move[os.path.basename(path)] = path
    return to_move


def do_experimental_move(to_move):
    if not DRY_RUN and to_move:
        # Clean rebuild in the GPU bundle: drop any existing experimental_ner/ in
        # gpu_bundle, then recreate it and move the newly-produced original-results
        # files in. (Guarded by `to_move` so a no-op re-run never wipes the bundle.)
        if os.path.isdir(EXPERIMENTAL_DIR):
            shutil.rmtree(EXPERIMENTAL_DIR)
        os.makedirs(EXPERIMENTAL_DIR, exist_ok=True)
    counts = Counter()
    for path in to_move.values():
        counts[safe_move(path, EXPERIMENTAL_DIR)] += 1
    return counts


# --------------------------------------------------------------------------- #
# 4. Determine article type of files in named_entity_xmls/
# --------------------------------------------------------------------------- #
def classify_destination(jats_type_map):
    """Type every file now in OUT_DIR.

    Returns (grobid_types, jats_types) as Counters. JATS files are read natively;
    GROBID files are typed by PMC-ID join to the JATS map (high_impact_xmls/).
    """
    grobid_types, jats_types = Counter(), Counter()
    for path in glob.glob(os.path.join(OUT_DIR, "*")):
        base = os.path.basename(path)
        pmc = pmcid_of(base)
        if base.endswith(GROBID_SUFFIX):
            grobid_types[jats_type_map.get(pmc, "unknown (no JATS twin)")] += 1
        elif base.endswith(JATS_SUFFIX):
            jats_types[jats_article_type(path) or "unknown"] += 1
    return grobid_types, jats_types


# --------------------------------------------------------------------------- #
# 5. Summary HTML
# --------------------------------------------------------------------------- #
CSS = """
 body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        margin: 2rem auto; max-width: 1100px; color: #222; line-height: 1.45; padding: 0 1rem; }
 h1 { margin-bottom: .25rem; } h2 { margin-top: 2.25rem; border-bottom: 1px solid #ddd; padding-bottom: .25rem; }
 h3 { margin: 1.2rem 0 .3rem; }
 .meta { color: #555; margin-bottom: 1rem; }
 .stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: .75rem; margin: 1rem 0 1.5rem; }
 .stat { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: .75rem 1rem; }
 .stat .v { font-size: 1.4rem; font-weight: 600; } .stat .k { color: #555; font-size: .85rem; }
 table { border-collapse: collapse; width: 100%; margin: .5rem 0 1rem; font-size: .9rem; }
 th, td { border: 1px solid #e1e4e8; padding: .35rem .55rem; text-align: left; vertical-align: middle; }
 th { background: #f6f8fa; }
 td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
 code { background: #f6f8fa; padding: 1px 4px; border-radius: 3px; font-size: .88em; }
 pre { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: .9rem 1rem; overflow-x: auto; font-size: .85rem; line-height: 1.4; white-space: pre-wrap; }
 .bar { display:inline-block; height:.72em; background:#3b7dd8; border-radius:2px; vertical-align: middle; }
 .bar.g { background:#2da44e; }
 .barwrap { display:flex; align-items:center; gap:.4rem; }
 .dim { color: #888; font-size: .85em; }
 .ok { color:#2da44e; font-weight:600; }
 .key  { background: #ddf4ff; border-left: 4px solid #0969da; padding: .6rem .9rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 .note { background: #fff8c5; border-left: 4px solid #d4a72c; padding: .5rem .75rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 ol.strategy > li { margin: .45rem 0; }
"""

STRATEGY_HTML = """
<ol class="strategy">
  <li><strong>Equivalence key = PMC ID.</strong> Both directories name files by PMC ID &mdash; JATS as
      <code>PMC&lt;digits&gt;.xml</code>, GROBID as <code>PMC&lt;digits&gt;.grobid.tei.xml</code>. A GROBID
      TEI carries no PMID/PMCID <em>inside</em> the document (only DOI/ISSN/ORCID), so the filename stem
      is the reliable shared identity. The two schemes are disjoint, so a JATS file and a GROBID file for
      the same paper never collide in the destination.</li>
  <li><strong>Unique publications, GROBID-preferred.</strong> The unique set is the union of the two
      sources by PMC ID. To guarantee that <em>all of</em> <code>grobid_xmls/</code> is included, GROBID
      wins every overlap: move <em>every</em> GROBID file, and move a JATS file <em>only if its PMC ID is
      absent from the GROBID set</em>. The overlapping JATS records are excluded in favour of the richer
      GROBID full text, so each publication appears exactly once.</li>
  <li><strong>Move &amp; verify.</strong> <code>shutil.move</code> each selected file into
      <code>named_entity_xmls/</code>; a file already in the destination is treated as done (resume-safe),
      and the move is verified by source/destination counts.</li>
  <li><strong>Determine article type.</strong> JATS files carry an authoritative
      <code>&lt;article article-type="&hellip;"&gt;</code> &mdash; read directly. GROBID TEI has
      <em>no</em> document-level type field (every <code>type="&hellip;"</code> in it is structural), so
      each GROBID file is typed by joining on its PMC ID to the matching JATS record &mdash; the
      overlapping JATS twins remain in <code>high_impact_xmls/</code> precisely because GROBID won the
      overlap. Same method as <code>ncbi_grobid_article_types.html</code>.</li>
  <li><strong>Classify which types carry original results.</strong> Each article type is mapped to one
      of three buckets &mdash; <em>original-results</em> (primary research with its own methods, data and
      findings: research articles, brief reports, case reports, clinical trials &mdash; clinical trials
      are tagged <code>research-article</code> by JATS, which has no dedicated type for them),
      <em>secondary-synthesis</em> (reviews / systematic reviews / meeting reports of prior results) and
      <em>non-research</em> (editorial, correspondence, news, corrections/retractions and other notices).</li>
  <li><strong>Isolate the original-results subset.</strong> Every file in the original-results bucket is
      moved a second time, out of <code>named_entity_xmls/</code> into <code>experimental_ner/</code>
      &mdash; collecting the primary-research papers for experimental NER and leaving the secondary /
      non-research files behind. GROBID files are categorised by PMC-ID join, JATS files by their own
      <code>article-type</code>; both schemes carry their filename through, so no paper is duplicated.</li>
  <li><strong>Summarise.</strong> This page is built from disk: what moved by source, the directory state
      before/after, the article-type distribution with an explicit <strong>GROBID vs non-GROBID
      breakdown</strong> per type, the original-results classification (also broken down GROBID vs
      non-GROBID), and the second move to <code>experimental_ner/</code>.</li>
</ol>
"""


def _fmt(x):
    return format(x, ",")


def _breakdown_table(grobid_types, jats_types):
    """Article-type table with an explicit GROBID vs non-GROBID (JATS) breakdown."""
    types = set(grobid_types) | set(jats_types)
    total = sum(grobid_types.values()) + sum(jats_types.values())
    rows = sorted(types, key=lambda t: (-(grobid_types.get(t, 0) + jats_types.get(t, 0)), t))
    maxc = max((grobid_types.get(t, 0) + jats_types.get(t, 0) for t in types), default=1)
    H = ["<table><thead><tr><th>article-type</th><th class='num'>GROBID</th>"
         "<th class='num'>non-GROBID<br><span class='dim'>(JATS)</span></th>"
         "<th class='num'>total</th><th class='num'>% of corpus</th><th>share</th>"
         "</tr></thead><tbody>"]
    gt = sum(grobid_types.values())
    jt = sum(jats_types.values())
    for t in rows:
        g, j = grobid_types.get(t, 0), jats_types.get(t, 0)
        tot = g + j
        w = max(1, int(round(300 * tot / maxc)))
        pct = (100 * tot / total) if total else 0
        H.append("<tr><td><code>%s</code></td><td class='num'>%s</td><td class='num'>%s</td>"
                 "<td class='num'>%s</td><td class='num'>%.1f%%</td>"
                 "<td><div class='barwrap'><span class='bar' style='width:%dpx'></span></div></td></tr>"
                 % (_html.escape(t), _fmt(g), _fmt(j), _fmt(tot), pct, w))
    H.append("<tfoot><tr><th>Total</th><th class='num'>%s</th><th class='num'>%s</th>"
             "<th class='num'>%s</th><th class='num'>100%%</th><th></th></tr></tfoot></table>"
             % (_fmt(gt), _fmt(jt), _fmt(total)))
    return "".join(H)


def _results_breakdown_table(grobid_types, jats_types):
    """Original-results buckets with GROBID vs non-GROBID breakdown + member types."""
    total = sum(grobid_types.values()) + sum(jats_types.values())
    g_by_bucket, j_by_bucket, members = Counter(), Counter(), {}
    for t, c in grobid_types.items():
        b = results_category(t)
        g_by_bucket[b] += c
        members.setdefault(b, set()).add(t)
    for t, c in jats_types.items():
        b = results_category(t)
        j_by_bucket[b] += c
        members.setdefault(b, set()).add(t)
    maxc = max(((g_by_bucket.get(b, 0) + j_by_bucket.get(b, 0)) for b, _ in RESULTS_BUCKETS), default=1)
    H = ["<table><thead><tr><th>results category</th><th>member article-types present</th>"
         "<th class='num'>GROBID</th><th class='num'>non-GROBID<br><span class='dim'>(JATS)</span></th>"
         "<th class='num'>total</th><th class='num'>% of corpus</th><th>share</th>"
         "</tr></thead><tbody>"]
    for b, _desc in RESULTS_BUCKETS:
        g, j = g_by_bucket.get(b, 0), j_by_bucket.get(b, 0)
        tot = g + j
        w = max(1, int(round(300 * tot / maxc))) if tot else 1
        pct = (100 * tot / total) if total else 0
        mem = ", ".join("<code>%s</code>" % _html.escape(m)
                        for m in sorted(members.get(b, ()))) or "&mdash;"
        cls = "bar g" if b == RESULTS_ORIGINAL else "bar"
        H.append("<tr><td><strong>%s</strong></td><td class='dim'>%s</td><td class='num'>%s</td>"
                 "<td class='num'>%s</td><td class='num'>%s</td><td class='num'>%.1f%%</td>"
                 "<td><div class='barwrap'><span class='%s' style='width:%dpx'></span></div></td></tr>"
                 % (_html.escape(b), mem, _fmt(g), _fmt(j), _fmt(tot), pct, cls, w))
    gt = sum(g_by_bucket.values())
    jt = sum(j_by_bucket.values())
    H.append("<tfoot><tr><th>Total</th><th></th><th class='num'>%s</th><th class='num'>%s</th>"
             "<th class='num'>%s</th><th class='num'>100%%</th><th></th></tr></tfoot></table>"
             % (_fmt(gt), _fmt(jt), _fmt(total)))
    return "".join(H)


def build_summary(jats_type_map):
    """Build summaries/named_entity_xmls.html from the state on disk."""
    # PMC IDs of each scheme in each location. A publication is counted wherever
    # it currently lives -- named_entity_xmls/ OR experimental_ner/ (the second
    # move) -- so the corpus stays coherent pre-, mid- and post-move.
    grobid_src_ids = {pmcid_of(p) for p in grobid_files(GROBID_DIR)}
    jats_src_ids   = {pmcid_of(p) for p in jats_files(HIGH_IMPACT_DIR)}
    dst_grobid_ids = {pmcid_of(p) for p in grobid_files(OUT_DIR)}
    dst_jats_ids   = {pmcid_of(p) for p in jats_files(OUT_DIR)}
    exp_grobid_ids = {pmcid_of(p) for p in grobid_files(EXPERIMENTAL_DIR)}
    exp_jats_ids   = {pmcid_of(p) for p in jats_files(EXPERIMENTAL_DIR)}

    # GROBID is the winner set (all of grobid_xmls, wherever it now lives); a JATS
    # publication belongs in the corpus only if its PMC ID is not a GROBID winner.
    all_grobid_ids = grobid_src_ids | dst_grobid_ids | exp_grobid_ids
    all_jats_ids   = (jats_src_ids | dst_jats_ids | exp_jats_ids) - all_grobid_ids
    moved_as_grobid_ids = all_grobid_ids
    moved_as_jats_ids   = all_jats_ids
    excluded_ids        = (jats_src_ids | dst_jats_ids | exp_jats_ids) & all_grobid_ids
    union_ids           = all_grobid_ids | all_jats_ids

    n_grobid   = len(moved_as_grobid_ids)
    n_jats     = len(moved_as_jats_ids)
    n_total    = n_grobid + n_jats
    n_excluded = len(excluded_ids)
    in_dest_total = len(dst_grobid_ids) + len(dst_jats_ids) + len(exp_grobid_ids) + len(exp_jats_ids)

    # Article types of the assembled corpus.
    grobid_types = Counter(jats_type_map.get(p, "unknown (no JATS twin)") for p in moved_as_grobid_ids)
    jats_types   = Counter(jats_type_map.get(p, "unknown") for p in moved_as_jats_ids)
    total_types  = Counter(); total_types.update(grobid_types); total_types.update(jats_types)
    n_original   = sum(c for t, c in total_types.items() if results_category(t) == RESULTS_ORIGINAL)

    progress = ("complete" if in_dest_total >= n_total and n_total
                else ("not started (projected from sources)" if in_dest_total == 0 else "partial"))

    # Second move: the original-results subset -> experimental_ner/. The breakdown
    # is the original-results slice of the per-type counters (so it stays correct
    # whether those files are still in named_entity_xmls/ or already in exp dir).
    exp_grobid_types = Counter({t: c for t, c in grobid_types.items()
                                if results_category(t) == RESULTS_ORIGINAL})
    exp_jats_types   = Counter({t: c for t, c in jats_types.items()
                                if results_category(t) == RESULTS_ORIGINAL})
    n_exp_grobid = sum(exp_grobid_types.values())
    n_exp_jats   = sum(exp_jats_types.values())
    n_exp_total  = n_exp_grobid + n_exp_jats               # == n_original
    n_in_exp     = len(exp_grobid_ids) + len(exp_jats_ids)  # already in experimental_ner/
    exp_progress = ("complete" if n_in_exp >= n_exp_total and n_exp_total
                    else ("not started (projected)" if n_in_exp == 0 else "partial"))

    H = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<title>Named-entity XML corpus &mdash; de-duplicated &amp; typed</title>",
         "<style>" + CSS + "</style></head><body>"]
    H.append("<h1>XMLs moved to <code>named_entity_xmls/</code> &mdash; "
             "de-duplicated across GROBID &amp; high-impact JATS</h1>")
    H.append("<p class='meta'>Generated by <code>named_entity_xml.py</code> (after "
             "<code>summaries/bacs/xmls_for_ner.html</code>) &middot; identity: <strong>PMC ID</strong> "
             "&middot; overlap resolved in favour of <strong>GROBID full text</strong> &middot; move "
             "state: <strong>%s</strong></p>" % _html.escape(progress))

    cards = [(_fmt(n_total), "files in named_entity_xmls/"),
             (_fmt(n_grobid), "GROBID (all of grobid_xmls)"),
             (_fmt(n_jats), "non-GROBID (high-impact JATS)"),
             (_fmt(n_excluded), "JATS excluded as duplicates"),
             (_fmt(len(total_types)), "distinct article types"),
             (_fmt(n_original), "original-results files (%.1f%%)"
              % (100 * n_original / n_total if n_total else 0))]
    H.append("<div class='stat-grid'>")
    for v, k in cards:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")

    H.append("<div class='key'>Every XML in <code>grobid_xmls/</code> (<strong>%s</strong>) is included. "
             "The unique-publication set is the union of the two sources by <strong>PMC ID</strong> "
             "(<strong>%s</strong> publications). On the <strong>%s</strong> papers present in both, the "
             "GROBID full-text version is kept and the JATS twin excluded, so each paper appears once: "
             "<strong>%s</strong> GROBID + <strong>%s</strong> JATS-only = <strong>%s</strong> files in "
             "<code>named_entity_xmls/</code>.</div>"
             % (_fmt(len(grobid_src_ids)), _fmt(len(union_ids)), _fmt(n_excluded),
                _fmt(n_grobid), _fmt(n_jats), _fmt(n_total)))

    H.append("<h2>1. Strategy</h2>" + STRATEGY_HTML)

    # 2. What moved, by source
    H.append("<h2>2. What moved, by source</h2><table>")
    H.append("<thead><tr><th>source</th><th>scheme</th><th class='num'>publications</th>"
             "<th>role</th></tr></thead><tbody>")
    H.append("<tr><td><code>grobid_xmls/</code></td><td><code>*.grobid.tei.xml</code></td>"
             "<td class='num'>%s</td><td>all moved &mdash; GROBID wins every overlap</td></tr>"
             % _fmt(n_grobid))
    H.append("<tr><td><code>high_impact_xmls/</code></td><td><code>*.xml</code></td>"
             "<td class='num'>%s</td><td>moved only where PMC ID absent from GROBID</td></tr>"
             % _fmt(n_jats))
    H.append("<tr><td><code>high_impact_xmls/</code> (overlap)</td><td><code>*.xml</code></td>"
             "<td class='num'>%s</td><td>excluded &mdash; GROBID twin taken instead</td></tr>"
             % _fmt(n_excluded))
    H.append("<tfoot><tr><th>moved total</th><th></th><th class='num'>%s</th><th></th></tr></tfoot></table>"
             % _fmt(n_total))

    # 4. Directory state before/after. Reconstructed from the publication universe
    # (source + destination), so the net before->after is correct whether this page
    # is built before, during, or after the move.
    before_grobid = len(all_grobid_ids)                         # every GROBID file
    after_grobid  = before_grobid - n_grobid                    # -> 0 (all move out)
    before_jats   = len(all_jats_ids)                           # every JATS file
    after_jats    = before_jats - n_jats                        # excluded twins stay
    H.append("<h2>3. Directory state (before &rarr; after the move)</h2><table>")
    H.append("<thead><tr><th>Directory</th><th class='num'>before</th><th class='num'>&Delta;</th>"
             "<th class='num'>after</th></tr></thead><tbody>")
    H.append("<tr><td><code>grobid_xmls/</code></td><td class='num'>%s</td><td class='num'>&minus;%s</td>"
             "<td class='num'>%s</td></tr>"
             % (_fmt(before_grobid), _fmt(n_grobid), _fmt(after_grobid)))
    H.append("<tr><td><code>high_impact_xmls/</code></td><td class='num'>%s</td>"
             "<td class='num'>&minus;%s</td><td class='num'>%s</td></tr>"
             % (_fmt(before_jats), _fmt(n_jats), _fmt(after_jats)))
    H.append("<tr><td><code>named_entity_xmls/</code> (new)</td><td class='num'>0</td>"
             "<td class='num'>+%s</td><td class='num'>%s</td></tr>"
             % (_fmt(n_total), _fmt(n_total)))
    H.append("</tbody></table>")
    H.append("<p class='dim'>After the move <code>high_impact_xmls/</code> retains the <strong>%s</strong> "
             "overlap JATS twins (deliberately excluded, and the source for GROBID article-type lookup). "
             "<code>grobid_xmls/</code> is emptied &mdash; all of it is in the corpus.</p>" % _fmt(n_excluded))

    # 4. Article types -- GROBID vs non-GROBID breakdown
    H.append("<h2>4. Article types of <code>named_entity_xmls/</code> "
             "&mdash; GROBID vs non-GROBID</h2>")
    H.append("<div class='note'><strong>GROBID TEI does not classify article type.</strong> A TEI document "
             "has no document-level type field, so the <strong>%s</strong> GROBID files were typed by "
             "joining on PMC ID to their high-impact JATS record (<code>high_impact_xmls/</code>). The "
             "<strong>%s</strong> non-GROBID files are JATS and carry an authoritative "
             "<code>article-type</code> read directly.</div>" % (_fmt(n_grobid), _fmt(n_jats)))
    H.append("<p class='dim'>Each row gives the GROBID (<code>*.grobid.tei.xml</code>) and non-GROBID "
             "(JATS <code>*.xml</code>) file counts for that type, their total and its share of the "
             "%s-file corpus.</p>" % _fmt(n_total))
    H.append(_breakdown_table(grobid_types, jats_types))

    # 5. Article types that contain original results
    H.append("<h2>5. Article types that contain original results</h2>")
    H.append("<p>Which of the article types above report <strong>original results</strong> &mdash; a study "
             "presenting its own methods, data and findings (original experiments, case reports, clinical "
             "trials, &hellip;) &mdash; versus those that synthesise prior work or carry no results. JATS "
             "has no dedicated <code>clinical-trial</code> type, so clinical trials are tagged "
             "<code>research-article</code> and counted in the original-results bucket.</p>")
    H.append(_results_breakdown_table(grobid_types, jats_types))
    H.append("<div class='key'><strong>%s of %s</strong> files (<strong>%.1f%%</strong>) in "
             "<code>named_entity_xmls/</code> are article types that contain original results "
             "(primary research). The remainder are secondary synthesis (reviews) or non-research "
             "(editorial / correspondence / notices).</div>"
             % (_fmt(n_original), _fmt(n_total), (100 * n_original / n_total if n_total else 0)))

    # 6. Move of the original-results subset -> experimental_ner/
    H.append("<h2>6. Move of original-results files to <code>experimental_ner/</code></h2>")
    H.append("<p>The <strong>%s</strong> original-results files from section 5 (every file whose article "
             "type is in the original-results bucket) are moved out of <code>named_entity_xmls/</code> into "
             "<code>experimental_ner/</code> &mdash; isolating the primary-research papers (their own "
             "methods, data and findings: research articles, brief reports, case reports, clinical trials, "
             "&hellip;) for experimental NER. Move state: <strong>%s</strong> "
             "(<strong>%s</strong> already in <code>experimental_ner/</code>).</p>"
             % (_fmt(n_exp_total), _html.escape(exp_progress), _fmt(n_in_exp)))
    H.append("<div class='stat-grid'>")
    for v, k in [(_fmt(n_exp_total), "original-results files moved"),
                 (_fmt(n_exp_grobid), "GROBID (*.grobid.tei.xml)"),
                 (_fmt(n_exp_jats), "non-GROBID JATS (*.xml)"),
                 (_fmt(n_total - n_exp_total), "left in named_entity_xmls/")]:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")

    H.append("<h3>6a. What moves to <code>experimental_ner/</code>, by type &amp; scheme</h3>")
    H.append(_breakdown_table(exp_grobid_types, exp_jats_types))

    H.append("<h3>6b. Directory state (before &rarr; after this move)</h3><table>")
    H.append("<thead><tr><th>Directory</th><th class='num'>before</th><th class='num'>&Delta;</th>"
             "<th class='num'>after</th></tr></thead><tbody>")
    H.append("<tr><td><code>named_entity_xmls/</code></td><td class='num'>%s</td>"
             "<td class='num'>&minus;%s</td><td class='num'>%s</td></tr>"
             % (_fmt(n_total), _fmt(n_exp_total), _fmt(n_total - n_exp_total)))
    H.append("<tr><td><code>experimental_ner/</code> (new)</td><td class='num'>0</td>"
             "<td class='num'>+%s</td><td class='num'>%s</td></tr>"
             % (_fmt(n_exp_total), _fmt(n_exp_total)))
    H.append("</tbody></table>")
    H.append("<p class='dim'>After this move <code>named_entity_xmls/</code> retains the <strong>%s</strong> "
             "secondary-synthesis and non-research files (reviews, editorials, correspondence, notices), "
             "while <code>experimental_ner/</code> holds only original-results papers. Each file keeps its "
             "scheme (<code>*.xml</code> / <code>*.grobid.tei.xml</code>), so no PMC ID is duplicated.</p>"
             % _fmt(n_total - n_exp_total))

    # 7. Notes
    H.append("<h2>7. Notes &amp; caveats</h2><ul>")
    H.append("<li><strong>One representation per paper.</strong> No PMC ID appears twice in "
             "<code>named_entity_xmls/</code>; overlapping papers are present only as GROBID TEI.</li>")
    H.append("<li><strong>Why GROBID wins the overlap.</strong> The overlapping JATS records are the "
             "high-impact eFetch XML; the GROBID export recovers the full body from the PDF, so it is the "
             "richer text for NER &mdash; and using it guarantees all of <code>grobid_xmls/</code> is "
             "included, as required.</li>")
    H.append("<li><strong>Mixed schemas downstream.</strong> The corpus holds NLM JATS "
             "(<code>&lt;sec&gt;</code>/<code>&lt;p&gt;</code>) and TEI "
             "(<code>&lt;div&gt;</code>/<code>&lt;head&gt;</code>/<code>&lt;p&gt;</code>) files.</li>")
    H.append("<li><strong>Original-results classification is a documented type-to-bucket mapping</strong> "
             "(<code>RESULTS_CATEGORY</code> in <code>named_entity_xml.py</code>), applied to the JATS "
             "<code>article-type</code> &mdash; not a read of each paper's body. It reflects what a type "
             "<em>conventionally</em> contains; edit the mapping to retune. Clinical trials are not a "
             "separate JATS type and so are inside <code>research-article</code>.</li>")
    H.append("<li><strong>Two moves, both <code>shutil.move</code>.</strong> First "
             "<code>grobid_xmls/</code> + <code>high_impact_xmls/</code> &rarr; "
             "<code>named_entity_xmls/</code>; then the original-results subset "
             "<code>named_entity_xmls/</code> &rarr; <code>experimental_ner/</code>. Originals are removed "
             "from their source dirs (set <code>DRY_RUN=1</code> to plan only). Both are resume-safe &mdash; "
             "a file already at its destination is left untouched &mdash; and this page is rewritten from "
             "disk on every run, so the corpus is reported wherever its files currently live.</li>")
    H.append("</ul>")
    H.append("<p class='dim'>Related: <code>summaries/bacs/xmls_for_ner.html</code>, "
             "<code>summaries/bacs/ncbi_grobid_article_types.html</code>, "
             "<code>summaries/grobid_xml_summary.html</code>.</p>")
    H.append("</body></html>")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(SUMMARY_HTML, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    return n_total


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    # Ensure the input corpora exist (empty is fine -- nothing to move).
    os.makedirs(GROBID_DIR, exist_ok=True)
    os.makedirs(HIGH_IMPACT_DIR, exist_ok=True)
    # Capture JATS types up front so GROBID files can be typed even post-move.
    jats_type_map = build_jats_type_map()
    grobid, jats_to_move, excluded = plan_move()

    print("grobid files to move (remaining): %d" % len(grobid))
    print("jats files to move (non-overlap): %d" % len(jats_to_move))
    print("jats excluded (overlap):          %d" % len(excluded))
    print("total to move this run:           %d" % (len(grobid) + len(jats_to_move)))
    if DRY_RUN:
        print("DRY_RUN=1 -> not moving; writing projected summary only.")

    try:
        counts = do_move(grobid, jats_to_move)
        print("move 1 outcome (-> named_entity_xmls):", dict(counts))

        # Second move: original-results subset -> experimental_ner/ (after move 1,
        # so the files are in named_entity_xmls/ to be selected from).
        exp = plan_experimental_move(jats_type_map)
        print("original-results files to move (-> experimental_ner): %d" % len(exp))
        exp_counts = do_experimental_move(exp)
        print("move 2 outcome (-> experimental_ner):", dict(exp_counts))
    finally:
        n = build_summary(jats_type_map)
        print("summary written: %s (%d files described)" % (SUMMARY_HTML, n))


if __name__ == "__main__":
    sys.exit(main())
