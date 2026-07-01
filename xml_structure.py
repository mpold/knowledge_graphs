#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xml_structure.py
================
Analyse the JATS XMLs in ``high_impact_xmls/`` for the sections typical of a
scientific publication -- title, abstract, introduction, methods, results,
discussion, conclusion/summary, references -- and report how completely each
article is structured. Writes ``summaries/ncbi_xml.html``.

Run it:  python xml_structure.py

WHY REGEX, NOT AN XML PARSER
----------------------------
~80% of these PMC records carry named character entities (e.g. ``&alpha;``) whose
definitions live in the external JATS DTD. A strict parser (xml.etree) raises
"undefined entity" on those files, and many scanned records are not well-formed at
all. So the analysis is done with tolerant, dependency-free string/regex scanning:
presence of a tag and the text of each ``<sec>``'s ``<title>`` is all we need, and
that survives entities, odd namespaces and OCR noise.

STRATEGY
--------
1. Per file, split the document into its JATS zones: ``<front>`` (metadata +
   ``<abstract>``), ``<body>`` (the narrative), ``<back>`` (references, funding,
   acknowledgements). Section harvesting is confined to ``<body>`` so abstract
   sub-sections don't inflate the counts.

2. Body shape. Classify each article's body as
     * structured     -- contains ``<sec>`` elements,
     * OCR scanned     -- only a ``<preformat pmc-ocr-text>`` block (no sections),
     * unstructured    -- a body with neither, or
     * front-only      -- no ``<body>`` at all (publisher front matter / abstract).

3. Section detection. Every ``<sec>`` inside the body is mapped to a canonical
   IMRaD bucket using BOTH its ``sec-type`` attribute (when present -- only ~a
   third carry one) and the text of its first ``<title>``:
       Introduction  <- intro / introduction / background
       Methods       <- methods / materials and methods / methodology /
                        experimental / statistical analysis / patients & methods
       Results       <- results / findings        (handles "Results and Discussion")
       Discussion    <- discussion
       Conclusion    <- conclusion(s) / summary / concluding remarks
   Document-level parts are detected directly: Title (``<article-title>``),
   Abstract (``<abstract>``), References (``<ref-list>`` / ``<ref>`` in ``<back>``).

4. Aggregate. Count, per section, how many articles contain it; tabulate body
   shapes; measure IMRaD completeness (Intro + Methods + Results + Discussion);
   and tally common ancillary parts (figures, tables, supplementary material,
   data-availability, funding, acknowledgements, COI, author contributions).

5. No-abstract subset. Pull out the articles that lack an <abstract> and analyse
   them separately (body presence, body shape, section coverage) -- a missing
   abstract often flags scanned/front-only records or editorials/letters.

6. Summarise in summaries/ncbi_xml.html (same look as brain_cancer_hits.html),
   embedding the strategy.
"""

import os
import re
import sys
import glob
import html as _html

BASE        = os.path.dirname(os.path.abspath(__file__))
XML_DIR     = os.path.join(BASE, "high_impact_xmls")
SUMMARY_DIR = os.path.join(BASE, "summaries")
SUMMARY_HTML= os.path.join(SUMMARY_DIR, "xml_structure.html")

# Canonical IMRaD-ish buckets reported in the coverage table, in reading order.
CORE_SECTIONS = ["Introduction", "Methods", "Results", "Discussion", "Conclusion"]

# Pre-compiled scanners --------------------------------------------------------
RE_BODY     = re.compile(r"<body\b.*?</body>", re.S)
RE_BACK     = re.compile(r"<back\b.*?</back>", re.S)
RE_ARTTITLE = re.compile(r"<article-title\b[^>]*>(.*?)</article-title>", re.S)
RE_ABSTRACT = re.compile(r"<abstract\b")
RE_PREFORMAT= re.compile(r"<preformat\b")
RE_REFLIST  = re.compile(r"<ref-list\b")
RE_REF      = re.compile(r"<ref\b")
RE_FIG      = re.compile(r"<fig\b")
RE_TABLEWRAP= re.compile(r"<table-wrap\b")
RE_SECSTART = re.compile(r"<sec\b")
RE_SECTYPE  = re.compile(r'sec-type="([^"]*)"')
RE_TITLE    = re.compile(r"<title\b[^>]*>(.*?)</title>", re.S)
RE_TAGS     = re.compile(r"<[^>]+>")
RE_WS       = re.compile(r"\s+")


def clean_text(s):
    s = RE_TAGS.sub(" ", s or "")
    s = _html.unescape(s)
    return RE_WS.sub(" ", s).strip().lower()


def classify_section(title_l, sectype):
    """Return the set of canonical buckets a <sec> belongs to (may be 0, 1 or 2)."""
    types = set(re.split(r"[^a-z]+", (sectype or "").lower()))
    b = set()

    def has(*kw):
        return any(k in title_l for k in kw)

    if (types & {"intro", "introduction", "background"}) or has("introduction") \
            or title_l == "background" or title_l.startswith("background"):
        b.add("Introduction")
    if (types & {"methods", "methodology", "materials"}) \
            or has("method", "materials and", "material and", "methodology",
                   "experimental", "statistical analys", "statistical method",
                   "patients and", "subjects and", "data collection", "study design"):
        b.add("Methods")
    if "results" in types or has("results", "findings"):
        b.add("Results")
    if "discussion" in types or has("discussion"):
        b.add("Discussion")
    if (types & {"conclusion", "conclusions"}) or has("conclusion", "concluding remarks") \
            or title_l == "summary" or (title_l.endswith(" summary") and "reporting" not in title_l):
        b.add("Conclusion")
    return b


def analyse_file(text):
    """Return a dict of structural facts for one article's XML text."""
    bm = RE_BODY.search(text)
    body = bm.group(0) if bm else ""
    km = RE_BACK.search(text)
    back = km.group(0) if km else ""

    # --- document-level parts ---
    tm = RE_ARTTITLE.search(text)
    has_title = bool(tm and clean_text(tm.group(1)))
    has_abstract = bool(RE_ABSTRACT.search(text))
    has_reflist = bool(RE_REFLIST.search(back) or RE_REFLIST.search(text))
    n_refs = len(RE_REF.findall(back)) if back else len(RE_REF.findall(text))

    # --- body shape ---
    has_body = bool(body)
    has_sec = bool(RE_SECSTART.search(body))
    has_ocr = bool(RE_PREFORMAT.search(body))
    if not has_body:
        shape = "front-only"
    elif has_sec:
        shape = "structured"
    elif has_ocr:
        shape = "OCR scanned"
    else:
        shape = "unstructured"

    # --- section harvesting (within body only) ---
    found = set()
    sec_titles = []
    if has_sec:
        starts = [m.start() for m in RE_SECSTART.finditer(body)]
        bounds = starts[1:] + [len(body)]
        for s, e in zip(starts, bounds):
            window = body[s:e]
            stype_m = RE_SECTYPE.search(body[s:s + 200])
            stype = stype_m.group(1) if stype_m else ""
            ttl_m = RE_TITLE.search(window)
            ttl = clean_text(ttl_m.group(1)) if ttl_m else ""
            if ttl:
                sec_titles.append(ttl)
            found |= classify_section(ttl, stype)

    # --- ancillary parts (presence) ---
    body_titles = " | ".join(sec_titles)
    back_l = back.lower()
    anc = {
        "Figures":            bool(RE_FIG.search(text)),
        "Tables":             bool(RE_TABLEWRAP.search(text)),
        "Supplementary":      ("supplementary" in body_titles) or ("supplementary-material" in text.lower()),
        "Data availability":  "data availab" in body_titles or "data availab" in back_l,
        "Funding":            "funding" in body_titles or "funding" in back_l,
        "Acknowledgements":   "acknowledg" in body_titles or "acknowledg" in back_l,
        "Conflict of interest": any(k in (body_titles + " " + back_l)
                                    for k in ("conflict of interest", "competing interest",
                                              "coi-statement", "competing financial")),
        "Author contributions": "author contribution" in (body_titles + " " + back_l),
    }

    return {
        "has_title": has_title, "has_abstract": has_abstract, "has_reflist": has_reflist,
        "has_body": has_body, "n_refs": n_refs, "shape": shape, "sections": found,
        "n_sec_titles": len(sec_titles), "anc": anc,
    }


# --------------------------------------------------------------------------- #
# Aggregate + HTML
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
  <li><strong>Tolerant scan, not a strict parser.</strong> ~80% of these records carry named entities
      (e.g. <code>&amp;alpha;</code>) defined only in the external JATS DTD, and many scanned records are
      not well-formed, so a strict XML parser fails. Each file is scanned as text: tag presence plus the
      text of each <code>&lt;sec&gt;</code>'s <code>&lt;title&gt;</code> is all the structure we need.</li>
  <li><strong>Zone the document.</strong> Split into <code>&lt;front&gt;</code> (metadata +
      <code>&lt;abstract&gt;</code>), <code>&lt;body&gt;</code> (narrative) and <code>&lt;back&gt;</code>
      (references, funding, acknowledgements). Sections are harvested only from <code>&lt;body&gt;</code>
      so abstract sub-sections do not inflate the counts.</li>
  <li><strong>Determine <code>&lt;body&gt;</code> presence/absence.</strong> A record either has a
      <code>&lt;body&gt;</code> (the full-text narrative) or it is front matter only &mdash; abstract +
      metadata with no <code>&lt;body&gt;</code>, which is what publisher-restricted and many older
      records return. This present/absent split is reported first, then refined into body
      <em>shape</em>: <em>structured</em> (has <code>&lt;sec&gt;</code>), <em>OCR scanned</em> (only a
      <code>&lt;preformat&gt;</code> OCR block), or <em>unstructured</em> (a body with neither).</li>
  <li><strong>Map each section to an IMRaD bucket</strong> using its <code>sec-type</code> (present on
      only ~1/3 of sections) <em>and</em> its title text: Introduction&larr;intro/background;
      Methods&larr;methods/materials/methodology/experimental/statistical analysis; Results&larr;results/
      findings; Discussion&larr;discussion; Conclusion&larr;conclusion(s)/summary. Title, Abstract and
      References are detected from <code>&lt;article-title&gt;</code>, <code>&lt;abstract&gt;</code> and
      <code>&lt;ref-list&gt;</code>/<code>&lt;ref&gt;</code> directly.</li>
  <li><strong>Aggregate</strong> per-section coverage, body shapes, IMRaD completeness (Intro+Methods+
      Results+Discussion) and ancillary parts (figures, tables, supplementary, data-availability,
      funding, acknowledgements, COI, author contributions), then render this page.</li>
  <li><strong>Separate the no-abstract subset.</strong> Articles with no <code>&lt;abstract&gt;</code>
      are pulled out and analysed on their own &mdash; <code>&lt;body&gt;</code> presence, body shape and
      section coverage &mdash; to see whether a missing abstract coincides with otherwise sparse
      structure (scanned / front-only records) or is just an editorial/letter that still carries full
      text.</li>
</ol>
"""


def _fmt(x):
    return format(x, ",")


def bar(v, vmax, px=260, cls="bar"):
    w = max(2, round((v / vmax) * px)) if vmax else 2
    return "<span class='%s' style='width:%dpx'></span>" % (cls, w)


def build_html(stats):
    n = stats["n"]
    sec_cov = stats["sec_cov"]          # canonical part -> count
    shapes = stats["shapes"]            # shape -> count
    anc = stats["anc"]                  # ancillary -> count
    imrad_full = stats["imrad_full"]
    body_present = stats["body_present"]
    body_absent = stats["body_absent"]
    abstract_n = stats["abstract_n"]
    title_n = stats["title_n"]
    reflist_n = stats["reflist_n"]
    refs_total = stats["refs_total"]
    refs_articles = stats["refs_articles"]
    secn_hist = stats["secn_hist"]

    pct = lambda c: (c / n * 100) if n else 0

    H = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<title>High-impact PMC XMLs &mdash; section structure</title>",
         "<style>" + CSS + "</style></head><body>"]
    H.append("<h1>High-impact PMC XMLs &mdash; section structure</h1>")
    H.append("<p class='meta'>Generated by <code>xml_structure.py</code> &middot; analysed "
             "<strong>%s</strong> XML files in <code>high_impact_xmls/</code> for the sections typical of "
             "a scientific publication.</p>" % _fmt(n))

    structured = shapes.get("structured", 0)
    cards = [
        (_fmt(n), "XML files analysed"),
        (_fmt(body_present), "&lt;body&gt; present (%.1f%%)" % pct(body_present)),
        (_fmt(body_absent), "&lt;body&gt; absent (%.1f%%)" % pct(body_absent)),
        (_fmt(structured), "structured body (%.1f%%)" % pct(structured)),
        (_fmt(reflist_n), "with reference list (%.1f%%)" % pct(reflist_n)),
        (_fmt(imrad_full), "full IMRaD (%.1f%%)" % pct(imrad_full)),
    ]
    H.append("<div class='stat-grid'>")
    for v, k in cards:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")

    H.append("<div class='key'>Of <strong>%s</strong> articles, <strong>%s (%.1f%%)</strong> contain a "
             "<code>&lt;body&gt;</code> (full-text narrative) and <strong>%s (%.1f%%)</strong> do not "
             "(front matter / abstract only). Among those with a body, <strong>%s (%.1f%%)</strong> are "
             "structured with <code>&lt;sec&gt;</code> sections and <strong>%s</strong> are scanned OCR "
             "text without sections. <strong>%s (%.1f%%)</strong> contain all four IMRaD sections "
             "(Introduction, Methods, Results, Discussion).</div>"
             % (_fmt(n), _fmt(body_present), pct(body_present), _fmt(body_absent), pct(body_absent),
                _fmt(structured), pct(structured), _fmt(shapes.get("OCR scanned", 0)),
                _fmt(imrad_full), pct(imrad_full)))

    H.append("<h2>1. Strategy</h2>" + STRATEGY_HTML)

    # --- 3. Section coverage ---
    H.append("<h2>2. Section coverage</h2>")
    H.append("<table><thead><tr><th>section</th><th class='num'>articles</th><th>coverage</th>"
             "<th class='num'>%</th></tr></thead><tbody>")
    rows = [("Title", title_n), ("Abstract", abstract_n)]
    rows += [(s, sec_cov.get(s, 0)) for s in CORE_SECTIONS]
    rows += [("References", reflist_n)]
    vmax = max(c for _, c in rows) or 1
    for name, c in rows:
        H.append("<tr><td>%s</td><td class='num'>%s</td><td>%s</td><td class='num dim'>%.1f%%</td></tr>"
                 % (name, _fmt(c), bar(c, vmax), pct(c)))
    H.append("</tbody></table>")
    H.append("<p class='dim'>Title/Abstract/References are document-level parts; Introduction&hellip;"
             "Conclusion are counted from <code>&lt;body&gt;</code> sections only.</p>")

    # --- 4. <body> presence / absence ---
    H.append("<h2>3. &lt;body&gt; section: presence / absence</h2>")
    H.append("<table><thead><tr><th>&lt;body&gt; section</th><th class='num'>articles</th><th>dist.</th>"
             "<th class='num'>%</th></tr></thead><tbody>")
    bmax = max(body_present, body_absent) or 1
    H.append("<tr><td>present <span class='dim'>&mdash; has full-text narrative</span></td>"
             "<td class='num'>%s</td><td>%s</td><td class='num dim'>%.1f%%</td></tr>"
             % (_fmt(body_present), bar(body_present, bmax, 260), pct(body_present)))
    H.append("<tr><td>absent <span class='dim'>&mdash; front matter / abstract only</span></td>"
             "<td class='num'>%s</td><td>%s</td><td class='num dim'>%.1f%%</td></tr>"
             % (_fmt(body_absent), bar(body_absent, bmax, 260, "bar o"), pct(body_absent)))
    H.append("</tbody></table>")
    H.append("<p class='dim'>Detected as the presence of a <code>&lt;body&gt;</code> element. Absence "
             "means PMC returned only <code>&lt;front&gt;</code> (publisher-restricted full text or an "
             "abstract-only deposit). The next table refines the present cases by shape.</p>")
    H.append("<div class='key'>The <strong>%s</strong> articles <em>without</em> a "
             "<code>&lt;body&gt;</code> have no full-text XML to mine. <code>ncbi_pdf.py</code> (a "
             "replica of <code>bacs/doi_openalex.py</code>: PMC proof-of-work &rarr; OpenAlex &rarr; "
             "CrossRef, with auto-recovery) targets exactly this set, exporting their PDFs into "
             "<code>ncbi_pdfs_grobid/</code> for downstream GROBID conversion to XML.</div>"
             % _fmt(body_absent))

    # --- 5. Body shape ---
    H.append("<h2>4. Body shape (articles with a &lt;body&gt;)</h2>")
    H.append("<table><thead><tr><th>shape</th><th class='num'>articles</th><th>dist.</th>"
             "<th class='num'>%</th></tr></thead><tbody>")
    present_shapes = {k: shapes.get(k, 0) for k in ("structured", "OCR scanned", "unstructured")}
    smax = max(present_shapes.values()) or 1
    explain = {"structured": "has &lt;sec&gt; sections",
               "OCR scanned": "only a &lt;preformat&gt; OCR block",
               "unstructured": "body present but neither sections nor OCR block"}
    for name in ("structured", "OCR scanned", "unstructured"):
        c = present_shapes[name]
        share = (c / body_present * 100) if body_present else 0
        H.append("<tr><td>%s <span class='dim'>&mdash; %s</span></td><td class='num'>%s</td>"
                 "<td>%s</td><td class='num dim'>%.1f%%</td></tr>"
                 % (name, explain[name], _fmt(c), bar(c, smax, 220, "bar o"), share))
    H.append("</tbody></table>")
    H.append("<p class='dim'>Shares are of the <strong>%s</strong> articles that have a "
             "<code>&lt;body&gt;</code>.</p>" % _fmt(body_present))

    # --- 6. No-abstract subset ---
    na = stats["na"]
    na_n = na["n"]
    na_pct = lambda c: (c / na_n * 100) if na_n else 0
    H.append("<h2>5. Articles lacking an &lt;abstract&gt; (analysed separately)</h2>")
    H.append("<p class='meta'><strong>%s</strong> of %s articles (%.1f%%) have no "
             "<code>&lt;abstract&gt;</code>. Their structure is analysed on its own below.</p>"
             % (_fmt(na_n), _fmt(n), (na_n / n * 100 if n else 0)))
    na_struct = na["shapes"].get("structured", 0)
    na_cards = [
        (_fmt(na_n), "no-abstract articles"),
        (_fmt(na["body_present"]), "&lt;body&gt; present (%.1f%%)" % na_pct(na["body_present"])),
        (_fmt(na["body_absent"]), "&lt;body&gt; absent (%.1f%%)" % na_pct(na["body_absent"])),
        (_fmt(na_struct), "structured (%.1f%%)" % na_pct(na_struct)),
        (_fmt(na["reflist"]), "with references (%.1f%%)" % na_pct(na["reflist"])),
        (_fmt(na["imrad_full"]), "full IMRaD (%.1f%%)" % na_pct(na["imrad_full"])),
    ]
    H.append("<div class='stat-grid'>")
    for v, k in na_cards:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")
    if na_n:
        H.append("<div class='cols'><div>")
        H.append("<h3 style='margin-top:0'>Section coverage (no-abstract set)</h3>")
        H.append("<table><thead><tr><th>section</th><th class='num'>articles</th><th>coverage</th>"
                 "<th class='num'>%</th></tr></thead><tbody>")
        na_rows = [("Title", na["title"])] + [(s, na["sec_cov"].get(s, 0)) for s in CORE_SECTIONS] \
                  + [("References", na["reflist"])]
        na_vmax = max(c for _, c in na_rows) or 1
        for name, c in na_rows:
            H.append("<tr><td>%s</td><td class='num'>%s</td><td>%s</td><td class='num dim'>%.1f%%</td></tr>"
                     % (name, _fmt(c), bar(c, na_vmax, 200), na_pct(c)))
        H.append("</tbody></table></div><div>")
        H.append("<h3 style='margin-top:0'>Body shape (no-abstract set)</h3>")
        H.append("<table><thead><tr><th>shape</th><th class='num'>articles</th><th class='num'>%</th>"
                 "</tr></thead><tbody>")
        for name in ("structured", "OCR scanned", "unstructured", "front-only"):
            c = na["shapes"].get(name, 0)
            if c:
                H.append("<tr><td>%s</td><td class='num'>%s</td><td class='num dim'>%.1f%%</td></tr>"
                         % (name, _fmt(c), na_pct(c)))
        H.append("</tbody></table></div></div>")
        H.append("<p class='dim'>'front-only' here = no <code>&lt;body&gt;</code> <em>and</em> no "
                 "<code>&lt;abstract&gt;</code> &mdash; bare metadata records.</p>")
    else:
        H.append("<p class='dim'>Every analysed article carries an <code>&lt;abstract&gt;</code>.</p>")

    # --- 7. IMRaD completeness ---
    H.append("<h2>6. IMRaD completeness (structured articles)</h2>")
    H.append("<div class='cols'><div><table><thead><tr><th>IMRaD sections present</th>"
             "<th class='num'>articles</th><th class='num'>%</th></tr></thead><tbody>")
    for k in range(4, -1, -1):
        c = stats["imrad_hist"].get(k, 0)
        H.append("<tr><td>%d of 4</td><td class='num'>%s</td><td class='num dim'>%.1f%%</td></tr>"
                 % (k, _fmt(c), pct(c)))
    H.append("</tbody></table></div><div>")
    H.append("<h3 style='margin-top:0'>Sections per article</h3><table><thead><tr>"
             "<th>section count</th><th class='num'>articles</th></tr></thead><tbody>")
    band = {"0": 0, "1-3": 0, "4-6": 0, "7-10": 0, "11-15": 0, "16+": 0}
    for k, c in secn_hist.items():
        if k == 0: band["0"] += c
        elif k <= 3: band["1-3"] += c
        elif k <= 6: band["4-6"] += c
        elif k <= 10: band["7-10"] += c
        elif k <= 15: band["11-15"] += c
        else: band["16+"] += c
    for name in ("0", "1-3", "4-6", "7-10", "11-15", "16+"):
        H.append("<tr><td>%s</td><td class='num'>%s</td></tr>" % (name, _fmt(band[name])))
    H.append("</tbody></table></div></div>")

    # --- 8. References ---
    H.append("<h2>7. References</h2>")
    avg = (refs_total / refs_articles) if refs_articles else 0
    H.append("<div class='stat-grid'>"
             "<div class='stat'><div class='v'>%s</div><div class='k'>articles with &lt;ref-list&gt;</div></div>"
             "<div class='stat'><div class='v'>%s</div><div class='k'>total &lt;ref&gt; entries</div></div>"
             "<div class='stat'><div class='v'>%.0f</div><div class='k'>mean refs / article</div></div>"
             "</div>" % (_fmt(reflist_n), _fmt(refs_total), avg))

    # --- 9. Ancillary parts ---
    H.append("<h2>8. Ancillary parts</h2>")
    H.append("<table><thead><tr><th>part</th><th class='num'>articles</th><th>coverage</th>"
             "<th class='num'>%</th></tr></thead><tbody>")
    amax = max(anc.values()) if anc else 1
    for name in ("Figures", "Tables", "Supplementary", "Data availability", "Funding",
                 "Acknowledgements", "Conflict of interest", "Author contributions"):
        c = anc.get(name, 0)
        H.append("<tr><td>%s</td><td class='num'>%s</td><td>%s</td><td class='num dim'>%.1f%%</td></tr>"
                 % (name, _fmt(c), bar(c, amax, 220, "bar g"), pct(c)))
    H.append("</tbody></table>")

    H.append("""<h2>9. Caveats</h2><ul>
  <li><strong>Section labels are heuristic.</strong> Buckets come from <code>sec-type</code> + title
      keywords; unconventional headings (e.g. "Results and Discussion" counts for both; a standalone
      "Summary" counts as Conclusion) are mapped by best effort.</li>
  <li><strong>Front-only and OCR records have no sections to find</strong> &mdash; they depress the
      IMRaD figures. Structured-body articles are the meaningful denominator for section coverage.</li>
  <li><strong>Presence, not quality.</strong> The analysis records that a section exists, not how
      complete or correctly OCR'd its text is.</li>
</ul>""")
    H.append("</body></html>")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(SUMMARY_HTML, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))


def main():
    os.makedirs(XML_DIR, exist_ok=True)      # create input dir if missing
    files = sorted(glob.glob(os.path.join(XML_DIR, "PMC*.xml")))
    if not files:
        raise SystemExit("error: no XML files in %s" % XML_DIR)

    n = 0
    title_n = abstract_n = reflist_n = imrad_full = 0
    body_present = body_absent = 0
    refs_total = refs_articles = 0
    sec_cov = {s: 0 for s in CORE_SECTIONS}
    shapes = {}
    anc = {}
    imrad_hist = {}
    secn_hist = {}

    # Separate accumulator for the subset of articles that LACK an <abstract>.
    na = {"n": 0, "title": 0, "reflist": 0, "body_present": 0, "body_absent": 0,
          "imrad_full": 0, "sec_cov": {s: 0 for s in CORE_SECTIONS}, "shapes": {}}

    for i, path in enumerate(files, 1):
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                text = fh.read()
        except OSError:
            continue
        r = analyse_file(text)
        n += 1
        title_n += r["has_title"]
        abstract_n += r["has_abstract"]
        reflist_n += r["has_reflist"]
        if r["has_body"]:
            body_present += 1
        else:
            body_absent += 1
        if r["n_refs"]:
            refs_total += r["n_refs"]
            refs_articles += 1
        shapes[r["shape"]] = shapes.get(r["shape"], 0) + 1
        for s in r["sections"]:
            if s in sec_cov:
                sec_cov[s] += 1
        imrad = sum(1 for s in ("Introduction", "Methods", "Results", "Discussion")
                    if s in r["sections"])
        imrad_hist[imrad] = imrad_hist.get(imrad, 0) + 1
        if imrad == 4:
            imrad_full += 1
        secn_hist[r["n_sec_titles"]] = secn_hist.get(r["n_sec_titles"], 0) + 1
        for k, v in r["anc"].items():
            if v:
                anc[k] = anc.get(k, 0) + 1

        if not r["has_abstract"]:
            na["n"] += 1
            na["title"] += r["has_title"]
            na["reflist"] += r["has_reflist"]
            if r["has_body"]:
                na["body_present"] += 1
            else:
                na["body_absent"] += 1
            na["shapes"][r["shape"]] = na["shapes"].get(r["shape"], 0) + 1
            for s in r["sections"]:
                if s in na["sec_cov"]:
                    na["sec_cov"][s] += 1
            if imrad == 4:
                na["imrad_full"] += 1

        if i % 2000 == 0:
            print("[scan] %d/%d" % (i, len(files)), file=sys.stderr)

    stats = {"n": n, "title_n": title_n, "abstract_n": abstract_n, "reflist_n": reflist_n,
             "imrad_full": imrad_full, "refs_total": refs_total, "refs_articles": refs_articles,
             "body_present": body_present, "body_absent": body_absent,
             "sec_cov": sec_cov, "shapes": shapes, "anc": anc, "imrad_hist": imrad_hist,
             "secn_hist": secn_hist, "na": na}
    build_html(stats)
    print("[done] analysed %d files -> %s" % (n, SUMMARY_HTML), file=sys.stderr)
    print("[done] <body> present=%d absent=%d | structured=%d ocr=%d full-IMRaD=%d"
          % (body_present, body_absent, shapes.get("structured", 0),
             shapes.get("OCR scanned", 0), imrad_full), file=sys.stderr)
    print("[done] no-abstract subset=%d (body present=%d, absent=%d)"
          % (na["n"], na["body_present"], na["body_absent"]), file=sys.stderr)


if __name__ == "__main__":
    main()
