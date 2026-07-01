#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""relationships.py -- high-confidence GENETIC-GENETIC relationships.

================================ STRATEGY ===================================
From the fully ID-normalized triples
(TRIPLES/triples_GENETIC_DISEASE_CHEMICAL_normalized.json, produced by triples.py)
keep only the triples that are gene-gene relations between two genuine,
HGNC-resolved genes -- the most trustworthy slice for a gene-interaction network.

A triple qualifies when BOTH its subject and object satisfy ALL of:
  * "type": "GENETIC"            (a gene/protein entity, not disease/chemical), AND
  (1) "hgnc_symbol" != null      (the surface was normalized to an HGNC gene), AND
  (2) "control": "no"            (the gene is NOT an experimental-control / reagent
                                  surface, per controls.py's annotation), AND its
                                  hgnc_symbol is not in CONTROL_SYMBOLS -- symbols
                                  treated as PURELY experimental controls (currently
                                  MKI67 / Ki-67, a proliferation / IHC-staining marker),
                                  excluded on top of the "control": "yes" genes,
AND, in addition,
  (3) the triple's sentence (keyed by pmid + section + sentence) contains at least one
      DISEASE or CHEMICAL entity somewhere in it -- i.e. some triple from that same
      sentence in TRIPLES/triples.json has a "type": "DISEASE" or "type": "CHEMICAL"
      subject or object. This keeps gene-gene relations stated in a disease/chemical
      context (the gene pair co-occurs with a disease or drug in the sentence).

Each qualifying triple written to TRIPLES/genetic_genetic.json is additionally
annotated with three top-level keys:
  * "genetic_subject" : the subject's "hgnc_symbol" value (str, or list if ambiguous),
  * "genetic_object"  : the object's  "hgnc_symbol" value (str, or list if ambiguous),
  * "pair_count"      : how many qualifying triples share this exact directional
                        subject->object hgnc_symbol pair (the pair's frequency),
  * "diseases"        : sorted MONDO disease labels of DISEASE entities in the triple's
                        sentence that were normalized via the "disease key" route
                        ("mondo_via" == "disease key"); drives the disease filter,
  * "chemicals"       : sorted ChEBI chemical labels of any CHEMICAL entity in the
                        triple's sentence (enables chemical-specific visualization),
  * "year"            : publication year of the source article (from the PMC accession,
                        via the databases/pmc_years.json cache built by pub_years.py;
                        null if the year is not in the cache).

DISEASE FILTER
--------------
The disease labels for the gene-gene network are taken from the intersection of
sentences shared by genetic_genetic.json and the normalized triples file (i.e. the
sentences each qualifying gene-gene triple came from). From those shared sentences we
keep only "mondo_label" values whose DISEASE entity has "mondo_via" == "disease key";
diseases reached by any other route, or left unnormalized, are excluded.
A summary (counts, the 100 most frequent unique subject->object hgnc_symbol pairs,
top predicates) is written to TRIPLES/relationships.html.
=============================================================================

Input  : TRIPLES/triples_GENETIC_DISEASE_CHEMICAL_normalized.json  (qualifying slice)
         TRIPLES/triples.json   (base triples -- to find DISEASE/CHEMICAL sentences)
Outputs: TRIPLES/genetic_genetic.json   the qualifying GENETIC-GENETIC triples
         TRIPLES/relationships.html      summary + strategy

Run from anywhere::  python relationships.py
"""

import html
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "TRIPLES"
IN = OUT_DIR / "triples_GENETIC_DISEASE_CHEMICAL_normalized.json"
TRIPLES_BASE = OUT_DIR / "triples.json"
GG_OUT = OUT_DIR / "genetic_genetic.json"
HTML_OUT = OUT_DIR / "relationships.html"
# PMC accession -> publication year, populated by pub_years.py (NCBI E-utilities)
YEARS_CACHE = ROOT / "databases" / "pmc_years.json"

DC_TYPES = ("DISEASE", "CHEMICAL")
# only DISEASE entities normalized through this MONDO route feed the disease filter
DISEASE_VIA = "disease key"


def sent_key(t):
    """Identity of the sentence a triple was extracted from."""
    return (t["pmid"], t["section"], t["sentence"])


def symstr(el):
    hs = el["hgnc_symbol"]
    return hs if isinstance(hs, str) else "|".join(hs)


# symbols treated as PURELY experimental controls regardless of the per-triple
# "control" annotation -- e.g. MKI67 (Ki-67), a proliferation / IHC-staining marker
CONTROL_SYMBOLS = {"MKI67"}


def qualifies(el):
    """A GENETIC endpoint that is HGNC-resolved and not an experimental control."""
    if (el.get("type") != "GENETIC" or el.get("hgnc_symbol") is None
            or el.get("control") != "no"):
        return False
    hs = el["hgnc_symbol"]
    syms = hs if isinstance(hs, list) else [hs]
    return not any(s in CONTROL_SYMBOLS for s in syms)


def main():
    data = json.loads(IN.read_text(encoding="utf-8"))

    # sentences (pmid+section+sentence) that contain >= 1 DISEASE or CHEMICAL entity,
    # determined from the base triples in TRIPLES/triples.json
    base = json.loads(TRIPLES_BASE.read_text(encoding="utf-8"))
    dc_sentences = {sent_key(t) for t in base
                    if t["subject"]["type"] in DC_TYPES or t["object"]["type"] in DC_TYPES}

    # PMC accession -> publication year (offline cache produced by pub_years.py)
    years = json.loads(YEARS_CACHE.read_text(encoding="utf-8")) if YEARS_CACHE.exists() else {}

    gg = [t for t in data if qualifies(t["subject"]) and qualifies(t["object"])]
    qual = [t for t in gg if sent_key(t) in dc_sentences]

    # DISEASE FILTER: the disease labels used for the gene-gene network come from the
    # *intersection of sentences* between genetic_genetic.json (the qualifying triples)
    # and the normalized triples file -- i.e. the sentences a qualifying gene-gene
    # triple was drawn from. From those shared sentences we keep only the "mondo_label"
    # of DISEASE entities normalized via the "disease key" route ("mondo_via" ==
    # DISEASE_VIA); diseases reached by any other route (or left unnormalized) are
    # excluded. ChEBI chemical labels are collected over the same shared sentences.
    gg_sents = {sent_key(t) for t in qual}
    sent_diseases = defaultdict(set)
    sent_chemicals = defaultdict(set)
    for t in data:
        sk = sent_key(t)
        if sk not in gg_sents:
            continue
        for el in (t["subject"], t["object"]):
            if (el.get("type") == "DISEASE" and el.get("mondo_via") == DISEASE_VIA
                    and el.get("mondo_label")):
                ml = el["mondo_label"]
                sent_diseases[sk].add(ml if isinstance(ml, str) else "|".join(ml))
            elif el.get("type") == "CHEMICAL" and el.get("chebi_label"):
                cl = el["chebi_label"]
                sent_chemicals[sk].add(cl if isinstance(cl, str) else "|".join(cl))

    # frequency of each unique subject->object hgnc_symbol pair (directional)
    pairs = Counter((symstr(t["subject"]), symstr(t["object"])) for t in qual)
    # annotate each qualifying triple with the endpoint hgnc_symbols, the number of
    # times its subject->object hgnc_symbol pair occurs, and the MONDO disease
    # labels present in its sentence
    for t in qual:
        t["genetic_subject"] = t["subject"]["hgnc_symbol"]
        t["genetic_object"] = t["object"]["hgnc_symbol"]
        t["pair_count"] = pairs[(symstr(t["subject"]), symstr(t["object"]))]
        t["diseases"] = sorted(sent_diseases.get(sent_key(t), ()))
        t["chemicals"] = sorted(sent_chemicals.get(sent_key(t), ()))
        t["year"] = years.get(t["pmid"].split(".")[0])

    GG_OUT.write_text(json.dumps(qual, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")

    preds = Counter(t["predicate"]["text"].lower() for t in qual)
    n_gg_all = sum(1 for t in data
                   if t["subject"]["type"] == "GENETIC" and t["object"]["type"] == "GENETIC")

    print(f"input triples: {len(data):,}  (GENETIC-GENETIC: {n_gg_all:,})")
    print(f"DISEASE/CHEMICAL-bearing sentences: {len(dc_sentences):,}")
    print(f"GENETIC-GENETIC (hgnc_symbol != null AND control == 'no'): {len(gg):,}")
    print(f"  ... of those, in a DISEASE/CHEMICAL sentence: {len(qual):,}  -> {GG_OUT.name}")
    print(f"  unique gene pairs (subject hgnc_symbol -> object hgnc_symbol): {len(pairs):,}")
    print("  top pairs: " + ", ".join(
        f"{a}->{b} {c}" for (a, b), c in pairs.most_common(6)))

    render_html(len(data), n_gg_all, len(gg), len(dc_sentences), qual, pairs, preds)
    print(f"Wrote {HTML_OUT}")


def render_html(n_in, n_gg_all, n_gg_qual, n_dc_sent, qual, pairs, preds):
    esc = html.escape
    n = len(qual)

    pred_rows = "".join(
        f'<tr><td><code>{esc(p)}</code></td><td class="num">{c:,}</td></tr>'
        for p, c in preds.most_common(40))
    pair_rows = "".join(
        f'<tr><td>{esc(a)}</td><td>{esc(b)}</td><td class="num">{c:,}</td></tr>'
        for (a, b), c in pairs.most_common(100))
    sample = "".join(
        f'<tr><td>{esc(t["subject"]["text"])} '
        f'<span class="hg">{esc(symstr(t["subject"]))}</span></td>'
        f'<td><code>{esc(t["predicate"]["text"])}</code></td>'
        f'<td>{esc(t["object"]["text"])} '
        f'<span class="hg">{esc(symstr(t["object"]))}</span></td>'
        f'<td><code>{esc(t["pmid"])}</code></td></tr>' for t in qual[:60])

    style = (
        " body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;"
        "margin:2rem auto;padding:0 1rem;color:#1a1a1a;}"
        " h1{font-size:1.45rem;} h2{font-size:1.15rem;margin-top:1.8rem;"
        "border-bottom:1px solid #ddd;padding-bottom:.3rem;}"
        " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;}"
        " th,td{border:1px solid #ccc;padding:.3rem .6rem;text-align:left;vertical-align:top;}"
        " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
        " .big{font-size:2rem;font-weight:700;}"
        " .headline{background:#f1faf3;border:1px solid #c7e7d2;border-radius:8px;"
        "padding:.8rem 1rem;margin:1rem 0;}"
        " code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.85em;}"
        " .hg{color:#1a7f37;font-size:.82em;} details{margin:.6rem 0;} summary{cursor:pointer;color:#357;}"
        " p.note{color:#444;font-size:.92em;} ol li{margin:.25rem 0;}")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>relationships &mdash; high-confidence GENETIC-GENETIC relations</title>
<style>{style}</style></head><body>
<h1>High-confidence GENETIC &rarr; GENETIC relationships</h1>
<p>Gene-gene relation triples between two genuine, HGNC-resolved genes, distilled from
<code>TRIPLES/triples_GENETIC_DISEASE_CHEMICAL_normalized.json</code>. Qualifying
triples are written to <code>TRIPLES/genetic_genetic.json</code>. Produced by
<code>relationships.py</code>.</p>
<div class="headline"><span class="big">{n:,}</span> qualifying GENETIC&ndash;GENETIC
triples &middot; <span class="big">{len(pairs):,}</span> unique gene pairs
&mdash; from {n_gg_qual:,} HGNC-resolved, non-control gene&ndash;gene triples restricted
to the {n_dc_sent:,} sentences carrying a DISEASE or CHEMICAL entity
(of {n_gg_all:,} GENETIC&ndash;GENETIC triples; {n_in:,} total).</div>

<h2>Strategy</h2>
<p>A triple qualifies when <strong>both</strong> its subject and object satisfy all of:</p>
<ol>
<li><code>"type": "GENETIC"</code> &mdash; a gene/protein entity (not disease/chemical).</li>
<li><code>"hgnc_symbol"</code> is <strong>not <code>null</code></strong> &mdash; the
surface was normalized to an HGNC gene (via roman key / greek key / greek_expanded).</li>
<li><code>"control": "no"</code> &mdash; the gene is <strong>not</strong> an
experimental-control / reagent surface (per <code>controls.py</code>'s annotation;
excludes loading controls, reporters, etc.). <strong>In addition</strong>, the symbol
must not be in a small denylist of genes viewed as <strong>purely experimental
controls</strong> &mdash; currently <code>MKI67</code> (Ki-67), a proliferation /
IHC-staining marker. So <code>MKI67</code> is excluded on top of the
<code>"control": "yes"</code> genes, and no <code>MKI67</code> gene&ndash;gene pairs
appear in the network.</li>
</ol>
<p>And, in addition, a context restriction on the <strong>sentence</strong> the triple
was extracted from:</p>
<ol start="4">
<li>The triple's sentence (keyed by <code>pmid</code> + <code>section</code> +
<code>sentence</code>) must contain <strong>at least one DISEASE or CHEMICAL entity</strong>
&mdash; i.e. some triple from that same sentence in <code>TRIPLES/triples.json</code> has a
<code>"type": "DISEASE"</code> or <code>"type": "CHEMICAL"</code> subject or object. This
keeps only gene&ndash;gene relations stated alongside a disease or drug
({n_dc_sent:,} such sentences).</li>
</ol>
<p>Both endpoints must pass criteria 1&ndash;3 and the triple must pass criterion 4;
the result is the high-confidence gene&ndash;gene slice stated in a disease/chemical
context, suitable for a gene-interaction network. Each qualifying triple in
<code>genetic_genetic.json</code> (which also carries the predicate, PMID, section,
sentence and the per-element <code>hgnc_symbol</code>/<code>control</code>) is annotated
with three extra top-level keys:</p>
<ul>
<li><code>"genetic_subject"</code> &mdash; the subject's <code>hgnc_symbol</code> value
(string, or a list when the symbol is ambiguous).</li>
<li><code>"genetic_object"</code> &mdash; the object's <code>hgnc_symbol</code> value
(string, or a list when the symbol is ambiguous).</li>
<li><code>"pair_count"</code> &mdash; how many qualifying triples share this exact
directional <code>genetic_subject</code> &rarr; <code>genetic_object</code> pair (the
frequency of that unique gene pair across <code>genetic_genetic.json</code>).</li>
<li><code>"diseases"</code> &mdash; the sorted MONDO disease labels driving the disease
filter (see below).</li>
<li><code>"chemicals"</code> &mdash; the sorted ChEBI chemical labels of any CHEMICAL
entity in the triple's sentence (enables chemical-specific filtering).</li>
<li><code>"year"</code> &mdash; the publication year of the source article, resolved
from its PMC accession via NCBI E-utilities and cached in
<code>databases/pmc_years.json</code> (built by <code>pub_years.py</code>);
<code>null</code> if not in the cache.</li>
</ul>

<h2>Disease filter</h2>
<p>The <code>"diseases"</code> annotation that powers disease-specific visualization is
built from the <strong>intersection of sentences</strong> shared by
<code>TRIPLES/genetic_genetic.json</code> and
<code>TRIPLES/triples_GENETIC_DISEASE_CHEMICAL_normalized.json</code> &mdash; that is,
the sentences each qualifying gene&ndash;gene triple was drawn from. From those shared
sentences we include only the <code>"mondo_label"</code> values of DISEASE entities
normalized through the <strong><code>"mondo_via": "disease key"</code></strong> route;
diseases reached by any other route, or left unnormalized
(<code>mondo_via</code> = <code>null</code>), are excluded. So a qualifying triple's
<code>"diseases"</code> list contains exactly the disease-key&ndash;normalized MONDO
labels mentioned in its own sentence.</p>

<h2>Top 100 gene pairs (subject &rarr; object, by # triples)</h2>
<details open><summary>show / hide table</summary>
<table><tr><th>subject hgnc_symbol</th><th>object hgnc_symbol</th>
<th class="num">triples</th></tr>{pair_rows}</table></details>

<h2>Top predicates</h2>
<details><summary>show / hide table</summary>
<table><tr><th>predicate</th><th class="num">triples</th></tr>{pred_rows}</table></details>

<h2>Sample qualifying triples (first 60)</h2>
<details><summary>show / hide table</summary>
<table><tr><th>subject</th><th>predicate</th><th>object</th><th>pmid</th></tr>
{sample}</table></details>
</body></html>
"""
    HTML_OUT.write_text(doc, encoding="utf-8")


if __name__ == "__main__":
    main()
