#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""nonchemical.py -- triage CHEMICAL entities that are NOT chemicals.

The CHEMICAL analog of controls.py (GENETIC) and phenotypes.py (DISEASE). Finds
the CHEMICAL-labelled surfaces (from sentences/ via chemical.py's
clean_chemical_ne.tsv) whose referent is NOT a chemical -- a gene/protein, a
nucleic-acid / genetic element, or a process/therapy concept -- rather than a
compound. Writes CHEMICAL/nonchemical.html, prints the analytical steps, and
annotates the chemical libraries with a nested ("non_chemical": "yes"/"no").

================================ STRATEGY ===================================
(0) UPSTREAM -- the aggregated CHEMICAL surfaces (CHEMICAL/clean_chemical_ne.tsv
    from chemical.py; or rebuilt from sentences/*.json with --rebuild), plus the
    ChEBI-matched set (chemical.json + chemical_ambiguous.json) and the HGNC
    approved-symbol set (hgnc_complete_set).

(1) CLASSIFICATION into three categories (first that applies wins):
      * nucleic acid / genetic element -- LINC#/circ*/MIR*/SNHG#/miR-*/let-7/
        sh|si|sgRNA/lncRNA/mRNA/cDNA/... (named genetic elements, not compounds).
      * process / therapy / concept -- a small curated set (ICB, ICI,
        chemotherapy, radiotherapy, immunotherapy, PDT, TACE, ...).
      * gene / protein -- a ChEBI-UNMATCHED surface that is an HGNC APPROVED symbol
        (high precision; aliases are skipped because chemical abbreviations
        coincide with gene aliases -- DDP=cisplatin, FDG, DAG, ADR), minus a
        DENYLIST of approved-symbol/chemical collisions (LPA, ADM=adriamycin,
        GEM=gemcitabine, ACR, CORT), plus a curated set of alias-only proteins
        (TRAIL, TDP-43, YAP, TAZ).

(2) ANNOTATION. Every entry in chemical.json, chemical_ambiguous.json and
    unmatched_chemical.json gets ("non_chemical": "yes") iff its key is a
    process/nucleic-acid/gene-protein surface, else "no". A surface that DID link
    to ChEBI is treated as a real chemical (gene/protein category requires
    ChEBI-unmatched), so only genetic-element/process strings can flag a matched
    entry. nonchemical.html lists every flagged surface by category.
=============================================================================

Run from anywhere (paths resolve relative to this file)::

    python nonchemical.py            # uses CHEMICAL/clean_chemical_ne.tsv if present
    python nonchemical.py --rebuild  # re-aggregate CHEMICAL from sentences/*.json
"""

import argparse
import glob
import html
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent
SENT_DIR = ROOT / "sentences"
OUT_DIR = ROOT / "CHEMICAL"
HGNC_PATH = ROOT / "databases" / "hgnc_complete_set_2026-05-01.json"
TSV = OUT_DIR / "clean_chemical_ne.tsv"
HTML_OUT = OUT_DIR / "nonchemical.html"
LABEL = "CHEMICAL"
LIB_FILES = ["chemical.json", "chemical_ambiguous.json", "unmatched_chemical.json"]
MATCHED_LIBS = {"chemical.json", "chemical_ambiguous.json"}

DASH_VARIANTS = "-‐‑‒–—―−⁃"
_TO_ASCII = str.maketrans({c: "-" for c in DASH_VARIANTS})
_DESPACE = re.compile(r"\s*-\s*")


def dash_normalize(text):
    return _DESPACE.sub("-", (text or "").translate(_TO_ASCII))


CAT_NUCLEIC = "nucleic acid / genetic element"
CAT_PROCESS = "process / therapy / concept"
CAT_GENE = "gene / protein"
CATEGORIES = [CAT_NUCLEIC, CAT_PROCESS, CAT_GENE]

_NUCLEIC = re.compile(
    r"(?i)(^linc\d|^circ[a-z0-9]|^snhg\d|^mir\d|^mir[- ]?\d|mir\w*hg$|"
    r"^(hsa-)?miR-?\d|^let-7|^lncrna|"
    r"(mrna|cdna|sirna|shrna|sgrna|mirna|lncrna|ncrna|pirna|snorna|"
    r"oligonucleotide|antisense|aptamer)\b)")
_PROCESS = {"icb", "ici", "icbs", "chemotherapy", "radiotherapy", "immunotherapy",
            "chemoradiation", "chemoradiotherapy", "photodynamic therapy", "pdt",
            "tace", "hyperthermia", "checkpoint blockade", "checkpoint inhibition"}
# approved HGNC symbols that are really chemicals in this corpus -> NOT gene/protein
DENYLIST = {"LPA", "ADM", "GEM", "ACR", "CORT"}
# alias-only proteins (not approved symbols) seen mislabelled CHEMICAL
CURATED_PROTEIN = {"TRAIL", "TDP-43", "YAP", "TAZ"}


def classify(v, is_matched, approved):
    """Return the non-chemical category for surface v, else None."""
    if _NUCLEIC.search(v):
        return CAT_NUCLEIC
    if v.casefold() in _PROCESS:
        return CAT_PROCESS
    if v in CURATED_PROTEIN:
        return CAT_GENE
    if (not is_matched) and v in approved and v not in DENYLIST:
        return CAT_GENE
    return None


# ============================================================ upstream
def aggregate_from_sentences():
    counts = Counter()
    nf = 0
    for fp in sorted(glob.glob(str(SENT_DIR / "*.json"))):
        try:
            rec = json.load(open(fp, encoding="utf-8"))
        except Exception:
            continue
        nf += 1
        for sent in rec.get("sentences", []):
            for ent in sent.get("entities", []):
                if ent.get("label") == LABEL:
                    counts[ent.get("text", "")] += 1
    agg = Counter()
    for form, n in counts.items():
        agg[dash_normalize(form)] += n
    rows = sorted(agg.items(), key=lambda kv: (-kv[1], kv[0]))
    return rows, {"source": f"{nf:,} sentences/*.json files",
                  "occ": sum(counts.values())}


def load_from_tsv():
    rows = []
    for ln in TSV.read_text(encoding="utf-8").splitlines()[1:]:
        p = ln.split("\t")
        if len(p) >= 2 and p[1].isdigit():
            rows.append((p[0], int(p[1])))
    return rows, {"source": str(TSV.relative_to(ROOT)),
                  "occ": sum(o for _, o in rows)}


def load_approved():
    docs = json.loads(HGNC_PATH.read_text(encoding="utf-8"))["response"]["docs"]
    return {d["symbol"] for d in docs if d.get("symbol")}


def matched_set():
    s = set()
    for f in MATCHED_LIBS:
        p = OUT_DIR / f
        if p.exists():
            s |= set(json.loads(p.read_text(encoding="utf-8")))
    return s


# ============================================================ annotation
def annotate_libraries(nonchem_values):
    counts = {}
    for f in LIB_FILES:
        path = OUT_DIR / f
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        yes = 0
        for k, e in data.items():
            if not isinstance(e, dict):
                continue
            flag = "yes" if k in nonchem_values else "no"
            e["non_chemical"] = flag
            yes += (flag == "yes")
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
        counts[f] = (yes, len(data) - yes)
    return counts


# ============================================================ run
def run(rows, info):
    approved = load_approved()
    matched = matched_set()
    cat = defaultdict(list)
    for v, o in rows:
        c = classify(v, v in matched, approved)
        if c:
            cat[c].append((v, o))
    for c in cat:
        cat[c].sort(key=lambda x: (-x[1], x[0].casefold()))

    nonchem = {v for items in cat.values() for v, _ in items}
    gocc = sum(o for items in cat.values() for _, o in items)
    corpus_n, corpus_occ = len(rows), info["occ"]

    print("=" * 70)
    print("NON-CHEMICAL TRIAGE (CHEMICAL entities)")
    print("=" * 70)
    print(f"(0) UPSTREAM: source = {info['source']}")
    print(f"      {corpus_n:,} distinct CHEMICAL surfaces, {corpus_occ:,} occurrences; "
          f"{len(matched):,} link to ChEBI; {len(approved):,} HGNC approved symbols")
    print("\n(1) CLASSIFICATION by category:")
    for c in CATEGORIES:
        print(f"      {c:34} {len(cat[c]):>4} surfaces, "
              f"{sum(o for _, o in cat[c]):>5} occ")
    print("\n    exclusion checks (real chemicals -- must NOT be flagged):")
    for ex in ["cisplatin", "TMZ", "glucose", "DOX", "LPA", "GEM", "ADM", "paclitaxel"]:
        print(f"      {ex:12} -> non-chemical? "
              f"{classify(ex, ex in matched, approved) is not None}")

    counts = annotate_libraries(nonchem)
    print("\n(2) ANNOTATION of chemical libraries (added \"non_chemical\": yes/no):")
    for f in LIB_FILES:
        if f in counts:
            y, n = counts[f]
            print(f"      {f:26} non_chemical=yes {y:>4}   non_chemical=no {n:>6}")

    print(f"\nTOTAL: {len(nonchem):,} non-chemical surfaces, {gocc:,} occ "
          f"({100*len(nonchem)/corpus_n:.1f}% of surfaces, "
          f"{100*gocc/corpus_occ:.1f}% of occ)")

    write_html(cat, nonchem, gocc, corpus_n, corpus_occ)
    print(f"\nWrote {HTML_OUT}")


def write_html(cat, nonchem, gocc, corpus_n, corpus_occ):
    esc = html.escape

    def section(c):
        items = cat[c]
        occ = sum(o for _, o in items)
        cells = "".join(f'<tr><td><code>{esc(v)}</code></td>'
                        f'<td class="num">{o:,}</td></tr>' for v, o in items)
        return (f'<h3>{esc(c)} &mdash; {len(items)} surfaces, {occ:,} occ</h3>'
                f'<table><tr><th>surface</th><th class="num">occ</th></tr>'
                f'{cells}</table>')

    summ = "".join(
        f'<tr><td>{esc(c)}</td><td class="num">{len(cat[c]):,}</td>'
        f'<td class="num">{sum(o for _, o in cat[c]):,}</td></tr>'
        for c in CATEGORIES if cat[c])

    style = (
        " body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:980px;"
        "margin:2rem auto;padding:0 1rem;color:#1a1a1a;}"
        " h1{font-size:1.45rem;} h2{font-size:1.15rem;margin-top:1.8rem;"
        "border-bottom:1px solid #ddd;padding-bottom:.3rem;} h3{font-size:1rem;margin-top:1.3rem;}"
        " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;}"
        " th,td{border:1px solid #ccc;padding:.28rem .6rem;text-align:left;}"
        " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
        " .big{font-size:2rem;font-weight:700;}"
        " .headline{background:#fbf3f6;border:1px solid #eed0e0;border-radius:8px;"
        "padding:.8rem 1rem;margin:1rem 0;}"
        " code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.9em;}"
        " p.note{color:#444;font-size:.92em;}")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>CHEMICAL entities that are not chemicals</title>
<style>{style}</style></head><body>
<h1>CHEMICAL entities that are not chemicals</h1>
<p>CHEMICAL-labelled surfaces (from <code>sentences/*.json</code>) whose referent is
<strong>not a chemical</strong> &mdash; a gene/protein, a nucleic-acid / genetic
element, or a process/therapy concept &mdash; i.e. NER mislabels. The CHEMICAL
analog of <code>controls.py</code> / <code>phenotypes.py</code>. Produced by
<code>nonchemical.py</code>.</p>
<div class="headline"><span class="big">{len(nonchem):,}</span> surfaces
&nbsp;&middot;&nbsp; <span class="big">{gocc:,}</span> occurrences &mdash;
{100*len(nonchem)/corpus_n:.1f}% of the {corpus_n:,} CHEMICAL surfaces
({100*gocc/corpus_occ:.1f}% of occurrences). All chemical libraries now carry a
nested <code>"non_chemical"</code> field (yes/no).</div>

<h2>Summary by category</h2>
<table><tr><th>category</th><th class="num">surfaces</th><th class="num">occ</th></tr>
{summ}
<tr><td><strong>Total</strong></td><td class="num"><strong>{len(nonchem):,}</strong></td><td class="num"><strong>{gocc:,}</strong></td></tr>
</table>

<h2>Method &amp; caveats</h2>
<p class="note"><strong>gene/protein</strong> = a ChEBI-unmatched surface equal to an
HGNC <em>approved</em> symbol (aliases skipped &mdash; chemical abbreviations
coincide with gene aliases: <code>DDP</code>=cisplatin, <code>FDG</code>,
<code>DAG</code>), minus a denylist of approved-symbol/chemical collisions
(<code>LPA</code>, <code>ADM</code>=adriamycin, <code>GEM</code>=gemcitabine,
<code>ACR</code>, <code>CORT</code>), plus curated alias-only proteins
(<code>TRAIL</code>, <code>TDP-43</code>, <code>YAP</code>, <code>TAZ</code>).
<strong>nucleic acid</strong> = named genetic elements (<code>LINC#</code>,
<code>circ*</code>, <code>MIR*</code>, <code>SNHG#</code>, miR/siRNA/lncRNA/...).
<strong>process/therapy</strong> = a curated concept set (<code>ICB</code>,
<code>chemotherapy</code>, <code>PDT</code>, ...). Real chemicals (incl. the
denylisted collisions and ChEBI-matched compounds) are not flagged.</p>

{''.join(section(c) for c in CATEGORIES if cat[c])}
</body></html>
"""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    HTML_OUT.write_text(doc, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rebuild", action="store_true")
    args = ap.parse_args()
    if args.rebuild or not TSV.exists():
        rows, info = aggregate_from_sentences()
    else:
        rows, info = load_from_tsv()
    run(rows, info)


if __name__ == "__main__":
    main()
