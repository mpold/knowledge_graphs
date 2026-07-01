#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""phenotypes.py -- triage DISEASE entities that are processes / phenotypes.

The DISEASE analog of controls.py (which flags GENETIC entities that are lab
controls/tools). Finds the DISEASE-labelled surfaces (from sentences/ via
disease.py's clean_disease_ne.tsv) whose referent is a biological PROCESS,
PHENOTYPE or pathological FINDING -- cytotoxicity, hypoxia, tumorigenesis,
metastasis, inflammation, fibrosis, ... -- rather than a disease entity. Writes
DISEASE/phenotypes.html, prints the analytical steps, and annotates the disease
libraries with a nested ("phenotype": "yes"/"no").

================================ STRATEGY ===================================
(0) UPSTREAM -- the aggregated DISEASE surfaces (DISEASE/clean_disease_ne.tsv from
    disease.py; or rebuilt from sentences/*.json with --rebuild).

(1) BROAD SCAN (illustrative). A naive "contains a process keyword" regex
    OVER-captures real diseases that merely contain the word (malignant glioma,
    radiation necrosis, lung injury). Printed only to show what refining removes.

(2) REFINED, ANCHORED CLASSIFICATION. A surface is a process/phenotype only if its
    core (casefold + British->American + plural-strip) EXACTLY equals a curated
    process/phenotype term, or ends with a productive process suffix
    (-toxicity / -genesis). Exact-core matching excludes disease phrases that
    merely contain the word ('malignant glioma' != 'glioma'-process; 'radiation
    necrosis' != 'necrosis').

(3) TWO CONFIDENCE TIERS.
      Tier A -- pure processes / phenomena that are NEVER a disease: cell death /
        cytotoxicity, stress / hypoxia / damage, oncogenic / growth processes
        (tumorigenesis, metastasis, invasion, senescence, ...).
      Tier B -- pathological STATES / findings that MONDO sometimes catalogs as
        diseases (inflammation, fibrosis, edema, hyperplasia, malignancy, ...):
        dual-use -- some DO link to a MONDO term, but here read as a finding.

(4) ANNOTATION. Every entry in disease.json, disease_ambiguous.json and
    unmatched_disease.json gets ("phenotype": "yes") iff its key is case-
    sensitively equal to a Tier A/B process/phenotype surface, else "no".
    phenotypes.html notes, per flagged surface, whether it links to MONDO
    (dual-use) or is unmatched, with the MONDO label where applicable.
=============================================================================

Run from anywhere (paths resolve relative to this file)::

    python phenotypes.py            # uses DISEASE/clean_disease_ne.tsv if present
    python phenotypes.py --rebuild  # re-aggregate DISEASE from sentences/*.json
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
OUT_DIR = ROOT / "DISEASE"
TSV = OUT_DIR / "clean_disease_ne.tsv"
HTML_OUT = OUT_DIR / "phenotypes.html"
LABEL = "DISEASE"
LIB_FILES = ["disease.json", "disease_ambiguous.json", "unmatched_disease.json"]

DASH_VARIANTS = "-‐‑‒–—―−⁃"
_TO_ASCII = str.maketrans({c: "-" for c in DASH_VARIANTS})
_DESPACE = re.compile(r"\s*-\s*")


def dash_normalize(text):
    return _DESPACE.sub("-", (text or "").translate(_TO_ASCII))


_UK_US = [(re.compile(p), r) for p, r in [
    (r"tumour", "tumor"), (r"oedema", "edema"), (r"oesophag", "esophag"),
    (r"haem", "hem"), (r"anaemi", "anemi"), (r"leukaemi", "leukemi"),
    (r"ischaemi", "ischemi"), (r"paediatr", "pediatr")]]


def britishize(s):
    for rx, rep in _UK_US:
        s = rx.sub(rep, s)
    return s


def singulars(s):
    out = []
    if s.endswith("ies") and len(s) > 4:
        out.append(s[:-3] + "y")          # malignancies -> malignancy
    if s.endswith("es") and len(s) > 4:
        out.append(s[:-2])
    if s.endswith("s") and not s.endswith("ss") and len(s) > 3:
        out.append(s[:-1])                # lesions -> lesion, tumorspheres -> tumorsphere
    return out


# ---- curated process / phenotype vocabulary (Tier A = pure process) ----
TIER_A = {
    "cell death / cytotoxicity": {
        "apoptosis", "apoptotic", "necrosis", "necrotic", "necroptosis",
        "pyroptosis", "ferroptosis", "cuproptosis", "autophagy", "autophagic",
        "cytotoxicity", "cytotoxic", "cell death", "viability", "cell viability",
        "anoikis", "oncosis", "paraptosis", "parthanatos"},
    "stress / hypoxia / damage": {
        "hypoxia", "hypoxic", "oxidative stress", "er stress",
        "endoplasmic reticulum stress", "dna damage", "replication stress",
        "genotoxic stress", "metabolic stress", "starvation"},
    "oncogenic / growth process": {
        "tumorigenesis", "tumorigenicity", "tumorigenic", "oncogenesis",
        "oncogenic", "carcinogenesis", "carcinogenic", "gliomagenesis",
        "leukemogenesis", "lymphomagenesis", "proliferation", "proliferative",
        "invasion", "invasiveness", "invasive", "migration", "metastasis",
        "metastases", "metastatic", "dissemination", "angiogenesis",
        "neovascularization", "vascularization", "emt",
        "epithelial-mesenchymal transition", "epithelial mesenchymal transition",
        "transformation", "malignant transformation", "immortalization",
        "senescence", "senescent", "stemness", "self-renewal", "sphere formation",
        "colony formation", "clonogenicity", "anchorage-independent growth",
        "tumorsphere", "gliomasphere", "neurosphere", "progression",
        "recurrence", "relapse", "regression", "remission", "growth retardation"},
}
# Tier B = pathological states / findings (dual-use; some are MONDO diseases)
TIER_B = {
    "pathological state / finding": {
        "inflammation", "inflammatory", "neuroinflammation", "fibrosis", "edema",
        "hyperplasia", "dysplasia", "metaplasia", "atrophy", "hypertrophy",
        "lesion", "injury", "damage", "dysfunction", "malignancy", "morbidity",
        "mortality", "resistance", "chemoresistance", "radioresistance",
        "hypersensitivity"},
}
_A_LOOKUP = {t: cat for cat, terms in TIER_A.items() for t in terms}
_B_LOOKUP = {t: cat for cat, terms in TIER_B.items() for t in terms}
_TOX = re.compile(r"(toxicity|toxicities)$")
_GENESIS = re.compile(r"genesis$")
TIER_A_CATS = list(TIER_A) + ["(suffix) -toxicity / -genesis"]
TIER_B_CATS = list(TIER_B)


def classify(v):
    """Return (category, tier 'A'/'B') if v is a process/phenotype, else None."""
    core = britishize(v.casefold())
    for c in [core] + singulars(core):
        if c in _A_LOOKUP:
            return _A_LOOKUP[c], "A"
        if c in _B_LOOKUP:
            return _B_LOOKUP[c], "B"
    if _TOX.search(core) or _GENESIS.search(core):   # neurotoxicity, gliomagenesis
        return "(suffix) -toxicity / -genesis", "A"
    return None


# ---- broad scan (illustrative over-capture) ----
_BROAD = re.compile(
    r"(?i)(apoptos|necros|ferroptos|cytotox|toxic|hypoxi|tumou?rigen|metasta|"
    r"senescen|invasi|inflammat|fibros|edema|oedema|hyperplas|malignan|sphere)")
# real diseases that the broad scan grabs but the refined pass MUST reject
EXCLUSION_CHECKS = ["malignant glioma", "radiation necrosis", "glioblastoma",
                    "breast cancer", "lung injury", "intestinal metaplasia"]


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


# ============================================================ annotation
def annotate_libraries(pheno_values):
    """Add ("phenotype": yes/no) to every entry in the disease libraries; return
    (matched_label: surface->MONDO label for matched libs, counts: f->(yes,no))."""
    matched_label, counts = {}, {}
    for f in LIB_FILES:
        path = OUT_DIR / f
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        yes = 0
        for k, e in data.items():
            if not isinstance(e, dict):
                continue
            flag = "yes" if k in pheno_values else "no"
            e["phenotype"] = flag
            yes += (flag == "yes")
            if flag == "yes" and "mondo_label" in e and f != "unmatched_disease.json":
                lbl = e["mondo_label"]
                matched_label[k] = lbl if isinstance(lbl, str) else "/".join(lbl)
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
        counts[f] = (yes, len(data) - yes)
    return matched_label, counts


# ============================================================ run
def run(rows, info):
    cat = defaultdict(list)         # category -> [(surface, occ, tier)]
    tier_of = {}
    for v, o in rows:
        r = classify(v)
        if r:
            c, t = r
            cat[c].append((v, o))
            tier_of[v] = t
    for c in cat:
        cat[c].sort(key=lambda x: (-x[1], x[0].casefold()))

    pheno_values = {v for items in cat.values() for v, _ in items}
    gocc = sum(o for items in cat.values() for _, o in items)
    nA = sum(len(cat[c]) for c in TIER_A_CATS)
    oA = sum(o for c in TIER_A_CATS for _, o in cat[c])
    nB = sum(len(cat[c]) for c in TIER_B_CATS)
    oB = sum(o for c in TIER_B_CATS for _, o in cat[c])
    corpus_n, corpus_occ = len(rows), info["occ"]

    print("=" * 70)
    print("PROCESS / PHENOTYPE TRIAGE (DISEASE entities)")
    print("=" * 70)
    print(f"(0) UPSTREAM: source = {info['source']}")
    print(f"      {corpus_n:,} distinct DISEASE surfaces, {corpus_occ:,} occurrences")
    nb = [(v, o) for v, o in rows if _BROAD.search(v)]
    print(f"\n(1) BROAD keyword scan (naive): {len(nb):,} surfaces "
          f"({sum(o for _, o in nb):,} occ) -- over-captures real diseases")
    print("\n(2/3) REFINED, tiered classification:")
    for c in TIER_A_CATS + TIER_B_CATS:
        tier = "A" if c in TIER_A_CATS else "B"
        print(f"      [{tier}] {c:34} {len(cat[c]):>4} surfaces, "
              f"{sum(o for _, o in cat[c]):>6} occ")
    print("\n    exclusion checks (real diseases -- must NOT be flagged):")
    for ex in EXCLUSION_CHECKS:
        print(f"      {ex:24} -> phenotype? {classify(ex) is not None}")

    matched_label, ann = annotate_libraries(pheno_values)
    n_dual = sum(1 for v in pheno_values if v in matched_label)
    print(f"\n(4) ANNOTATION of disease libraries (added \"phenotype\": yes/no):")
    for f in LIB_FILES:
        if f in ann:
            y, n = ann[f]
            print(f"      {f:26} phenotype=yes {y:>4}   phenotype=no {n:>6}")

    print(f"\nTOTAL: {len(pheno_values):,} process/phenotype surfaces, {gocc:,} occ"
          f"  ({100*len(pheno_values)/corpus_n:.1f}% of surfaces, "
          f"{100*gocc/corpus_occ:.1f}% of occ)")
    print(f"      Tier A (pure process/phenomenon): {nA:,} surfaces, {oA:,} occ")
    print(f"      Tier B (pathological state, dual-use): {nB:,} surfaces, {oB:,} occ"
          f"  ({n_dual} also link to a MONDO term)")

    write_html(cat, tier_of, matched_label, pheno_values, gocc, nA, oA, nB, oB,
               corpus_n, corpus_occ, n_dual)
    print(f"\nWrote {HTML_OUT}")


def write_html(cat, tier_of, matched_label, pheno_values, gocc, nA, oA, nB, oB,
               corpus_n, corpus_occ, n_dual):
    esc = html.escape

    def section(c):
        items = cat[c]
        occ = sum(o for _, o in items)
        crows = []
        for v, o in items:
            ml = matched_label.get(v)
            link = (f'<span class="dual">MONDO: {esc(ml)}</span>' if ml
                    else '<span class="un">unmatched</span>')
            crows.append(f'<tr><td><code>{esc(v)}</code></td>'
                         f'<td class="num">{o:,}</td><td>{link}</td></tr>')
        return (f'<h3>{esc(c)} &mdash; {len(items)} surfaces, {occ:,} occ</h3>'
                f'<table><tr><th>surface</th><th class="num">occ</th>'
                f'<th>MONDO link</th></tr>{"".join(crows)}</table>')

    summ = "".join(
        f'<tr><td>{esc(c)}</td><td>{("A" if c in TIER_A_CATS else "B")}</td>'
        f'<td class="num">{len(cat[c]):,}</td>'
        f'<td class="num">{sum(o for _, o in cat[c]):,}</td></tr>'
        for c in TIER_A_CATS + TIER_B_CATS if cat[c])

    style = (
        " body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:980px;"
        "margin:2rem auto;padding:0 1rem;color:#1a1a1a;}"
        " h1{font-size:1.45rem;} h2{font-size:1.15rem;margin-top:1.8rem;"
        "border-bottom:1px solid #ddd;padding-bottom:.3rem;} h3{font-size:1rem;margin-top:1.3rem;}"
        " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;}"
        " th,td{border:1px solid #ccc;padding:.28rem .6rem;text-align:left;}"
        " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
        " .big{font-size:2rem;font-weight:700;}"
        " .headline{background:#f6f3fb;border:1px solid #ddd0ee;border-radius:8px;"
        "padding:.8rem 1rem;margin:1rem 0;}"
        " code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.9em;}"
        " .dual{color:#9a6700;} .un{color:#999;} p.note{color:#444;font-size:.92em;}")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>DISEASE entities that are processes / phenotypes</title>
<style>{style}</style></head><body>
<h1>DISEASE entities that are processes / phenotypes</h1>
<p>DISEASE-labelled surfaces (from <code>sentences/*.json</code>) whose referent is
a biological <strong>process, phenotype or pathological finding</strong> &mdash;
cytotoxicity, hypoxia, tumorigenesis, metastasis, inflammation, &hellip; &mdash;
<em>not</em> a disease entity. The DISEASE analog of <code>controls.py</code>.
Produced by <code>phenotypes.py</code>.</p>
<div class="headline"><span class="big">{len(pheno_values):,}</span> surfaces
&nbsp;&middot;&nbsp; <span class="big">{gocc:,}</span> occurrences &mdash;
{100*len(pheno_values)/corpus_n:.1f}% of the {corpus_n:,} DISEASE surfaces
({100*gocc/corpus_occ:.1f}% of occurrences). {n_dual} also link to a MONDO term
(dual-use Tier&nbsp;B). All disease libraries now carry a nested
<code>"phenotype"</code> field (yes/no).</div>

<h2>Summary by category</h2>
<table><tr><th>category</th><th>tier</th><th class="num">surfaces</th><th class="num">occ</th></tr>
{summ}
<tr><td><strong>Tier A &mdash; pure process/phenomenon</strong></td><td>A</td><td class="num"><strong>{nA:,}</strong></td><td class="num"><strong>{oA:,}</strong></td></tr>
<tr><td><strong>Tier B &mdash; pathological state (dual-use)</strong></td><td>B</td><td class="num"><strong>{nB:,}</strong></td><td class="num"><strong>{oB:,}</strong></td></tr>
<tr><td><strong>Total</strong></td><td></td><td class="num"><strong>{len(pheno_values):,}</strong></td><td class="num"><strong>{gocc:,}</strong></td></tr>
</table>

<h2>Method &amp; caveats</h2>
<p class="note">A surface is flagged only if its core (casefold + British&rarr;American
+ plural-strip) <strong>exactly equals</strong> a curated process/phenotype term, or
ends with a productive process suffix (<code>-toxicity</code>/<code>-genesis</code>).
Exact-core matching keeps real diseases that merely contain the word out
(<code>malignant glioma</code>, <code>radiation necrosis</code>,
<code>lung injury</code> are <em>not</em> flagged). <strong>Tier A</strong> terms
are never diseases; <strong>Tier B</strong> are pathological states MONDO sometimes
catalogs (<code>fibrosis</code>, <code>inflammation</code>, &hellip;) &mdash; the
<em>MONDO link</em> column shows which actually resolved (dual-use).</p>

<h2>Tier A &mdash; pure processes / phenomena ({nA} surfaces, {oA:,} occ)</h2>
{''.join(section(c) for c in TIER_A_CATS if cat[c])}

<h2>Tier B &mdash; pathological states / findings ({nB} surfaces, {oB:,} occ)</h2>
{''.join(section(c) for c in TIER_B_CATS if cat[c])}
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
