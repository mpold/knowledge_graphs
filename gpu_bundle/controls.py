#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""controls.py -- identify GENETIC entities that are experimental controls/tools.

Finds the GENETIC named entities (from the BioBERT output under sentences/) whose
referent is a *laboratory control or tool* -- a reporter, an epitope tag, a
housekeeping/loading normalizer, the Cas9 nuclease, or a non-targeting/scramble
control -- rather than a gene/protein that participates in the biology or disease
under study. Writes GENETIC/controls.html and prints the analytical steps.

================================ STRATEGY ===================================
(0) UPSTREAM -- corpus of GENETIC entities. The same aggregation roman.py STAGE 1
    performs: collect every entity labelled GENETIC from sentences/*.json,
    dash-normalize (unify dash glyphs to '-', strip whitespace around '-'), and
    sum occurrences per normalized surface. (Fast path: read the already-built
    GENETIC/clean_genetic_ne.tsv; otherwise rebuild from sentences/.)

(1) BROAD SCAN (illustrative). A naive "does the value contain a marker keyword"
    regex pass over the six classes. This OVER-captures: it also matches tagged /
    fusion constructs of studied genes (Flag-AKT1, PTEN-GFP), gene knockdowns whose
    target merely starts with a control token (shNCAM, shNCL -> "NC..."), isotype
    genes that are real study subjects (TUBB3), and unrelated words (tubacin,
    Cfp1). Printed only to show what the refinement removes.

(2) REFINED, ANCHORED CLASSIFICATION (the result). An entity counts as a control
    ONLY if its referent IS the control itself:
      * reporter / tag / housekeeping: the value, after stripping generic
        decorations (leading anti-/vector prefixes Ad-/AAV-/LV-/MV-/oHSV-/Lenti-/
        NP-; trailing ' protein/mRNA/gene/reporter/antibody/staining/complex';
        surrounding '+'/'-'), must EQUAL a known marker. Exact-equality is what
        excludes fusions: 'flag-akt1' != 'flag', 'pten-gfp' != 'gfp'.
      * Cas9 tool: any value containing 'cas9' (the nuclease is always a reagent).
      * non-targeting / scramble control: a sh/si/sg prefix + a control token
        (control|ctrl|scramble(d)|scr|nc|ntc) that is ANCHORED (not a prefix of a
        gene name -- so 'shNC' qualifies but 'shNCAM'/'shNCL' do not).
      * reporter-targeting control: sh/si/sg-GFP / -Luc / -luciferase (knock-down
        of a reporter, used as a negative control).

(3) TWO CONFIDENCE TIERS.
      Tier A -- reporters, tags, Cas9, scramble/NT controls, reporter-targeting
        controls: essentially NEVER the biological subject.
      Tier B -- housekeeping/loading controls (GAPDH, ACTB/beta-actin, tubulin,
        vinculin, lamin B, B2M, HPRT1, TBP, ...): these map to real HGNC genes and
        are dual-use, but in this corpus serve overwhelmingly as normalizers. A
        few strings the text alone cannot fully disambiguate (HA = the tag vs viral
        hemagglutinin; Cas9 in Cas9-engineering papers) are noted as caveats.

(4) ANNOTATION of the HGNC-linkage libraries. Every entry (key) in the seven
    libraries -- greek.json, greek_ambiguous.json, greek_complex.json,
    greek_cosine.json, roman_cosine.json, roman.json, roman_ambiguous.json -- is
    tagged with a nested ("control": "yes") when its key is CASE-SENSITIVELY equal
    to a Tier A/B control entity value, else ("control": "no"). controls.html in
    turn notes, for each Tier A/B entity, "yes"/"no" (whether it occurs as a key in
    any library) and references the specific library file(s) it was found in.
    (This re-annotates the libraries in place; re-run controls.py after roman.py /
    greek.py, which rewrite those files.)

Outputs: GENETIC/controls.html (full enumerated list by class, with occurrences,
the per-entity library match + file references, the by-class summary, the tier
roll-up, and the strategy/caveats); plus the in-place ("control") annotation of the
seven libraries above.
Run from anywhere (paths resolve relative to this file)::

    python controls.py            # uses clean_genetic_ne.tsv if present
    python controls.py --rebuild  # re-aggregate from sentences/*.json
=============================================================================
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
SUM = ROOT / "GENETIC"
TSV = SUM / "clean_genetic_ne.tsv"
HTML_OUT = SUM / "controls.html"
LABEL = "GENETIC"

# HGNC-linkage libraries to cross-reference and annotate with ("control": yes/no)
LIB_FILES = ["greek.json", "greek_ambiguous.json", "greek_complex.json",
             "greek_cosine.json", "roman_cosine.json", "roman.json",
             "roman_ambiguous.json"]

# dash glyph variants -> ASCII '-', then whitespace around '-' removed (roman.py STAGE 1)
DASH_VARIANTS = "-‐‑‒–—―−⁃"
_TO_ASCII = str.maketrans({c: "-" for c in DASH_VARIANTS})
_DESPACE = re.compile(r"\s*-\s*")


def dash_normalize(text):
    return _DESPACE.sub("-", (text or "").translate(_TO_ASCII))


# ============================================================ (0) upstream
def aggregate_from_sentences():
    """Collect GENETIC entities from sentences/*.json, dash-normalize, sum occ.
    Returns (rows: [(value, occ)], info)."""
    counts = Counter()
    n_files = 0
    for fp in sorted(glob.glob(str(SENT_DIR / "*.json"))):
        try:
            rec = json.load(open(fp, encoding="utf-8"))
        except Exception:
            continue
        n_files += 1
        for sent in rec.get("sentences", []):
            for ent in sent.get("entities", []):
                if ent.get("label") == LABEL:
                    counts[ent.get("text", "")] += 1
    agg = Counter()
    for form, n in counts.items():
        agg[dash_normalize(form)] += n
    rows = sorted(agg.items(), key=lambda kv: (-kv[1], kv[0]))
    info = {"source": f"{n_files:,} sentences/*.json files",
            "unique_before": len(counts), "unique_after": len(agg),
            "occ": sum(counts.values())}
    return rows, info


def load_from_tsv():
    rows = []
    for ln in TSV.read_text(encoding="utf-8").splitlines()[1:]:
        p = ln.split("\t")
        if len(p) >= 2 and p[1].isdigit():
            rows.append((p[0], int(p[1])))
    info = {"source": str(TSV.relative_to(ROOT)), "unique_before": None,
            "unique_after": len(rows), "occ": sum(o for _, o in rows)}
    return rows, info


# ============================================================ markers / rules
LEAD = re.compile(r"^(anti-|ad-|aav-|lv-|mv-|ohsv-|lenti-|np-)+", re.I)
TRAIL = re.compile(r"(\s*antibody|\s*reporter gene|\s*reporter|\s*gene|\s*mrna|"
                   r"\s*protein|\s*staining|\s*complex)+$", re.I)


def core(v):
    """Strip generic decorations, casefold -> the bare marker candidate."""
    s = LEAD.sub("", v.strip())
    s = TRAIL.sub("", s)
    return s.strip(" +-").strip().casefold()


REPORTER = {"gfp", "egfp", "eyfp", "yfp", "rfp", "mrfp", "mcherry", "tdtomato",
            "dsred", "venus", "mneongreen", "mscarlet", "mkate", "cerulean",
            "citrine", "luciferase", "luc", "firefly luciferase",
            "renilla luciferase", "renilla", "nanoluc", "lacz", "gus",
            "β-gal", "beta-gal", "β-galactosidase", "beta-galactosidase",
            "sa-β-gal", "sen-β-gal"}
TAG = {"flag", "3xflag", "ha", "his", "6xhis", "his6", "v5", "myc-tag",
       "strep-tag", "t7-tag", "sbp-tag"}
HOUSE = {"gapdh", "actin", "β-actin", "beta-actin", "actb", "α-tubulin",
         "alpha-tubulin", "β-tubulin", "beta-tubulin", "tubulin", "vinculin",
         "vcl", "lamin b", "lamin b1", "laminb1", "lmnb1", "b2m",
         "β-2-microglobulin", "beta-2-microglobulin", "hprt", "hprt1", "18s",
         "18s rrna", "28s", "tbp", "rplp0", "rpl13a", "cyclophilin a", "ppia",
         "36b4"}
# control sh/si/sg suffix, ANCHORED so 'shNC' matches but 'shNCAM'/'shNCL' do not
CTRL = re.compile(r"(?i)(?:^|[^A-Za-z])(sh|si|sg)-?"
                  r"(control|ctrl|scrambled|scramble|scr|ntc|nc)(?![A-Za-z])")
REP_KD = re.compile(r"(?i)^(sh|si|sg)-?(gfp|luc|luciferase)(-\w+)?$")

REPORTER_CLASS = "Reporter / visualization protein"
TAG_CLASS = "Epitope / purification tag"
HOUSE_CLASS = "Housekeeping / loading control"
CAS9_CLASS = "Cas9 nuclease (editing tool)"
NTC_CLASS = "Non-targeting / scramble control"
REPKD_CLASS = "Reporter-targeting control (sh/si/sg-GFP/Luc)"
TIER_A = [REPORTER_CLASS, TAG_CLASS, CAS9_CLASS, NTC_CLASS, REPKD_CLASS]
TIER_B = [HOUSE_CLASS]
ORDER = TIER_A + TIER_B


def classify(v):
    """Return the control class for entity `v`, or None (refined, anchored)."""
    s = v.casefold()
    if "cas9" in s:
        return CAS9_CLASS
    if (CTRL.search(s) or s in ("scramble", "scrambled", "ntc")
            or re.fullmatch(r"(?i)(ntc|nc)\s*sirna", s)):
        return NTC_CLASS
    if REP_KD.match(s):
        return REPKD_CLASS
    c = core(v)
    if c in REPORTER:
        return REPORTER_CLASS
    if c in TAG:
        return TAG_CLASS
    if c in HOUSE:
        return HOUSE_CLASS
    return None


# ---- (1) broad/naive scan, only to demonstrate over-capture in the console ----
_BROAD = {
    REPORTER_CLASS: re.compile(r"(?i)(?<![A-Za-z])(e?GFP|YFP|[ce]FP|m?RFP|mCherry|"
        r"tdTomato|ds[- ]?red|luciferase|\bluc\b|firefly|renilla|lacZ|β-?gal|"
        r"beta-?gal|galactosidase|venus|mscarlet|mneon\w*)(?![A-Za-z])"),
    TAG_CLASS: re.compile(r"(?i)(?<![A-Za-z])(3?x?FLAG|HA|6?x?His|V5|Myc-?tag|"
        r"Strep-?tag|T7-?tag)(?![A-Za-z])"),
    HOUSE_CLASS: re.compile(r"(?i)(GAPDH|β-?actin|beta-?actin|ACTB|actin|tubulin|"
        r"TUBB\w*|TUBA\w*|vinculin|VCL|lamin\s?B\d?|LMNB\d?|B2M|HPRT1?|18S|TBP|"
        r"cyclophilin|PPIA)"),
    CAS9_CLASS: re.compile(r"(?i)cas9"),
    NTC_CLASS: re.compile(r"(?i)(sh|si|sg)-?(control|ctrl|nc|ntc|scr|scramble)"),
}


def broad_counts(rows):
    out = {}
    for cls, rx in _BROAD.items():
        hits = [(v, o) for v, o in rows if rx.search(v)]
        out[cls] = (len(hits), sum(o for _, o in hits))
    return out


# entities that the broad scan wrongly grabs but the refined pass MUST exclude
EXCLUSION_CHECKS = ["Flag-AKT1-WT", "PTENT277A-GFP", "TUBB3-V5", "Ha-ras",
                    "Ccr2RFP", "Il34LacZ", "shNCAM", "shNCL", "TUBB3", "tubacin",
                    "Cfp1", "MOVCAR-shNCAM"]


# ============================================================ (4) annotation
def annotate_libraries(control_values):
    """Add a nested ("control": "yes"/"no") to every entry in each library:
    "yes" iff the entry's key is case-sensitively equal to a Tier A/B control
    entity value, else "no". Rewrites each file in place (indent=2, ensure_ascii=
    False) and returns (file_keys: f->set(keys), counts: f->(n_yes, n_no))."""
    file_keys, counts = {}, {}
    for f in LIB_FILES:
        path = SUM / f
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        yes = 0
        for k, entry in data.items():
            if isinstance(entry, dict):
                flag = "yes" if k in control_values else "no"
                entry["control"] = flag
                yes += (flag == "yes")
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
        file_keys[f] = set(data)
        counts[f] = (yes, len(data) - yes)
    return file_keys, counts


# ============================================================ reporting
def run(rows, info):
    cat = defaultdict(list)
    for v, o in rows:
        k = classify(v)
        if k:
            cat[k].append((v, o))
    for k in cat:
        cat[k].sort(key=lambda x: (-x[1], x[0].casefold()))

    grand = {v for k in ORDER for v, _ in cat[k]}
    gocc = sum(o for k in ORDER for _, o in cat[k])
    nA = sum(len(cat[k]) for k in TIER_A)
    oA = sum(o for k in TIER_A for _, o in cat[k])
    nB = sum(len(cat[k]) for k in TIER_B)
    oB = sum(o for k in TIER_B for _, o in cat[k])
    corpus_n = len(rows)
    corpus_occ = info["occ"]

    # ---- console: reproduce the analytical steps -------------------------
    print("=" * 70)
    print("CONTROL-MARKER ANALYSIS")
    print("=" * 70)
    print(f"(0) UPSTREAM corpus: source = {info['source']}")
    if info["unique_before"]:
        print(f"      unique surfaces {info['unique_before']:,} -> "
              f"{info['unique_after']:,} after dash-normalization")
    print(f"      {corpus_n:,} distinct GENETIC entities, {corpus_occ:,} occurrences")

    print("\n(1) BROAD keyword scan (naive -- OVER-captures, shown for contrast):")
    for cls in _BROAD:
        n, o = broad_counts(rows)[cls]
        print(f"      {cls:42} {n:>5} entities, {o:>6} occ")

    print("\n(2) REFINED anchored classification (the result):")
    for k in ORDER:
        print(f"      {k:42} {len(cat[k]):>5} entities, "
              f"{sum(o for _, o in cat[k]):>6} occ")

    print("\n    exclusion checks (must NOT be classified as controls):")
    for ex in EXCLUSION_CHECKS:
        print(f"      {ex:18} -> control? {classify(ex) is not None}")

    print(f"\n(3) TOTALS: {len(grand):,} distinct control entities, {gocc:,} occ"
          f"  ({100*len(grand)/corpus_n:.1f}% of entities, "
          f"{100*gocc/corpus_occ:.1f}% of occ)")
    print(f"      Tier A (pure tools/reporters/controls): {nA:,} entities, {oA:,} occ")
    print(f"      Tier B (dual-use housekeeping)        : {nB:,} entities, {oB:,} occ")

    # ---- (4) annotate the libraries + cross-reference each control entity ----
    control_values = grand
    file_keys, ann = annotate_libraries(control_values)
    refs = {v: [f for f in LIB_FILES if f in file_keys and v in file_keys[f]]
            for v in control_values}
    n_yes = sum(1 for v in control_values if refs[v])
    print("\n(4) ANNOTATION of the seven libraries (added key \"control\": yes/no):")
    for f in LIB_FILES:
        if f in ann:
            y, n = ann[f]
            print(f"      {f:22} control=yes {y:>4}   control=no {n:>6}")
    print(f"      control entities present as a library key: "
          f"yes {n_yes}/{len(control_values)}, no {len(control_values) - n_yes}")

    write_html(cat, grand, gocc, nA, oA, nB, oB, corpus_n, corpus_occ, refs, n_yes)
    print(f"\nWrote {HTML_OUT}")


def write_html(cat, grand, gocc, nA, oA, nB, oB, corpus_n, corpus_occ, refs, n_yes):
    esc = html.escape

    def section(k, dual=False):
        items = cat[k]
        occ = sum(o for _, o in items)
        crows = []
        for v, o in items:
            rf = refs.get(v, [])
            flag = "yes" if rf else "no"
            cls = "ctl-yes" if rf else "ctl-no"
            files = (" ".join(f'<code>{esc(f)}</code>' for f in rf) if rf else "&mdash;")
            crows.append(f'<tr><td><code>{esc(v)}</code></td>'
                         f'<td class="num">{o:,}</td>'
                         f'<td class="{cls}">{flag}</td><td>{files}</td></tr>')
        cells = "".join(crows)
        note = (' <span class="dual">dual-use &mdash; map to real HGNC genes; '
                'flagged because in this corpus they serve overwhelmingly as '
                'normalizers</span>' if dual else "")
        return (f'<h3>{esc(k)} &mdash; {len(items)} entities, {occ:,} occ{note}</h3>'
                f'<table><tr><th>entity</th><th class="num">occ</th>'
                f'<th>control&nbsp;key?</th><th>library file(s)</th></tr>'
                f'{cells}</table>')

    summary_rows = "".join(
        f'<tr><td>{esc(k)}</td><td class="num">{len(cat[k]):,}</td>'
        f'<td class="num">{sum(o for _, o in cat[k]):,}</td></tr>' for k in ORDER)

    style = (
        " body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:960px;"
        "margin:2rem auto;padding:0 1rem;color:#1a1a1a;}"
        " h1{font-size:1.45rem;} h2{font-size:1.15rem;margin-top:1.8rem;"
        "border-bottom:1px solid #ddd;padding-bottom:.3rem;}"
        " h3{font-size:1rem;margin-top:1.3rem;}"
        " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;}"
        " th,td{border:1px solid #ccc;padding:.28rem .6rem;text-align:left;}"
        " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
        " .big{font-size:2rem;font-weight:700;}"
        " .headline{background:#fff7ed;border:1px solid #f0d9b5;border-radius:8px;"
        "padding:.8rem 1rem;margin:1rem 0;}"
        " code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.9em;}"
        " .dual{color:#9a6700;font-size:.85em;} p.note{color:#444;font-size:.92em;}"
        " .ctl-yes{color:#1a7f37;font-weight:700;text-align:center;}"
        " .ctl-no{color:#999;text-align:center;}")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>GENETIC entities that are experimental controls / tools</title>
<style>{style}</style></head><body>
<h1>GENETIC entities used as experimental controls / tools</h1>
<p>Entities from <code>sentences/*.json</code> whose referent is a
<strong>control or laboratory tool</strong> &mdash; a reporter, an epitope tag, a
housekeeping/loading normalizer, the Cas9 nuclease, or a non-targeting/scramble
control &mdash; <em>not</em> a gene/protein participating in the biology or disease
under study. Produced by <code>controls.py</code>.</p>
<div class="headline"><span class="big">{len(grand):,}</span> distinct entities
&nbsp;&middot;&nbsp; <span class="big">{gocc:,}</span> occurrences &mdash;
{100*len(grand)/corpus_n:.1f}% of the {corpus_n:,} distinct GENETIC entities
({100*gocc/corpus_occ:.1f}% of {corpus_occ:,} occurrences). Of these,
<strong>{n_yes}</strong> occur as a key in one of the seven HGNC-linkage libraries
(marked <span class="ctl-yes">yes</span> below and annotated
<code>"control":"yes"</code> in that file); the remaining
{len(grand) - n_yes} are <span class="ctl-no">no</span> (unmatched tools/reagents
absent from the libraries).</div>
<p class="note"><strong>Annotation.</strong> Every entry in
<code>greek.json</code>, <code>greek_ambiguous.json</code>,
<code>greek_complex.json</code>, <code>greek_cosine.json</code>,
<code>roman_cosine.json</code>, <code>roman.json</code> and
<code>roman_ambiguous.json</code> now carries a nested
<code>"control"</code> field &mdash; <code>"yes"</code> when its key is
case-sensitively equal to a Tier&nbsp;A/B control entity below, else
<code>"no"</code>. The <em>control&nbsp;key?</em> and <em>library file(s)</em>
columns below report the reverse lookup for each control entity.</p>

<h2>Summary by class</h2>
<table><tr><th>class</th><th class="num">entities</th><th class="num">occ</th></tr>
{summary_rows}
<tr><td><strong>Tier A &mdash; pure tools/reporters/controls</strong></td><td class="num"><strong>{nA:,}</strong></td><td class="num"><strong>{oA:,}</strong></td></tr>
<tr><td><strong>Tier B &mdash; dual-use housekeeping</strong></td><td class="num"><strong>{nB:,}</strong></td><td class="num"><strong>{oB:,}</strong></td></tr>
<tr><td><strong>Total</strong></td><td class="num"><strong>{len(grand):,}</strong></td><td class="num"><strong>{gocc:,}</strong></td></tr>
</table>

<h2>Method &amp; caveats</h2>
<p class="note">Each entity is anchored to a known marker after stripping generic
decorations (<code>anti-</code>, vector prefixes <code>Ad-/AAV-/LV-/MV-/oHSV-/Lenti-/NP-</code>,
trailing <code> protein/mRNA/gene/reporter/antibody</code> and <code>+/-</code>), so
<strong>tagged / fusion constructs of studied genes are excluded</strong>
(<code>Flag-AKT1</code>, <code>PTEN-GFP</code>, <code>TUBB3-V5</code>,
<code>Ha-ras</code>=HRAS, reporter knock-ins <code>Ccr2RFP</code>/<code>Il34LacZ</code>);
<strong>gene knockdowns are excluded</strong> (<code>shNCAM</code>, <code>shNCL</code>);
and isotype genes that are genuine study subjects (<code>TUBB3</code>) are excluded.
<strong>Tier A</strong> markers are essentially never the biological subject;
<strong>Tier B</strong> (housekeeping/loading) are real genes used here as
normalizers &mdash; a judgment call, and a few strings (<code>HA</code> = viral
hemagglutinin; Cas9-engineering studies) the entity text alone cannot fully
disambiguate.</p>

<h2>Tier A &mdash; pure tools, reporters &amp; controls ({nA} entities, {oA:,} occ)</h2>
{''.join(section(k) for k in TIER_A)}

<h2>Tier B &mdash; housekeeping / loading controls ({nB} entities, {oB:,} occ)</h2>
{''.join(section(k, dual=True) for k in TIER_B)}
</body></html>
"""
    SUM.mkdir(parents=True, exist_ok=True)
    HTML_OUT.write_text(doc, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rebuild", action="store_true",
                    help="re-aggregate GENETIC entities from sentences/*.json "
                         "instead of reading GENETIC/clean_genetic_ne.tsv")
    args = ap.parse_args()
    if args.rebuild or not TSV.exists():
        rows, info = aggregate_from_sentences()
    else:
        rows, info = load_from_tsv()
    run(rows, info)


if __name__ == "__main__":
    main()
