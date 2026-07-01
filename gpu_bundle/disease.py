#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""disease.py -- normalize DISEASE named entities to the MONDO disease ontology.

The DISEASE analog of the GENETIC pipeline (roman.py/greek.py link gene surfaces
to HGNC; this links disease surfaces to MONDO). Reads the BioBERT DISEASE spans
from sentences/*.json and maps each to a MONDO term by exact equality of a
transformed key against MONDO labels + synonyms.

Inputs
------
  sentences/*.json                  BioBERT NER output (DISEASE entities)
  databases/mondo-clingen.json      MONDO ontology (OBO-Graph JSON: graphs[].nodes[])

Outputs (under DISEASE/)
-------
  clean_disease_ne.tsv              aggregated DISEASE surfaces + occurrences
  disease.json                      surfaces linked to exactly ONE MONDO term
  disease_ambiguous.json            surfaces linked to >=2 MONDO terms
  unmatched_disease.json            surfaces linked to no MONDO term
  disease.html                      self-documenting strategy + results report,
                                    incl. a per-file nested-key summary (disease.json
                                    / disease_ambiguous.json) and the distinct MONDO
                                    terms (mondo_id -> mondo_label + summed
                                    occurrences) for disease.json alone AND
                                    single + ambiguous

================================ STRATEGY ===================================
(0) UPSTREAM. Collect every entity labelled DISEASE from sentences/*.json,
    dash-normalize (unify dash glyphs to '-', strip whitespace around '-'), and
    sum occurrences per normalized surface (the same STAGE 1 roman.py performs).

(0b) MONDO INDEXING. For each non-deprecated MONDO_* node, index its label and
    every synonym (hasExact/Related/Broad/NarrowSynonym) under four key
    normalizations: literal, casefold, hyphen<->space fold, separator-deletion.
    Each surface key -> the set of (MONDO id, label) it denotes; the matching
    field is recorded with priority label > exact > related > narrow > broad.

(1) MATCH CASCADE, tried in order; the first pass that hits wins (match_mode):
      1. curated synonym       hand map for corpus abbreviations MONDO lacks
                               (GBM->glioblastoma, GSC-> ... ). str->single term.
      2. case-sensitive        exact equality against label/synonym surfaces
      3. case-insensitive      casefold equality
      4. hyphen/whitespace     casefold + '-'<->space
      5. normalized            strip plural (tumors->tumor) and British->American
                               spelling (tumour->tumor, leukaemia->leukemia), then
                               re-match case-insensitively
      6. separator-deletion    casefold + drop hyphens/spaces

(2) ROUTING by distinct MONDO-term count of the winning pass:
      1 term  -> disease.json (single)        2+ -> disease_ambiguous.json
      0       -> unmatched_disease.json
=============================================================================

Run from anywhere (paths resolve relative to this file)::

    python disease.py
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
MONDO_PATH = ROOT / "databases" / "mondo-clingen.json"
LABEL = "DISEASE"

CLEAN_TSV = OUT_DIR / "clean_disease_ne.tsv"
SINGLE_OUT = OUT_DIR / "disease.json"
AMBIG_OUT = OUT_DIR / "disease_ambiguous.json"
UNMATCHED_OUT = OUT_DIR / "unmatched_disease.json"
HTML_OUT = OUT_DIR / "disease.html"

# synonym predicate -> short field name + priority (lower = preferred)
_SYN_FIELD = {"hasExactSynonym": "exact", "hasRelatedSynonym": "related",
              "hasNarrowSynonym": "narrow", "hasBroadSynonym": "broad"}
FIELD_RANK = {"label": 0, "exact": 1, "related": 2, "narrow": 3, "broad": 4}

DASH_VARIANTS = "-‐‑‒–—―−⁃"
_TO_ASCII = str.maketrans({c: "-" for c in DASH_VARIANTS})
_DESPACE = re.compile(r"\s*-\s*")


def dash_normalize(text):
    return _DESPACE.sub("-", (text or "").translate(_TO_ASCII))


def hsfold(s):
    return s.casefold().replace("-", " ")


def delsep(s):
    return s.casefold().replace("-", "").replace(" ", "")


# British -> American spelling folds (word-internal), applied before re-matching
_UK_US = [(re.compile(p), r) for p, r in [
    (r"tumour", "tumor"), (r"oesophag", "esophag"), (r"oedema", "edema"),
    (r"haemorrh", "hemorrh"), (r"haemat", "hemat"), (r"haemangi", "hemangi"),
    (r"haem", "hem"), (r"anaemi", "anemi"), (r"leukaemi", "leukemi"),
    (r"ischaemi", "ischemi"), (r"paediatr", "pediatr"), (r"coeliac", "celiac"),
    (r"oestrogen", "estrogen"), (r"caemia", "cemia"), (r"aemia", "emia")]]
_IRREGULAR = {"metastases": "metastasis", "foci": "focus",
              "carcinomata": "carcinoma", "adenomata": "adenoma",
              "lymphomata": "lymphoma", "sarcomata": "sarcoma",
              "neoplasiae": "neoplasia"}


def britishize(s):
    for rx, rep in _UK_US:
        s = rx.sub(rep, s)
    return s


def singular_candidates(s):
    """Candidate singular forms of a (casefolded) plural surface."""
    out = []
    if s in _IRREGULAR:
        out.append(_IRREGULAR[s])
    if s.endswith("ies") and len(s) > 4:
        out.append(s[:-3] + "y")
    if s.endswith("es") and len(s) > 4:
        out.append(s[:-2])
    if s.endswith("s") and not s.endswith("ss") and len(s) > 3:
        out.append(s[:-1])
    return out


# ---- corpus abbreviations MONDO does not carry as a synonym ----
# surface -> either a MONDO id ("MONDO:xxxxxxx", pinned to that single term) or a
# MONDO label/synonym string (matched case-insensitively). Most disease
# abbreviations (GBM, HCC, NSCLC, CRC, ...) are already MONDO synonyms and link via
# the case passes; these are high-frequency ones MONDO lacks, read in this
# (neuro-oncology-heavy) corpus's sense. Mapped by id to pin a single term
# (e.g. the bare label 'medulloblastoma' is also a synonym of the adult/childhood
# variants, so a label lookup would be ambiguous).
CURATED = {
    "MB":   "MONDO:0007959",   # medulloblastoma (corpus is neuro-oncology)
    "GB":   "MONDO:0018177",   # glioblastoma (cf. GBM)
    "LGG":  "MONDO:0021637",   # low grade glioma
    "TNBC": "MONDO:0005494",   # triple-negative breast carcinoma
}


# ============================================================ (0) upstream
def stage1():
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
    forms = defaultdict(set)
    for form, n in counts.items():
        key = dash_normalize(form)
        agg[key] += n
        forms[key].add(form)
    rows = sorted(agg.items(), key=lambda kv: (-kv[1], kv[0]))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"{k}\t{n}\t{len(forms[k])}" for k, n in rows)
    CLEAN_TSV.write_text("clean_disease_ne\toccurrences\tn_source_forms\n"
                         + body + "\n", encoding="utf-8")
    info = {"n_files": n_files, "occ": sum(counts.values()),
            "unique_before": len(counts), "unique_after": len(agg)}
    return rows, info


# ============================================================ (0b) MONDO index
def mondo_curie(node_id):
    return node_id.rsplit("/", 1)[-1].replace("_", ":")


def build_indexes(mondo_path):
    docs = json.loads(Path(mondo_path).read_text(encoding="utf-8"))
    nodes = docs["graphs"][0]["nodes"]
    idx, idx_ci, idx_hs, idx_del = ({}, {}, {}, {})
    label_of = {}
    n_terms = 0

    def add(index, key, mid, field):
        slot = index.setdefault(key, {})
        # keep the highest-priority field per (key, term)
        if mid not in slot or FIELD_RANK[field] < FIELD_RANK[slot[mid]]:
            slot[mid] = field

    for nd in nodes:
        nid = nd.get("id", "")
        if "MONDO_" not in nid:
            continue
        meta = nd.get("meta", {})
        if meta.get("deprecated"):
            continue
        lbl = nd.get("lbl")
        if not lbl:
            continue
        mid = mondo_curie(nid)
        label_of[mid] = lbl
        n_terms += 1
        surfaces = [(lbl, "label")]
        for syn in meta.get("synonyms", []):
            val = syn.get("val")
            fld = _SYN_FIELD.get(syn.get("pred"))
            if val and fld:
                surfaces.append((val, fld))
        for s, field in surfaces:
            add(idx, s, mid, field)
            add(idx_ci, s.casefold(), mid, field)
            add(idx_hs, hsfold(s), mid, field)
            add(idx_del, delsep(s), mid, field)
    return {"idx": idx, "ci": idx_ci, "hs": idx_hs, "del": idx_del,
            "label_of": label_of, "n_terms": n_terms}


def best_field(slot):
    return min((f for f in slot.values()), key=lambda f: FIELD_RANK[f])


# ============================================================ (1) cascade
def match_cascade(value, IX):
    def hit(key, index):
        return index.get(key)

    cur = CURATED.get(value)
    if cur is not None:
        if cur.startswith("MONDO:"):              # pinned directly to one term
            if cur in IX["label_of"]:
                return {cur}, "curated", "curated synonym"
        else:                                     # a label/synonym string
            slot = IX["ci"].get(cur.casefold())
            if slot:
                return set(slot), best_field(slot), "curated synonym"
    slot = hit(value, IX["idx"])
    if slot:
        return set(slot), best_field(slot), "case-sensitive"
    slot = hit(value.casefold(), IX["ci"])
    if slot:
        return set(slot), best_field(slot), "case-insensitive"
    slot = hit(hsfold(value), IX["hs"])
    if slot:
        return set(slot), best_field(slot), "hyphen/whitespace"
    # normalized: British->American + plural-strip, matched case-insensitively
    base = britishize(value.casefold())
    cands = [base] + singular_candidates(base)
    if base != value.casefold():
        for c in [value.casefold()] + singular_candidates(value.casefold()):
            if c not in cands:
                cands.append(c)
    for c in cands:
        slot = hit(c, IX["ci"])
        if slot:
            return set(slot), best_field(slot), "normalized"
    slot = hit(delsep(value), IX["del"])
    if slot:
        return set(slot), best_field(slot), "separator-deletion"
    return set(), None, None


_MODES = ["curated synonym", "case-sensitive", "case-insensitive",
          "hyphen/whitespace", "normalized", "separator-deletion"]


def mondo_term_stats(single, ambiguous=None):
    """Aggregate the linked libraries by MONDO TERM: distinct mondo_id -> label and
    summed occurrences (several surfaces can share one term, e.g. GBM and
    glioblastoma). A single entry credits its one term; an ambiguous entry credits
    EACH of its candidate terms the entry's full occurrence count. Returns
    (ranked [(mondo_id, label, occ)] desc, total_occ)."""
    occ_by_id = defaultdict(int)
    label_of = {}
    for e in single.values():
        occ_by_id[e["mondo_id"]] += e["occurrences"]
        label_of[e["mondo_id"]] = e["mondo_label"]
    for e in (ambiguous or {}).values():
        for mid, lbl in zip(e["mondo_id"], e["mondo_label"]):
            occ_by_id[mid] += e["occurrences"]
            label_of[mid] = lbl
    ranked = sorted(((mid, label_of[mid], occ) for mid, occ in occ_by_id.items()),
                    key=lambda r: (-r[2], r[0]))
    return ranked, sum(occ_by_id.values())


# ============================================================ driver
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mondo", default=str(MONDO_PATH))
    ap.parse_args()

    rows, info = stage1()
    print(f"STAGE 1  files={info['n_files']:,}  DISEASE occ={info['occ']:,}  "
          f"unique {info['unique_before']:,}->{info['unique_after']:,} "
          f"-> {CLEAN_TSV}")

    IX = build_indexes(MONDO_PATH)
    print(f"MONDO    {IX['n_terms']:,} non-deprecated terms indexed "
          f"(label + synonyms, 4 key folds)")

    single, ambiguous, unmatched = {}, {}, {}
    mode_counts = Counter()
    occ = {"single": 0, "ambiguous": 0, "unmatched": 0}
    for value, n in rows:
        ids, field, mode = match_cascade(value, IX)
        if not ids:
            unmatched[value] = {"occurrences": n}
            occ["unmatched"] += n
            continue
        mids = sorted(ids)
        labels = [IX["label_of"][m] for m in mids]
        base = {"occurrences": n, "match_field": field, "match_mode": mode}
        if len(mids) == 1:
            single[value] = {**base, "mondo_id": mids[0], "mondo_label": labels[0]}
            occ["single"] += n
        else:
            ambiguous[value] = {**base, "mondo_id": mids, "mondo_label": labels,
                                "n_terms": len(mids)}
            occ["ambiguous"] += n
        mode_counts[mode] += 1

    for path, lib in ((SINGLE_OUT, single), (AMBIG_OUT, ambiguous),
                      (UNMATCHED_OUT, unmatched)):
        path.write_text(json.dumps(lib, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")

    tot = len(single) + len(ambiguous) + len(unmatched)
    print(f"\nMATCH    single -> {SINGLE_OUT.name}: {len(single):,} "
          f"({occ['single']:,} occ)")
    print(f"         ambiguous -> {AMBIG_OUT.name}: {len(ambiguous):,} "
          f"({occ['ambiguous']:,} occ)")
    print(f"         unmatched -> {UNMATCHED_OUT.name}: {len(unmatched):,} "
          f"({occ['unmatched']:,} occ)")
    print(f"         partition: {len(single)} + {len(ambiguous)} + "
          f"{len(unmatched)} = {tot:,} (of {len(rows):,} surfaces)")
    print(f"         by mode: " + ", ".join(f"{m} {mode_counts[m]}" for m in _MODES))
    linked_occ = occ["single"] + occ["ambiguous"]
    print(f"         linked occ: {linked_occ:,} / {info['occ']:,} "
          f"({100*linked_occ/info['occ']:.1f}%)")
    sterm_rank, sterm_occ = mondo_term_stats(single)
    print(f"         distinct MONDO terms (disease.json, single): {len(sterm_rank):,}; "
          f"summed occurrences: {sterm_occ:,}")
    term_rank, term_occ = mondo_term_stats(single, ambiguous)
    print(f"         distinct MONDO terms (single+ambiguous): {len(term_rank):,}; "
          f"summed term occurrences: {term_occ:,}")

    render_html(rows, info, IX, single, ambiguous, unmatched, occ, mode_counts)
    print(f"\nWrote {HTML_OUT}")


def _typename(v):
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, int):
        return "int"
    if isinstance(v, float):
        return "float"
    if isinstance(v, str):
        return "str"
    if isinstance(v, list):
        inner = sorted({_typename(x) for x in v}) or ["empty"]
        return f"list[{'|'.join(inner)}]"
    if isinstance(v, dict):
        return "object"
    if v is None:
        return "null"
    return type(v).__name__


def _short(v, n=60):
    s = json.dumps(v, ensure_ascii=False)
    return s if len(s) <= n else s[:n - 1] + "…"


def render_nested_keys(named_libs):
    """HTML summary of the nested keys per library: a cross-file presence matrix
    plus, for each file, every key with its coverage, value type(s) and an example.
    named_libs = [(filename, lib_dict), ...]."""
    esc = html.escape
    perfile = []
    allkeys = []
    for fname, d in named_libs:
        present, types, example, order = Counter(), defaultdict(set), {}, []
        for e in d.values():
            if not isinstance(e, dict):
                continue
            for k, v in e.items():
                if k not in present:
                    order.append(k)
                present[k] += 1
                types[k].add(_typename(v))
                example.setdefault(k, v)
        perfile.append((fname, len(d), order, present, types, example))
        for k in order:
            if k not in allkeys:
                allkeys.append(k)
    hdr = "".join(f'<th>{esc(f)}</th>' for f, _, _, _, _, _ in perfile)
    matrix = "".join(
        '<tr><td><code>' + esc(k) + '</code></td>' + "".join(
            '<td class="yes">&#10003;</td>' if k in pf[2] else '<td class="no"></td>'
            for pf in perfile) + '</tr>' for k in allkeys)
    parts = ["<h2>Nested keys (per file)</h2>",
             "<p>The nested keys inside each entry of <code>disease.json</code> and "
             "<code>disease_ambiguous.json</code> &mdash; coverage, value type(s) and "
             "an example. (<code>phenotypes.py</code> later adds a "
             "<code>\"phenotype\"</code> key to both files.)</p>",
             f'<table><tr><th>nested key</th>{hdr}</tr>{matrix}</table>']
    for fname, n, order, present, types, example in perfile:
        rows = "".join(
            f'<tr><td><code>{esc(k)}</code></td>'
            f'<td class="num">{present[k]:,}/{n:,}</td>'
            f'<td><code>{esc(" | ".join(sorted(types[k])))}</code></td>'
            f'<td><code>{esc(_short(example[k]))}</code></td></tr>' for k in order)
        parts.append(
            f'<h3><code>{esc(fname)}</code> &mdash; {n:,} entries, {len(order)} '
            f'nested key(s)</h3><table><tr><th>nested key</th>'
            f'<th class="num">present</th><th>type(s)</th><th>example</th></tr>'
            f'{rows}</table>')
    return "\n".join(parts)


# ============================================================ report
def render_html(rows, info, IX, single, ambiguous, unmatched, occ, mode_counts):
    esc = html.escape
    n_surf = len(rows)
    linked = len(single) + len(ambiguous)
    linked_occ = occ["single"] + occ["ambiguous"]

    def occ_of(d):
        return sorted(d.items(), key=lambda kv: -kv[1]["occurrences"])

    top_single = occ_of(single)[:25]
    top_amb = occ_of(ambiguous)[:20]
    top_un = occ_of(unmatched)[:30]

    srows = "".join(
        f'<tr><td><code>{esc(v)}</code></td><td class="num">{e["occurrences"]:,}</td>'
        f'<td><code>{esc(e["mondo_id"])}</code></td><td>{esc(e["mondo_label"])}</td>'
        f'<td><code>{esc(e["match_mode"])}</code>/{esc(e["match_field"])}</td></tr>'
        for v, e in top_single)
    arows = "".join(
        f'<tr><td><code>{esc(v)}</code></td><td class="num">{e["occurrences"]:,}</td>'
        f'<td class="num">{e["n_terms"]}</td>'
        f'<td>{esc(", ".join(e["mondo_label"][:4]))}'
        f'{"&hellip;" if len(e["mondo_label"])>4 else ""}</td></tr>'
        for v, e in top_amb)
    urows = "".join(
        f'<tr><td><code>{esc(v)}</code></td><td class="num">{e["occurrences"]:,}</td></tr>'
        for v, e in top_un)
    mrows = "".join(
        f'<tr><td><code>{esc(m)}</code></td><td class="num">{mode_counts[m]:,}</td></tr>'
        for m in _MODES if mode_counts[m])
    # distinct MONDO terms in disease.json (single matches): id -> label, summed occ
    sterm_rank, sterm_occ = mondo_term_stats(single)
    strows = "".join(
        f'<tr><td><code>{esc(mid)}</code></td><td>{esc(lbl)}</td>'
        f'<td class="num">{o:,}</td></tr>' for mid, lbl, o in sterm_rank)
    # distinct MONDO terms across single + ambiguous
    term_rank, term_occ = mondo_term_stats(single, ambiguous)
    trows = "".join(
        f'<tr><td><code>{esc(mid)}</code></td><td>{esc(lbl)}</td>'
        f'<td class="num">{o:,}</td></tr>' for mid, lbl, o in term_rank)
    nested_section = render_nested_keys([("disease.json", single),
                                         ("disease_ambiguous.json", ambiguous)])

    style = (
        " body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;"
        "margin:2rem auto;padding:0 1rem;color:#1a1a1a;}"
        " h1{font-size:1.45rem;} h2{font-size:1.15rem;margin-top:1.8rem;"
        "border-bottom:1px solid #ddd;padding-bottom:.3rem;}"
        " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;}"
        " th,td{border:1px solid #ccc;padding:.28rem .6rem;text-align:left;vertical-align:top;}"
        " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
        " .big{font-size:2rem;font-weight:700;}"
        " .headline{background:#f1f6fb;border:1px solid #cdddf0;border-radius:8px;"
        "padding:.8rem 1rem;margin:1rem 0;}"
        " code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.9em;}"
        " .yes{text-align:center;color:#1a7f37;font-weight:700;} .no{text-align:center;background:#fafafa;}"
        " details{margin:.6rem 0;} summary{cursor:pointer;color:#357;}")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>disease &mdash; DISEASE entity normalization to MONDO</title>
<style>{style}</style></head><body>
<h1>disease &mdash; DISEASE entity &rarr; MONDO normalization</h1>
<p>The DISEASE analog of the GENETIC pipeline: BioBERT DISEASE spans from
<code>sentences/*.json</code> linked to the <strong>MONDO</strong> disease
ontology by exact equality of a transformed key against MONDO labels and
synonyms. Produced by <code>disease.py</code>.</p>
<div class="headline"><span class="big">{linked:,}</span> of {n_surf:,} distinct
DISEASE surfaces link to MONDO ({100*linked/n_surf:.1f}%);
<span class="big">{100*linked_occ/info['occ']:.1f}%</span> of the
{info['occ']:,} DISEASE occurrences are linked.
Indexed against {IX['n_terms']:,} MONDO terms.</div>

<h2>Buckets</h2>
<table><tr><th>bucket</th><th class="num">surfaces</th><th class="num">occurrences</th></tr>
<tr><td>single MONDO term &rarr; <code>disease.json</code></td><td class="num">{len(single):,}</td><td class="num">{occ['single']:,}</td></tr>
<tr><td>ambiguous (&ge;2 terms) &rarr; <code>disease_ambiguous.json</code></td><td class="num">{len(ambiguous):,}</td><td class="num">{occ['ambiguous']:,}</td></tr>
<tr><td>unmatched &rarr; <code>unmatched_disease.json</code></td><td class="num">{len(unmatched):,}</td><td class="num">{occ['unmatched']:,}</td></tr>
<tr><td><strong>total</strong></td><td class="num"><strong>{n_surf:,}</strong></td><td class="num"><strong>{info['occ']:,}</strong></td></tr>
</table>

<h2>Yield by matching pass</h2>
<table><tr><th>match_mode</th><th class="num">surfaces</th></tr>{mrows}</table>

<h2>Top single-term matches</h2>
<table><tr><th>surface</th><th class="num">occ</th><th>MONDO id</th><th>MONDO label</th><th>via</th></tr>
{srows}</table>

<h2>Top ambiguous matches</h2>
<table><tr><th>surface</th><th class="num">occ</th><th class="num">#terms</th><th>candidate labels</th></tr>
{arows}</table>

<h2>Top unmatched surfaces</h2>
<p>High-frequency DISEASE surfaces with no MONDO label/synonym match (generic
mentions like <code>tumor</code>-class terms, abbreviations MONDO lacks, or
phrases) &mdash; candidates for curation.</p>
<table><tr><th>surface</th><th class="num">occ</th></tr>{urows}</table>

{nested_section}

<h2>Distinct MONDO terms in <code>disease.json</code> (single matches)</h2>
<p>The distinct <code>mondo_id</code> values in <code>disease.json</code> (the
single-term matches), each mapped to its <code>mondo_label</code> and the sum of
the <code>occurrences</code> of every surface that resolves to it (several surfaces
can share one term, e.g. <code>GBM</code> and <code>glioblastoma</code>).</p>
<div class="headline"><span class="big">{len(sterm_rank):,}</span> distinct MONDO
terms &nbsp;&middot;&nbsp; <span class="big">{sterm_occ:,}</span> summed occurrences.</div>
<details><summary>show / hide all {len(sterm_rank):,} MONDO terms (by &Sigma; occ)</summary>
<table><tr><th>mondo_id</th><th>mondo_label</th><th class="num">&Sigma; occ</th></tr>
{strows}</table></details>

<h2>Distinct MONDO terms linked (single + ambiguous)</h2>
<p>The distinct <code>mondo_id</code> values across <code>disease.json</code> and
<code>disease_ambiguous.json</code>, each mapped to its <code>mondo_label</code> and
the sum of the <code>occurrences</code> of every surface that resolves to it (an
ambiguous surface credits its full count to <em>each</em> candidate term).</p>
<div class="headline"><span class="big">{len(term_rank):,}</span> distinct MONDO
terms &nbsp;&middot;&nbsp; <span class="big">{term_occ:,}</span> summed occurrences
(term-attributed).</div>
<details><summary>show / hide all {len(term_rank):,} MONDO terms (by &Sigma; occ)</summary>
<table><tr><th>mondo_id</th><th>mondo_label</th><th class="num">&Sigma; occ</th></tr>
{trows}</table></details>

<h2>Reproduce</h2><pre><code>python disease.py</code></pre>
</body></html>
"""
    HTML_OUT.write_text(doc, encoding="utf-8")


if __name__ == "__main__":
    main()
