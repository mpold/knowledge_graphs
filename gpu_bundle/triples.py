#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""triples.py -- extract relation triples from the BioBERT NER sentences.

The sentences/*.json files are BioBERT output (NER models
biobert_genetic/disease/chemical): each sentence carries typed entity spans
(GENETIC / DISEASE / CHEMICAL) with character offsets and the sentence text.
BioBERT did the entity recognition; this script derives (subject, predicate,
object) TRIPLES from the co-occurring entities: for each pair of CONSECUTIVE
entities in a sentence (sorted by position, so nothing else lies between them) the
connecting text is the candidate predicate, and the pair is kept as a triple only
when that text carries a relation cue (a verb / relational stem -- induced,
inhibited, expression, associated, -mediated, ...), which filters out pure
coordination ("and", "or", "in"). No separate relation-extraction model is run
(none is available here); the predicate is the surface text between two entities.
GENETIC entity surfaces are hyphen-normalized in the output (dash-glyph variants ->
ASCII '-', whitespace around '-' stripped); DISEASE/CHEMICAL text is kept verbatim.

Input  : sentences/*.json
Output : TRIPLES/triples.json   list of triples, each carrying the typed elements,
           PMID and sentence:
             {subject:{text,type}, predicate:{text,type:"relation"},
              object:{text,type}, pmid, section, sentence}
         TRIPLES/triples_GENETIC_normalized.json  the triples with each GENETIC
              subject/object annotated with an hgnc_symbol (case-sensitive match of
              its text to a roman*.json key, else a greek*.json key, else a
              greek*.json greek_expanded value) plus the matched entry's control
              ("yes"/"no") flag
         TRIPLES/triples_GENETIC_DISEASE_normalized.json  the above, plus each
              DISEASE subject/object annotated with a mondo_label (case-sensitive
              match of its text to a disease*.json key)
         TRIPLES/triples_GENETIC_DISEASE_CHEMICAL_normalized.json  the above, plus
              each CHEMICAL subject/object annotated with a chebi_label (case-
              sensitive match of its text to a chemical*.json key) -- fully
              ID-normalized triples
         TRIPLES/triples.html   summary (counts, type-pair matrix, top predicates,
              by section, GENETIC->HGNC normalization strategy + coverage, samples)

Run from anywhere::  python triples.py
"""

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
OUT_DIR = ROOT / "TRIPLES"
GEN_DIR = ROOT / "GENETIC"           # roman*.json + greek*.json libraries
DIS_DIR = ROOT / "DISEASE"           # disease*.json libraries
CHE_DIR = ROOT / "CHEMICAL"          # chemical*.json libraries
JSON_OUT = OUT_DIR / "triples.json"
NORM_OUT = OUT_DIR / "triples_GENETIC_normalized.json"
NORM2_OUT = OUT_DIR / "triples_GENETIC_DISEASE_normalized.json"
NORM3_OUT = OUT_DIR / "triples_GENETIC_DISEASE_CHEMICAL_normalized.json"
HTML_OUT = OUT_DIR / "triples.html"

_WS = re.compile(r"\s+")
MAX_PRED = 120                      # cap on connecting-text length (chars)

# hyphen normalization for GENETIC entities (roman.py STAGE 1): unify dash-glyph
# variants to ASCII '-' and strip whitespace around '-' (e.g. 'Pisd - ps1' -> 'Pisd-ps1')
DASH_VARIANTS = "-‐‑‒–—―−⁃"
_TO_ASCII = str.maketrans({c: "-" for c in DASH_VARIANTS})
_DESPACE = re.compile(r"\s*-\s*")


def dash_normalize(text):
    return _DESPACE.sub("-", (text or "").translate(_TO_ASCII))


def entity_text(ent):
    """GENETIC entity surfaces are hyphen-normalized; other types kept verbatim."""
    return dash_normalize(ent["text"]) if ent.get("label") == "GENETIC" else ent["text"]


def build_genetic_maps():
    """Three case-sensitive lookup maps for GENETIC normalization to HGNC; each maps
    to a {hgnc_symbol, control} record (the 'control' yes/no flag is carried over
    from the roman/greek libraries, where controls.py added it):
      roman_map      : roman*.json KEY            -> record
      greek_key_map  : greek*.json KEY (original, glyph-bearing surface) -> record
      greek_exp_map  : greek*.json greek_expanded -> record
    (the roman/greek libraries live under GENETIC/.)"""
    roman_map, greek_key_map, greek_exp_map = {}, {}, {}
    for f in ("roman.json", "roman_ambiguous.json", "roman_cosine.json"):
        p = GEN_DIR / f
        if p.exists():
            for key, e in json.loads(p.read_text(encoding="utf-8")).items():
                roman_map[key] = {"hgnc_symbol": e.get("hgnc_symbol"),
                                  "control": e.get("control")}
    for f in ("greek.json", "greek_ambiguous.json", "greek_complex.json",
              "greek_cosine.json"):
        p = GEN_DIR / f
        if p.exists():
            for key, e in json.loads(p.read_text(encoding="utf-8")).items():
                if e.get("hgnc_symbol") is None:
                    continue
                rec = {"hgnc_symbol": e.get("hgnc_symbol"), "control": e.get("control")}
                greek_key_map.setdefault(key, rec)
                ge = e.get("greek_expanded")
                if ge is not None:
                    greek_exp_map.setdefault(ge, rec)
    return roman_map, greek_key_map, greek_exp_map


def normalize_genetic(triple, roman_map, greek_key_map, greek_exp_map):
    """Return a copy of `triple` whose GENETIC subject/object carry hgnc_symbol +
    control:
    (1) text == a roman*.json key (case-sensitive) -> that entry; else
    (2) text == a greek*.json key (original glyph surface) -> that entry; else
    (3) text == a greek*.json greek_expanded value -> that entry; else null.
    hgnc_via records the route; control is the 'yes'/'no' flag from the matched
    roman/greek entry (null if unmatched)."""
    out = {**triple, "subject": dict(triple["subject"]), "object": dict(triple["object"])}
    for role in ("subject", "object"):
        el = out[role]
        if el["type"] != "GENETIC":
            continue
        txt = el["text"]
        if txt in roman_map:
            rec, via = roman_map[txt], "roman key"
        elif txt in greek_key_map:
            rec, via = greek_key_map[txt], "greek key"
        elif txt in greek_exp_map:
            rec, via = greek_exp_map[txt], "greek_expanded"
        else:
            rec, via = None, None
        el["hgnc_symbol"] = rec["hgnc_symbol"] if rec else None
        el["hgnc_via"] = via
        el["control"] = rec["control"] if rec else None
    return out


def build_disease_map():
    """Case-sensitive map disease*.json KEY -> mondo_label (DISEASE/ libraries)."""
    dmap = {}
    for f in ("disease.json", "disease_ambiguous.json"):
        p = DIS_DIR / f
        if p.exists():
            for key, e in json.loads(p.read_text(encoding="utf-8")).items():
                dmap[key] = e.get("mondo_label")
    return dmap


def normalize_disease(triple, disease_map):
    """Return a copy of `triple` whose DISEASE subject/object carry a mondo_label:
    text == a disease*.json key (case-sensitive) -> that key's mondo_label, else
    null. mondo_via records the match. GENETIC/CHEMICAL elements pass through."""
    out = {**triple, "subject": dict(triple["subject"]), "object": dict(triple["object"])}
    for role in ("subject", "object"):
        el = out[role]
        if el["type"] != "DISEASE":
            continue
        txt = el["text"]
        if txt in disease_map:
            el["mondo_label"], el["mondo_via"] = disease_map[txt], "disease key"
        else:
            el["mondo_label"], el["mondo_via"] = None, None
    return out


def build_chemical_map():
    """Case-sensitive map chemical*.json KEY -> chebi_label (CHEMICAL/ libraries)."""
    cmap = {}
    for f in ("chemical.json", "chemical_ambiguous.json"):
        p = CHE_DIR / f
        if p.exists():
            for key, e in json.loads(p.read_text(encoding="utf-8")).items():
                cmap[key] = e.get("chebi_label")
    return cmap


def normalize_chemical(triple, chemical_map):
    """Return a copy of `triple` whose CHEMICAL subject/object carry a chebi_label:
    text == a chemical*.json key (case-sensitive) -> that key's chebi_label, else
    null. chebi_via records the match. GENETIC/DISEASE elements pass through."""
    out = {**triple, "subject": dict(triple["subject"]), "object": dict(triple["object"])}
    for role in ("subject", "object"):
        el = out[role]
        if el["type"] != "CHEMICAL":
            continue
        txt = el["text"]
        if txt in chemical_map:
            el["chebi_label"], el["chebi_via"] = chemical_map[txt], "chemical key"
        else:
            el["chebi_label"], el["chebi_via"] = None, None
    return out


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


def nested_key_stats(triples):
    """Summarize the nested keys at each level of the (fully-normalized) triples:
    the triple object, the subject/object element objects (combined), and the
    predicate object. Returns scope -> (key_order, present, types, example, total)."""
    def scan(objs):
        present, types, example, order = Counter(), defaultdict(set), {}, []
        for o in objs:
            for k, v in o.items():
                if k not in present:
                    order.append(k)
                present[k] += 1
                types[k].add(_typename(v))
                example.setdefault(k, v)
        return order, present, types, example, len(objs)
    elements = [t["subject"] for t in triples] + [t["object"] for t in triples]
    return {
        "triple object": scan(triples),
        "subject / object element": scan(elements),
        "predicate object": scan([t["predicate"] for t in triples]),
    }

# relation cue: the connecting text must carry a verb / relational stem to count as
# a predicate (otherwise it is coordination / a list / a preposition, not a relation)
RELCUE = re.compile(
    r"induc|inhibit|activat|regulat|express|\bbind|bound|target|associat|interact|"
    r"phosphorylat|suppress|promot|mediat|overexpress|knock\s?down|knock\s?out|"
    r"silenc|deplet|mutat|mutant|treat|encod|\bcaus|block|modulat|stimulat|increas|"
    r"decreas|reduc|elevat|attenuat|amelior|enhanc|repress|sensiti|resist|cleav|"
    r"degrad|secret|recruit|antagon|agoni|deficien|depend|correlat|involv|signal|"
    r"\brole\b|\beffect|abolish|abrogat|disrupt|restor|rescu|prevent|trigger|driv|"
    r"confer|\bloss\b|deletion|amplif|fusion|translocat|methylat|acetylat|ubiquitin|"
    r"glycosylat|inhibitor|activator|agonist|antagonist", re.I)


def extract(sentences, pmid):
    """Yield triples for one document's sentences."""
    for s in sentences:
        text = s.get("text", "")
        section = s.get("section_type", "")
        ents = [e for e in s.get("entities", []) if e.get("start") is not None]
        if len(ents) < 2:
            continue
        ents.sort(key=lambda e: e["start"])
        for a, b in zip(ents, ents[1:]):
            raw = text[a["end"]:b["start"]]
            pred = _WS.sub(" ", raw).strip()
            if not pred or len(pred) > MAX_PRED or not re.search("[A-Za-z]", pred):
                continue
            if not RELCUE.search(pred):
                continue
            yield {
                "subject": {"text": entity_text(a), "type": a["label"]},
                "predicate": {"text": pred, "type": "relation"},
                "object": {"text": entity_text(b), "type": b["label"]},
                "pmid": pmid,
                "section": section,
                "sentence": text,
            }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    triples = []
    n_files = n_sent = 0
    for fp in sorted(glob.glob(str(SENT_DIR / "*.json"))):
        try:
            rec = json.load(open(fp, encoding="utf-8"))
        except Exception:
            continue
        n_files += 1
        pmid = Path(rec.get("source_file", Path(fp).name)).stem
        sents = rec.get("sentences", [])
        n_sent += len(sents)
        triples.extend(extract(sents, pmid))

    JSON_OUT.write_text(json.dumps(triples, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")

    # ---- GENETIC -> HGNC normalization (separate output) ----
    roman_map, greek_key_map, greek_exp_map = build_genetic_maps()
    normalized = [normalize_genetic(t, roman_map, greek_key_map, greek_exp_map)
                  for t in triples]
    NORM_OUT.write_text(json.dumps(normalized, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
    via = Counter()
    symcount = Counter()
    ctrl = Counter()
    ctrl_sym = Counter()
    for t in normalized:
        for role in ("subject", "object"):
            el = t[role]
            if el["type"] != "GENETIC":
                continue
            via[el["hgnc_via"]] += 1
            ctrl[el["control"]] += 1
            if el["hgnc_symbol"] is not None:
                s = el["hgnc_symbol"]
                s = s if isinstance(s, str) else "/".join(s)
                symcount[s] += 1
                if el["control"] == "yes":
                    ctrl_sym[s] += 1
    gnorm = {"elems": sum(via.values()), "roman": via["roman key"],
             "greek_key": via["greek key"], "greek": via["greek_expanded"],
             "none": via[None], "top": symcount,
             "ctrl_yes": ctrl["yes"], "ctrl_no": ctrl["no"], "ctrl_null": ctrl[None],
             "ctrl_top": ctrl_sym}

    # ---- DISEASE -> MONDO normalization (chained onto the GENETIC-normalized set) ----
    disease_map = build_disease_map()
    normalized2 = [normalize_disease(t, disease_map) for t in normalized]
    NORM2_OUT.write_text(json.dumps(normalized2, ensure_ascii=False, indent=2) + "\n",
                         encoding="utf-8")
    dvia = Counter()
    mlcount = Counter()
    for t in normalized2:
        for role in ("subject", "object"):
            el = t[role]
            if el["type"] != "DISEASE":
                continue
            dvia[el["mondo_via"]] += 1
            if el["mondo_label"] is not None:
                ml = el["mondo_label"]
                mlcount[ml if isinstance(ml, str) else "/".join(ml)] += 1
    dnorm = {"elems": sum(dvia.values()), "matched": dvia["disease key"],
             "none": dvia[None], "top": mlcount}

    # ---- CHEMICAL -> ChEBI normalization (chained onto GENETIC+DISEASE-normalized) ----
    chemical_map = build_chemical_map()
    normalized3 = [normalize_chemical(t, chemical_map) for t in normalized2]
    NORM3_OUT.write_text(json.dumps(normalized3, ensure_ascii=False, indent=2) + "\n",
                         encoding="utf-8")
    cvia = Counter()
    clcount = Counter()
    for t in normalized3:
        for role in ("subject", "object"):
            el = t[role]
            if el["type"] != "CHEMICAL":
                continue
            cvia[el["chebi_via"]] += 1
            if el["chebi_label"] is not None:
                cl = el["chebi_label"]
                clcount[cl if isinstance(cl, str) else "/".join(cl)] += 1
    cnorm = {"elems": sum(cvia.values()), "matched": cvia["chemical key"],
             "none": cvia[None], "top": clcount}

    typepair = Counter((t["subject"]["type"], t["object"]["type"]) for t in triples)
    preds = Counter(t["predicate"]["text"].lower() for t in triples)
    by_section = Counter(t["section"] for t in triples)
    papers = {t["pmid"] for t in triples}
    print(f"files={n_files:,}  sentences={n_sent:,}  triples={len(triples):,} "
          f"(from {len(papers):,} papers) -> {JSON_OUT}")
    print("  top type pairs: " + ", ".join(
        f"{a}->{b} {c}" for (a, b), c in typepair.most_common(6)))
    print("  top predicates: " + ", ".join(
        f"{p!r} {c}" for p, c in preds.most_common(8)))
    print(f"  GENETIC normalization: {gnorm['elems']:,} GENETIC elements -> "
          f"hgnc_symbol via roman key {gnorm['roman']:,}, greek key "
          f"{gnorm['greek_key']:,}, greek_expanded {gnorm['greek']:,}, "
          f"unmatched {gnorm['none']:,} -> {NORM_OUT.name}")
    print(f"    control flag: yes {gnorm['ctrl_yes']:,}, no {gnorm['ctrl_no']:,}, "
          f"unannotated {gnorm['ctrl_null']:,}")
    print(f"  DISEASE normalization: {dnorm['elems']:,} DISEASE elements -> "
          f"mondo_label via disease key {dnorm['matched']:,}, "
          f"unmatched {dnorm['none']:,} -> {NORM2_OUT.name}")
    print(f"  CHEMICAL normalization: {cnorm['elems']:,} CHEMICAL elements -> "
          f"chebi_label via chemical key {cnorm['matched']:,}, "
          f"unmatched {cnorm['none']:,} -> {NORM3_OUT.name}")

    nkeys = nested_key_stats(normalized3)
    render_html(triples, n_files, n_sent, len(papers), typepair, preds, by_section,
                gnorm, dnorm, cnorm, nkeys)
    print(f"Wrote {HTML_OUT}")


def render_html(triples, n_files, n_sent, n_papers, typepair, preds, by_section,
                gnorm, dnorm, cnorm, nkeys):
    esc = html.escape
    types = ["GENETIC", "DISEASE", "CHEMICAL"]
    nk_html = ""
    for scope, (order, present, ktypes, example, total) in nkeys.items():
        rows = "".join(
            f'<tr><td><code>{esc(k)}</code></td><td class="num">{present[k]:,}/{total:,}</td>'
            f'<td><code>{esc(" | ".join(sorted(ktypes[k])))}</code></td>'
            f'<td><code>{esc(_short(example[k]))}</code></td></tr>' for k in order)
        nk_html += (f'<h3>{esc(scope)} &mdash; {len(order)} keys (n={total:,})</h3>'
                    f'<details><summary>show / hide table</summary>'
                    f'<table><tr><th>nested key</th><th class="num">present</th>'
                    f'<th>type(s)</th><th>example</th></tr>{rows}</table></details>')
    gtop = "".join(
        f'<tr><td>{esc(s)}</td><td class="num">{c:,}</td></tr>'
        for s, c in gnorm["top"].most_common(25))
    gctrl = "".join(
        f'<tr><td>{esc(s)}</td><td class="num">{c:,}</td></tr>'
        for s, c in gnorm["ctrl_top"].most_common(25))
    dtop = "".join(
        f'<tr><td>{esc(s)}</td><td class="num">{c:,}</td></tr>'
        for s, c in dnorm["top"].most_common(25))
    ctop = "".join(
        f'<tr><td>{esc(s)}</td><td class="num">{c:,}</td></tr>'
        for s, c in cnorm["top"].most_common(25))

    head = "".join(f'<th>{t}</th>' for t in types)
    matrix = "".join(
        f'<tr><th>{a}</th>' + "".join(
            f'<td class="num">{typepair.get((a, b), 0):,}</td>' for b in types)
        + '</tr>' for a in types)
    pred_rows = "".join(
        f'<tr><td><code>{esc(p)}</code></td><td class="num">{c:,}</td></tr>'
        for p, c in preds.most_common(40))
    sec_rows = "".join(
        f'<tr><td>{esc(s or "&mdash;")}</td><td class="num">{c:,}</td></tr>'
        for s, c in by_section.most_common())

    style = (
        " body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1040px;"
        "margin:2rem auto;padding:0 1rem;color:#1a1a1a;}"
        " h1{font-size:1.45rem;} h2{font-size:1.15rem;margin-top:1.8rem;"
        "border-bottom:1px solid #ddd;padding-bottom:.3rem;} h3{font-size:1rem;margin-top:1.2rem;}"
        " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;}"
        " th,td{border:1px solid #ccc;padding:.3rem .6rem;text-align:left;vertical-align:top;}"
        " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
        " .big{font-size:2rem;font-weight:700;}"
        " .headline{background:#eef4fb;border:1px solid #cdddf0;border-radius:8px;"
        "padding:.8rem 1rem;margin:1rem 0;}"
        " code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.9em;}"
        " p.note{color:#444;font-size:.92em;}")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>triples &mdash; BioBERT entity relation triples</title>
<style>{style}</style></head><body>
<h1>Relation triples from BioBERT NER sentences</h1>
<p>(subject, predicate, object) triples derived from co-occurring BioBERT entities:
consecutive typed entities in a sentence joined by connecting text that carries a
relation cue. Full data in <code>TRIPLES/triples.json</code> (each triple keeps the
typed elements, PMID, section and sentence). Produced by <code>triples.py</code>.</p>
<div class="headline"><span class="big">{len(triples):,}</span> triples
&mdash; from {n_sent:,} sentences across {n_papers:,} papers ({n_files:,} files).</div>

<h2>Triples by entity-type pair (subject &rarr; object)</h2>
<details open><summary>show / hide table</summary>
<table><tr><th>subj \\ obj</th>{head}</tr>{matrix}</table></details>

<h2>Top predicates</h2>
<details><summary>show / hide table</summary>
<table><tr><th>predicate (connecting text)</th><th class="num">triples</th></tr>{pred_rows}</table></details>

<h2>Triples by section</h2>
<details><summary>show / hide table</summary>
<table><tr><th>section</th><th class="num">triples</th></tr>{sec_rows}</table></details>

<h2>GENETIC normalization to HGNC</h2>
<p>Each <strong>GENETIC</strong> subject / object is normalized to an HGNC symbol and
written to <code>TRIPLES/triples_GENETIC_normalized.json</code> (as
<code>hgnc_symbol</code> + <code>hgnc_via</code> on the element). The match is
<strong>case-sensitive string equality</strong>, in two routes:</p>
<ol>
<li><strong>roman key.</strong> If the (hyphen-normalized) GENETIC text equals a
<em>key</em> in <code>GENETIC/roman.json</code> / <code>roman_ambiguous.json</code> /
<code>roman_cosine.json</code>, take that entry's <code>hgnc_symbol</code>.</li>
<li><strong>greek key.</strong> Otherwise, if the text equals a <em>key</em> (the
original, glyph-bearing surface, e.g. <code>IFN-&gamma;</code>) in
<code>GENETIC/greek.json</code> / <code>greek_ambiguous.json</code> /
<code>greek_complex.json</code> / <code>greek_cosine.json</code>, take that entry's
<code>hgnc_symbol</code>.</li>
<li><strong>greek_expanded.</strong> Otherwise, if the text equals a
<code>greek_expanded</code> value (the glyph-spelled form, e.g.
<code>IFN-gamma</code>) in the same greek libraries, take that entry's
<code>hgnc_symbol</code>.</li>
</ol>
<p>The symbol is a single value (single match), a list (ambiguous / complex), or
<code>null</code> (no match). The matched roman/greek entry's
<strong><code>control</code></strong> flag (<code>"yes"</code>/<code>"no"</code>,
added by <code>controls.py</code>: whether the gene is really an experimental
control / reagent rather than a studied gene) is copied onto the GENETIC element
too. DISEASE/CHEMICAL elements are unchanged.</p>
<div class="headline">{gnorm["elems"]:,} GENETIC subject/object elements &mdash;
normalized to an HGNC symbol via <strong>roman key {gnorm["roman"]:,}</strong>,
<strong>greek key {gnorm["greek_key"]:,}</strong>,
<strong>greek_expanded {gnorm["greek"]:,}</strong>; unmatched {gnorm["none"]:,}
({100*(gnorm["roman"]+gnorm["greek_key"]+gnorm["greek"])/max(1,gnorm["elems"]):.1f}% normalized).</div>
<h3>Most frequent normalized HGNC symbols (GENETIC elements)</h3>
<details><summary>show / hide table</summary>
<table><tr><th>hgnc_symbol</th><th class="num">elements</th></tr>{gtop}</table></details>
<h3>control annotation (GENETIC elements)</h3>
<p>The <code>control</code> flag carried onto each GENETIC element from the
roman/greek libraries &mdash; <code>"yes"</code> = the gene surface is an
experimental control / reagent (per <code>controls.py</code>), <code>"no"</code> =
a genuine gene, <code>null</code> = no library match.</p>
<div class="headline">{gnorm["ctrl_yes"]:,} GENETIC elements flagged
<strong>control = yes</strong>; {gnorm["ctrl_no"]:,} control = no;
{gnorm["ctrl_null"]:,} unannotated (no match).</div>
<details><summary>show / hide control = yes symbols</summary>
<table><tr><th>hgnc_symbol (control = yes)</th><th class="num">elements</th></tr>{gctrl}</table></details>

<h2>DISEASE normalization to MONDO</h2>
<p>Each <strong>DISEASE</strong> subject / object is normalized to a MONDO label and
written to <code>TRIPLES/triples_GENETIC_DISEASE_normalized.json</code> (as
<code>mondo_label</code> + <code>mondo_via</code> on the element), by
<strong>case-sensitive string equality</strong> of its text to a <em>key</em> in
<code>DISEASE/disease.json</code> / <code>disease_ambiguous.json</code>; the label is
a single value (single match), a list (ambiguous), or <code>null</code> (no match).
GENETIC (already annotated with <code>hgnc_symbol</code>) and CHEMICAL elements pass
through unchanged.</p>
<div class="headline">{dnorm["elems"]:,} DISEASE subject/object elements &mdash;
normalized to a MONDO label via <strong>disease key {dnorm["matched"]:,}</strong>;
unmatched {dnorm["none"]:,}
({100*dnorm["matched"]/max(1,dnorm["elems"]):.1f}% normalized).</div>
<h3>Most frequent normalized MONDO labels (DISEASE elements)</h3>
<details><summary>show / hide table</summary>
<table><tr><th>mondo_label</th><th class="num">elements</th></tr>{dtop}</table></details>

<h2>CHEMICAL normalization to ChEBI</h2>
<p>Each <strong>CHEMICAL</strong> subject / object is normalized to a ChEBI label and
written to <code>TRIPLES/triples_GENETIC_DISEASE_CHEMICAL_normalized.json</code> (as
<code>chebi_label</code> + <code>chebi_via</code> on the element), by
<strong>case-sensitive string equality</strong> of its text to a <em>key</em> in
<code>CHEMICAL/chemical.json</code> / <code>chemical_ambiguous.json</code>; the label
is a single value (single match), a list (ambiguous), or <code>null</code> (no
match). GENETIC (<code>hgnc_symbol</code>) and DISEASE (<code>mondo_label</code>)
annotations carry through. This file is the fully ID-normalized triple set.</p>
<div class="headline">{cnorm["elems"]:,} CHEMICAL subject/object elements &mdash;
normalized to a ChEBI label via <strong>chemical key {cnorm["matched"]:,}</strong>;
unmatched {cnorm["none"]:,}
({100*cnorm["matched"]/max(1,cnorm["elems"]):.1f}% normalized).</div>
<h3>Most frequent normalized ChEBI labels (CHEMICAL elements)</h3>
<details><summary>show / hide table</summary>
<table><tr><th>chebi_label</th><th class="num">elements</th></tr>{ctop}</table></details>

<h2>Nested keys of <code>triples_GENETIC_DISEASE_CHEMICAL_normalized.json</code></h2>
<p>The keys at each level of the fully ID-normalized triples &mdash; the triple
object, the subject/object element objects (combined), and the predicate object
&mdash; with coverage (present / total), value type(s) and an example. Type-specific
element keys (<code>hgnc_symbol</code>/<code>mondo_label</code>/<code>chebi_label</code>
and their <code>*_via</code>) appear only on elements of the matching type, so their
coverage is partial.</p>
{nk_html}

<h2>Method &amp; caveats</h2>
<p class="note">Entities are BioBERT NER spans (<code>biobert_genetic/disease/chemical</code>);
this is not a trained relation-extraction model. A triple is a pair of
<em>consecutive</em> entities (nothing else between them) whose connecting text
(&le;{MAX_PRED} chars) contains a relation cue (verb / relational stem: induced,
inhibited, expression, associated, -mediated, &hellip;); pure coordination
(&ldquo;and&rdquo;, &ldquo;or&rdquo;, &ldquo;in&rdquo;) is excluded. The predicate
is raw surface text; treat these as candidate relations. <strong>GENETIC</strong>
entity surfaces are hyphen-normalized (dash-glyph variants &rarr; <code>-</code>,
whitespace around <code>-</code> stripped, e.g. <code>Pisd - ps1</code>&rarr;
<code>Pisd-ps1</code>); DISEASE/CHEMICAL text is kept verbatim.</p>
</body></html>
"""
    HTML_OUT.write_text(doc, encoding="utf-8")


if __name__ == "__main__":
    main()
