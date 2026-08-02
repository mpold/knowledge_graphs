#!/usr/bin/env python3
"""
high_confidence_g.py -- extract the high-confidence Gene-context relation triples from
the scored RE output, summarize the threshold statistics, and draw the brain-cancer
gene-gene relationship graph.

This is the "G" (gene-only) variant of high_confidence.py -- and, since high_confidence.py
was DEPRECATED, the maintained step-3 entry point. Where high_confidence.py applies the
"G_D_C" filter (gene-gene IN a disease/chemical context), this step applies the "G" filter
(gene-gene, context-agnostic). Its "_G" output names never clobber the other's, so both can
still target the same data root.

TYPED, SIGNED EDGES + MULTI-MODEL MERGE (the two things high_confidence.py does not do)
  * The graph categorizes each edge by the RELATION the model predicted -- activates /
    inhibits / binds / interacts / associated, prefixed "not " when the statement is
    negated -- instead of by polarity alone. With the BioRED checkpoint in the routing
    that label is SIGNED, so `A inhibits B` and `A activates B` are no longer the same
    edge. Colour encodes the relation, negated edges are dashed, and the per-sentence
    tooltip carries each sentence's own label (an edge shows its dominant relation, so a
    minority reading stays visible there).
  * relation_extraction.py --route-mode additive scores a pair with EVERY applicable
    checkpoint, so the input can hold two triples per pair -- the binary PPI verdict and
    the typed BioRED one, sharing a pair_id. --merge collapses them to one triple per
    pair; the default `union` keeps every pair EITHER model kept -- so BioRED-only edges
    (relations BioInfer missed, plus everything BioRED types) are in the graph -- and
    prefers the TYPED label wherever BioRED fired, so a pair both models claim appears
    once, signed, scored by the higher of the two. `--merge gate` is the stricter variant
    (keep only what the binary model also kept) for when BioRED-only edges prove noisy;
    see also typed / none, and `score_by_model` / `models` / `corroborated` on each
    merged triple.

  The "G" step -- the ONE logic difference from high_confidence.py:
  the DISEASE-or-CHEMICAL sentence-context requirement (G_D_C rule 5) is DROPPED.
  A triple QUALIFIES (the "G" filter) when ALL of these hold:
    1. score >= SCORE                                   (export default 0.8)
    2. annotated control:no   -- >=1 GENETIC endpoint has control == "no"
    3. NOT annotated control:yes -- no endpoint has control == "yes"
    4. NOT hgnc_symbol == "MKI67" on either endpoint
  (The G_D_C filter's 5th condition -- the sentence must also carry a DISEASE or
  CHEMICAL entity -- is intentionally NOT applied in this step, so the "G" universe is
  strictly larger: every high-confidence gene-gene relation, disease/chemical context or not.)

All inputs are read from the pipeline output tree (the writable run dir gpu.py produced),
which defaults to ``kaggle_working/`` next to this script; override with ``--data-root``.

Inputs  : <data-root>/TRIPLES/triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json  (scored + normalized)
          <data-root>/{CHEMICAL,databases,sentences}/...                             (drug targets + year + corpus size)
          <data-root>/DISEASE/disease.json, <data-root>/CHEMICAL/chemical.json  (--nodes all: phenotype / non_chemical flags)
Outputs : <data-root>/TRIPLES/high_confidence_G.json     qualifying triples at --score (default 0.8)
          (--nodes all writes the same three outputs with an "_M" suffix instead)
          <data-root>/summaries/high_confidence_G.html    gene-gene graph with a 0.5..0.99 in-browser score slider
                                                          and a "match text in sentence" filter (substring, or /regex/)
          <root>/<current_dir>_YYYY_MM_DD_G.html          a copy of that graph, named after the directory
                                                           holding this script plus today's date;
                                                           e.g. lung_large_2026_07_19_G.html

The output filenames all differ from those written by high_confidence.py (which uses the
"_G_D_C" JSON, "high_confidence.html" graph, and unsuffixed "<dir>_<date>.html" copy) so the
two scripts can be run against the same data root without clobbering each other's outputs.

  * `--nodes all` widens the graph beyond gene-gene. The gene-only view draws an edge only
    when BOTH endpoints carry a single HGNC symbol, which discards everything BioRED adds
    beyond gene-gene -- on the reference run, 649 of its 655 solo pairs. With --nodes all,
    DISEASE and CHEMICAL endpoints become nodes too, identified by their normalized id
    (mondo_label / chebi_label), giving gene-disease / chemical-gene / chemical-disease /
    chemical-chemical edges. Node hygiene reuses the pipeline's own annotations: control
    (genes), phenotype (DISEASE/disease.json) and non_chemical (CHEMICAL/chemical.json),
    plus DISEASE_IGNORE for labels too generic to be a useful node. Outputs take an "_M"
    suffix so they never clobber the "_G" ones, node shape encodes the type (gene dot /
    disease diamond / chemical square), a "Node type" filter appears in the panel, and the
    score slider opens at >=0.8 (most gene-disease edges sit in 0.5-0.8, so the gene-only
    default of >=0.99 would show an almost empty canvas).

Run::  python high_confidence_g.py [--data-root kaggle_working] [--score 0.8] [--thresholds 0.8,0.95,0.99] [--no-graph]
       python high_confidence_g.py --nodes all          # + disease and chemical nodes ("_M" outputs)
       python high_confidence_g.py --merge gate         # drop pairs only BioRED claimed
       python high_confidence_g.py --merge none         # pre-merge behaviour (duplicates survive)
"""
import argparse
import collections
import datetime
import html
import json
import re
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# The pipeline outputs (gpu.py's writable run dir) live under kaggle_working/ by default;
# every input below is resolved beneath this data root. Override with --data-root.
DATA_ROOT = ROOT / "kaggle_working"

# module-level paths; (re)bound to DATA_ROOT by set_data_root() so --data-root can retarget them
OUT_DIR = XML_DIR = SENT_DIR = None
RE_FILE = PMC_YEARS = TARGET_FILE = None
DISEASE_LIB = CHEM_LIB = None
JSON_OUT = GRAPH_OUT = None


def set_data_root(data_root):
    """Point every input/output path at `data_root` (the pipeline's output tree)."""
    global DATA_ROOT, OUT_DIR, XML_DIR, SENT_DIR, RE_FILE, PMC_YEARS
    global TARGET_FILE, DISEASE_LIB, CHEM_LIB, JSON_OUT, GRAPH_OUT
    DATA_ROOT = Path(data_root).resolve()
    OUT_DIR = DATA_ROOT / "TRIPLES"
    XML_DIR = DATA_ROOT / "experimental_ner"   # input XML corpus (may be empty in the bundle)
    SENT_DIR = DATA_ROOT / "sentences"         # one JSON per source document (corpus-size fallback)
    RE_FILE = OUT_DIR / "triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json"
    PMC_YEARS = DATA_ROOT / "databases" / "pmc_years.json"
    TARGET_FILE = DATA_ROOT / "CHEMICAL" / "chemical_to_target.json"   # gene -> corpus chemicals (in_corpus_GENETIC flag)
    # normalization libraries carrying the in-place phenotype / non_chemical flags (--nodes all)
    DISEASE_LIB = DATA_ROOT / "DISEASE" / "disease.json"
    CHEM_LIB = DATA_ROOT / "CHEMICAL" / "chemical.json"
    # "_G" (gene-only) output names, distinct from high_confidence.py's "_G_D_C"/"high_confidence.html".
    JSON_OUT = OUT_DIR / "high_confidence_G.json"
    GRAPH_OUT = DATA_ROOT / "summaries" / "high_confidence_G.html"


set_data_root(DATA_ROOT)

GRAPH_BASE = 0.5           # graph universe = qualifying triples at this score (lowest in-browser slider stop)
VIS_URL = "https://unpkg.com/vis-network@9.1.9/standalone/umd/vis-network.min.js"

# predicate.text -> graph edge category. The typed BioRED labels carry a SIGN, which is
# the point of training that model at all; an unrecognized label falls through as itself.
REL_CAT = {"upregulator/activator": "activates", "downregulator/inhibitor": "inhibits",
           "binds": "binds", "interacts": "interacts", "associated": "associated",
           # the rare chemical-chemical BioRED types add no sign -> unsigned bucket
           "cotreatment": "associated", "comparison": "associated",
           "drug-interaction": "associated", "conversion": "associated"}
RBASE = {"activates": "#2e9e5b", "inhibits": "#e0533d", "binds": "#3b7dd8",
         "interacts": "#6b3fa0", "associated": "#8a8f98"}
# a negated statement ("X does not inhibit Y") is a fact worth keeping and worth seeing as
# distinct: same category name prefixed with "not ", drawn in one warning colour.
RCOLOR = {**RBASE, **{f"not {k}": "#d59a2e" for k in RBASE}}


# ----- multi-type nodes (--nodes all) ----------------------------------------
# The gene-only graph draws an edge only when BOTH endpoints carry a single HGNC symbol,
# which silently discards everything the BioRED checkpoint adds beyond gene-gene:
# gene-disease, chemical-gene, chemical-disease, chemical-chemical. --nodes all keeps
# them, with each endpoint identified by its NORMALIZED id (hgnc_symbol / mondo_label /
# chebi_label) so the same entity under different surfaces is one node.
NODE_STYLE = {                                   # kind -> vis shape + default colours
    "gene": {"shape": "dot", "bg": "#cfe3ff", "border": "#2b6cb0"},
    "disease": {"shape": "diamond", "bg": "#ffe0e0", "border": "#b3243b"},
    "chemical": {"shape": "square", "bg": "#ece0f8", "border": "#6b3fa0"},
}
# Disease terms too generic to be a useful node -- they attach to everything and turn the
# view into a hairball around one hub. Matched case-insensitively against BOTH the surface
# text and the normalized MONDO label, because the two cases differ:
#   neoplasm  a LABEL: "tumor"/"tumors"/"tumour" all normalize onto it. On the reference run
#             this is the entry that does the work -- 167 endpoint mentions, ~40% of all
#             disease mentions. Dropping them is a deliberate trade: "gene X associated with
#             tumor" carries no information in a cancer corpus. Remove this entry to keep them.
#   cancer    a LABEL, 18 mentions. Note the match is exact, so "breast cancer" and
#             "lung neoplasm" are NOT caught by these two.
#   os        a SURFACE: "OS" is overall survival, a metric the NER mislabels as DISEASE.
#             It never appears as a MONDO label, so a label-only check would never fire.
DISEASE_IGNORE = {"neoplasm", "cancer", "os"}

# Ambiguous disease surfaces whose normalization is wrong FOR THIS CORPUS: the mention is
# real, only the mapping is off, so overriding beats ignoring. Keys lowercase surfaces; the
# value replaces whatever mondo_label the normalizer produced.
#   lcc  disease.py maps it to "leukoencephalopathy with calcifications and cysts" (6
#        mentions); in a lung-large-cell corpus LCC is large cell carcinoma.
DISEASE_SURFACE = {"lcc": "lung large cell carcinoma"}

# Distinct MONDO terms that name the SAME disease at the same granularity. Splitting one
# disease across several nodes is an artefact of ontology structure, not a finding: MONDO
# carries "lung cancer" (MONDO:0008903), "lung carcinoma" (MONDO:0005138) and
# "lung neoplasm" as separate terms, and which one a mention lands on depends only on the
# surface the author happened to write ("lung cancer" vs "lung tumor"). Keys are matched
# lowercased against the resolved label; the value is the label the merged node takes.
# SUBTYPES ARE DELIBERATELY NOT MERGED -- non-small cell lung carcinoma, small cell lung
# carcinoma, lung adenocarcinoma, lung large cell carcinoma and pulmonary LCNEC are
# different diseases and stay different nodes.
DISEASE_MERGE = {
    "lung carcinoma": "lung cancer",     # MONDO:0005138 -> MONDO:0008903's plainer label
    "lung neoplasm": "lung cancer",      # surfaces "lung tumor(s)" / "lung tumours"
    "breast carcinoma": "breast cancer",
}


def model_roles(t):
    """{'ppi','biored'} -- which training corpus stands behind a (merged) triple. A
    BioRED-style checkpoint is recognized by name (as in compare_re.py and the merge);
    every other checkpoint counts as the binary/base model. Lets the graph answer 'which
    edges does each training set actually contribute?'"""
    names = t.get("models") or [((t.get("predicate") or {}).get("model") or "")]
    return {("biored" if "biored" in (n or "").lower() else "ppi") for n in names if n}


def src_tag(roles):
    """'ppi' | 'biored' | 'both' for a set of roles."""
    return "both" if {"ppi", "biored"} <= set(roles) else ("biored" if "biored" in roles else "ppi")


def rel_cat(t):
    """Edge category for a triple: its relation label, prefixed 'not ' when negated."""
    lab = ((t.get("predicate") or {}).get("text") or "").strip().lower()
    cat = REL_CAT.get(lab, lab or "interacts")
    return f"not {cat}" if t.get("polarity") == "negated" else cat


# ----- multi-model merge (BioRED alongside PPI) -------------------------------
# relation_extraction.py --route-mode additive scores each pair with EVERY applicable
# checkpoint, so triples_re.json can carry two triples per entity pair: the binary PPI
# verdict ("interacts") and the typed, signed BioRED one ("downregulator/inhibitor").
# They share a pair_id. This step collapses them to one triple per pair:
#
#   union (default) keep every pair EITHER model kept -- including the ones only BioRED
#                  found (gene-gene relations BioInfer missed) -- and prefer the TYPED
#                  label wherever BioRED fired, so a pair both models claim enters the
#                  graph once, signed. Score is the max of the two.
#   gate           the stricter variant: keep a pair only if the BINARY model also kept
#                  it. BioInfer is sentence-scoped and cleanly labelled, so it is the
#                  more conservative judge of WHETHER an edge exists; use this if
#                  BioRED-only edges prove noisy (its annotation is document-level).
#   typed          the typed model alone
#   none           no merge -- duplicates survive into graph_payload(), which then
#                  collapses them per (pair, sentence) by max score
#
# With only one checkpoint in the file every policy is a no-op.
MERGE_POLICIES = ("gate", "union", "typed", "none")


def _pair_key(t):
    """Stable identity of one entity pair in one sentence. relation_extraction.py writes
    pair_id; older files fall back to the text-level key."""
    return t.get("pair_id") or (t.get("pmid"), t.get("sentence"),
                                (t.get("subject") or {}).get("text"),
                                (t.get("object") or {}).get("text"))


def classify_models(triples, typed_name=None, gate_name=None):
    """(all models, typed model, gate model) present in the file. The typed one is
    identified by name (matching compare_re.py); override with --typed-model/--gate-model."""
    models = sorted({(t.get("predicate") or {}).get("model") for t in triples
                     if (t.get("predicate") or {}).get("model")})
    typed = typed_name or next((m for m in models if "biored" in m.lower()), None)
    gate = gate_name or next((m for m in models if m != typed), None)
    return models, typed, gate


def merge_models(triples, policy="gate", typed_name=None, gate_name=None):
    """One triple per entity pair under `policy`. Returns (triples, summary_string)."""
    models, typed, gate = classify_models(triples, typed_name, gate_name)
    if policy == "none":
        return triples, f"merge: skipped (--merge none); {len(models)} model(s): {', '.join(models) or 'none'}"
    if len(models) < 2 or typed is None:
        why = "only one checkpoint in the file" if len(models) < 2 else \
              f"no typed (BioRED-style) checkpoint among {', '.join(models)}"
        return triples, f"merge: nothing to merge -- {why}"

    groups, order = collections.defaultdict(dict), {}
    for i, t in enumerate(triples):
        k, m = _pair_key(t), t["predicate"].get("model")
        cur = groups[k].get(m)
        if cur is None or float(t.get("score") or 0.0) > float(cur.get("score") or 0.0):
            groups[k][m] = t
        order.setdefault(k, i)

    out = []
    n_drop = n_typed_lab = n_corr = 0
    for k in sorted(groups, key=lambda key: order[key]):
        by = groups[k]
        g, ty = by.get(gate), by.get(typed)
        if policy == "gate" and g is None:          # typed model alone claimed this pair
            n_drop += 1
            continue
        if policy == "typed":
            if ty is None:
                n_drop += 1
                continue
            by = {typed: ty}
        base = ty or g or next(iter(by.values()))   # typed label wins where it exists
        t = dict(base)
        t["subject"], t["object"] = dict(base["subject"]), dict(base["object"])
        t["predicate"] = dict(base["predicate"])
        scores = {m: float(v.get("score") or 0.0) for m, v in by.items()}
        t["score"] = max(scores.values())
        t["score_by_model"] = {m: round(s, 4) for m, s in sorted(scores.items())}
        t["models"] = sorted(by)
        t["corroborated"] = len(by) > 1
        n_corr += bool(t["corroborated"])
        n_typed_lab += bool(ty is not None)          # label came from the typed model
        out.append(t)
    dropped = f", {n_drop:,} dropped by the gate" if policy == "gate" else \
              (f", {n_drop:,} without a typed verdict dropped" if policy == "typed" else "")
    return out, (f"merge[{policy}]: {len(triples):,} triples ({', '.join(models)}) -> {len(out):,} pairs; "
                 f"gate={gate} typed={typed}; {n_corr:,} corroborated by both, "
                 f"{n_typed_lab:,} took the typed label{dropped}")


# ----- qualifying filter -----------------------------------------------------
def syms(e):
    v = e.get("hgnc_symbol")
    return (set(v) if isinstance(v, list) else {v}) if v is not None else set()


def has_mki67(t):
    return "MKI67" in (syms(t["subject"]) | syms(t["object"]))


def ctrl_no(t):
    return t["subject"].get("control") == "no" or t["object"].get("control") == "no"


def ctrl_yes(t):
    return t["subject"].get("control") == "yes" or t["object"].get("control") == "yes"


def qualifies(t, T):
    # gene-only filter: no DISEASE/CHEMICAL sentence-context requirement (that was the "_G_D_C" filter)
    sc = t.get("score")
    return (isinstance(sc, (int, float)) and sc >= T
            and ctrl_no(t) and not ctrl_yes(t) and not has_mki67(t))


def load_type_flags():
    """{'phenotype': {surface: yes/no}, 'non_chemical': {surface: yes/no}} from the
    normalization libraries. phenotypes.py and nonchemical.py annotate those files in
    place; the flags never reach the RE triples, so they are joined here on the surface
    text. Missing library -> empty map (nothing filtered)."""
    out = {"phenotype": {}, "non_chemical": {}}
    for path, key in ((DISEASE_LIB, "phenotype"), (CHEM_LIB, "non_chemical")):
        try:
            lib = json.loads(Path(path).read_text(encoding="utf-8"))
        except Exception:
            continue
        for surface, v in lib.items():
            if isinstance(v, dict) and v.get(key):
                out[key][str(surface).strip().lower()] = v[key]
    return out


def node_of(e, flags):
    """(kind, normalized_id) for one endpoint, or None when it cannot be a node.

    Dropped: endpoints with no normalized id (an un-normalized surface is not an entity we
    can merge across documents), lab controls / MKI67 (genes), process-or-phenotype
    surfaces (diseases), gene-or-process surfaces mislabelled CHEMICAL, and over-generic
    disease labels."""
    t, txt = e.get("type"), (e.get("text") or "").strip().lower()
    if t == "GENETIC":
        s = single(e.get("hgnc_symbol"))
        if not s or s == "MKI67" or e.get("control") == "yes":
            return None
        return ("gene", s)
    if t == "DISEASE":
        s = DISEASE_SURFACE.get(txt) or single(e.get("mondo_label"))   # override a wrong mapping
        if not s or txt in DISEASE_IGNORE or s.strip().lower() in DISEASE_IGNORE:
            return None
        if flags["phenotype"].get(txt) == "yes":
            return None
        s = s.strip()
        return ("disease", DISEASE_MERGE.get(s.lower(), s))   # one disease, one node
    if t == "CHEMICAL":
        s = single(e.get("chebi_label"))
        if not s or flags["non_chemical"].get(txt) == "yes":
            return None
        return ("chemical", s.strip())
    return None


def node_drop_report(triples, flags, top=8):
    """Lines accounting for every endpoint that did NOT become a node, so the node filters
    are visible rather than silent. The 'no id' buckets are the actionable ones: they are
    normalization gaps upstream (disease.py / chemical.py), not deliberate exclusions."""
    kept, drops = collections.Counter(), collections.Counter()
    unnorm, overridden = collections.Counter(), 0
    merged = collections.Counter()
    for t in triples:
        for side in ("subject", "object"):
            e = t[side]
            typ, txt = e.get("type"), (e.get("text") or "").strip()
            key = txt.lower()
            if typ == "DISEASE":
                if key in DISEASE_SURFACE:
                    overridden += 1
                lbl0 = DISEASE_SURFACE.get(key) or single(e.get("mondo_label"))
                if lbl0 and lbl0.strip().lower() in DISEASE_MERGE:
                    merged[f"{lbl0.strip()} -> {DISEASE_MERGE[lbl0.strip().lower()]}"] += 1
            n = node_of(e, flags)
            if n:
                kept[n[0]] += 1
                continue
            if typ == "GENETIC":
                drops["gene: no HGNC symbol" if not single(e.get("hgnc_symbol"))
                      else "gene: control / MKI67"] += 1
            elif typ == "DISEASE":
                lbl = DISEASE_SURFACE.get(key) or single(e.get("mondo_label"))
                if not lbl:
                    drops["disease: no MONDO id"] += 1
                    unnorm[f"DISEASE {txt}"] += 1
                elif key in DISEASE_IGNORE or lbl.strip().lower() in DISEASE_IGNORE:
                    drops["disease: generic (DISEASE_IGNORE)"] += 1
                else:
                    drops["disease: phenotype flag"] += 1
            elif typ == "CHEMICAL":
                if not single(e.get("chebi_label")):
                    drops["chemical: no ChEBI id"] += 1
                    unnorm[f"CHEMICAL {txt}"] += 1
                else:
                    drops["chemical: non_chemical flag"] += 1
    out = [f"  nodes kept: {', '.join(f'{k} {c:,}' for k, c in kept.most_common())}"]
    if overridden:
        out.append(f"  DISEASE_SURFACE overrides applied: {overridden:,} mention(s) "
                   f"({', '.join(sorted(DISEASE_SURFACE))})")
    if merged:
        out.append("  DISEASE_MERGE (one disease, one node): "
                   + ", ".join(f"{k} x{c}" for k, c in merged.most_common()))
    if drops:
        out.append("  endpoints dropped: " + ", ".join(f"{k} {c:,}" for k, c in drops.most_common()))
    if unnorm:
        out.append("  top un-normalized surfaces (fix upstream, not here): "
                   + ", ".join(f"{t} {s!r} x{c}" for (t, s), c in
                               ((tuple(k.split(" ", 1)), c) for k, c in unnorm.most_common(top))))
    return out


def qualifies_multi(t, T, flags):
    """Multi-type filter: score, both endpoints resolvable to distinct nodes, and the
    control rules applied ONLY where a gene is involved (a chemical-disease edge has no
    GENETIC endpoint to carry a control flag, so requiring one would drop every one)."""
    sc = t.get("score")
    if not (isinstance(sc, (int, float)) and sc >= T) or ctrl_yes(t) or has_mki67(t):
        return False
    a, b = node_of(t["subject"], flags), node_of(t["object"], flags)
    if not a or not b or a == b:
        return False
    if "GENETIC" in (t["subject"].get("type"), t["object"].get("type")) and not ctrl_no(t):
        return False        # same "annotated control:no" requirement as the gene-only filter
    return True


def stats(d, T):
    sub = [t for t in d if isinstance(t.get("score"), (int, float)) and t["score"] >= T]
    f = [t for t in sub if qualifies(t, T)]
    return {"T": T, "triples": len(sub), "sentences": len({t.get("sentence") for t in sub}),
            "control_no": sum(1 for t in sub if ctrl_no(t)),
            "control_yes": sum(1 for t in sub if ctrl_yes(t)),
            "mki67": sum(1 for t in sub if has_mki67(t)),
            "filt": len(f), "filt_sentences": len({t.get("sentence") for t in f})}


# ----- brain-cancer gene-gene graph ------------------------------------------
def single(v):
    return v.strip() if isinstance(v, str) and v.strip() else None


def _drug_targets():
    """(n_chemicals, chemicals, source, colour-class) per gene, from chemical_to_target.json:
    corpus genes flagged in_corpus_GENETIC, used to shade drug-target nodes."""
    try:
        c2t = json.loads(TARGET_FILE.read_text(encoding="utf-8"))
    except Exception:
        c2t = {}
    tgt, chems_by_gene, tsrc, tcat = {}, {}, {}, {}
    for g, v in c2t.items():
        if not v.get("in_corpus_GENETIC"):
            continue
        chems = v.get("chemicals") or []
        tgt[g] = v.get("n_chemicals") or len(chems)
        chems_by_gene[g] = sorted({c.get("chebi_label") for c in chems if c.get("chebi_label")})
        vrs = [str(r) for c in chems for r in (c.get("via_roles") or [])]
        has_db = any(r.startswith("DGIdb") for r in vrs)
        has_ch = any(not r.startswith("DGIdb") for r in vrs)
        tsrc[g] = "dgidb" if (has_db and not has_ch) else ("chebi+dgidb" if has_db else "chebi")
        # colour class: green = approved anti-neoplastic, amber = approved (other), pink = other target
        if any(r.startswith("DGIdb-antineoplastic") for r in vrs):
            tcat[g] = "green"
        elif any(r.startswith("DGIdb-approved") for r in vrs):
            tcat[g] = "amber"
        else:
            tcat[g] = "other"
    return tgt, chems_by_gene, tsrc, tcat


def graph_payload(triples):
    """Gene-gene (single HGNC symbol, non-self) graph; per-sentence scores so the
    in-browser confidence toggle can re-filter to >=0.99."""
    try:
        years = json.loads(PMC_YEARS.read_text(encoding="utf-8"))
    except Exception:
        years = {}
    tgt, chems_by_gene, tsrc, tcat = _drug_targets()
    # Edges are categorized by the RELATION the model predicted (activates / inhibits /
    # binds / interacts / associated, "not X" when negated), not merely by polarity: with
    # the BioRED checkpoint in the routing that label is signed, and dropping it here would
    # throw away the entire reason for training it.
    dir_sent = collections.defaultdict(set)           # (s,o,cat) -> sentences (direction + relation)
    pair_sent = collections.defaultdict(dict)          # pair -> {sentence: [maxscore, pmid, cat, spec, src]}
    pair_src = collections.defaultdict(set)            # pair -> model roles behind ANY of its triples
    node_sent = collections.defaultdict(dict)          # node -> {sentence: maxscore}
    for t in triples:
        s, o = single(t["subject"].get("hgnc_symbol")), single(t["object"].get("hgnc_symbol"))
        if not (s and o) or s == o:
            continue
        sc = float(t.get("score") or 0.0)
        cat = rel_cat(t)
        spec = t.get("modality") == "speculated"
        roles = model_roles(t)
        sent = t.get("sentence", "")
        pm = (t.get("pmid") or "?").replace(".grobid.tei", "")
        dir_sent[(s, o, cat)].add(sent)
        pair_src[frozenset((s, o))] |= roles
        cur = pair_sent[frozenset((s, o))].get(sent)
        if cur is None or sc > cur[0]:
            pair_sent[frozenset((s, o))][sent] = [sc, pm, cat, spec, src_tag(roles)]
        for nd in (s, o):
            if node_sent[nd].get(sent, -1) < sc:
                node_sent[nd][sent] = sc
    edges = []
    for pr, sd in pair_sent.items():
        cands = [(len(ss), f, to, cat) for (f, to, cat), ss in dir_sent.items() if frozenset((f, to)) == pr]
        cands.sort(reverse=True)
        _, ff, ft, fcat = cands[0]           # best-supported (direction, relation) wins the edge
        sents = [{"pmid": pm, "text": sent[:300], "sc": round(sc, 4), "yr": years.get(pm),
                  "rc": cat, "sp": sp, "sr": sr}
                 for sent, (sc, pm, cat, sp, sr) in sd.items()]
        sents.sort(key=lambda z: (-z["sc"], z["pmid"]))
        edges.append({"from": ff, "to": ft, "cat": fcat, "color": RCOLOR.get(fcat, "#888"),
                      "neg": fcat.startswith("not "), "src": src_tag(pair_src[pr]), "sents": sents})
    nodes = [{"id": nd, "label": nd, "kind": "gene", "shape": NODE_STYLE["gene"]["shape"],
              "bg": NODE_STYLE["gene"]["bg"], "border": NODE_STYLE["gene"]["border"],
              "sent95": len(sd), "sent99": sum(1 for v in sd.values() if v >= 0.99),
              "target": tgt.get(nd, 0), "chems": chems_by_gene.get(nd, []), "tsource": tsrc.get(nd, ""), "tcat": tcat.get(nd, "other")}
             for nd, sd in node_sent.items()]
    return {"nodes": nodes, "edges": edges}


def graph_payload_multi(triples, flags):
    """Gene + disease + chemical graph (--nodes all). Same edge model as the gene-only
    payload -- one edge per unordered node pair, best-supported (direction, relation) wins,
    one record per sentence -- but nodes are typed and identified by their normalized id."""
    try:
        years = json.loads(PMC_YEARS.read_text(encoding="utf-8"))
    except Exception:
        years = {}
    tgt, chems_by_gene, tsrc, tcat = _drug_targets()
    dir_sent = collections.defaultdict(set)            # (a,b,cat) -> sentences
    pair_sent = collections.defaultdict(dict)          # pair -> {sentence: [score, pmid, cat, spec, src]}
    pair_src = collections.defaultdict(set)            # pair -> model roles behind ANY of its triples
    node_sent = collections.defaultdict(dict)          # node -> {sentence: maxscore}
    kind_of = {}
    for t in triples:
        a, b = node_of(t["subject"], flags), node_of(t["object"], flags)
        if not a or not b or a == b:
            continue
        sc, cat = float(t.get("score") or 0.0), rel_cat(t)
        spec = t.get("modality") == "speculated"
        roles = model_roles(t)
        sent = t.get("sentence", "")
        pm = (t.get("pmid") or "?").replace(".grobid.tei", "")
        ka, kb = a[1], b[1]
        kind_of[ka], kind_of[kb] = a[0], b[0]
        dir_sent[(ka, kb, cat)].add(sent)
        pair_src[frozenset((ka, kb))] |= roles
        cur = pair_sent[frozenset((ka, kb))].get(sent)
        if cur is None or sc > cur[0]:
            pair_sent[frozenset((ka, kb))][sent] = [sc, pm, cat, spec, src_tag(roles)]
        for nd in (ka, kb):
            if node_sent[nd].get(sent, -1) < sc:
                node_sent[nd][sent] = sc
    edges = []
    for pr, sd in pair_sent.items():
        cands = [(len(ss), f, to, cat) for (f, to, cat), ss in dir_sent.items() if frozenset((f, to)) == pr]
        cands.sort(reverse=True)
        _, ff, ft, fcat = cands[0]
        sents = [{"pmid": pm, "text": sent[:300], "sc": round(sc, 4), "yr": years.get(pm),
                  "rc": cat, "sp": sp, "sr": sr}
                 for sent, (sc, pm, cat, sp, sr) in sd.items()]
        sents.sort(key=lambda z: (-z["sc"], z["pmid"]))
        edges.append({"from": ff, "to": ft, "cat": fcat, "color": RCOLOR.get(fcat, "#888"),
                      "neg": fcat.startswith("not "), "src": src_tag(pair_src[pr]), "sents": sents})
    nodes = []
    for nd, sd in node_sent.items():
        k = kind_of[nd]
        st = NODE_STYLE[k]
        nodes.append({"id": nd, "label": nd, "kind": k, "shape": st["shape"],
                      "bg": st["bg"], "border": st["border"],
                      "sent95": len(sd), "sent99": sum(1 for v in sd.values() if v >= 0.99),
                      # drug-target shading stays a GENE property (a disease has no targets)
                      "target": tgt.get(nd, 0) if k == "gene" else 0,
                      "chems": chems_by_gene.get(nd, []) if k == "gene" else [],
                      "tsource": tsrc.get(nd, "") if k == "gene" else "",
                      "tcat": tcat.get(nd, "other")})
    return {"nodes": nodes, "edges": edges}


def read_pubmed_query():
    """Pull the PubMed query out of summaries/pubmed_query.html (the step-1 publications
    summary), i.e. the value after 'Query (read strictly from <STDIN>):'. Looks under the
    data root first, then next to this script; returns '' if the file/line is absent."""
    for cand in (DATA_ROOT / "summaries" / "pubmed_query.html",
                 ROOT / "summaries" / "pubmed_query.html"):
        try:
            txt = cand.read_text(encoding="utf-8")
        except Exception:
            continue
        m = re.search(r"Query \(read strictly from &lt;STDIN&gt;\):\s*<code>(.*?)</code>", txt, re.S)
        if m:
            return html.unescape(m.group(1)).strip()
    return ""


def get_vis_lib():
    try:
        with urllib.request.urlopen(VIS_URL, timeout=30) as r:
            return r.read().decode("utf-8")
    except Exception as e:
        print(f"  [graph] could not fetch vis-network ({type(e).__name__}); HTML will use the CDN (needs internet to view)")
        return None


GRAPH_TEMPLATE = r"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
__LIBTAG__
<style>
 html,body{margin:0;height:100%;background:#ffffff;color:#1c2330;font:14px/1.5 Segoe UI,Arial,sans-serif}
 #net{position:absolute;top:0;left:0;right:0;bottom:0;background:#ffffff}
 /* max-height + overflow keep the panel inside the viewport: without them a tall control list runs off the
    bottom edge and that overflowing tail covers the graph with no way to scroll it back into reach. */
 #panel{position:absolute;top:12px;right:12px;z-index:5;background:rgba(255,255,255,.97);border:1px solid #cdd5e0;border-radius:10px;padding:14px 16px;max-width:320px;max-height:calc(100vh - 24px);overflow-y:auto;overscroll-behavior:contain;box-shadow:0 2px 12px rgba(0,0,0,.18);color:#1c2330}
 #panel h1{font-size:13px;margin:0 0 8px;color:#1c2330;font-variant:small-caps;letter-spacing:.4px}
 .row{margin:8px 0}
 input[type=range]{width:150px;max-width:100%;vertical-align:middle}
 .legend{display:flex;flex-wrap:wrap;align-items:center;gap:4px 12px}
 .legend b{display:inline-block;width:11px;height:11px;border-radius:2px;margin-right:5px;vertical-align:-1px;border:1px solid #999}
 .sw{display:inline-block;width:10px;height:10px;border-radius:2px;vertical-align:-1px}
 .mut{color:#5b6677;font-size:12px} b{color:#2b6cb0}
 #conf{width:190px;cursor:pointer}
 select,#search,#genefilter,#drugsearch,#textfilter{max-width:100%;background:#fff;border:1px solid #cdd5e0;color:#1c2330;border-radius:5px;padding:3px 6px;font-size:13px}
 #search,#genefilter,#drugsearch,#textfilter{width:200px}
 mark{background:#ffe680;color:inherit;border-radius:2px;padding:0 1px}
 #catfilters label{display:block;cursor:pointer;white-space:nowrap;font-size:12px;margin:1px 0}
 #catfilters{border:1px solid #cdd5e0;border-radius:6px;padding:4px 6px;max-height:140px;overflow:auto}
 #catfilters .cnt{color:#5b6677;font-size:11px}
 #zoom button,#srcbtns button{background:#eef2f7;color:#1c2330;border:1px solid #cdd5e0;border-radius:6px;padding:4px 10px;cursor:pointer;margin-right:6px;font-size:13px}
 #zoom button:hover,#srcbtns button:hover{background:#dde4ee}
 #srcbtns button{margin-bottom:4px}
 #srcbtns button.on{background:#0969da;border-color:#0969da;color:#fff;font-weight:600}
 #srcbtns button:disabled{opacity:.45;cursor:default}
 .vis-tooltip{max-width:480px!important;white-space:normal!important;background:#fff!important;color:#1a1a1a!important;border:1px solid #999!important;border-radius:8px!important;padding:8px 10px!important;box-shadow:0 4px 16px rgba(0,0,0,.35)!important;font:12px/1.45 Segoe UI,Arial,sans-serif!important}
 .eth{font-size:13px;margin-bottom:6px} .stip{padding:3px 0;border-top:1px solid #e3e3e3}
 .pm{display:inline-block;background:#eef3fb;color:#2b6cb0;border-radius:4px;padding:0 5px;margin-right:5px;font-weight:600;font-size:11px;text-decoration:none}
 a.pm:hover{background:#d6e6fb;text-decoration:underline} .more{margin-top:5px;color:#888;font-style:italic}
 #info{max-height:240px;overflow:auto} #info .stip{border-top:1px solid #e3e3e3}
 #toggle{position:absolute;top:12px;left:12px;z-index:6;background:#fff;color:#1c2330;border:1px solid #cdd5e0;border-radius:8px;padding:5px 11px;font-size:18px;line-height:1;cursor:pointer;box-shadow:0 2px 8px rgba(0,0,0,.18)}
 #panel.collapsed{display:none}
 @media (max-width:700px){
  #panel{left:12px;right:12px;max-width:none;max-height:62vh;overflow:auto;top:56px}
  .vis-tooltip{max-width:88vw!important}
 }
</style></head><body>
<button id="toggle" aria-label="Toggle controls">&#9776;</button>
<div id="panel">
 __PUBMED_QUERY__
 <h1>Gene&ndash;gene interactions</h1>
 <div class="row mut" id="pubinfo"></div>
 <div class="row legend"><b style="background:#cfe3ff;border-color:#2b6cb0"></b>gene <b style="background:#1b7837;border-color:#145a28"></b>approved anti-neoplastic <b style="background:#e08600;border-color:#9a6700"></b>approved (other) <b style="background:#c2185b;border-color:#7a0f3a"></b>ChEBI</div>
 <div class="row mut">Drug-target genes (corpus chemicals): deeper colour = more chemicals. <b style="color:#1b7837">Green</b> = DGIdb approved anti-neoplastic, <b style="color:#e08600">amber</b> = DGIdb approved (non-anti-neoplastic), <b style="color:#c2185b">pink</b> = ChEBI, absent in DGIdb.</div>
 <div class="row legend"><span class="sw" style="background:#2e9e5b"></span>activates <span class="sw" style="background:#e0533d"></span>inhibits <span class="sw" style="background:#3b7dd8"></span>binds <span class="sw" style="background:#6b3fa0"></span>interacts <span class="sw" style="background:#8a8f98"></span>associated <span class="sw" style="background:#d59a2e"></span>negated</div>
 <div class="row mut">Edge colour = the relation the model predicted. <b>activates</b>/<b>inhibits</b> are signed and come from the BioRED checkpoint; <b>interacts</b> is the unsigned PPI verdict. An edge takes its best-supported direction and relation; hover for the per-sentence labels.</div>
 <div class="row">Min unique sentences/edge: <b id="thv">1</b><br><input id="thr" type="range" min="1" max="10" value="1"></div>
 <div class="row">Min cluster size: <select id="mincluster"><option>1</option><option selected>2</option><option>3</option><option>4</option><option>5</option><option>6</option><option>7</option><option>8</option><option>9</option><option>10</option><option>11</option><option>12</option></select></div>
 <div class="row">Min connections: <select id="mindeg"><option selected>1</option><option>2</option><option>3</option><option>4</option></select><div class="mut">Hides genes linked to fewer than this many others; thins the hairball's single-link fringe.</div></div>
 <div class="row">Year: <b id="yrlab"></b><br><input id="yrlo" type="range" style="width:74px"> <input id="yrhi" type="range" style="width:74px"></div>
 <div class="row">Search gene: <input id="search" placeholder="e.g. EGFR" autocomplete="off"></div>
 <div class="row">Filter to gene:<br><input id="genefilter" placeholder="e.g. EGFR (+neighbors)" autocomplete="off"> <select id="hops"><option value="1">1 hop</option><option value="2">2 hops</option></select></div>
 <div class="row">Search drug: <input id="drugsearch" placeholder="e.g. nivolumab" autocomplete="off"></div>
 <div class="row">Filter to drug:<br><select id="chemfilter"><option value="">(all drugs)</option></select></div>
 <div class="row">Match text in sentence:<br><input id="textfilter" placeholder="e.g. phosphorylat or /inhibit(s|ed)?/" autocomplete="off">
  <div class="mut">Case-insensitive substring; wrap in / / for a regex. Keeps only edges with a matching sentence.</div></div>
__KINDROW__
 <div class="row mut">Relation type <span class="mut">(as predicted by the RE model; &ldquo;not X&rdquo; = negated statement, drawn dashed)</span>:</div><div id="catfilters"></div>
 <div class="row mut">Counts read <em>total &middot; in view</em>: the total is every edge of that type in the file, &ldquo;in view&rdquo; is how many survive the current score, year, text, min-connections and min-cluster settings. <span style="color:#b3243b">A red 0</span> means the type is ticked but everything of it is pruned &mdash; usually its edges sit in components smaller than <b>Min cluster size</b>, so lower that (or the score) to see them.</div>
 <div class="row mut">Training set behind the edge:</div>
 <div class="row" id="srcbtns"></div>
 <div class="row mut" id="srchint">Which corpus the relation was learned from &mdash; <b>PPI-only</b> = found by the BioInfer/PPI model alone, <b>BioRED-only</b> = by the BioRED model alone (typed and often signed), <b>both</b> = the two agreed a relation is there. Edges here can carry several sentences from different models; an edge counts as &ldquo;both&rdquo; if any of its support is corroborated.</div>
 <div class="row" id="zoom"><button id="zin">+ Zoom in</button><button id="zout">&minus; Zoom out</button><button id="zfit">Fit</button></div>
 <div class="row">Relationship score: <b id="scval">&ge;0.99</b><br><input id="conf" type="range" min="0" max="13" step="1" value="__CONFDEF__" aria-label="Minimum relationship score"></div>
 <div class="row mut" id="stats"></div>
 <div class="row mut" id="info">Click a node or edge for details.</div>
</div>
<div id="net"></div>
<script>
const DATA=__PAYLOAD__;
const CCOLOR=__CCOLOR__;
const MINY=__MINY__, MAXY=__MAXY__;
const MAXTGT=Math.max(1,...DATA.nodes.map(n=>n.target||0));
// Node colour: genes keep the drug-target shading (deeper = more corpus chemicals); disease
// and chemical nodes take their type colour, and every node its type SHAPE, so the three
// kinds stay distinguishable without relying on colour alone.
function nodeColor(n){if(n.kind&&n.kind!=='gene')return {background:n.bg,border:n.border};if(!n.target)return {background:n.bg||'#cfe3ff',border:n.border||'#2b6cb0'};const t=n.target/MAXTGT,L=(a,b)=>Math.round(a+(b-a)*t);if(n.tcat==='green')return {background:'rgb('+L(200,27)+','+L(230,120)+','+L(201,55)+')',border:'#145a28'};if(n.tcat==='amber')return {background:'rgb('+L(255,224)+','+L(231,134)+','+L(179,0)+')',border:'#9a6700'};return {background:'rgb('+L(255,194)+','+L(217,24)+','+L(232,91)+')',border:'#7a0f3a'};}
const KIND={};DATA.nodes.forEach(n=>{KIND[n.id]=n.kind||'gene';});
// --- training-set provenance ---------------------------------------------------------
// Each edge records which corpus produced it: 'ppi' (BioInfer only), 'biored' (BioRED only,
// i.e. relations the binary model never claimed) or 'both'. The buttons isolate each, which
// is how you SEE what a training set contributes rather than inferring it from counts.
let SRC_MODE='all';
const SRC_BTN=[['all','All'],['ppi','PPI-only'],['biored','BioRED-only'],['both','Both agreed']];
function buildSrcButtons(){
 const n={all:DATA.edges.length,ppi:0,biored:0,both:0};
 DATA.edges.forEach(e=>{n[e.src||'ppi']=(n[e.src||'ppi']||0)+1;});
 document.getElementById('srcbtns').innerHTML=SRC_BTN.map(([k,lab])=>
   '<button class="srcb'+(k===SRC_MODE?' on':'')+'" data-src="'+k+'"'+(n[k]?'':' disabled')+'>'
   +lab+' ('+(n[k]||0)+')</button>').join('');
 document.querySelectorAll('.srcb').forEach(b=>b.addEventListener('click',()=>{
   SRC_MODE=b.getAttribute('data-src');
   document.querySelectorAll('.srcb').forEach(x=>x.classList.toggle('on',x===b));
   build(+thr.value);}));
 // one checkpoint in the run -> nothing to separate; keep the row but say so
 if(!n.biored&&!n.both)document.getElementById('srchint').innerHTML=
   'Only one RE checkpoint stands behind this graph, so every edge is <b>PPI-only</b>; '
   +'run step&nbsp;2 with both models (<code>--route-mode additive</code>) to split them.';
}
function activeKinds(){const b=[...document.querySelectorAll('.kindf')];return b.length?new Set(b.filter(c=>c.checked).map(c=>c.value)):null;}
const net=document.getElementById('net'); let network=null;
// Labels fade in as you zoom: small graphs always show every symbol; dense views reveal labels as the
// zoom scale climbs from LABEL_LO to LABEL_HI. Opacity is driven through the shared node-font colour, so
// one setOptions call recolours all labels (per-node font carries only size, inheriting this colour).
const LABEL_SMALL=60, LABEL_LO=0.05, LABEL_HI=0.25;
let LABEL_N=0, LABEL_A=-1;
function labelOpacity(sc){ if(LABEL_N<=LABEL_SMALL)return 1; return Math.max(0,Math.min(1,(sc-LABEL_LO)/(LABEL_HI-LABEL_LO))); }
function updateLabels(){ if(!network)return; const a=Math.round(labelOpacity(network.getScale())*20)/20; if(a===LABEL_A)return; LABEL_A=a; network.setOptions({nodes:{font:{color:'rgba(26,26,26,'+a+')'}}}); }
const labelById={};DATA.nodes.forEach(n=>{labelById[n.id]=n.label;});
function esc(s){return (s||'').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');}
function pmA(p){return '<a class=pm target=_blank rel=noopener href="https://www.ncbi.nlm.nih.gov/pmc/articles/'+p+'/">'+p+'</a>';}
const SCORE_STEPS=[0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99]; // 0.05 up to 0.95, then 0.01
function activeConf(){const s=document.getElementById('conf');let i=s?parseInt(s.value):SCORE_STEPS.length-1;if(isNaN(i))i=SCORE_STEPS.length-1;return SCORE_STEPS[Math.max(0,Math.min(SCORE_STEPS.length-1,i))];}
function activeCats(){return new Set(Array.from(document.querySelectorAll(".catf:checked")).map(c=>c.value));}
function activeYears(){const a=parseInt(document.getElementById('yrlo').value),b=parseInt(document.getElementById('yrhi').value);return [Math.min(a,b),Math.max(a,b)];}
function passYear(yr,lo,hi){return (yr!=null&&yr>=lo&&yr<=hi)||(yr==null&&lo<=MINY&&hi>=MAXY);}
// free-text sentence filter: "foo bar" = case-insensitive substring, "/foo(bar)?/" = regex
// (an unparseable regex falls back to a literal substring match, so typing is never an error)
function activeText(){return (document.getElementById('textfilter').value||'').trim();}
function reEsc(s){return s.replace(/[.*+?^${}()|[\]\\]/g,'\\$&');}
function textMatcher(q){
 if(!q)return null;
 const m=/^\/(.*)\/([a-z]*)$/.exec(q);
 let src=null,flags='i';
 if(m){try{new RegExp(m[1],m[2]);src=m[1];flags=(m[2].indexOf('i')>=0?m[2]:m[2]+'i').replace(/g/g,'');}catch(err){src=null;}}
 if(src===null)src=reEsc(q);
 const re=new RegExp(src,flags);
 return {test:t=>re.test(t||''),hlre:new RegExp(src,flags+'g')};
}
let TM=null;   // matcher in force for the current view; used to highlight hits in sentence text
// mark hits on the RAW text (so a query containing <, > or & still highlights), escaping each piece as we go
function hl(t){
 t=t||'';
 if(!TM)return esc(t);
 const re=TM.hlre;re.lastIndex=0;
 let out='',last=0,m;
 while((m=re.exec(t))!==null){
  if(!m[0].length){re.lastIndex++;continue;}   // zero-length match (e.g. /x*/): skip, never loop
  out+=esc(t.slice(last,m.index))+'<mark>'+esc(m[0])+'</mark>';
  last=m.index+m[0].length;
 }
 return out+esc(t.slice(last));
}
function visSents(e,conf,lo,hi,tm){return e.sents.filter(s=>s.sc>=conf&&passYear(s.yr,lo,hi)&&(!tm||tm.test(s.text)));}
const SRCLAB={ppi:'PPI',biored:'BioRED',both:'PPI+BioRED'};
function edgeHead(e,vis){const np=new Set(vis.map(s=>s.pmid)).size;return '<div class=eth><b>'+esc(labelById[e.from])+' &rarr; '+esc(labelById[e.to])+'</b> ('+vis.length+' sentences &middot; '+np+' PMIDs &middot; '+e.cat+' &middot; '+(SRCLAB[e.src]||'PPI')+')</div>';}
// per-sentence relation tag: the label the model gave THIS sentence (an edge shows its
// dominant relation, so a minority reading -- e.g. one "inhibits" under an "interacts"
// edge -- would otherwise be invisible). "?" marks a speculated statement.
function relTag(s){return s.rc?' <span class=mut style="color:'+(CCOLOR[s.rc]||'#888')+'">['+esc(s.rc)+(s.sp?' ?':'')+']</span>'+(s.sr?' <span class=mut>'+esc(SRCLAB[s.sr]||s.sr)+'</span>':''):'';}
function edgeTip(e,vis){const d=document.createElement('div');let h=edgeHead(e,vis);const lim=20;vis.slice(0,lim).forEach(s=>{h+='<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span>'+relTag(s)+' '+hl(s.text)+'</div>';});if(vis.length>lim)h+='<div class=more>+'+(vis.length-lim)+' more</div>';d.innerHTML=h;return d;}
function scaleNode(s){return 6+Math.sqrt(s)*3.4;}
function fontSize(c){return c<5?13:2*Math.max(13,Math.min(Math.round(c*2.2),48));}
// Which node kinds have their NAME DRAWN on the canvas. Disease and chemical names are long,
// repeat across many edges and out-shout the gene symbols simply by being wordy, so those
// nodes are drawn unlabelled -- shape and colour say what they are, and the name is one hover
// (or click) away in the tooltip and the info panel. Add 'disease'/'chemical' here to get the
// drawn names back.
const LABEL_KINDS=new Set(['gene']);
function nodeLabel(n){return LABEL_KINDS.has(n.kind||'gene')?n.label:'';}
function activeMinCluster(){const v=parseInt((document.getElementById('mincluster')||{}).value);return isNaN(v)?2:v;}
function activeMinDegree(){const v=parseInt((document.getElementById('mindeg')||{}).value);return isNaN(v)?1:v;}
function build(thr){
 const conf=activeConf(), cats=activeCats(); const [ylo,yhi]=activeYears(); const mc=activeMinCluster(); const md=activeMinDegree();
 const txt=activeText(), tm=textMatcher(txt); TM=tm;
 const kinds=activeKinds();
 let edges=[];
 DATA.edges.forEach(e=>{ if(!cats.has(e.cat))return; if(SRC_MODE!=='all'&&(e.src||'ppi')!==SRC_MODE)return; if(kinds&&!(kinds.has(KIND[e.from])&&kinds.has(KIND[e.to])))return; const vis=visSents(e,conf,ylo,yhi,tm); if(vis.length>=thr) edges.push({e:e,vis:vis,w:vis.length}); });
 const gf=(document.getElementById('genefilter').value||'').trim().toLowerCase();
 const chemSel=document.getElementById('chemfilter').value;
 let focusActive=false, focusLabel='';
 if(gf||chemSel){
   focusActive=true;
   const seeds=new Set(), labs=[];
   if(gf){const fn=DATA.nodes.find(n=>n.label.toLowerCase()===gf)||DATA.nodes.find(n=>n.label.toLowerCase().indexOf(gf)===0);if(fn){seeds.add(fn.id);labs.push(fn.label);}else labs.push('(no gene: '+gf+')');}
   if(chemSel){DATA.nodes.forEach(n=>{if((n.chems||[]).indexOf(chemSel)>=0)seeds.add(n.id);});labs.push('chem: '+chemSel);}
   const hops=parseInt(document.getElementById('hops').value)||1;
   const ag={};edges.forEach(o=>{(ag[o.e.from]=ag[o.e.from]||[]).push(o.e.to);(ag[o.e.to]=ag[o.e.to]||[]).push(o.e.from);});
   const seen=new Set(seeds);let fr=[...seeds];
   for(let h=0;h<hops;h++){const nf=[];fr.forEach(x=>{(ag[x]||[]).forEach(y=>{if(!seen.has(y)){seen.add(y);nf.push(y);}});});fr=nf;}
   edges=edges.filter(o=>seen.has(o.e.from)&&seen.has(o.e.to));
   focusLabel=labs.join(', ');
 }
 // min-cluster prunes small connected components; it stays live under a text/year/etc. filter (set it
 // to 1 to see every match), and is skipped only under gene/chem focus where you want the neighborhood
 if(!focusActive){
   if(md>1){ // single-pass degree filter: measure each gene's links once, drop the ones below md (peels the fringe)
     const deg={};edges.forEach(o=>{deg[o.e.from]=(deg[o.e.from]||0)+1;deg[o.e.to]=(deg[o.e.to]||0)+1;});
     edges=edges.filter(o=>deg[o.e.from]>=md&&deg[o.e.to]>=md);
   }
   const adj={};
   edges.forEach(o=>{(adj[o.e.from]=adj[o.e.from]||[]).push(o.e.to);(adj[o.e.to]=adj[o.e.to]||[]).push(o.e.from);});
   const comp={};let cid=0;
   for(const n in adj){if(comp[n]!==undefined)continue;const sk=[n];comp[n]=cid;while(sk.length){const x=sk.pop();(adj[x]||[]).forEach(y=>{if(comp[y]===undefined){comp[y]=cid;sk.push(y);}});}cid++;}
   const csz={};for(const n in comp)csz[comp[n]]=(csz[comp[n]]||0)+1;
   edges=edges.filter(o=>csz[comp[o.e.from]]>=mc);
 }
 updateCatCounts(edges,cats);   // edges is final here (category, score, year, text, degree, cluster)
 const keep=new Set();edges.forEach(o=>{keep.add(o.e.from);keep.add(o.e.to);});
 const nss={};edges.forEach(o=>{o.vis.forEach(s=>{(nss[o.e.from]=nss[o.e.from]||new Set()).add(s.text);(nss[o.e.to]=nss[o.e.to]||new Set()).add(s.text);});});
 const nsz=id=>(nss[id]?nss[id].size:0);
 const allCatsSel=[...new Set(DATA.edges.map(e=>e.cat))].every(c=>cats.has(c));
 const nodes=DATA.nodes.filter(n=>keep.has(n.id)).map(n=>({id:n.id,label:nodeLabel(n),value:nsz(n.id),size:scaleNode(nsz(n.id)),shape:n.shape||'dot',title:n.label+((n.kind&&n.kind!=='gene')?'  ['+n.kind+']':'')+' — '+nsz(n.id)+' unique sentences (in view)'+(n.target?' · drug target: '+n.target+' chemicals'+(n.tcat==='green'?' (approved anti-neoplastic)':(n.tcat==='amber'?' (approved)':' (ChEBI)')):''),color:nodeColor(n),font:{size:fontSize(nsz(n.id))}}));
 const eds=edges.map((o,i)=>({id:i,from:o.e.from,to:o.e.to,value:o.w,width:Math.min(1+o.w*0.7,10),color:{color:o.e.color,opacity:0.6},dashes:!!o.e.neg,title:edgeTip(o.e,o.vis)}));
 const vpub=new Set();edges.forEach(o=>o.vis.forEach(s=>vpub.add(s.pmid)));
 const nkinds=new Set(nodes.map(n=>KIND[n.id]));
 document.getElementById('stats').innerHTML='Showing <b>'+nodes.length+'</b> '+(nkinds.size>1?'nodes':'genes')+', <b>'+eds.length+'</b> edges, <b>'+vpub.size+'</b> publications (&ge;'+conf+')'+(txt?' &middot; text: <b>'+esc(txt)+'</b>':'')+(focusActive?' &middot; focus: <b>'+esc(focusLabel)+'</b>':'');
 const data={nodes:new vis.DataSet(nodes),edges:new vis.DataSet(eds)};
 const options={layout:{improvedLayout:false},physics:{stabilization:{iterations:200},barnesHut:{gravitationalConstant:-14000,springLength:130,springConstant:0.02,avoidOverlap:0.3}},interaction:{hover:true,tooltipDelay:120},nodes:{shape:'dot',scaling:{min:6,max:60},font:{color:'rgba(26,26,26,0)'}},edges:{smooth:false,arrowStrikethrough:false,hoverWidth:0,selectionWidth:0,arrows:{to:{enabled:true,scaleFactor:0.6}}}};
 if(network)network.destroy();
 network=new vis.Network(net,data,options);
 LABEL_N=nodes.length; LABEL_A=-1;
 network.on('stabilizationIterationsDone',()=>{network.setOptions({physics:false});network.fit({animation:false});LABEL_A=-1;updateLabels();network.redraw();});
 network.on('zoom',updateLabels);
 network.on('animationFinished',updateLabels);
 const _e=edges;
 network.on('click',p=>{const info=document.getElementById('info');
   if(p.nodes.length){const n=DATA.nodes.find(x=>x.id===p.nodes[0]);info.innerHTML='<b>'+n.label+'</b>: '+nsz(n.id)+' unique sentences (in view)';}
   else if(p.edges.length){const o=_e[p.edges[0]];info.innerHTML=edgeHead(o.e,o.vis)+o.vis.map(s=>'<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span> '+hl(s.text)+'</div>').join('');}});
}
const thr=document.getElementById('thr');
let CATTOT={};
function buildCatFilters(){CATTOT={};DATA.edges.forEach(e=>CATTOT[e.cat]=(CATTOT[e.cat]||0)+1);const cats=Object.keys(CATTOT).sort((a,b)=>CATTOT[b]-CATTOT[a]);document.getElementById('catfilters').innerHTML=cats.map(c=>'<label><input type=checkbox class=catf value="'+esc(c)+'" checked> <span class=sw style="background:'+(CCOLOR[c]||'#888')+'"></span> '+esc(c)+' <span class=cnt data-cat="'+esc(c)+'">('+CATTOT[c]+')</span></label>').join('');document.querySelectorAll('.catf').forEach(c=>c.addEventListener('change',()=>build(+thr.value)));}
// Counts are LIVE: "(total · N in view)" is recomputed from the edges actually drawn. The
// total is the whole payload (everything qualifying at >=0.5), which is NOT what you see --
// the score slider, year range, text filter, min-connections and min-cluster size all prune
// afterwards. A category whose edges are all pruned now reads 0 (in red) instead of looking
// available: that is the case where ticking it alone leaves the canvas blank, typically
// because its edges form components smaller than "Min cluster size".
function updateCatCounts(edges,cats){const seen={};edges.forEach(o=>seen[o.e.cat]=(seen[o.e.cat]||0)+1);
 document.querySelectorAll('#catfilters .cnt').forEach(el=>{const c=el.getAttribute('data-cat');
  el.textContent=cats.has(c)?('('+CATTOT[c]+' · '+(seen[c]||0)+' in view)'):('('+CATTOT[c]+' · off)');
  el.style.color=(cats.has(c)&&!seen[c])?'#b3243b':'';});}
thr.addEventListener('input',()=>{document.getElementById('thv').textContent=thr.value;build(+thr.value);});
const mcl=document.getElementById('mincluster');
mcl.addEventListener('change',()=>build(+thr.value));
document.getElementById('mindeg').addEventListener('change',()=>build(+thr.value));
const confEl=document.getElementById('conf');function updScore(){document.getElementById('scval').innerHTML='&ge;'+activeConf().toFixed(2);}updScore();confEl.addEventListener('input',()=>{updScore();build(+thr.value);});
document.getElementById('genefilter').addEventListener('change',()=>build(+thr.value));
document.getElementById('genefilter').addEventListener('keydown',ev=>{if(ev.key==='Enter')build(+thr.value);});
document.getElementById('hops').addEventListener('change',()=>build(+thr.value));
const txtEl=document.getElementById('textfilter');let txtTimer=null;   // debounce: each keystroke would otherwise rebuild the whole network
txtEl.addEventListener('input',()=>{clearTimeout(txtTimer);txtTimer=setTimeout(()=>build(+thr.value),350);});
txtEl.addEventListener('keydown',ev=>{if(ev.key==='Enter'){clearTimeout(txtTimer);build(+thr.value);}});
const searchBox=document.getElementById('search');
function doSearch(q){q=(q||'').trim();const info=document.getElementById('info');if(!q)return;const hit=DATA.nodes.find(n=>n.label.toLowerCase()===q.toLowerCase())||DATA.nodes.find(n=>n.label.toLowerCase().indexOf(q.toLowerCase())===0);if(!hit){info.innerHTML='No gene matching "'+q+'"';return;}try{network.selectNodes([hit.id]);network.focus(hit.id,{scale:1.3,animation:true});info.innerHTML='<b>'+hit.label+'</b>';}catch(e){info.innerHTML='<b>'+hit.label+'</b> not in current view';}}
searchBox.addEventListener('keydown',ev=>{if(ev.key==='Enter')doSearch(searchBox.value);});
searchBox.addEventListener('change',()=>doSearch(searchBox.value));
function zoomBy(f){if(!network)return;const s=network.getScale();network.moveTo({scale:s*f,animation:{duration:200}});}
document.getElementById('zin').addEventListener('click',()=>zoomBy(1.25));
document.getElementById('zout').addEventListener('click',()=>zoomBy(0.8));
document.getElementById('zfit').addEventListener('click',()=>{if(network)network.fit({animation:true});});
const yl=document.getElementById('yrlo'),yh=document.getElementById('yrhi');
[yl,yh].forEach(el=>{el.min=MINY;el.max=MAXY;});yl.value=MINY;yh.value=MAXY;
function updYr(){const a=activeYears();document.getElementById('yrlab').textContent=a[0]+'–'+a[1];}
updYr();
yl.addEventListener('input',()=>{updYr();build(+thr.value);});
yh.addEventListener('input',()=>{updYr();build(+thr.value);});
const chemGenes={};DATA.nodes.forEach(n=>(n.chems||[]).forEach(c=>{(chemGenes[c]=chemGenes[c]||[]).push(n.label);}));const csel=document.getElementById('chemfilter');Object.keys(chemGenes).sort().forEach(c=>{const g=chemGenes[c].slice().sort();const o=document.createElement('option');o.value=c;o.textContent=c+' → '+g.join(', ');csel.appendChild(o);});csel.addEventListener('change',()=>build(+thr.value));
const drugBox=document.getElementById('drugsearch');function findDrug(q){q=(q||'').trim().toLowerCase();if(!q)return;const info=document.getElementById('info');const opts=[...csel.options].filter(o=>o.value);const m=opts.find(o=>o.value.toLowerCase()===q)||opts.find(o=>o.value.toLowerCase().indexOf(q)===0)||opts.find(o=>o.value.toLowerCase().indexOf(q)>=0);if(m){csel.value=m.value;build(+thr.value);info.innerHTML='Drug filter: <b>'+esc(m.value)+'</b>';}else{info.innerHTML='No drug matching "'+esc(q)+'"';}}drugBox.addEventListener('keydown',ev=>{if(ev.key==='Enter')findDrug(drugBox.value);});drugBox.addEventListener('change',()=>findDrug(drugBox.value));
document.getElementById('toggle').addEventListener('click',()=>document.getElementById('panel').classList.toggle('collapsed'));
if(window.innerWidth<=700)document.getElementById('panel').classList.add('collapsed');
window.addEventListener('resize',()=>{if(network)network.redraw();});
(function(){const pm=new Set();DATA.edges.forEach(e=>e.sents.forEach(s=>pm.add(s.pmid)));document.getElementById('pubinfo').innerHTML='<b>'+pm.size+'</b> out of <b>'+__NXML__+'</b> produced high-score gene-gene interactions';})();
document.querySelectorAll('.kindf').forEach(c=>c.addEventListener('change',()=>build(+thr.value)));
buildCatFilters();buildSrcButtons();build(1);
</script></body></html>"""


KIND_ROW = (' <div class="row mut">Node type:</div>\n'
            ' <div class="row legend" id="kindfilters">'
            '<label><input type=checkbox class=kindf value="gene" checked> '
            '<span class="sw" style="background:#cfe3ff;border:1px solid #2b6cb0"></span> gene</label> '
            '<label><input type=checkbox class=kindf value="disease" checked> '
            '<span class="sw" style="background:#ffe0e0;border:1px solid #b3243b;'
            'transform:rotate(45deg)"></span> disease</label> '
            '<label><input type=checkbox class=kindf value="chemical" checked> '
            '<span class="sw" style="background:#ece0f8;border:1px solid #6b3fa0"></span> chemical</label>'
            '</div>\n'
            ' <div class="row mut">Shapes: gene &#9679; &middot; disease &#9670; &middot; chemical &#9632;. '
            'Only <b>gene</b> names are drawn on the canvas &mdash; disease and chemical names are long '
            'and would bury the symbols, so <b>hover</b> (or click) those nodes to read them. '
            'An edge is shown only when BOTH its endpoint types are ticked.</div>')


def render_graph(payload, lib, miny, maxy, nxml, pubmed_query="", multi=False):
    if lib:
        libtag = "<script>\n" + lib.replace("</script>", "<\\/script>") + "\n</script>"
    else:
        libtag = f'<script src="{VIS_URL}"></script>'
    # legend line above the title: the PubMed query this corpus came from (empty -> nothing shown)
    qrow = (f'<div class="row" id="pubmedq" style="font-size:12px;line-height:1.35;font-weight:600;'
            f'color:#2b6cb0;word-break:break-word;margin-bottom:6px">PubMed query = '
            f'{html.escape(pubmed_query)}</div>'
            if pubmed_query else "")
    title = html.escape(pubmed_query) if pubmed_query else (
        "High-confidence gene / disease / chemical relations" if multi
        else "High-confidence brain-cancer gene-gene interactions (gene-only filter)")
    # the wider view's edges sit mostly in 0.5-0.8 (BioRED gene-disease relations), so opening
    # it at the gene-only default of >=0.99 would show an almost empty canvas
    return (GRAPH_TEMPLATE.replace("__LIBTAG__", libtag)
            .replace("__TITLE__", title)
            .replace("__KINDROW__", KIND_ROW if multi else "")
            .replace("__CONFDEF__", "6" if multi else "13")
            .replace("__PUBMED_QUERY__", qrow)
            .replace("__PAYLOAD__", json.dumps(payload, ensure_ascii=False))
            .replace("__CCOLOR__", json.dumps(RCOLOR))
            .replace("__NXML__", str(nxml))
            .replace("__MINY__", str(miny)).replace("__MAXY__", str(maxy)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=str(DATA_ROOT),
                    help="pipeline output tree to read inputs from (default: kaggle_working/ next to this script)")
    ap.add_argument("--score", type=float, default=0.8, help="threshold for the exported JSON (default 0.8)")
    ap.add_argument("--thresholds", default="0.8,0.95,0.99", help="thresholds summarized in the console stats")
    ap.add_argument("--no-graph", action="store_true", help="skip the brain-cancer graph HTML")
    ap.add_argument("--merge", choices=MERGE_POLICIES, default="union",
                    help="how to combine several RE checkpoints scoring the same pair "
                         "(relation_extraction.py --route-mode additive): union (default) keeps "
                         "every pair either model kept and prefers the typed label; gate is the "
                         "stricter variant that keeps only pairs the binary model also kept; "
                         "typed uses the BioRED-style model alone; none disables the merge")
    ap.add_argument("--typed-model", default=None,
                    help="checkpoint name supplying the typed/signed labels (default: the one "
                         "whose name contains 'biored')")
    ap.add_argument("--gate-model", default=None,
                    help="checkpoint name used as the existence gate (default: the other one)")
    ap.add_argument("--nodes", choices=["gene", "all"], default="gene",
                    help="gene (default): the gene-gene graph, '_G' outputs. all: also draw "
                         "DISEASE and CHEMICAL nodes -- the gene-disease / chemical-gene / "
                         "chemical-disease / chemical-chemical edges the BioRED model adds, "
                         "which the gene-only view discards; writes '_M' outputs instead")
    args = ap.parse_args()

    set_data_root(args.data_root)
    multi = args.nodes == "all"
    if multi:                        # '_M' names: never clobber the gene-only outputs
        globals()["JSON_OUT"] = OUT_DIR / "high_confidence_M.json"
        globals()["GRAPH_OUT"] = DATA_ROOT / "summaries" / "high_confidence_M.html"
    flags = load_type_flags() if multi else None
    keep_fn = (lambda t, T: qualifies_multi(t, T, flags)) if multi else qualifies
    if not RE_FILE.exists():
        raise SystemExit(f"ERROR: {RE_FILE} not found under data root {DATA_ROOT} "
                         f"(run the gpu_bundle pipeline first, or pass --data-root).")

    d = json.loads(RE_FILE.read_text(encoding="utf-8"))
    d, merge_note = merge_models(d, args.merge, args.typed_model, args.gate_model)
    print(merge_note)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) JSON export at --score
    kept = [t for t in d if keep_fn(t, args.score)]
    JSON_OUT.write_text(json.dumps(kept, ensure_ascii=False), encoding="utf-8")
    ksents = len({t.get("sentence") for t in kept})
    print(f"exported {len(kept):,} triples ({ksents:,} sentences) at score>={args.score} -> {JSON_OUT}")

    # 2) threshold stats (for the console summary below)
    rows = [stats(d, float(x)) for x in args.thresholds.split(",")]

    # 3) the graph (universe = qualifying at GRAPH_BASE; in-browser score slider 0.5..0.99)
    if not args.no_graph:
        universe = kept if args.score <= GRAPH_BASE else [t for t in d if keep_fn(t, GRAPH_BASE)]
        payload = graph_payload_multi(universe, flags) if multi else graph_payload(universe)
        yrs = [s["yr"] for e in payload["edges"] for s in e["sents"] if s.get("yr")]
        miny, maxy = (min(yrs), max(yrs)) if yrs else (2000, 2026)
        lib = get_vis_lib()
        # corpus size = # source documents; experimental_ner/ is empty in the bundle
        # (it was a runtime symlink), so fall back to the per-document sentences/ files.
        nxml = len(list(XML_DIR.glob("*.xml"))) or len(list(SENT_DIR.glob("*.json")))
        GRAPH_OUT.parent.mkdir(parents=True, exist_ok=True)
        pubmed_query = read_pubmed_query()
        GRAPH_OUT.write_text(render_graph(payload, lib, miny, maxy, nxml, pubmed_query, multi),
                             encoding="utf-8")
        n99 = sum(1 for e in payload["edges"] if any(s["sc"] >= 0.99 for s in e["sents"]))
        npubs = len({s["pmid"] for e in payload["edges"] for s in e["sents"]})
        kinds = collections.Counter(n.get("kind", "gene") for n in payload["nodes"])
        what = ", ".join(f"{c:,} {k}" for k, c in kinds.most_common()) if multi else \
               f"{len(payload['nodes']):,} genes"
        print(f"wrote graph -> {GRAPH_OUT}  ({what}; {len(payload['edges']):,} edges; "
              f"{n99:,} edges have a >=0.99 sentence; {npubs:,} of {nxml:,} input XMLs produced triples)")
        cats = collections.Counter(e["cat"] for e in payload["edges"])
        signed = sum(c for k, c in cats.items()
                     if (k[4:] if k.startswith("not ") else k) in ("activates", "inhibits"))
        print(f"  edge relations: {', '.join(f'{k} {c:,}' for k, c in cats.most_common())}"
              f"  ({signed:,} signed)")
        if multi:
            # report over everything at the graph threshold, NOT `universe`: the latter only
            # contains triples whose endpoints already resolved, so it can never show a drop
            for line in node_drop_report(
                    [t for t in d if isinstance(t.get("score"), (int, float))
                     and t["score"] >= GRAPH_BASE], flags):
                print(line)
            kind_of = {n["id"]: n.get("kind", "gene") for n in payload["nodes"]}
            tp = collections.Counter(tuple(sorted((kind_of[e["from"]], kind_of[e["to"]])))
                                     for e in payload["edges"])
            print(f"  edges by node-type pair: "
                  f"{', '.join(f'{a}-{b} {c:,}' for (a, b), c in tp.most_common())}")

        # also drop a copy in the root dir, named "<current directory>_YYYY_MM_DD_<G|M>.html"
        # (illegal filename characters underscored so the name is always valid).
        dirname = re.sub(r'[\\/:*?"<>|\s]+', "_", ROOT.name)
        today = datetime.date.today().strftime("%Y_%m_%d")
        dest = ROOT / f"{dirname}_{today}_{'M' if multi else 'G'}.html"
        shutil.copy2(GRAPH_OUT, dest)
        print(f"copied graph -> {dest}")

    for r in rows:
        print(f"  score>={r['T']}: {r['triples']:,} triples / {r['sentences']:,} sent; "
              f"G {r['filt']:,} triples / {r['filt_sentences']:,} sent"
              + (f"; M {sum(1 for t in d if keep_fn(t, r['T'])):,} triples" if multi else ""))


if __name__ == "__main__":
    main()
