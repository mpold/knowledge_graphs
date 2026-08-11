#!/usr/bin/env python3
"""
high_confidence_g.py -- extract the high-confidence Gene-context relation triples from
the scored RE output, summarize the threshold statistics, and draw the gene / disease /
chemical relationship graph.

This is the "G" (gene-only) variant of high_confidence.py -- and, since high_confidence.py
was DEPRECATED, the maintained step-3 entry point. Where high_confidence.py applies the
"G_D_C" filter (gene-gene IN a disease/chemical context), this step applies the "G" filter
(gene-gene, context-agnostic). Its "_M" output names never clobber the other's, so both can
still target the same data root.

The "G" filter still defines the gene-gene slice reported per threshold in the console
summary, but it no longer has a graph of its own: the gene-only view ("_G" outputs) was
retired, and the multi-type graph below is the only one drawn.

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
          <data-root>/DISEASE/disease.json, <data-root>/CHEMICAL/chemical.json  (phenotype / non_chemical flags)
Outputs : <data-root>/TRIPLES/high_confidence_M.json     qualifying triples at --score (default 0.8)
          <data-root>/summaries/high_confidence_M.html    the graph, with a 0.5..0.99 in-browser score slider
                                                          and a "match text in sentence" filter (substring, or /regex/)
          <root>/<current_dir>_YYYY_MM_DD_M.html          a copy of that graph, named after the directory
                                                           holding this script plus today's date;
                                                           e.g. lung_large_2026_07_19_M.html

The output filenames all differ from those written by high_confidence.py (which uses the
"_G_D_C" JSON, "high_confidence.html" graph, and unsuffixed "<dir>_<date>.html" copy) so the
two scripts can be run against the same data root without clobbering each other's outputs.

  * The graph spans gene, DISEASE and CHEMICAL nodes. A gene-only view draws an edge only
    when BOTH endpoints carry a single HGNC symbol, which discards everything BioRED adds
    beyond gene-gene -- on the reference run, 649 of its 655 solo pairs -- and that is why
    it was retired. DISEASE and CHEMICAL endpoints are nodes too, identified by their
    normalized id (mondo_label / chebi_label), giving gene-disease / chemical-gene /
    chemical-disease / chemical-chemical edges. Node hygiene reuses the pipeline's own
    annotations: control (genes), phenotype (DISEASE/disease.json) and non_chemical
    (CHEMICAL/chemical.json), plus DISEASE_IGNORE for labels too generic to be a useful
    node and CHEMICAL_IGNORE for chemicals that are not compounds under study. Outputs
    take an "_M" suffix so they never clobber high_confidence.py's, node shape encodes the
    type (gene dot / disease diamond / chemical square), a "Node type" filter appears in
    the panel, and the score slider opens at >=0.8 (most gene-disease edges sit in 0.5-0.8,
    so the old gene-only default of >=0.99 would show an almost empty canvas).

Run::  python high_confidence_g.py [--data-root kaggle_working] [--score 0.8] [--thresholds 0.8,0.95,0.99] [--no-graph]
       python high_confidence_g.py --merge gate         # drop pairs only BioRED claimed
       python high_confidence_g.py --merge none         # pre-merge behaviour (duplicates survive)
"""
import argparse
import bisect
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
# every lung_* project caches its slice of the shared corpus under this name
CONTRIB_NAME = "corpus_contrib.json"
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
    # normalization libraries carrying the in-place phenotype / non_chemical flags
    DISEASE_LIB = DATA_ROOT / "DISEASE" / "disease.json"
    CHEM_LIB = DATA_ROOT / "CHEMICAL" / "chemical.json"
    # "_M" output names, distinct from high_confidence.py's "_G_D_C"/"high_confidence.html".
    JSON_OUT = OUT_DIR / "high_confidence_M.json"
    GRAPH_OUT = DATA_ROOT / "summaries" / "high_confidence_M.html"


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


# ----- multi-type nodes -------------------------------------------------------
# A gene-only graph draws an edge only when BOTH endpoints carry a single HGNC symbol,
# which silently discards everything the BioRED checkpoint adds beyond gene-gene:
# gene-disease, chemical-gene, chemical-disease, chemical-chemical. This graph keeps
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

# Chemical surfaces/labels that ChEBI resolves correctly but that are not COMPOUNDS UNDER STUDY,
# and so do not belong in a list the panel calls "drugs".
#   tyrosine  an amino acid, and in this corpus never a treatment: its 43 edges come from
#             "tyrosine phosphorylation", "phosphorylation at tyrosine 759", "tyrosine kinase"
#             -- the residue being named, which the NER hands over as a CHEMICAL mention.
#   glucose   87 edges of culture condition and metabolism: "under low glucose", "glucose
#             starvation", "glucose-induced up-regulation", "gene sets of glucose metabolism".
#             The cells' medium, not an agent anyone administered.
# The match is exact on the surface or the ChEBI label, so the compounds that merely contain
# these names survive: "tyrosine kinase inhibitor", "O(4)-phosphotyrosine", "2-deoxy-D-glucose"
# (a real glycolysis inhibitor) all stay nodes.
CHEMICAL_IGNORE = {"tyrosine", "glucose"}

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
#   none           no merge -- duplicates survive into graph_payload_multi(), which then
#                  collapses them per (pair, sentence) by max score
#
# With only one checkpoint in the file every policy is a no-op.
MERGE_POLICIES = ("gate", "union", "typed", "none")

# How to decide WHICH model's verdict wins when several scored the same pair. The policies
# above choose which pairs survive; this chooses whose number represents the survivor.
#
#   rank (default)  compare each model's score to that model's OWN distribution and let the
#                   higher percentile win. Two checkpoints calibrated against different
#                   corpora do not share a scale even after calibration -- on the reference
#                   run a raw max() handed PPI 65 of 69 corroborated pairs, because its
#                   isotonic map saturated at 1.0 while BioRED's Platt fit tops out at 0.90.
#                   That is an artefact of the calibrators, not evidence about the pair.
#   calibrated      the previous behaviour: plain max() of the calibrated probabilities.
#                   Use it to reproduce pre-2026-08 runs.
#
# Either way the surviving `score` is the winner's own calibrated probability, so it stays a
# probability and every downstream threshold keeps its meaning. Only the CHOICE changes.
MERGE_SCALES = ("rank", "calibrated")


def model_percentiles(triples):
    """(pct(model, score), n_by_model) -- empirical percentile of a score within the
    distribution of scores that model wrote across the whole file."""
    dist = collections.defaultdict(list)
    for t in triples:
        m = (t.get("predicate") or {}).get("model")
        s = t.get("score")
        if m and isinstance(s, (int, float)):
            dist[m].append(float(s))
    for v in dist.values():
        v.sort()

    def pct(model, score):
        v = dist.get(model)
        return bisect.bisect_right(v, float(score)) / len(v) if v else 0.0

    return pct, {m: len(v) for m, v in dist.items()}


def score_ceilings(triples):
    """model -> highest composite score it achieved anywhere in the file. A threshold above
    a model's ceiling is unreachable for it: the cutoff silently becomes other-model-only."""
    out = {}
    for t in triples:
        m = (t.get("predicate") or {}).get("model")
        s = t.get("score")
        if m and isinstance(s, (int, float)):
            out[m] = max(out.get(m, 0.0), float(s))
    return out


def ceiling_report(triples, thresholds):
    """Lines naming any threshold that one model's scores can never reach."""
    ceil = score_ceilings(triples)
    if not ceil:
        return []
    lines = ["  score ceilings (max composite reached per model): "
             + ", ".join(f"{m.split('/')[-1]} {c:.4f}" for m, c in sorted(ceil.items()))]
    for T in sorted({float(t) for t in thresholds}):
        blocked = sorted(m for m, c in ceil.items() if c < T)
        if blocked and len(blocked) < len(ceil):
            reachable = sorted(m for m in ceil if m not in blocked)
            lines.append(f"  !! score>={T:g} is unreachable for "
                         f"{', '.join(m.split('/')[-1] for m in blocked)} -- that cutoff is "
                         f"effectively {'/'.join(m.split('/')[-1] for m in reachable)}-only")
        elif blocked:
            lines.append(f"  !! score>={T:g} is above EVERY model's ceiling -- it selects nothing")
    return lines


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


def merge_models(triples, policy="gate", typed_name=None, gate_name=None, scale="rank"):
    """One triple per entity pair under `policy`. Returns (triples, summary_string).

    `scale` decides whose score represents a pair both models claimed -- see MERGE_SCALES."""
    models, typed, gate = classify_models(triples, typed_name, gate_name)
    pct, _ = model_percentiles(triples)
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
    won = collections.Counter()          # who supplied the score on corroborated pairs
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
        # which model's number represents this pair. Ranks are compared within each model's
        # own distribution (score breaks ties), so a saturating calibrator cannot win by
        # construction; the kept value is still the winner's calibrated probability.
        if scale == "rank":
            win = max(scores, key=lambda m: (pct(m, scores[m]), scores[m]))
            t["score_pct_by_model"] = {m: round(pct(m, s), 4) for m, s in sorted(scores.items())}
        else:
            win = max(scores, key=lambda m: scores[m])
        t["score"] = scores[win]
        t["score_by_model"] = {m: round(s, 4) for m, s in sorted(scores.items())}
        t["score_from"] = win
        t["models"] = sorted(by)
        t["corroborated"] = len(by) > 1
        if t["corroborated"]:
            won[win] += 1
        n_corr += bool(t["corroborated"])
        n_typed_lab += bool(ty is not None)          # label came from the typed model
        out.append(t)
    dropped = f", {n_drop:,} dropped by the gate" if policy == "gate" else \
              (f", {n_drop:,} without a typed verdict dropped" if policy == "typed" else "")
    scored = (f"\n  score[{scale}]: on the {n_corr:,} corroborated pairs the kept score came from "
              + ", ".join(f"{m.split('/')[-1]} {c:,}" for m, c in won.most_common())) if n_corr else ""
    return out, (f"merge[{policy}]: {len(triples):,} triples ({', '.join(models)}) -> {len(out):,} pairs; "
                 f"gate={gate} typed={typed}; {n_corr:,} corroborated by both, "
                 f"{n_typed_lab:,} took the typed label{dropped}{scored}")


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
        # ChEBI first, then NCIt: ChEBI is a small-molecule ontology, so every antibody in the
        # corpus would otherwise be dropped for lacking an id it can never have. chemical.py
        # only writes ncit_label where the ChEBI cascade found nothing, so the two cannot fight.
        s = single(e.get("chebi_label")) or single(e.get("ncit_label"))
        if not s or flags["non_chemical"].get(txt) == "yes":
            return None
        if txt in CHEMICAL_IGNORE or s.strip().lower() in CHEMICAL_IGNORE:
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
                clbl = single(e.get("chebi_label")) or single(e.get("ncit_label"))
                if not clbl:
                    drops["chemical: no ChEBI or NCIt id"] += 1
                    unnorm[f"CHEMICAL {txt}"] += 1
                elif key in CHEMICAL_IGNORE or clbl.strip().lower() in CHEMICAL_IGNORE:
                    drops["chemical: not a compound under study (CHEMICAL_IGNORE)"] += 1
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


# ----- the graph -------------------------------------------------------------
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


def graph_payload_multi(triples, flags):
    """Gene + disease + chemical graph. One edge per unordered node pair, best-supported
    (direction, relation) wins, one record per sentence; nodes are typed and identified by
    their normalized id.

    Edges are categorized by the RELATION the model predicted (activates / inhibits /
    binds / interacts / associated, "not X" when negated), not merely by polarity: with
    the BioRED checkpoint in the routing that label is signed, and dropping it here would
    throw away the entire reason for training it."""
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


def corpus_contrib(triples, flags):
    """This run's contribution to the shared lung corpus: which publications and partners stand
    behind each entity, as SETS rather than counts.

    Sets, because the sibling corpora overlap -- only by one paper today, but summing counts
    would make the denominator wrong the moment two queries pull the same article, and a
    denominator quietly off by a few is worse than one that is obviously off."""
    flags = flags or {"phenotype": {}, "non_chemical": {}}
    pm, ents = set(), set()
    doc = {"gene": collections.defaultdict(set), "chemical": collections.defaultdict(set)}
    par = {"gene": collections.defaultdict(set), "chemical": collections.defaultdict(set)}
    for t in triples:
        p = (t.get("pmid") or "?").split(".")[0]
        pm.add(p)
        a, b = node_of(t["subject"], flags), node_of(t["object"], flags)
        for n in (a, b):
            if n:
                ents.add(n[1])
                if n[0] in doc:
                    doc[n[0]][n[1]].add(p)
        if a and b and a[1] != b[1]:
            if a[0] in par:
                par[a[0]][a[1]].add(b[1])
            if b[0] in par:
                par[b[0]][b[1]].add(a[1])
    return {"pm": sorted(pm), "ents": sorted(ents),
            "doc": {k: {i: sorted(v) for i, v in d.items()} for k, d in doc.items()},
            "par": {k: {i: sorted(v) for i, v in d.items()} for k, d in par.items()}}


def shared_background(data_root, contrib, src_file):
    """Merge this run's contribution with every sibling lung_* corpus into one background.

    Every lung analysis then divides by the same denominator, so an odds ratio from one
    directory means the same thing as one from another. Each directory caches its own
    contribution next to its databases (keyed to the source file's size and mtime, so an
    updated pipeline invalidates it); a sibling that has never been run is simply absent, and
    the page says which corpora it actually covers rather than pretending to totality."""
    root = Path(data_root).resolve()
    cache = root / "databases" / CONTRIB_NAME
    cache.parent.mkdir(parents=True, exist_ok=True)
    st = src_file.stat()
    payload = dict(contrib, _src={"size": st.st_size, "mtime": int(st.st_mtime)})
    cache.write_text(json.dumps(payload), encoding="utf-8")

    family = root.parent.parent            # <...>/relationship_graphs/<project>/<data_root>
    here = root.parent.name
    parts, missing = {here: contrib}, []
    for sib in sorted(family.glob("lung_*")):
        if not sib.is_dir() or sib.name == here:
            continue
        f = sib / root.name / "databases" / CONTRIB_NAME
        if not f.exists():
            missing.append(sib.name)
            continue
        try:
            parts[sib.name] = json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            missing.append(sib.name)

    pm, ents = set(), set()
    doc = {"gene": collections.defaultdict(set), "chemical": collections.defaultdict(set)}
    par = {"gene": collections.defaultdict(set), "chemical": collections.defaultdict(set)}
    for c in parts.values():
        pm.update(c.get("pm", ()))
        ents.update(c.get("ents", ()))
        for kind in doc:
            for i, v in (c.get("doc", {}).get(kind) or {}).items():
                doc[kind][i].update(v)
            for i, v in (c.get("par", {}).get(kind) or {}).items():
                par[kind][i].update(v)
    return {"n": len(pm), "ne": len(ents),
            "doc": {k: {i: len(v) for i, v in d.items()} for k, d in doc.items()},
            "par": {k: {i: len(v) for i, v in d.items()} for k, d in par.items()},
            "src": sorted(parts), "missing": missing}


def background_index(triples, flags):
    """Per-entity denominators over ALL normalized triples, for the in-browser Fisher test.

    Two backgrounds, because the panel ranks on two different units and a test must use the
    unit it is ranking on:
      doc  publications mentioning the entity, and `n`, the corpus publication count
      par  distinct partners the entity is related to, and `ne`, the corpus entity count

    Deliberately score-unfiltered: the page tests a filtered view against the corpus, so the
    corpus must not already be filtered by the cutoff being tested. Keyed through node_of, the
    same identity the graph nodes use -- including its ignore lists -- so a count joins to a
    node without any string matching.

    The population is "documents with an extracted candidate pair", not "documents mentioning
    the entity": a gene named only in a methods section never enters either column of the
    table. That is what a p-value from this index is about."""
    flags = flags or {"phenotype": {}, "non_chemical": {}}   # tolerate a caller with no flags
    pm, ents = set(), set()
    doc = {"gene": collections.defaultdict(set), "chemical": collections.defaultdict(set)}
    par = {"gene": collections.defaultdict(set), "chemical": collections.defaultdict(set)}
    for t in triples:
        p = (t.get("pmid") or "?").split(".")[0]
        pm.add(p)
        a, b = node_of(t["subject"], flags), node_of(t["object"], flags)
        for n in (a, b):
            if n:
                ents.add(n[1])
                if n[0] in doc:
                    doc[n[0]][n[1]].add(p)
        if a and b and a[1] != b[1]:
            if a[0] in par:
                par[a[0]][a[1]].add(b[1])
            if b[0] in par:
                par[b[0]][b[1]].add(a[1])
    return {"n": len(pm), "ne": len(ents),
            "doc": {k: {i: len(v) for i, v in d.items()} for k, d in doc.items()},
            "par": {k: {i: len(v) for i, v in d.items()} for k, d in par.items()}}


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
 #panel,#lpanel{position:absolute;top:12px;z-index:5;background:rgba(255,255,255,.97);border:1px solid #cdd5e0;border-radius:10px;padding:14px 16px;max-width:320px;max-height:calc(100vh - 24px);overflow-y:auto;overscroll-behavior:contain;box-shadow:0 2px 12px rgba(0,0,0,.18);color:#1c2330}
 #panel{right:12px}
 /* left column: what the picture LOOKS like (names, spacing, grouping) -- the right one keeps
    what is IN it (thresholds, filters, provenance), so the two columns split by question */
 #lpanel{left:12px}
 #panel h1{font-size:13px;margin:0 0 8px;color:#1c2330;font-variant:small-caps;letter-spacing:.4px}
 .row{margin:8px 0}
 input[type=range]{width:150px;max-width:100%;vertical-align:middle}
 .legend{display:flex;flex-wrap:wrap;align-items:center;gap:4px 12px}
 /* round: the node legend stands for nodes, which are drawn as circles */
 .legend b{display:inline-block;width:11px;height:11px;border-radius:50%;margin-right:5px;vertical-align:-1px;border:1px solid #999}
 .sw{display:inline-block;width:10px;height:10px;border-radius:2px;vertical-align:-1px}
 .mut{color:#5b6677;font-size:12px} b{color:#2b6cb0}
 /* the left column's prose is longer than its controls, so it hides behind an "i" per section */
 .ihelp{display:inline-flex;align-items:center;justify-content:center;width:16px;height:16px;margin-left:6px;padding:0;border:1px solid #cdd5e0;border-radius:50%;background:#eef2f7;color:#2b6cb0;font:600 11px/1 Segoe UI,Arial,sans-serif;cursor:pointer;vertical-align:1px}
 .ihelp:hover{background:#dde4ee}
 .ihelp.on{background:#0969da;border-color:#0969da;color:#fff}
 .help{display:none;margin-top:4px} .help.open{display:block}
 #conf{width:190px;cursor:pointer}
 select,#search,#genefilter,#drugsearch,#textfilter{max-width:100%;background:#fff;border:1px solid #cdd5e0;color:#1c2330;border-radius:5px;padding:3px 6px;font-size:13px}
 #search,#genefilter,#drugsearch,#textfilter{width:200px}
 mark{background:#ffe680;color:inherit;border-radius:2px;padding:0 1px}
 #catfilters label{display:block;cursor:pointer;white-space:nowrap;font-size:12px;margin:1px 0}
 /* no inner scroller: the relation list is a dozen rows at most, and a box that scrolls inside
    a panel that also scrolls hides ticked types from anyone who does not think to scroll it */
 #catfilters{border:1px solid #cdd5e0;border-radius:6px;padding:4px 6px}
 #catfilters .cnt{color:#5b6677;font-size:11px}
 #siglist{margin-top:5px}
 #siglist .sig{display:flex;align-items:center;gap:6px;padding:1px 2px;font-size:12px;cursor:pointer;border-radius:3px}
 #siglist .sig:hover{background:#eef2f7}
 #siglist .nm{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
 #siglist .vl{width:52px;text-align:right;color:#5b6677;font-variant-numeric:tabular-nums;font-size:11px}
 #siglist .vl b{color:#1c2330}
 #siglist .pc{width:40px;text-align:right;color:#2b6cb0;font-size:11px;font-variant-numeric:tabular-nums}
 #sighdr{font-size:11px;color:#5b6677;text-align:right;margin-top:5px}
 #sigtab,#sigboot{background:#eef2f7;color:#1c2330;border:1px solid #cdd5e0;border-radius:6px;padding:2px 8px;cursor:pointer;font-size:12px}
 #sigtab:hover,#sigboot:hover{background:#dde4ee}
 #sigboot{margin-left:14px}
 #sigboot:disabled{opacity:.6;cursor:default}
 /* The table replaces the CANVAS, not the controls: it covers the graph and the left column, but
    stops short of the right panel (layoutTable measures it), which stays live above the overlay.
    Filtering while reading the table is the point -- every change re-runs the ranking under it. */
 #sigtable{display:none;position:absolute;top:0;left:0;right:0;bottom:0;z-index:20;background:#fff;overflow:auto;padding:16px 20px}
 #sigtable.open+#panel,#panel.overtable{z-index:25}
 #sigtable.open{display:block}
 /* inline-block so the "i" sits beside the title rather than under it */
 #sigtable h2{font-size:15px;margin:0;color:#1c2330;display:inline-block}
 #sigtable table{border-collapse:collapse;font-size:13px;margin-top:10px;min-width:620px}
 #sigtable th,#sigtable td{padding:3px 8px;border-bottom:1px solid #e9edf3;text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
 #sigtable th{position:sticky;top:0;background:#fff;color:#2b6cb0;cursor:pointer;white-space:nowrap;border-bottom:1px solid #cdd5e0}
 #sigtable th.on{font-weight:700;text-decoration:underline}
 #sigtable td.nm,#sigtable th.nm{text-align:left;position:sticky;left:0;background:#fff}
 #sigtable tbody tr:hover td.nm{background:#eef2f7}
 #sigtbody{overflow-x:auto}   /* the table scrolls sideways inside the page, name column pinned */
 #sigtable tbody tr:hover{background:#eef2f7}
 #sigyr{margin:10px 0 2px;font-size:13px}
 /* an undefined test is not a missing feature: say so where the columns would have been */
 #signotest{margin:8px 0 2px;padding:6px 10px;border:1px solid #e5cf9a;border-left:4px solid #d59a2e;
            border-radius:6px;background:#fdf7e8;font-size:12px;color:#5b4a1f;max-width:760px}
 /* the bars live here now: a table row has the width for them, a 320px panel row does not */
 #sigtable td.bar,#sigtable th.bar{width:130px;padding-right:0}
 #sigtable td.bar span{display:block;height:9px;background:#2b6cb0;border-radius:3px;min-width:1px}
 #sigclose{margin-left:14px;background:#eef2f7;border:1px solid #cdd5e0;border-radius:6px;padding:3px 12px;cursor:pointer;font-size:13px}
 #sigclose:hover{background:#dde4ee}
 #sighelp2{max-width:640px}
 #zoom button,#srcbtns button,#labelbtns button,#orbtns button,#orbtns2 button{background:#eef2f7;color:#1c2330;border:1px solid #cdd5e0;border-radius:6px;padding:4px 10px;cursor:pointer;margin-right:6px;font-size:13px}
 #zoom button:hover,#srcbtns button:hover,#labelbtns button:hover,#orbtns button:hover,#orbtns2 button:hover{background:#dde4ee}
 #srcbtns button,#labelbtns button,#orbtns button{margin-bottom:4px}
 #orbtns{margin-top:5px}
 #disfilters{max-height:148px;overflow-y:auto;margin-top:4px}
 #disfilters label{display:block}
 #disclear:hover{color:#0969da}
 #orbtns2 button{margin-left:6px;margin-right:0}
 /* the zoom row sets the left panel's width (see fitLeftPanel), so the last button must not
    carry a trailing margin -- 6px of it would push Fit onto a second line */
 #zoom button:last-child{margin-right:0}
 #comention{width:100%}
 #srcbtns button.on,#labelbtns button.on,#orbtns button.on,#orbtns2 button.on{background:#0969da;border-color:#0969da;color:#fff;font-weight:600}
 #srcbtns button:disabled{opacity:.45;cursor:default}
 .vis-tooltip{max-width:480px!important;white-space:normal!important;background:#fff!important;color:#1a1a1a!important;border:1px solid #999!important;border-radius:8px!important;padding:8px 10px!important;box-shadow:0 4px 16px rgba(0,0,0,.35)!important;font:12px/1.45 Segoe UI,Arial,sans-serif!important}
 .eth{font-size:13px;margin-bottom:6px} .stip{padding:3px 0;border-top:1px solid #e3e3e3}
 .pm{display:inline-block;background:#eef3fb;color:#2b6cb0;border-radius:4px;padding:0 5px;margin-right:5px;font-weight:600;font-size:11px;text-decoration:none}
 a.pm:hover{background:#d6e6fb;text-decoration:underline} .more{margin-top:5px;color:#888;font-style:italic}
 #info{max-height:240px;overflow:auto} #info .stip{border-top:1px solid #e3e3e3}
 /* the details box is 12px prose; its heading takes the panel's normal size so the block reads
    as a section rather than as more small print */
 #info .ihead{font-size:14px;color:#1c2330;margin-bottom:4px}
 /* one screen, two columns is a desktop luxury: on a phone they stack, one at each edge */
 @media (max-width:700px){
  #panel,#lpanel{left:12px;right:12px;max-width:none;max-height:42vh;overflow:auto}
  #lpanel{top:12px}
  #panel{top:auto;bottom:12px}
  .vis-tooltip{max-width:88vw!important}
 }
</style></head><body>
<div id="lpanel">
 __PUBMED_QUERY__
 <div class="row">Year: <b id="yrlab"></b><br><input id="yrlo" type="range" style="width:74px"> <input id="yrhi" type="range" style="width:74px"></div>
 <div class="row" id="lblrow">Draw names for: <button class="ihelp" aria-label="About drawn names" aria-expanded="false">i</button><br><span id="labelbtns"><button class="lblb" data-kind="disease">Disease names</button><button class="lblb" data-kind="chemical">Drug names</button></span>
  <div class="mut help">Gene symbols are always drawn. Disease and drug names are long, repeat across many edges and bury the symbols by sheer wordiness, so they start off &mdash; switch them on to read a neighbourhood, off to see its shape. When on they are set in a <b>condensed</b> face at <b>80%</b> of the size the same weight of gene would take, and drug names are <i>italic</i> as well: the extra cues let a long name sit beside a symbol without out-shouting it, and let you tell the three kinds apart in one glance without reading them. Relative size still tracks evidence within each kind, and the font slider scales everything together. Applies in place; the layout is not recomputed.</div></div>
 <div class="row">Shrink periphery: <b id="shv">0%</b> <button class="ihelp" aria-label="About shrink periphery" aria-expanded="false">i</button><br><input id="shrink" type="range" min="0" max="90" step="5" value="0" aria-label="Shrink periphery">
  <div class="mut help">Pulls everything past the middle distance inward, so the small clusters physics flung to the edges come back where you can read them without zooming out. The centre is left alone.</div></div>
 <div class="row">Expand center: <b id="exv">0%</b> <button class="ihelp" aria-label="About expand center" aria-expanded="false">i</button><br><input id="expand" type="range" min="0" max="100" step="5" value="0" aria-label="Expand center">
  <div class="mut help">Blows the crowded core outward and carries the rest along, thinning the hairball without moving anything past its neighbours. Both reshape the drawn layout only &mdash; radially, from the stabilized positions, so nothing changes order and returning a slider to 0% restores the layout exactly.</div></div>
 <div class="row" id="tissuerow"><label><input type=checkbox id="tissuestack" checked> Stack same-tissue diseases</label> <button class="ihelp" aria-label="About tissue stacking" aria-expanded="false">i</button>
  <div class="mut help">Drops the disease nodes naming one tissue onto a single spot, overlapping, so <em>lung cancer</em>, <em>lung adenocarcinoma</em> and <em>non-small cell lung carcinoma</em> read as one place on the canvas instead of three. They stay separate nodes with their own edges and tooltips &mdash; only their positions are pooled, after the layout settles. Tissue is read from the name (<span id="tissuen"></span>).</div></div>
 <div class="row" id="cmrow">Co-mention links: <button class="ihelp" aria-label="About co-mention links" aria-expanded="false">i</button><br><select id="comention"><option value="">(off)</option></select>
  <div class="mut help">Draws a dashed grey link from every node whose <em>visible</em> sentences name that disease &mdash; its full name or its acronym &mdash; even where no model predicted a relation. Nodes already wired to it by a drawn relation keep that edge and get no second one, so a dashed link reads &ldquo;co-mentioned, nothing predicted&rdquo;. Co-occurrence only, never a claim; added after all filtering, so it changes nothing the thresholds keep.</div></div>
 <div class="row" id="disrow">Keep diseases: <button class="ihelp" aria-label="About the disease keep-list" aria-expanded="false">i</button>
  <span id="disclear" class="mut" style="cursor:pointer;text-decoration:underline;float:right">clear</span><br>
  <div class="legend" id="disfilters"></div>
  <div class="mut help">Ticking a disease <b>keeps it</b>: the canvas then shows those diseases only, and drops the rest along with the edges that led to them. Genes and drugs are untouched &mdash; this filters one node type, not the picture. <b>Tick nothing and nothing is filtered</b>, which is the default; untick the last one (or press <em>clear</em>) to go back. The list is rebuilt on every redraw from the diseases actually in front of you, with the number of drawn edges each one carries &mdash; so it follows the score, year, text and every other control. It is built from the view <em>before</em> this filter is applied, the same rule the significance cut follows: otherwise ticking one disease would empty the menu you are ticking from. A disease you have ticked stays listed even when the current filters leave it nothing, shown <span style="color:#b3243b">(0)</span>, so you can always untick it.</div></div>
 <div class="row" id="zoom"><button id="zin">+ Zoom in</button><button id="zout">&minus; Zoom out</button><button id="zfit">Fit</button></div>
 <div class="row" id="qrow">Significance: <b id="qsv">off (show all)</b> <button class="ihelp" aria-label="About the significance cutoff" aria-expanded="false">i</button><br><input id="qsig" type="range" min="0" max="9" step="1" value="0" aria-label="Significance cutoff">
  <div id="orbtns"><button class="orb" data-or="gt1">OR&gt;1</button><button class="orb" data-or="lt1">OR&lt;1</button></div>
  <div class="mut help" id="qhelp">Hides genes and drugs whose enrichment was <em>tested and missed</em> the cutoff. Entities the test could not reach &mdash; fewer than five corpus papers, or nothing outside the view to contrast against &mdash; are <b>kept</b>: they were never judged, so they cannot have failed, and on this corpus they are over half the ranked genes. The buttons ask a different question &mdash; which side of 1 &mdash; so they <em>do</em> drop the untested, which sit on neither: over-represented in this view (OR&gt;1) or under-represented (OR&lt;1). Each toggles on its own and press a lit one again to clear it. <b>Both lit is not the same as neither</b>: it keeps every entity the test placed on a side and drops the ones it could not place &mdash; the &ldquo;only what was actually measured&rdquo; view. The depleted side is a finding too: in an adenocarcinoma view the small-cell markers DLL3 and ASCL1 land there. Diseases stay throughout: they carry no ranking of their own. <b>A high q is not a small node.</b> The cut asks whether an entity is over-represented <em>against the pooled six-corpus lung background</em>, not whether it matters here: EGFR carries this view at rank 2 by publications, yet sits at OR&nbsp;1.05, q&nbsp;0.88, because it is just as common in the squamous and small-cell papers. Read the ranking for what the view is made of, and this for what is distinctive about it.</div></div>
 <div class="row" id="sigrow">Significance in view <button class="ihelp" aria-label="About significance in view" aria-expanded="false">i</button><br>
  <select id="sigkind"><option value="gene">genes</option><option value="chemical">drugs</option></select>
  <select id="sigmeasure"><option value="pub">by publications</option><option value="deg">by partners</option><option value="sent">by sentences</option></select>
  <div class="mut" id="signote" style="display:none"></div>
  <div id="sighdr">publications &nbsp;&nbsp;z</div>
  <div id="siglist"></div>
  <button id="sigtab">Full table view</button>
  <div class="mut help" id="sighelp">The top six of whatever the graph is <em>currently drawing</em>: every control reshapes this too &mdash; score, year, text, node and relation type, training set, the support thresholds and the structural pruning &mdash; and it is recomputed on every redraw, so the ranking and the picture can never disagree. Widen the filters to read it as the corpus; narrow them to ask the same question of a slice. Co-mention links are excluded, since this counts relations. <b>Publications</b> counts the distinct papers behind a node's relations &mdash; the measure the edge thicknesses use; <b>partners</b> counts the distinct entities it is related to (breadth, not weight); <b>sentences</b> counts the unique sentences supporting them. The second column is <b>z</b>: standard deviations above the mean on a log&#8321;&#8320; scale, <em>among its own kind</em>, since a gene is only remarkable among genes. The percentile is there too (hover a row), but it saturates &mdash; every one of a top six reads 99.9%, while z still separates them. In the panel list, click a row to select and centre that node; if the current filters have removed it, the details box says so rather than moving the view. Table rows do not navigate &mdash; they are there to be read and sorted. <b>Bootstrap CIs</b>, in the table, resamples the view's <em>publications</em> with replacement 300 times and reports 95% intervals for each count and each rank &mdash; papers, because that is the independent unit; resampling sentences would give intervals several times too tight. Ranks at the top are firm (1&ndash;2) and the tail is not (a gene ranked 500th may belong anywhere from 264th to 1447th), which is the honest width of &ldquo;top ten&rdquo;. <b>OR vs corpus</b> and <b>q</b> ask a different question: is this entity over-represented in the current view compared with the whole normalized corpus? A 2&times;2 over publications &mdash; in view or not, mentions it or not &mdash; by Fisher's exact test, with a Haldane-corrected odds ratio, a Woolf interval and Benjamini-Hochberg q-values over the entities with at least 5 corpus papers. Type &ldquo;immunotherapy&rdquo; and PDCD1, CD274 and CTLA4 come out at OR 5&ndash;9; EGFR comes out <em>depleted</em>. It cannot say an entity is specific to lung adenocarcinoma &mdash; every paper here is a lung paper, so the contrast is view-against-corpus, never corpus-against-literature. <b>Full table view</b> opens the whole ranking &mdash; every node of that kind, all three counts, both statistics and a bar, sortable by any column, with its own year handles.</div></div>
 <div class="row mut" id="info">Click a node or edge for details.</div>
</div>
<div id="panel">
 <div class="row legend"><b style="background:#cfe3ff;border-color:#2b6cb0"></b>gene <b style="background:#1b7837;border-color:#145a28"></b>approved anti-neoplastic <b style="background:#e08600;border-color:#9a6700"></b>approved (other) <b style="background:#c2185b;border-color:#7a0f3a"></b>ChEBI <button class="ihelp" data-help="col" aria-label="About node colours" aria-expanded="false">i</button></div>
 <div class="row mut help" data-help="col">Drug-target genes (corpus chemicals): deeper colour = more chemicals. <b style="color:#1b7837">Green</b> = DGIdb approved anti-neoplastic, <b style="color:#e08600">amber</b> = DGIdb approved (non-anti-neoplastic), <b style="color:#c2185b">pink</b> = ChEBI, absent in DGIdb.</div>
 <div class="row">Min unique sentences/edge: <b id="thv">1</b><br><input id="thr" type="range" min="1" max="10" value="1"></div>
 <div class="row">Min unique publications: <b id="mpv">1</b> <button class="ihelp" aria-label="About min unique publications" aria-expanded="false">i</button><br><input id="minpub" type="range" min="1" max="10" value="1">
  <div class="mut help">Distinct PMIDs behind an edge; raise it to drop relations that rest on one paper repeating itself.</div></div>
 <div class="row">Font size: <b id="fsv">50%</b><br><input id="fscale" type="range" min="10" max="100" step="5" value="50" aria-label="Label font size"></div>
 <div class="row">Min cluster size: <select id="mincluster"><option selected>2</option><option>3</option><option>4</option><option>5</option><option>6</option><option>7</option><option>8</option><option>9</option><option>10</option><option>11</option><option>12</option></select></div>
 <div class="row">Min connections: <select id="mindeg"><option selected>2</option><option>3</option><option>4</option><option>5</option><option>6</option></select> <button class="ihelp" aria-label="About min connections" aria-expanded="false">i</button><div class="mut help">Hides genes linked to fewer than this many others; thins the hairball's single-link fringe. A single pass: nodes that lose links in it can finish below the bar.</div></div>
 <div class="row">Search gene: <input id="search" placeholder="e.g. EGFR" autocomplete="off"></div>
 <div class="row">Filter to gene:<br><input id="genefilter" placeholder="e.g. EGFR (+neighbors)" autocomplete="off"> <select id="hops"><option value="1">1 hop</option><option value="2">2 hops</option></select></div>
 <div class="row">Search drug: <input id="drugsearch" placeholder="e.g. nivolumab" autocomplete="off"></div>
 <div class="row">Filter to drug:<br><select id="chemfilter"><option value="">(all drugs)</option></select></div>
 <div class="row">Match text in sentence: <button class="ihelp" aria-label="About the text filter" aria-expanded="false">i</button><br><input id="textfilter" placeholder="e.g. phosphorylat or /inhibit(s|ed)?/" autocomplete="off">
  <div class="mut help">Case-insensitive substring; wrap in / / for a regex. Keeps only edges with a matching sentence, and shows just those sentences. The thresholds above weigh an edge's <em>full</em> support, so a match is never dropped for evidence the query happened to hide &mdash; min-publications judges all of an edge's papers, not just the matching ones. <b>Min connections</b> and <b>Min cluster size</b> are the exception: they describe the picture, so they are re-applied to what the query leaves.</div></div>
__KINDROW__
 <div class="row">Relation type <button class="ihelp" data-help="rel" aria-label="About relation types" aria-expanded="false">i</button>
  <div class="mut help" data-help="rel">As predicted by the RE model; &ldquo;not X&rdquo; = negated statement, drawn dashed. Unticking one hides <em>sentences</em> with that label, and any edge left without support.</div></div><div id="catfilters"></div>
 <div class="row mut help" data-help="rel">Edge colour = the relation the model predicted. <b>activates</b>/<b>inhibits</b> are signed and come from the BioRED checkpoint; <b>interacts</b> is the unsigned PPI verdict. An edge takes its best-supported direction, and is drawn as the relation most of its sentences <em>in view</em> carry &mdash; so narrowing the filters can recolour an edge. Hover for the per-sentence labels. Thickness and arrowhead size follow the number of <b>independent publications</b> behind the edge, not its sentence count &mdash; one paper repeating itself never thickens a line.</div>
 <div class="row mut help" data-help="rel">Counts read <em>total &middot; in view</em>: the total is every edge in the file carrying at least one sentence of that type (an edge with mixed readings counts under each, so the totals exceed the edge count), &ldquo;in view&rdquo; is how many survive the current score, year, text, min-publications, min-connections and min-cluster settings. <span style="color:#b3243b">A red 0</span> means the type is ticked but everything of it is pruned &mdash; usually its edges sit in components smaller than <b>Min cluster size</b>, so lower that (or the score) to see them.</div>
 <div class="row">Training set behind the edge <button class="ihelp" data-help="src" aria-label="About training sets" aria-expanded="false">i</button></div>
 <div class="row" id="srcbtns"></div>
 <div class="row mut help" data-help="src" id="srchint">Which corpus the relation was learned from &mdash; <b>PPI-only</b> = found by the BioInfer/PPI model alone, <b>BioRED-only</b> = by the BioRED model alone (typed and often signed), <b>both</b> = the two agreed a relation is there. These cut at the <em>sentence</em>, like the relation types: an edge whose support is split between the models appears under each button with that model's sentences only, so the counts sum past the edge total. Since <b>interacts</b> is the binary model's only positive label, <b>PPI-only</b> and &ldquo;All with just interacts ticked&rdquo; are the same view.</div>
 <div class="row">Relationship score: <b id="scval">&ge;0.99</b><br><input id="conf" type="range" min="0" max="13" step="1" value="__CONFDEF__" aria-label="Minimum relationship score"></div>
 <div class="row mut" id="stats"></div>
</div>
<div id="net"></div>
<div id="sigtable"><h2 id="sigttl"></h2><button class="ihelp" data-help="sig" aria-label="About significance in view" aria-expanded="false">i</button>
 <div class="mut help" data-help="sig" id="sighelp2"></div>
 <div class="mut" id="sigsub"></div>
 <div class="row" id="qrow2">Significance: <b id="qsv2">off (show all)</b>
  <input id="qsig2" type="range" min="0" max="9" step="1" value="0" style="width:160px;vertical-align:middle" aria-label="Significance cutoff">
  <span id="orbtns2"><button class="orb2" data-or="gt1">OR&gt;1</button><button class="orb2" data-or="lt1">OR&lt;1</button></span>
  <button class="ihelp" aria-label="About the significance cutoff" aria-expanded="false">i</button>
  <div class="mut help" id="qhelp2"></div></div>
 <div class="mut" id="signote2" style="display:none"></div>
 <div id="signotest" style="display:none"></div>
 <div id="sigyr" class="row"><select id="sigkind2"><option value="gene">genes</option><option value="chemical">drugs</option></select>
  <select id="sigmeasure2"><option value="pub">by publications</option><option value="deg">by partners</option><option value="sent">by sentences</option></select>
  &nbsp; Year: <b id="yrlab2"></b> <input id="yrlo2" type="range" style="width:120px"> <input id="yrhi2" type="range" style="width:120px">
  <button id="sigboot">Bootstrap CIs</button>
  <button id="sigclose">Graph view</button></div>
 <div id="sigtbody"></div></div>
<script>
const DATA=__PAYLOAD__;
const BG=__BACKGROUND__;   // per-entity document counts over the whole normalized corpus
const CCOLOR=__CCOLOR__;
const MINY=__MINY__, MAXY=__MAXY__;
const MAXTGT=Math.max(1,...DATA.nodes.map(n=>n.target||0));
// Node colour: genes keep the drug-target shading (deeper = more corpus chemicals); disease
// and chemical nodes take their type colour, and every node its type SHAPE, so the three
// kinds stay distinguishable without relying on colour alone.
function nodeColor(n){if(n.kind&&n.kind!=='gene')return {background:n.bg,border:n.border};if(!n.target)return {background:n.bg||'#cfe3ff',border:n.border||'#2b6cb0'};const t=n.target/MAXTGT,L=(a,b)=>Math.round(a+(b-a)*t);if(n.tcat==='green')return {background:'rgb('+L(200,27)+','+L(230,120)+','+L(201,55)+')',border:'#145a28'};if(n.tcat==='amber')return {background:'rgb('+L(255,224)+','+L(231,134)+','+L(179,0)+')',border:'#9a6700'};return {background:'rgb('+L(255,194)+','+L(217,24)+','+L(232,91)+')',border:'#7a0f3a'};}
const KIND={};DATA.nodes.forEach(n=>{KIND[n.id]=n.kind||'gene';});
// --- training-set provenance ---------------------------------------------------------
// Every SENTENCE records which corpus produced it: 'ppi' (BioInfer only), 'biored' (BioRED
// only, i.e. a reading the binary model never claimed) or 'both' (the merge corroborated it).
// The buttons isolate each, which is how you SEE what a training set contributes rather than
// inferring it from counts -- and they cut at the sentence, so an edge the two models split
// between them shows up under each, carrying only that model's half of the evidence.
let SRC_MODE='all';
const SRC_BTN=[['all','All'],['ppi','PPI-only'],['biored','BioRED-only'],['both','Both agreed']];
function buildSrcButtons(){
 // counted like the relation types: edges holding at least one sentence from that source, so a
 // pair the two models split between them counts under each (the parts sum past the edge total)
 const n={all:DATA.edges.length,ppi:0,biored:0,both:0};
 DATA.edges.forEach(e=>{const r={};e.sents.forEach(s=>r[sentSrc(e,s)]=1);for(const k in r)n[k]=(n[k]||0)+1;});
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
const net=document.getElementById('net'); let network=null, NODEDS=null;
// the layout exactly as physics left it, plus the tissue groups of what is on screen; every
// position control replays from these, so they compose and none of them accumulates drift
let BASEPOS=null, TGROUPS=[];
function applyLayoutShape(){
 if(!network||!BASEPOS)return;
 reshapeRadial(network,BASEPOS,activeShrink(),activeExpand());
 if(activeTissueStack()&&TGROUPS.length)stackTissues(network,TGROUPS);
 network.redraw();
}
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
// The relation-type boxes filter SENTENCES, not whole edges. An edge's payload "cat" is the
// label that won over the whole file, so once a score/year/text filter narrows the support the
// sentences still in view can all carry a different label -- SOCS1-MALAT1 is 'interacts' over
// its nine sentences, but the three that mention "lung" all say 'activates'. Filtering,
// counting and colouring on the per-sentence label is what keeps the panel honest: "activates
// 0 in view" then means no visible sentence says activates, which is what the reader sees.
function sentCat(e,s){return s.rc||e.cat;}
// Same story for the training set: e.src is the union of model roles over ALL the pair's
// triples, so a pair where PPI said "interacts" in one sentence and BioRED said "activates" in
// another is tagged 'both' -- and an edge-level filter then hid its PPI-only sentences from the
// PPI-only view. Filtering per sentence makes the buttons agree with the relation boxes: since
// 'interacts' is the PPI checkpoint's only positive label, "PPI-only" and "all sources with
// just interacts ticked" now describe the same 14125 edges instead of 10395 and 14125.
function sentSrc(e,s){return s.sr||e.src||'ppi';}
function visSents(e,conf,lo,hi,tm,cats,mode){return e.sents.filter(s=>s.sc>=conf&&passYear(s.yr,lo,hi)&&(!tm||tm.test(s.text))&&(!cats||cats.has(sentCat(e,s)))&&(!mode||mode==='all'||sentSrc(e,s)===mode));}
// The relation an edge is DRAWN as: whichever label carries the most sentences in view. A tie
// keeps the payload's own cat if it is among the leaders, else the last alphabetically --
// mirroring the reverse sort that picked cat in the first place.
function viewCat(e,vis){
 if(!vis.length)return e.cat;
 const n={};vis.forEach(s=>{const c=sentCat(e,s);n[c]=(n[c]||0)+1;});
 let best=null;
 for(const c in n){if(best===null||n[c]>n[best]||(n[c]===n[best]&&(c===e.cat||(best!==e.cat&&c>best))))best=c;}
 return best;
}
const SRCLAB={ppi:'PPI',biored:'BioRED',both:'PPI+BioRED'};
// the provenance of what is ON SCREEN: roles pooled over the visible sentences, so an edge
// reads PPI+BioRED only while sentences from both are actually in view
function viewSrc(e,vis){
 if(!vis.length)return e.src||'ppi';
 let p=false,b=false;
 vis.forEach(s=>{const r=sentSrc(e,s);if(r==='both'){p=true;b=true;}else if(r==='biored')b=true;else p=true;});
 return p&&b?'both':(b?'biored':'ppi');
}
function edgeHead(e,vis,cat){const np=new Set(vis.map(s=>s.pmid)).size;return '<div class=eth><b>'+esc(labelById[e.from])+' &rarr; '+esc(labelById[e.to])+'</b> ('+vis.length+' sentences &middot; '+np+' PMIDs &middot; '+esc(cat||viewCat(e,vis))+' &middot; '+(SRCLAB[viewSrc(e,vis)]||'PPI')+')</div>';}
// per-sentence relation tag: the label the model gave THIS sentence (an edge shows its
// dominant relation, so a minority reading -- e.g. one "inhibits" under an "interacts"
// edge -- would otherwise be invisible). "?" marks a speculated statement.
function relTag(s){return s.rc?' <span class=mut style="color:'+(CCOLOR[s.rc]||'#888')+'">['+esc(s.rc)+(s.sp?' ?':'')+']</span>'+(s.sr?' <span class=mut>'+esc(SRCLAB[s.sr]||s.sr)+'</span>':''):'';}
function edgeTip(e,vis,cat){const d=document.createElement('div');let h=edgeHead(e,vis,cat);const lim=20;vis.slice(0,lim).forEach(s=>{h+='<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span>'+relTag(s)+' '+hl(s.text)+'</div>';});if(vis.length>lim)h+='<div class=more>+'+(vis.length-lim)+' more</div>';d.innerHTML=h;return d;}
// --- co-mention links -----------------------------------------------------------------
// Sentences name a disease far more often than the models emit a relation for it: 1195 edges
// in this corpus mention NSCLC, yet 70 of their endpoint nodes carry no NSCLC edge at all --
// the triple went to "lung cancer" while the sentence said "the lung cancer of NSCLC". These
// links show that co-occurrence for what it is: dashed, grey, undirected, never given a
// relation colour, so nothing here can be misread as something a checkpoint predicted.
// Aliases come from the label itself -- words separated by any run of spaces/hyphens, an
// interchangeable cancer/carcinoma/tumour/neoplasm tail, and the initials when they spell an
// acronym of three or more letters ("non-small cell lung carcinoma" -> NSCLC).
const CM_STOP=new Set(['of','the','and','with','in','a','to']);
const CM_TAIL='(?:carcinomas?|cancers?|tumou?rs?|neoplasms?)';
function diseaseAliases(label){
 const words=(label||'').toLowerCase().split(/[^a-z0-9]+/).filter(Boolean);
 if(!words.length)return null;
 const alts=[words.map((w,i)=>(i===words.length-1&&new RegExp('^'+CM_TAIL+'$').test(w))?CM_TAIL:reEsc(w)).join('[-\\s]+')];
 const ac=words.filter(w=>!CM_STOP.has(w)).map(w=>w[0]).join('').toUpperCase();
 if(ac.length>=3)alts.push('\\b'+ac+'s?\\b');
 return new RegExp(alts.join('|'),'i');
}
function cmTip(link,dis){
 const d=document.createElement('div');const lim=20;
 let h='<div class=eth><b>'+esc(labelById[link.nd]||link.nd)+' &middot;&middot;&middot; '+esc(labelById[dis]||dis)+'</b> ('
  +link.sents.length+' sentences &middot; co-mentioned, no relation predicted)</div>';
 link.sents.slice(0,lim).forEach(s=>{h+='<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span> '+hl(s.text)+'</div>';});
 if(link.sents.length>lim)h+='<div class=more>+'+(link.sents.length-lim)+' more</div>';
 d.innerHTML=h;return d;
}
function activeComention(){return (document.getElementById('comention')||{}).value||'';}
// --- the disease keep-list ------------------------------------------------------------------
// A WHITELIST over disease nodes: empty means no constraint (the default), and any tick means
// "these diseases only". It touches nothing else -- a gene-gene edge has no disease endpoint and
// cannot be cut by it -- so this narrows which diseases the picture is allowed to contain,
// rather than narrowing the picture.
const DIS_KEEP=new Set();
function disPass(id){return !DIS_KEEP.size||KIND[id]!=='disease'||DIS_KEEP.has(id);}
// Rebuilt on every redraw from the edges the rest of the controls left, so the menu is always
// the diseases actually in front of you. Two things it must not do: it must not be built from
// its OWN leavings, or ticking one disease would empty the list you are ticking from (so the
// caller passes the pre-filter edges, the same order the q cut uses); and it must not drop a
// ticked disease that the current filters have starved, or that tick could never be undone --
// those stay listed at (0).
function buildDiseaseFilter(edges){
 const box=document.getElementById('disfilters'); if(!box)return;
 const seen={};
 edges.forEach(o=>[o.e.from,o.e.to].forEach(nd=>{if(KIND[nd]==='disease')seen[nd]=(seen[nd]||0)+1;}));
 const ids=Object.keys(seen);
 DIS_KEEP.forEach(id=>{if(seen[id]===undefined)ids.push(id);});
 const lab=id=>String(labelById[id]||id);
 ids.sort((a,b)=>(seen[b]||0)-(seen[a]||0)||lab(a).localeCompare(lab(b)));
 box.innerHTML=ids.length
   ? ids.map(id=>{const n=seen[id]||0;
       return '<label><input type=checkbox class=disf value="'+esc(id)+'"'+(DIS_KEEP.has(id)?' checked':'')
         +'> '+esc(lab(id))+' <span class=cnt'+(n?'':' style="color:#b3243b"')+'>('+n+')</span></label>';}).join('')
   : '<span class="mut">no diseases in view</span>';
 box.querySelectorAll('.disf').forEach(c=>c.addEventListener('change',()=>{
   if(c.checked)DIS_KEEP.add(c.value);else DIS_KEEP.delete(c.value);
   build(+document.getElementById('thr').value);
 }));
}
// --- same-tissue disease stacking ------------------------------------------------------
// One tissue is spread over many disease nodes -- 22 of them name the lung here, carrying 4196
// connections between them -- because the corpus says "lung cancer" where it means NSCLC and
// vice versa. Merging them would be a claim about ontology this script has no business making,
// so instead their POSITIONS are pooled: same tissue, same spot, overlapping. Each keeps its own
// edges, size and tooltip; only the layout is touched, and only after it has settled.
// Tissue is read off the label. \b in front of "renal" is what keeps adrenal out of the kidney.
const TISSUE_RULES=[['lung',/\b(lung|pulmonary|nsclc|sclc|bronch|pleural)/i],
 ['breast',/\b(breast|mammary)/i],['liver',/\b(liver|hepat)/i],['stomach',/\b(gastric|stomach)/i],
 ['colon',/\b(colorect|colon|rectal|bowel|intestin)/i],['prostate',/\bprostat/i],
 ['brain',/\b(brain|glio|astrocyt|neuroblast|medulloblast|cerebr|meningi)/i],['pancreas',/\bpancrea/i],
 ['ovary',/\bovari/i],['kidney',/\b(renal|kidney|nephro)/i],['thyroid',/\bthyroid/i],
 ['skin',/\b(melanom|skin|cutaneous)/i],['esophagus',/\b(esophag|oesophag)/i],
 ['bladder',/\b(bladder|urothelial)/i],['cervix',/\bcervi/i],['uterus',/\b(uterine|endometri)/i],
 ['blood',/\b(leukemi|lymphom|myelom|myeloid|myelodysplas|marrow)/i],   // marrow before bone: it is haematopoietic
 ['head/neck',/\b(head and neck|nasopharyng|laryn|oral|tongue)/i],['bone',/\b(osteo|bone)/i],
 ['heart',/\b(cardiac|cardio|myocard|heart)/i]];
function tissueOf(label){for(const [t,re] of TISSUE_RULES)if(re.test(label||''))return t;return '';}
function activeTissueStack(){const el=document.getElementById('tissuestack');return el?!!el.checked:true;}
// drawn disease nodes that share a tissue, two or more of them -- nothing else is worth moving
function tissueGroups(nodes){
 const by={};
 nodes.forEach(n=>{if(KIND[n.id]!=='disease')return;const t=tissueOf(labelById[n.id]||n.id);if(t)(by[t]=by[t]||[]).push(n.id);});
 return Object.keys(by).filter(t=>by[t].length>1).map(t=>({tissue:t,ids:by[t]}));
}
// --- radial reshaping ------------------------------------------------------------------
// Two knobs on the DRAWN layout, both radial around the graph's centre of mass, both monotone
// in the radius: no node ever passes another, so the picture stays the one physics produced --
// only the spacing changes. The boundary between "centre" and "periphery" is the median radius,
// which puts half the nodes on each side whatever the graph's shape.
//   expand: r <= R -> r*(1+e)  (the core swells)   r > R -> r + R*e  (the rest rides out rigidly)
//   shrink: past the expanded boundary b, distances compress by (1-s), pulling the far clusters in
// Both read BASEPOS, the untouched stabilized positions, so the sliders never compound: every
// move is computed from the original layout and 0% restores it exactly.
function activeShrink(){const v=parseInt((document.getElementById('shrink')||{}).value);return isNaN(v)?0:v/100;}
function activeExpand(){const v=parseInt((document.getElementById('expand')||{}).value);return isNaN(v)?0:v/100;}
function reshapeRadial(net,pos,s,e){
 const ids=Object.keys(pos);
 if(!ids.length)return;
 let cx=0,cy=0;
 ids.forEach(id=>{cx+=pos[id].x;cy+=pos[id].y;});
 cx/=ids.length;cy/=ids.length;
 const rs=ids.map(id=>Math.hypot(pos[id].x-cx,pos[id].y-cy)).sort((a,b)=>a-b);
 const R=rs[Math.floor(rs.length/2)]||1, b=R*(1+e);
 ids.forEach(id=>{
  const dx=pos[id].x-cx,dy=pos[id].y-cy,r=Math.hypot(dx,dy);
  if(r<1e-6)return;
  let nr=e?(r<=R?r*(1+e):r+R*e):r;
  if(s&&nr>b)nr=b+(nr-b)*(1-s);
  const k=nr/r;
  net.moveNode(id,cx+dx*k,cy+dy*k);
 });
}
// pool each group onto its own centre: a ring tight enough that the discs overlap, wide enough
// that every node stays individually clickable
function stackTissues(net,groups){
 groups.forEach(g=>{
  const pos=net.getPositions(g.ids);let x=0,y=0,n=0;
  g.ids.forEach(id=>{const p=pos[id];if(p){x+=p.x;y+=p.y;n++;}});
  if(!n)return;
  x/=n;y/=n;
  const R=Math.min(30,8+g.ids.length*0.6);
  g.ids.forEach((id,k)=>{const a=2*Math.PI*k/g.ids.length;net.moveNode(id,x+R*Math.cos(a),y+R*Math.sin(a));});
 });
}
// Edge weight = INDEPENDENT PUBLICATIONS, not sentences: one paper restating a finding five
// times must not draw a line five times heavier than five papers agreeing once each. 11,025 of
// this corpus's 31,149 edges carry more sentences than papers, so the two measures differ for a
// third of the graph. Square-rooted because the distribution is brutally skewed -- 89% of edges
// rest on a single publication and the heaviest on 269 -- and capped so one hub cannot draw a
// band across the canvas. The arrowhead grows with the same measure, so direction stays legible
// on the thick edges instead of being swallowed by them.
function edgeWidth(np){return Math.min(1+2.0*Math.sqrt(Math.max(0,np-1)),12);}
function arrowScale(np){return Math.min(0.5+0.18*Math.sqrt(Math.max(0,np-1)),1.6);}
// heads the details box whenever it lists sentences -- above the "A -> B" line, not in the
// hover tooltip, which is transient and already framed by the edge you are pointing at
const INFO_HEAD='<div class="ihead">Sentences of interest</div>';
function scaleNode(s){return 6+Math.sqrt(s)*3.4;}
function fontSize(c){return c<5?13:2*Math.max(13,Math.min(Math.round(c*2.2),48));}
// Label size is a per-node property (a busy gene is drawn larger), so the slider is a single
// multiplier on top of fontSize(): every label shrinks by the same factor and the size ranking
// survives. The base size rides along on each node as _fs, which is what lets the slider relabel
// the existing DataSet instead of rebuilding -- a rebuild would re-run the layout.
let FSCALE=1;
// FS_BASE: the slider reading that means "the size fontSize() computed". Everything above it
// enlarges, so the old ceiling (the built-in size) now sits mid-scale at 50% and the slider
// runs to twice that. The divisor -- not the range -- is what fixes where 1.0x lands.
const FS_BASE=50;
function activeFontScale(){const el=document.getElementById('fscale');const v=el?parseInt(el.value):FS_BASE;return (isNaN(v)?FS_BASE:v)/FS_BASE;}
function scaledFont(b){return Math.max(4,Math.round(b*FSCALE));}
// Which node kinds have their NAME DRAWN on the canvas. Disease and chemical names are long,
// repeat across many edges and out-shout the gene symbols simply by being wordy, so those
// nodes start unlabelled -- shape and colour say what they are, and the name is one hover
// (or click) away. The two buttons put them back: the set is live, so toggling relabels the
// nodes already on screen (a rebuild would re-run the layout for a question about text).
const LABEL_KINDS=new Set(['gene']);
// Gene symbols are short and are the thing you read positionally; disease and chemical names are
// long prose that competes with them for the same canvas. So when those names are switched on
// they are set in a CONDENSED face at 80% of the size a gene of the same weight would take, and
// drugs are additionally italic -- three cues (shape, width, slope) doing what size alone had to.
// The ratio is applied to the BASE size, not the rendered one, so the font slider and the
// per-node "busier is bigger" ranking both survive untouched.
const NARROW_FACE='Arial Narrow,Liberation Sans Narrow,Segoe UI Semilight,Arial,sans-serif';
const KIND_FS={gene:1,disease:0.8,chemical:0.8};
// Italic has no plain switch in vis-network: the slanted face lives in the multi-font "ital"
// slot, which only applies to text the label marks up. So a drug label carries <i>...</i> and
// the node opts into multi:'html'. Anything that builds a label has to go through labelFor(),
// or a drug would render its own markup as literal text.
function labelFor(id,kind,text){
 if(!LABEL_KINDS.has(kind||'gene'))return '';
 const t=text===undefined?(labelById[id]||''):text;
 return (kind==='chemical')?'<i>'+esc(t)+'</i>':t;
}
function nodeFont(kind,base){
 const f={size:scaledFont(base)};
 if(kind==='disease')f.face=NARROW_FACE;
 if(kind==='chemical'){f.face=NARROW_FACE;f.multi='html';f.ital={face:NARROW_FACE,mod:'italic'};}
 return f;
}
function nodeLabel(n){return labelFor(n.id,n.kind||'gene',n.label);}
function applyLabelKinds(){
 if(!NODEDS)return;
 NODEDS.update(NODEDS.get().map(n=>{const k=KIND[n.id]||'gene';return {id:n.id,label:labelFor(n.id,k)};}));
}
function activeMinCluster(){const v=parseInt((document.getElementById('mincluster')||{}).value);return isNaN(v)?2:v;}
function activeMinDegree(){const v=parseInt((document.getElementById('mindeg')||{}).value);return isNaN(v)?2:v;}
function activeMinPub(){const v=parseInt((document.getElementById('minpub')||{}).value);return isNaN(v)?1:v;}
// Gene/chem focus: expand `hops` steps from the seeds over the CURRENT edge set and keep the
// edges whose both endpoints are in reach. It has to run again after the text lens -- expanding
// text-blind and then letting the query delete edges leaves the neighbourhood the query erased:
// PIK3CA-AKT1 (two MALAT1 neighbours) drawn while every MALAT1 edge is gone. Re-expanding over
// what survived means the seed either appears with its links or the canvas is empty.
function focusKeep(edges,seeds,hops){
 const ag={};edges.forEach(o=>{(ag[o.e.from]=ag[o.e.from]||[]).push(o.e.to);(ag[o.e.to]=ag[o.e.to]||[]).push(o.e.from);});
 const seen=new Set(seeds);let fr=[...seeds];
 for(let h=0;h<hops;h++){const nf=[];fr.forEach(x=>{(ag[x]||[]).forEach(y=>{if(!seen.has(y)){seen.add(y);nf.push(y);}});});fr=nf;}
 return edges.filter(o=>seen.has(o.e.from)&&seen.has(o.e.to));
}
function build(thr){
 const conf=activeConf(), cats=activeCats(); const [ylo,yhi]=activeYears(); const mc=activeMinCluster(); const md=activeMinDegree(); const mp=activeMinPub(); FSCALE=activeFontScale();
 const txt=activeText(), tm=textMatcher(txt); TM=tm;
 const kinds=activeKinds();
 let edges=[];
 // The text query is a LENS, not a threshold input: an edge's support is what survives the
 // score, year, relation and source settings, and the thresholds below judge THAT. Feeding the
 // matched-only sentences to them instead made typing a word collapse the graph -- "lung" cut
 // MALAT1-non-small cell lung carcinoma from 8 publications to 2 (under min-publications 7),
 // and thinning the edge set dropped whole genes under min-connections, whose components then
 // fell under min-cluster: an empty canvas from a query with thousands of hits. So support is
 // measured text-blind here; the lens comes further down, and only min-connections/min-cluster
 // are re-judged on what it leaves (they describe the picture, so they have to).
 // Some node types have no relations OF THEIR OWN: the RE routing pairs a disease with a gene
 // or a chemical, never with another disease, so ticking disease alone can never draw an edge
 // and the canvas went blank. When the ticked types admit no edge anywhere in the file, the
 // nodes are drawn on their own instead -- each one that still has evidence through a partner
 // of an unticked type -- so a disease-only view shows the disease landscape (sized by
 // evidence, pooled by tissue) rather than nothing. `orphan` collects that evidence as we go.
 let kindOK=false; const orphan={};
 DATA.edges.forEach(e=>{
   const kf=!kinds||(kinds.has(KIND[e.from])&&kinds.has(KIND[e.to]));
   if(kf)kindOK=true;
   const sup=visSents(e,conf,ylo,yhi,null,cats,SRC_MODE); if(sup.length<thr)return;
   const np=new Set(sup.map(s=>s.pmid)).size;
   if(mp>1&&np<mp)return;
   if(kf){edges.push({e:e,sup:sup,vis:sup,w:sup.length,np:np,cat:viewCat(e,sup)});return;}
   if(kinds)[e.from,e.to].forEach(nd=>{if(kinds.has(KIND[nd]))(orphan[nd]=orphan[nd]||[]).push(...sup);});
 });
 const gf=(document.getElementById('genefilter').value||'').trim().toLowerCase();
 const chemSel=document.getElementById('chemfilter').value;
 let focusActive=false, focusLabel='', focusSeeds=null, focusHops=1;
 if(gf||chemSel){
   focusActive=true;
   const seeds=new Set(), labs=[];
   if(gf){const fn=DATA.nodes.find(n=>n.label.toLowerCase()===gf)||DATA.nodes.find(n=>n.label.toLowerCase().indexOf(gf)===0);if(fn){seeds.add(fn.id);labs.push(fn.label);}else labs.push('(no gene: '+gf+')');}
   if(chemSel){DATA.nodes.forEach(n=>{if((n.chems||[]).indexOf(chemSel)>=0)seeds.add(n.id);});labs.push('chem: '+chemSel);}
   focusSeeds=seeds;focusHops=parseInt(document.getElementById('hops').value)||1;
   edges=focusKeep(edges,focusSeeds,focusHops);
   focusLabel=labs.join(', ');
 }
 // the lens: keep the edges that still say the word, and show those sentences only
 if(tm)edges=edges.filter(o=>{const v=o.sup.filter(s=>tm.test(s.text));if(!v.length)return false;o.vis=v;o.w=v.length;o.np=new Set(v.map(s=>s.pmid)).size;o.cat=viewCat(o.e,v);return true;});
 // and re-cut the neighbourhood on what the lens left, so a focus view never shows the seed's
 // neighbours to each other with the seed itself missing
 if(tm&&focusActive)edges=focusKeep(edges,focusSeeds,focusHops);
 // Structure comes LAST, on the edges that are actually drawn: min-connections peels the
 // single-link fringe, then min-cluster drops the small components -- so "min cluster size 6"
 // is a statement about the picture, not about some graph behind it. (It once ran before the
 // lens, and typing a word then shattered a 6-node floor into pairs and triples.) The support
 // thresholds above stay text-blind, which is what keeps a query from collapsing the view:
 // min-publications judges an edge's whole evidence, only these two judge the query's leavings.
 // Skipped under gene/chem focus, where you asked for a neighborhood and want all of it.
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
 // --- the ranking, then the significance cut -------------------------------------------------
 // Order matters and is deliberate. The enrichment is computed from the view BEFORE the q cut,
 // and the cut only hides what it judged: if the ranking were recomputed on its own survivors,
 // dragging the slider would move the very numbers it filters on, and each notch would be
 // answering a different question than the one it displays. So the q shown next to a node is
 // always the q it was judged by.
 DRAWN_EDGES=edges;BOOT=null;          // a new view invalidates any bootstrap taken of the old one
 (function(){const pm=new Set(), nd=new Set();
  edges.forEach(o=>{nd.add(o.e.from);nd.add(o.e.to);o.vis.forEach(s=>pm.add(s.pmid));});
  VIEW_PUBS=pm.size;VIEW_ENTS=nd.size;})();
 sigCompute(edges);
 ENRICH=enrichCompute(VIEW_PUBS,VIEW_ENTS);
 QMAP=qMapFor(VIEW_PUBS,VIEW_ENTS);    // q for genes AND drugs, so the cut can judge both
 const qcut=activeQ();
 if(qcut!==null||OR_MODE.size){
  // a node survives if it was tested, reached the threshold and points the way the buttons ask;
  // diseases carry no ranking of their own, so they stay as context rather than being cut on
  // evidence they never had
  const ok=id=>KIND[id]==='disease'||(qPass(id)&&orPass(id));
  edges=edges.filter(o=>ok(o.e.from)&&ok(o.e.to));
 }
 // the disease keep-list, built from what the controls above left and applied after it, so the
 // menu never filters itself out of existence (see buildDiseaseFilter)
 buildDiseaseFilter(edges);
 if(DIS_KEEP.size)edges=edges.filter(o=>disPass(o.e.from)&&disPass(o.e.to));
 sigRender();
 if(SIGTAB_OPEN)sigTable();
 updateCatCounts(edges,cats);   // edges is final here (category, score, year, degree, cluster, text, q)
 // co-mentions ride on the edges that survived: every endpoint whose visible sentences name the
 // chosen disease gets one dashed link to it, with those sentences (deduped) as its evidence
 const cmDis=activeComention();
 let cmLinks=[];
 if(cmDis){
   const re=diseaseAliases(labelById[cmDis]||cmDis), by={};
   // a node already wired to the disease by a drawn relation needs no second, weaker link:
   // the dashed one means "co-mentioned, and nothing predicted it" for the view you are in
   const linked=new Set();edges.forEach(o=>{if(o.e.from===cmDis)linked.add(o.e.to);if(o.e.to===cmDis)linked.add(o.e.from);});
   edges.forEach(o=>{const ms=o.vis.filter(s=>re.test(s.text));if(!ms.length)return;
     [o.e.from,o.e.to].forEach(nd=>{if(nd!==cmDis&&!linked.has(nd))(by[nd]=by[nd]||[]).push(...ms);});});
   cmLinks=Object.keys(by).map(nd=>{const seen=new Set(),ss=[];
     by[nd].forEach(s=>{if(!seen.has(s.text)){seen.add(s.text);ss.push(s);}});
     ss.sort((a,b)=>b.sc-a.sc);return {nd:nd,sents:ss};});
 }
 const keep=new Set();edges.forEach(o=>{keep.add(o.e.from);keep.add(o.e.to);});
 if(cmLinks.length)keep.add(cmDis);   // the disease itself may have no surviving relation edge
 const nss={};edges.forEach(o=>{o.vis.forEach(s=>{(nss[o.e.from]=nss[o.e.from]||new Set()).add(s.text);(nss[o.e.to]=nss[o.e.to]||new Set()).add(s.text);});});
 // the edgeless fallback, sized by the same measure as everything else: unique sentences in view
 const isoSz={}, isoPub=new Set();
 if(kinds&&!kindOK)Object.keys(orphan).forEach(nd=>{
   const t=new Set(),p=new Set();
   orphan[nd].forEach(s=>{if(!tm||tm.test(s.text)){t.add(s.text);p.add(s.pmid);}});
   if(t.size){isoSz[nd]=t.size;keep.add(nd);p.forEach(x=>isoPub.add(x));}});
 const nIso=Object.keys(isoSz).length;
 const nsz=id=>(nss[id]?nss[id].size:(isoSz[id]||0));
 const nodes=DATA.nodes.filter(n=>keep.has(n.id)).map(n=>{const k=n.kind||'gene';
  const fs=fontSize(nsz(n.id))*(KIND_FS[k]||1);
  return {id:n.id,label:nodeLabel(n),value:nsz(n.id),size:scaleNode(nsz(n.id)),shape:n.shape||'dot',title:n.label+((n.kind&&n.kind!=='gene')?'  ['+n.kind+']':'')+' — '+nsz(n.id)+' unique sentences (in view)'+(n.target?' · drug target: '+n.target+' chemicals'+(n.tcat==='green'?' (approved anti-neoplastic)':(n.tcat==='amber'?' (approved)':' (ChEBI)')):''),color:nodeColor(n),_fs:fs,font:nodeFont(k,fs)};});
 // no `value`: vis would then scale the width itself and ignore edgeWidth()
 const eds=edges.map((o,i)=>({id:i,from:o.e.from,to:o.e.to,width:edgeWidth(o.np),color:{color:CCOLOR[o.cat]||o.e.color,opacity:0.6},dashes:o.cat.indexOf('not ')===0,arrows:{to:{enabled:true,scaleFactor:arrowScale(o.np)}},title:edgeTip(o.e,o.vis,o.cat)}));
 // undirected and unarrowed: a shared sentence has no subject and object
 const _cm={};
 cmLinks.forEach((L,k)=>{const id='cm'+k;_cm[id]=L;
   eds.push({id:id,from:L.nd,to:cmDis,width:1,dashes:[3,4],color:{color:'#9aa4b2',opacity:0.45},
             arrows:{to:{enabled:false}},title:cmTip(L,cmDis)});});
 const vpub=new Set();edges.forEach(o=>o.vis.forEach(s=>vpub.add(s.pmid)));isoPub.forEach(p=>vpub.add(p));
 const nkinds=new Set(nodes.map(n=>KIND[n.id]));
 // name what is actually on screen: "378 diseases" beats "378 genes" in a disease-only view
 const KLAB={gene:'genes',disease:'diseases',chemical:'chemicals'};
 document.getElementById('stats').innerHTML='Showing <b>'+nodes.length+'</b> '+(nkinds.size===1?(KLAB[[...nkinds][0]]||'nodes'):'nodes')+', <b>'+edges.length+'</b> edges, <b>'+vpub.size+'</b> publications (&ge;'+conf+')'+(txt?' &middot; text: <b>'+esc(txt)+'</b>':'')+(focusActive?' &middot; focus: <b>'+esc(focusLabel)+'</b>':'')+(nIso?' &middot; <b>'+nIso+'</b> drawn unconnected: nothing in the corpus relates these node types to each other':'')+(cmLinks.length?' &middot; <b>'+cmLinks.length+'</b> co-mention links to <b>'+esc(labelById[cmDis]||cmDis)+'</b>':'');
 const data={nodes:new vis.DataSet(nodes),edges:new vis.DataSet(eds)};
 NODEDS=data.nodes;
 const options={layout:{improvedLayout:false},physics:{stabilization:{iterations:200},barnesHut:{gravitationalConstant:-14000,springLength:130,springConstant:0.02,avoidOverlap:0.3}},interaction:{hover:true,tooltipDelay:120},nodes:{shape:'dot',scaling:{min:6,max:60},font:{color:'rgba(26,26,26,0)'}},edges:{smooth:false,arrowStrikethrough:false,hoverWidth:0,selectionWidth:0,arrows:{to:{enabled:true,scaleFactor:0.6}}}};
 if(network)network.destroy();
 network=new vis.Network(net,data,options);
 LABEL_N=nodes.length; LABEL_A=-1;
 // reshape once physics is off, so nothing pulls the moved nodes back, and fit afterwards so the
 // new positions are inside the view. BASEPOS is captured first: it is the layout to replay from.
 TGROUPS=tissueGroups(nodes);
 network.on('stabilizationIterationsDone',()=>{network.setOptions({physics:false});
   BASEPOS=network.getPositions();
   applyLayoutShape();
   network.fit({animation:false});LABEL_A=-1;updateLabels();network.redraw();});
 network.on('zoom',updateLabels);
 network.on('animationFinished',updateLabels);
 fitLeftPanel();    // self-correcting: a panel narrowed by anything recovers on the next redraw
 const _e=edges;
 network.on('click',p=>{const info=document.getElementById('info');
   if(p.nodes.length){const n=DATA.nodes.find(x=>x.id===p.nodes[0]);info.innerHTML='<b>'+n.label+'</b>: '+nsz(n.id)+' unique sentences (in view)';}
   else if(p.edges.length&&_cm[p.edges[0]]){const L=_cm[p.edges[0]];info.innerHTML=INFO_HEAD+cmTip(L,cmDis).innerHTML;}
   else if(p.edges.length){const o=_e[p.edges[0]];info.innerHTML=INFO_HEAD+edgeHead(o.e,o.vis,o.cat)+o.vis.map(s=>'<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span> '+hl(s.text)+'</div>').join('');}});
}
const thr=document.getElementById('thr');
let CATTOT={};
// One edge counts under EVERY relation its sentences carry (a nine-sentence edge holding both
// 'interacts' and 'activates' readings counts once under each), so the totals add up to more
// than the edge count -- they answer "how many edges can show me this label", which is the
// question the tick boxes and the "in view" half answer too.
function catsOf(e){const c={};e.sents.forEach(s=>c[sentCat(e,s)]=1);return Object.keys(c);}
function buildCatFilters(){CATTOT={};DATA.edges.forEach(e=>catsOf(e).forEach(c=>CATTOT[c]=(CATTOT[c]||0)+1));const cats=Object.keys(CATTOT).sort((a,b)=>CATTOT[b]-CATTOT[a]);document.getElementById('catfilters').innerHTML=cats.map(c=>'<label><input type=checkbox class=catf value="'+esc(c)+'" checked> <span class=sw style="background:'+(CCOLOR[c]||'#888')+'"></span> '+esc(c)+' <span class=cnt data-cat="'+esc(c)+'">('+CATTOT[c]+')</span></label>').join('');document.querySelectorAll('.catf').forEach(c=>c.addEventListener('change',()=>build(+thr.value)));}
// Counts are LIVE: "(total · N in view)" is recomputed from the sentences actually drawn, so
// "N in view" is exactly the number of drawn edges that can show you a sentence tagged with
// that relation -- 0 means no visible sentence carries the label, never "0 but there it is in
// the tooltip". The total is the whole payload (everything qualifying at >=0.5), which is NOT
// what you see: the score slider, year range, text filter, min-publications, min-connections
// and min-cluster size all prune afterwards. A category pruned to nothing reads 0 (in red)
// instead of looking available: ticking it alone would leave the canvas blank, typically
// because its edges form components smaller than "Min cluster size".
function updateCatCounts(edges,cats){const seen={};edges.forEach(o=>{const c={};o.vis.forEach(s=>c[sentCat(o.e,s)]=1);for(const k in c)seen[k]=(seen[k]||0)+1;});
 document.querySelectorAll('#catfilters .cnt').forEach(el=>{const c=el.getAttribute('data-cat');
  el.textContent=cats.has(c)?('('+CATTOT[c]+' · '+(seen[c]||0)+' in view)'):('('+CATTOT[c]+' · off)');
  el.style.color=(cats.has(c)&&!seen[c])?'#b3243b':'';});}
thr.addEventListener('input',()=>{document.getElementById('thv').textContent=thr.value;build(+thr.value);});
const mpb=document.getElementById('minpub');
mpb.addEventListener('input',()=>{document.getElementById('mpv').textContent=mpb.value;build(+thr.value);});
const fsc=document.getElementById('fscale');
fsc.addEventListener('input',()=>{FSCALE=activeFontScale();document.getElementById('fsv').textContent=fsc.value+'%';
 // rebuild the WHOLE font object per node, not just its size: passing {size} alone would drop
 // the condensed face and the italic slot the disease and drug labels carry
 if(NODEDS)NODEDS.update(NODEDS.get().map(n=>({id:n.id,font:nodeFont(KIND[n.id]||'gene',n._fs||13)})));});
// Each "i" reveals its own section's prose. Most sections keep their help inside the row, but
// a few explain a control from BELOW it -- the relation counts, the training-set hint -- and
// those are tagged data-help="<group>" so one button can open every part of its section at once.
document.querySelectorAll('.ihelp').forEach(b=>b.addEventListener('click',()=>{
 const g=b.getAttribute('data-help');
 const row=b.closest?b.closest('.row'):null;
 const hs=g?[...document.querySelectorAll('.help[data-help="'+g+'"]')]:(row?[row.querySelector('.help')]:[]);
 const targets=hs.filter(Boolean);
 if(!targets.length)return;
 const open=!targets[0].classList.contains('open');
 targets.forEach(h=>h.classList.toggle('open',open));
 b.classList.toggle('on',open);
 b.setAttribute('aria-expanded',open?'true':'false');}));
// the position controls all replay from BASEPOS, so they reshape in place -- no rebuild, no relayout
// two handles, one cutoff: the panel's and the table's mirror each other and rebuild the view
function qChanged(fromTable){
 const a=document.getElementById('qsig'), b=document.getElementById('qsig2');
 if(fromTable&&b)a.value=b.value; else if(b)b.value=a.value;
 const t=qLabel();
 ['qsv','qsv2'].forEach(id=>{const e=document.getElementById(id);if(e)e.innerHTML=t;});
 build(+thr.value);
}
['qsig','qsig2'].forEach((id,i)=>{const el=document.getElementById(id);
 if(el)el.addEventListener('input',()=>qChanged(i===1));});
function orPaint(){document.querySelectorAll('.orb,.orb2').forEach(b=>
  b.classList.toggle('on',OR_MODE.has(b.getAttribute('data-or'))));}
document.querySelectorAll('.orb,.orb2').forEach(b=>b.addEventListener('click',()=>{
 const m=b.getAttribute('data-or');
 if(OR_MODE.has(m))OR_MODE.delete(m);else OR_MODE.add(m);   // each button toggles on its own
 orPaint();
 build(+thr.value);
}));
(function(){const c=document.getElementById('disclear');
 if(c)c.addEventListener('click',()=>{if(!DIS_KEEP.size)return;DIS_KEEP.clear();build(+thr.value);});})();
const shr=document.getElementById('shrink'), exp=document.getElementById('expand');
shr.addEventListener('input',()=>{document.getElementById('shv').textContent=shr.value+'%';applyLayoutShape();});
exp.addEventListener('input',()=>{document.getElementById('exv').textContent=exp.value+'%';applyLayoutShape();});
// name buttons: gene symbols are always drawn, these two add the wordy kinds back
document.querySelectorAll('.lblb').forEach(b=>b.addEventListener('click',()=>{
 const k=b.getAttribute('data-kind');
 if(LABEL_KINDS.has(k))LABEL_KINDS.delete(k);else LABEL_KINDS.add(k);
 b.classList.toggle('on',LABEL_KINDS.has(k));
 applyLabelKinds();}));
// the gene-only graph has nothing else to name
(function(){if(!DATA.nodes.some(n=>(n.kind||'gene')!=='gene')){const r=document.getElementById('lblrow');if(r)r.style.display='none';}})();
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
// --- significance in view ----------------------------------------------------------------
// The graph shows which relations survive your filters; this ranks the NODES behind them, so
// you can ask "and who carries this view" without counting edges by eye. It is computed from
// the drawn edge list on every redraw -- open every filter and it reads as the corpus, narrow
// them and it answers the same question of a slice.
//   publications = distinct papers behind the node's relations (what edge thickness uses)
//   partners     = distinct entities it is related to -- breadth, not weight
//   sentences    = unique sentences supporting those relations
// Publications leads because it is the one that cannot be inflated by a single talkative paper.
// Computed from the edges the graph is actually DRAWING, so the ranking is always an answer
// about the picture in front of you: every control that shapes the graph -- score, year, text,
// relation type, training set, node type, the support thresholds and the structural pruning --
// reshapes this too, and build() recomputes it on every pass. Co-mention links are excluded:
// they are co-occurrence, and this counts relations.
let SIG=[];
function sigCompute(drawn){
 const acc={};
 const touch=id=>acc[id]||(acc[id]={id:id,label:labelById[id]||id,kind:KIND[id]||'gene',
                                    pubs:new Set(),partners:new Set(),sents:new Set()});
 (drawn||[]).forEach(o=>{
  const a=touch(o.e.from),b=touch(o.e.to);
  a.partners.add(o.e.to);b.partners.add(o.e.from);
  o.vis.forEach(s=>{a.pubs.add(s.pmid);b.pubs.add(s.pmid);a.sents.add(s.text);b.sents.add(s.text);});
 });
 SIG=Object.keys(acc).map(id=>{const n=acc[id];
   return {id:id,label:n.label,kind:n.kind,pub:n.pubs.size,deg:n.partners.size,sent:n.sents.size};});
 sigStats();
}
// Two statistics per node, per measure, computed within its own kind -- a gene is only
// remarkable among genes. PERCENTILE is a mid-rank over the whole corpus (ties share their
// rank), which needs no assumption about the distribution; this one is brutally heavy-tailed,
// so a mean would be meaningless. Z is on log10 counts, where the tail is closer to symmetric,
// and answers the other question: not "how many below you" but "how far out".
function sigStats(){
 ['gene','disease','chemical'].forEach(k=>{
  const rows=SIG.filter(s=>s.kind===k), n=rows.length;
  if(!n)return;
  ['pub','deg','sent'].forEach(m=>{
   const vals=rows.map(s=>s[m]).sort((a,b)=>a-b);
   const lg=rows.map(s=>Math.log10(Math.max(1,s[m])));
   const mu=lg.reduce((a,b)=>a+b,0)/n;
   const sd=Math.sqrt(lg.reduce((a,b)=>a+(b-mu)*(b-mu),0)/n)||1;
   const bound=(v,inc)=>{let lo=0,hi=n;while(lo<hi){const mid=(lo+hi)>>1;
     if(inc?vals[mid]<=v:vals[mid]<v)lo=mid+1;else hi=mid;}return lo;};
   rows.forEach(s=>{const b=bound(s[m],false),e=bound(s[m],true);
    (s.pct=s.pct||{})[m]=100*(b+(e-b)/2)/n;
    (s.z=s.z||{})[m]=(Math.log10(Math.max(1,s[m]))-mu)/sd;});
  });
 });
}
function pctStr(p){return (p>=99.95?'99.9':p.toFixed(1))+'%';}
function zStr(z){return (z>=0?'+':'')+z.toFixed(1);}
// --- bootstrap over publications ---------------------------------------------------------
// The counts are one draw from the literature; these say how much of the ranking would survive
// another. Resample the PAPERS of the current view with replacement -- papers, because that is
// the independent unit: sentences inside one article are the same authors saying the same thing
// twice, and resampling them would give intervals several times too tight. Recompute the measure
// and the whole ordering per replicate, then take percentile intervals of each node's count and
// of its rank. A rank interval spanning 3-40 is the honest reading of "top ten": the leaderboard
// is a sample statistic, and this is the width of the uncertainty behind it.
// --- enrichment against the corpus --------------------------------------------------------
// The ranking says who carries this view; this says whether they carry it MORE than they carry
// the corpus. Per entity, a 2x2 over publications -- in view / not, mentions it / not -- tested
// with Fisher's exact test against BG, the document counts over every normalized triple in the
// file, score-unfiltered on purpose: a corpus already filtered by the cutoff you are testing is
// no denominator at all. The odds ratio carries a Woolf interval on the 0.5-corrected log OR
// (an approximation, not the conditional MLE), and q is Benjamini-Hochberg over the entities
// actually testable in this view.
//
// What it cannot say: whether a gene is specific to lung adenocarcinoma. Every paper here is a
// lung paper, so the contrast is view-against-corpus, never corpus-against-literature.
// log-factorials, grown on demand: sized to the corpus is enough for the tables this page
// builds (the four cells always sum to N), but a table that silently runs off its end returns
// exp(NaN) = 0 -- a p-value of zero that looks like a discovery. So it extends instead.
let LGF=(function(){const N=(BG&&BG.n?BG.n:0)+8, a=new Float64Array(N+1);
 for(let i=2;i<=N;i++)a[i]=a[i-1]+Math.log(i);return a;})();
function lgf(i){
 if(i>=LGF.length){const old=LGF, n=Math.max(i+1,old.length*2), a=new Float64Array(n);
  a.set(old);for(let k=old.length;k<n;k++)a[k]=a[k-1]+Math.log(k);LGF=a;}
 return LGF[i];
}
function lhyper(a,r1,r2,n){return lgf(r1)-lgf(a)-lgf(r1-a)+lgf(n-r1)-lgf(r2-a)-lgf(n-r1-r2+a)+lgf(r2)+lgf(n-r2)-lgf(n);}
function fisherP(a,b,c,d){                      // two-sided, by summing tables no likelier than seen
 const n=a+b+c+d, r1=a+b, r2=a+c;
 if(!n||!r1||!r2||r1===n||r2===n)return 1;
 const lo=Math.max(0,r1+r2-n), hi=Math.min(r1,r2);
 const obs=lhyper(a,r1,r2,n), eps=1e-9;
 let p=0;
 for(let k=lo;k<=hi;k++){const l=lhyper(k,r1,r2,n);if(l<=obs+eps)p+=Math.exp(l);}
 return Math.min(1,p);
}
function oddsRatio(a,b,c,d){                    // Haldane-Anscombe: 0.5 keeps a zero cell finite
 const A=a+0.5,B=b+0.5,C=c+0.5,D=d+0.5, or=(A*D)/(B*C);
 const se=Math.sqrt(1/A+1/B+1/C+1/D), l=Math.log(or);
 return {or:or,lo:Math.exp(l-1.96*se),hi:Math.exp(l+1.96*se)};
}
const ENRICH_MIN=5;                             // below ~5 in the corpus Fisher cannot reach significance
// The unit follows the ranking. By publications the table is papers-with-the-entity out of the
// corpus' papers; by PARTNERS it is entities-related-to-it out of the corpus' entities, against
// a separate partner background -- testing breadth against a paper denominator would be
// comparing two different things. Sentences fall back to the publication test on purpose: they
// are pseudo-replicates, and a Fisher test on them would be anticonservative by a large factor.
function enrichUnit(){return sigMeas()==='deg'?'par':'doc';}
// The corpus count is known for every node whether or not its row can be TESTED -- it comes
// straight from BG. Tying it to the test result printed a blank next to a positive view count,
// which reads as "0 corpus papers behind 37 in view", an impossibility: the corpus contains the
// view. Untestable now means an empty OR and q, never an empty denominator.
function bgCount(s){const m=(BG&&BG[enrichUnit()]&&BG[enrichUnit()][s.kind])||{};
 const v=m[s.id];return v===undefined?null:v;}
function enrichCompute(viewPubs,viewEnts){
 const kind=sigKind(), unit=enrichUnit();
 const bg=(BG&&BG[unit]&&BG[unit][kind])||{};
 const N=(unit==='par'?(BG&&BG.ne):(BG&&BG.n))||0;
 const V=unit==='par'?viewEnts:viewPubs;
 const val=s=>unit==='par'?s.deg:s.pub;
 if(!N||!V)return null;
 const out={}, tested=[];
 let nocontrast=0, belowfloor=0;
 SIG.forEach(s=>{
  if(s.kind!==kind)return;
  const c=bg[s.id];                             // corpus documents (or partners) for it
  if(!c||c<ENRICH_MIN){belowfloor++;return;}
  const a=Math.min(val(s),c), b=Math.max(0,V-a), cc=Math.max(0,c-a), d=Math.max(0,N-V-cc);
  // Nothing outside the view: every paper (or partner) this entity has is already in front of
  // you, so there is no contrast to test. The Haldane 0.5 would still hand back a huge odds
  // ratio and a tiny q -- an artefact of the correction, not a finding -- so the row says
  // nothing instead. On an unfiltered view that is every row, which is the honest answer.
  if(cc<1){nocontrast++;return;}
  const r=oddsRatio(a,b,cc,d);
  out[s.id]={a:a,c:c,or:r.or,lo:r.lo,hi:r.hi,p:fisherP(a,b,cc,d)};
  tested.push(s.id);
 });
 tested.sort((x,y)=>out[x].p-out[y].p);         // Benjamini-Hochberg over what was testable here
 const m=tested.length;
 let prev=1;
 for(let i=m-1;i>=0;i--){const id=tested[i];
  prev=Math.min(prev,out[id].p*m/(i+1));out[id].q=prev;}
 return {kind:kind,unit:unit,n:N,view:V,tested:m,nocontrast:nocontrast,belowfloor:belowfloor,rows:out};
}
function orStr(e){return e?e.or.toFixed(2)+' <span class=mut>('+e.lo.toFixed(2)+'&ndash;'+e.hi.toFixed(2)+')</span>':'&mdash;';}
function qStr(e){return e?(e.q<1e-4?e.q.toExponential(1):e.q.toFixed(4)):'&mdash;';}
let DRAWN_EDGES=[], BOOT=null, ENRICH=null, VIEW_PUBS=0, VIEW_ENTS=0, QMAP={}, ORMAP={};
// --- direction of enrichment --------------------------------------------------------------
// The q cut says how sure; these two say which way. An odds ratio above 1 means the entity is
// over-represented in this view against the shared corpus, below 1 under-represented -- and the
// depleted side is a finding, not a leftover: in an adenocarcinoma view the small-cell markers
// DLL3 and ASCL1 land there at q well under 0.05. Untested entities have no OR and are hidden by
// either button, the same rule the q slider follows; diseases carry no ranking and are exempt.
const OR_TESTS={gt1:{fn:v=>v>1},lt1:{fn:v=>v<1}};
// A SET, not one mode: the two directions are independent filters, so both can be lit at once.
// Both on is not the same as neither on -- it keeps every entity the test placed on a side and
// drops the ones it could not place, which is the "show me only what was actually measured"
// view. (An odds ratio of exactly 1 is on neither side and falls out too; with the Haldane
// correction that is a near-empty case, but it is the honest reading of "over" and "under".)
const OR_MODE=new Set();
function orPass(id){
 if(!OR_MODE.size)return true;
 const v=ORMAP[id];
 if(v===undefined)return false;      // untested: on neither side of 1, so no button keeps it
 for(const m of OR_MODE)if(OR_TESTS[m].fn(v))return true;
 return false;
}
// The q cut hides what the test JUDGED and rejected. An entity that could not be tested at all --
// corpus count below ENRICH_MIN, or nothing outside the view to contrast against -- was never
// judged, so it is not "below the cutoff": it has no cutoff to be below. Treating the two alike
// deleted 1,502 of this corpus' 2,904 ranked genes at q<=0.05 -- more than half the graph, and
// the whole long tail of it -- on evidence that was never gathered. So the slider keeps the
// untested and cuts only the judged.
// The OR buttons ask a DIFFERENT question -- which side of 1 the entity sits on -- and an
// untested entity sits on neither, so orPass() above still hides it. That asymmetry is the point:
// "not shown to be enriched" and "shown not to be enriched" are not the same claim.
function qPass(id){const q=activeQ();return q===null||QMAP[id]===undefined||QMAP[id]<=q;}
// q for BOTH rankable kinds, so the slider can cut genes and drugs in one pass. The table shows
// one kind at a time; the graph has to judge whatever it draws.
function qMapFor(viewPubs,viewEnts){
 const keep=sigKind(), out={};
 ORMAP={};
 ['gene','chemical'].forEach(k=>{
  const sel=document.getElementById('sigkind');
  if(sel)sel.value=k;                  // enrichCompute reads the selector; borrow it, then restore
  const e=enrichCompute(viewPubs,viewEnts);
  if(e)Object.keys(e.rows).forEach(id=>{out[id]=e.rows[id].q;ORMAP[id]=e.rows[id].or;});
 });
 const sel=document.getElementById('sigkind');
 if(sel)sel.value=keep;
 return out;
}
// Stops, not a linear range: q is read on a log scale and the conventional cutoffs are what a
// reader wants to land on. Index 0 is off, so the graph opens unfiltered.
const Q_STEPS=[null,0.5,0.2,0.1,0.05,0.01,0.001,1e-4,1e-5,1e-6];
function activeQ(){const el=document.getElementById('qsig');
 let i=el?parseInt(el.value):0;if(isNaN(i))i=0;
 return Q_STEPS[Math.max(0,Math.min(Q_STEPS.length-1,i))];}
function qLabel(){const q=activeQ();
 return q===null?'off (show all)':('q &le; '+(q>=0.001?q:q.toExponential(0)));}
function bootstrapCIs(B){
 const kind=sigKind(), meas=sigMeas();
 const kidx=new Map(), ids=[];
 SIG.forEach(s=>{if(s.kind===kind){kidx.set(s.id,ids.length);ids.push(s.id);}});
 const K=ids.length;
 if(!K)return null;
 const pidx=new Map(), inc=[], eList=[];
 const pi=pm=>{let k=pidx.get(pm);if(k===undefined){k=inc.length;pidx.set(pm,k);inc.push(null);}return k;};
 // per paper: which nodes it supports, and with what weight (1 per paper, or its sentence count)
 const acc=[];
 DRAWN_EDGES.forEach(o=>{
  const a=kidx.has(o.e.from)?kidx.get(o.e.from):-1, b=kidx.has(o.e.to)?kidx.get(o.e.to):-1;
  const pset=new Set();
  o.vis.forEach(s=>{
   const p=pi(s.pmid);pset.add(p);
   if(!acc[p])acc[p]=new Map();
   [a,b].forEach(n=>{if(n<0)return;let t=acc[p].get(n);if(!t){t=new Set();acc[p].set(n,t);}t.add(s.text);});
  });
  if(a>=0||b>=0)eList.push({a:a,b:b,p:[...pset]});
 });
 const P=inc.length;
 for(let p=0;p<P;p++){
  const m=acc[p]||new Map(), n=[], w=[];
  m.forEach((texts,node)=>{n.push(node);w.push(meas==='sent'?texts.size:1);});
  inc[p]={n:Int32Array.from(n),w:Float64Array.from(w)};
 }
 const counts=[], ranks=[];
 for(let k=0;k<K;k++){counts.push(new Float64Array(B));ranks.push(new Int32Array(B));}
 const val=new Float64Array(K), draw=new Int32Array(P), order=new Int32Array(K);
 for(let k=0;k<K;k++)order[k]=k;
 for(let r=0;r<B;r++){
  draw.fill(0);
  for(let i=0;i<P;i++)draw[(Math.random()*P)|0]++;
  val.fill(0);
  if(meas==='deg'){                     // distinct partners: an edge counts once if it survives
   for(let i=0;i<eList.length;i++){const e=eList[i];
    let live=false;
    for(let q=0;q<e.p.length;q++)if(draw[e.p[q]]>0){live=true;break;}
    if(!live)continue;
    if(e.a>=0)val[e.a]++;
    if(e.b>=0)val[e.b]++;}
  }else{                                // publications / sentences: additive over sampled papers
   for(let p=0;p<P;p++){const c=draw[p];if(!c)continue;
    const {n,w}=inc[p];
    for(let q=0;q<n.length;q++)val[n[q]]+=c*w[q];}
  }
  const ord=Array.prototype.slice.call(order).sort((x,y)=>val[y]-val[x]);
  for(let k=0;k<K;k++){const node=ord[k];ranks[node][r]=k+1;}
  for(let k=0;k<K;k++)counts[k][r]=val[k];
 }
 const q=(arr,p)=>{const a=Array.prototype.slice.call(arr).sort((x,y)=>x-y);
   return a[Math.min(a.length-1,Math.max(0,Math.round(p*(a.length-1))))];};
 const out={kind:kind,meas:meas,B:B,count:{},rank:{}};
 ids.forEach((id,k)=>{out.count[id]=[q(counts[k],0.025),q(counts[k],0.975)];
                      out.rank[id]=[q(ranks[k],0.025),q(ranks[k],0.975)];});
 return out;
}
const SIG_TOP=6;                       // a teaser in the panel; the full ranking is the table
const SIG_LAB={pub:'publications',deg:'partners',sent:'sentences'};
// The chemical nodes are whatever ChEBI recognised in the text, which mixes therapeutics with
// reagents -- LY294002 ranks high here and has never been in a patient. Say so wherever the
// drug list is on screen, rather than leaving the reader to spot it.
const SIG_DRUGNOTE='The list contains both (1) clinically used drugs and (2) lab chemicals used only in experiments.';
function sigNote(id,kind){const el=document.getElementById(id);if(!el)return;
 el.innerHTML=kind==='chemical'?SIG_DRUGNOTE:'';
 el.style.display=kind==='chemical'?'':'none';}
let SIG_ROWS=[];
function sigRender(){
 const kind=(document.getElementById('sigkind')||{}).value||'gene';
 const meas=(document.getElementById('sigmeasure')||{}).value||'pub';
 SIG_ROWS=SIG.filter(s=>s.kind===kind).sort((a,b)=>(b[meas]-a[meas])||a.label.localeCompare(b.label)).slice(0,SIG_TOP);
 document.getElementById('siglist').innerHTML=SIG_ROWS.map((s,i)=>
   '<div class="sig" data-i="'+i+'" title="'+esc(s.label)+' &mdash; '+s.pub+' publications &middot; '
   +s.deg+' partners &middot; '+s.sent+' sentences (corpus-wide) &middot; percentile '+pctStr(s.pct[meas])
   +' &middot; z '+zStr(s.z[meas])+' among '+s.kind+'s">'
   +'<span class=nm>'+(i+1)+'. '+esc(s.label)+'</span>'
   // one number, the ranked measure -- the other two are a column click away in the table
   +'<span class=vl>'+s[meas]+'</span>'
   // z, not the percentile: every one of a top-6 sits at 99.9%, while z still separates them
   +'<span class=pc>'+zStr(s.z[meas])+'</span></div>').join('')
   ||'<div class=mut>no '+esc(kind)+' nodes in this graph</div>';
 const hdr=document.getElementById('sighdr');
 if(hdr)hdr.innerHTML=SIG_LAB[meas]+' &nbsp;&nbsp;z';
 sigNote('signote',kind);
 document.querySelectorAll('#siglist .sig').forEach(el=>
   el.addEventListener('click',()=>sigFocus(SIG_ROWS[parseInt(el.getAttribute('data-i'))])));
}
// select and centre the node; a node the current filters removed reports that instead of
// silently doing nothing (or worse, moving the view to where it is not)
function sigFocus(s){
 if(!s)return;
 const info=document.getElementById('info');
 const head='<b>'+esc(s.label)+'</b> <span class=mut>(corpus-wide: '+s.pub+' publications &middot; '
   +s.deg+' partners &middot; '+s.sent+' sentences)</span>';
 try{network.selectNodes([s.id]);network.focus(s.id,{scale:1.3,animation:true});info.innerHTML=head;}
 catch(err){info.innerHTML=head+'<div class=mut>Not in the current view &mdash; loosen the filters to see it.</div>';}
}
// Two pairs of selects, one range: the panel's are the source of truth and the table's mirror
// them, so switching to drugs inside the table leaves the panel showing drugs when you close it.
function sigSync(fromTable){
 const k=document.getElementById('sigkind'), m=document.getElementById('sigmeasure');
 const k2=document.getElementById('sigkind2'), m2=document.getElementById('sigmeasure2');
 if(fromTable&&k2){k.value=k2.value;m.value=m2.value;}
 else if(k2){k2.value=k.value;m2.value=m.value;}
 SIGTAB_SORT=null;                      // the ranking changed; drop a column sort from before it
 ENRICH=enrichCompute(VIEW_PUBS,VIEW_ENTS);   // by-partners tests partners, not papers
 sigRender();
 if(SIGTAB_OPEN)sigTable();
}
['sigkind','sigmeasure'].forEach(id=>{const el=document.getElementById(id);if(el)el.addEventListener('change',()=>sigSync(false));});
['sigkind2','sigmeasure2'].forEach(id=>{const el=document.getElementById(id);if(el)el.addEventListener('change',()=>sigSync(true));});
// NB: the first sigCompute() runs with the year sliders, below -- before they carry values
// activeYears() is NaN and every sentence falls outside the window
// --- table view -------------------------------------------------------------------------
// The panel list is a top-15 teaser; this is the whole ranking, every node of the kind with
// all three counts and both statistics, sortable by any column. It replaces the canvas rather
// than floating over it -- reading a table and reading a graph are different jobs.
let SIGTAB_OPEN=false, SIGTAB_SORT=null, SIGTAB_ROWS=[];
const SIGCOLS=[['nm','node',s=>esc(s.label)],['pub','publications',s=>s.pub],['deg','partners',s=>s.deg],
               ['sent','sentences',s=>s.sent],['pct','percentile',s=>pctStr(s.pct[sigMeas()])],
               ['z','z (log₁₀)',s=>zStr(s.z[sigMeas()])],
               ['cic','count 95% CI',s=>bootStr(BOOT&&BOOT.count[s.id])],
               ['cir','rank 95% CI',s=>bootStr(BOOT&&BOOT.rank[s.id])],
               ['corpus',()=>enrichUnit()==='par'?'corpus partners':'corpus papers',s=>{const v=bgCount(s);return v===null?'&mdash;':v;}],
               ['or','OR vs corpus',s=>orStr(ENRICH&&ENRICH.rows[s.id])],
               ['q','q (BH)',s=>qStr(ENRICH&&ENRICH.rows[s.id])],
               // the bar is drawn from whichever measure the table is sorted by, so it says which
               ['bar',()=>'relative '+SIG_LAB[barMeas()],null]];
function barMeas(){const m=SIGTAB_SORT||sigMeas();return ['pub','deg','sent'].indexOf(m)>=0?m:sigMeas();}
// Columns only exist once they carry something. Twelve of them ran off the right edge of the
// space left beside the panel, which is how a computed enrichment can look like a missing one:
// the bootstrap pair appears when you press the button, the enrichment trio when there is a
// contrast to test, and the table stays narrow enough to read until then.
function sigCols(){return SIGCOLS.filter(([k])=>
  (k==='cic'||k==='cir')?!!BOOT
  :(k==='corpus')?!!(BG&&BG.n)
  :((k==='or'||k==='q')?!!(ENRICH&&ENRICH.tested):true));}
const SIGSORTABLE=new Set(['nm','pub','deg','sent','pct','z','corpus','or','q']);   // intervals are not a sort key
function sigSortVal(s,k){                       // enrichment columns sort on their own numbers
 const e=ENRICH&&ENRICH.rows[s.id];
 if(k==='corpus'){const v=bgCount(s);return v===null?-1:v;}
 if(k==='or')return e?e.or:-1;
 if(k==='q')return e?-e.q:-Infinity;            // smallest q first
 return s[k];
}
function bootStr(iv){return iv?(iv[0]===iv[1]?String(iv[0]):iv[0]+'&ndash;'+iv[1]):'&mdash;';}
function sigKind(){return (document.getElementById('sigkind')||{}).value||'gene';}
function sigMeas(){return (document.getElementById('sigmeasure')||{}).value||'pub';}
function sigTable(){
 const kind=sigKind(), meas=SIGTAB_SORT||sigMeas();
 // the same two predicates the graph applies, so the views can never disagree about what survives
 SIGTAB_ROWS=SIG.filter(s=>s.kind===kind)
   .filter(s=>qPass(s.id)&&orPass(s.id))
   .sort((a,b)=>{
   if(meas==='nm')return a.label.localeCompare(b.label);
   if(meas==='pct'||meas==='z')return (b[meas][sigMeas()]-a[meas][sigMeas()])||a.label.localeCompare(b.label);
   return (sigSortVal(b,meas)-sigSortVal(a,meas))||a.label.localeCompare(b.label);});
 const KL={gene:'genes',chemical:'drugs',disease:'diseases'};
 document.getElementById('sigttl').innerHTML='Significance in view &mdash; '+(KL[kind]||kind);
 const yr=activeYears();
 document.getElementById('sigsub').innerHTML=SIGTAB_ROWS.length+' '+(KL[kind]||kind)+' in the graph as currently drawn &mdash; <b>'
  +yr[0]+'&ndash;'+yr[1]+'</b>, score &ge;'+activeConf()+', and every other filter in force. '
  +'Percentile and z are computed among those '+(KL[kind]||kind)+' for <b>'+
  ({pub:'publications',deg:'partners',sent:'sentences'}[sigMeas()])+'</b>; click a heading to sort. The controls on the right stay live &mdash; filter while you read and the table follows.'
  +(ENRICH?(' &middot; enrichment counts '+(ENRICH.unit==='par'?'<b>partners</b>':'<b>publications</b>')
    +': this view’s <b>'+ENRICH.view+'</b> against the shared lung corpus'
    +((BG&&BG.src&&BG.src.length>1)?(' ('+BG.src.join(' + ')+')'):'')+'’ <b>'+ENRICH.n
    +'</b>, Fisher exact on <b>'+ENRICH.tested+'</b> testable '+(KL[kind]||kind)+' (≥'+ENRICH_MIN
    +' in the corpus), q = Benjamini-Hochberg.'
    +(ENRICH.belowfloor?(' <b>'+ENRICH.belowfloor+'</b> fall below that floor'):'')
    +(ENRICH.nocontrast?((ENRICH.belowfloor?' and ':' ')+'<b>'+ENRICH.nocontrast+'</b> have nothing outside the view to contrast against'):'')
    +((ENRICH.belowfloor||ENRICH.nocontrast)?(', so their OR and q are blank &mdash; their corpus counts are still shown'
      +(ENRICH.tested?'':'; narrow the view to make the comparison mean something')+'.'):'')):'');
 sigNote('signote2',kind);
 // when nothing can be tested the OR/q columns are absent; without this they would just look missing
 const nt=document.getElementById('signotest');
 if(nt){
  const dead=ENRICH&&!ENRICH.tested;
  nt.style.display=dead?'':'none';
  if(dead)nt.innerHTML='<b>No enrichment to test at this view.</b> All <b>'+ENRICH.nocontrast+'</b> '
    +(KL[kind]||kind)+' with enough corpus '+(ENRICH.unit==='par'?'partners':'papers')
    +' already have every one of them on screen, so there is nothing outside the view to compare against and '
    +'<em>OR vs corpus</em> and <em>q</em> are undefined rather than missing. Narrow the view to create the contrast '
    +'&mdash; raising <b>Min unique publications</b> to 2 or <b>Min connections</b> to 2 is usually enough.';
 }
 // the bar column tracks whatever the table is sorted by, scaled to the leading row
 const bmeas=barMeas();
 const bmax=Math.max(1,...SIGTAB_ROWS.map(s=>s[bmeas]));
 const COLS=sigCols();
 const head='<tr>'+COLS.map(([k,lab])=>'<th class="'+(k==='nm'?'nm':k)+(k===meas?' on':'')+'" data-k="'+k+'">'
   +(typeof lab==='function'?lab():lab)+'</th>').join('')+'</tr>';
 const body=SIGTAB_ROWS.map((s,i)=>'<tr data-i="'+i+'">'+COLS.map(([k,lab,f])=>
   k==='bar'?'<td class=bar><span style="width:'+Math.max(1,Math.round(100*s[bmeas]/bmax))+'%"></span></td>'
   :'<td class="'+(k==='nm'?'nm':'')+'">'+(k==='nm'?(i+1)+'. '+f(s):f(s))+'</td>').join('')+'</tr>').join('');
 document.getElementById('sigtbody').innerHTML='<table><thead>'+head+'</thead><tbody>'+body+'</tbody></table>';
 document.querySelectorAll('#sigtbody th').forEach(th=>{const k=th.getAttribute('data-k');
   if(SIGSORTABLE.has(k))th.addEventListener('click',()=>{SIGTAB_SORT=k;sigTable();});
   else th.style.cursor='default';});
 // rows are not links: a click used to close the table and jump the graph to that node, which
 // is the wrong default while you are reading a table. Headings still sort; the panel's top-six
 // list still navigates.
}
// Leave room for the right panel and lift it above the overlay, so every filter stays usable
// while the table is open -- build() re-runs the ranking and redraws the table under it. On a
// narrow screen there is no room to sit side by side, so the table covers everything as before.
function layoutTable(){
 const st=document.getElementById('sigtable'), pn=document.getElementById('panel');
 if(!st||!pn)return;
 if(!SIGTAB_OPEN||window.innerWidth<=700){st.style.right='';pn.classList.remove('overtable');return;}
 const w=Math.ceil(pn.getBoundingClientRect().width);
 st.style.right=(w>0?w+24:0)+'px';
 pn.classList.add('overtable');
}
function sigOpen(){SIGTAB_OPEN=true;SIGTAB_SORT=null;
 const k2=document.getElementById('sigkind2'),m2=document.getElementById('sigmeasure2');
 if(k2){k2.value=sigKind();m2.value=sigMeas();}   // open showing what the panel was showing
 sigTable();document.getElementById('sigtable').classList.add('open');layoutTable();}
function sigClose(){SIGTAB_OPEN=false;document.getElementById('sigtable').classList.remove('open');layoutTable();}
// one source for the explanation: the panel's copy is authored in the template, the table's is
// filled from it at load, so the two can never drift apart
(function(){[['sighelp','sighelp2'],['qhelp','qhelp2']].forEach(([x,y])=>{
 const a=document.getElementById(x),b=document.getElementById(y);
 if(a&&b)b.innerHTML=a.innerHTML;});})();
// run on demand, not on every redraw: it is the one thing here that costs real work, and it is
// only meaningful once you have settled on the view you want to read
const BOOT_B=300;
document.getElementById('sigboot').addEventListener('click',()=>{
 const b=document.getElementById('sigboot');
 b.textContent='resampling…';b.disabled=true;
 setTimeout(()=>{                       // let the label paint before the loop blocks the thread
  try{BOOT=bootstrapCIs(BOOT_B);}catch(err){BOOT=null;}
  b.textContent=BOOT?(BOOT_B+'× bootstrap'):'Bootstrap CIs';b.disabled=false;
  sigTable();
 },20);
});
document.getElementById('sigtab').addEventListener('click',sigOpen);
document.getElementById('sigclose').addEventListener('click',sigClose);
document.addEventListener('keydown',ev=>{if(ev.key==='Escape'&&SIGTAB_OPEN)sigClose();});
const searchBox=document.getElementById('search');
function doSearch(q){q=(q||'').trim();const info=document.getElementById('info');if(!q)return;const hit=DATA.nodes.find(n=>n.label.toLowerCase()===q.toLowerCase())||DATA.nodes.find(n=>n.label.toLowerCase().indexOf(q.toLowerCase())===0);if(!hit){info.innerHTML='No gene matching "'+q+'"';return;}try{network.selectNodes([hit.id]);network.focus(hit.id,{scale:1.3,animation:true});info.innerHTML='<b>'+hit.label+'</b>';}catch(e){info.innerHTML='<b>'+hit.label+'</b> not in current view';}}
searchBox.addEventListener('keydown',ev=>{if(ev.key==='Enter')doSearch(searchBox.value);});
searchBox.addEventListener('change',()=>doSearch(searchBox.value));
function zoomBy(f){if(!network)return;const s=network.getScale();network.moveTo({scale:s*f,animation:{duration:200}});}
document.getElementById('zin').addEventListener('click',()=>zoomBy(1.25));
document.getElementById('zout').addEventListener('click',()=>zoomBy(0.8));
document.getElementById('zfit').addEventListener('click',()=>{if(network)network.fit({animation:true});});
// Two pairs of year handles -- one in the panel, one in the table view -- driving one range.
// Whichever you drag, the other follows, the ranking is recomputed for the new window, and the
// graph rebuilds; the table redraws only while it is open.
const yl=document.getElementById('yrlo'),yh=document.getElementById('yrhi');
const yl2=document.getElementById('yrlo2'),yh2=document.getElementById('yrhi2');
[yl,yh,yl2,yh2].forEach(el=>{if(el){el.min=MINY;el.max=MAXY;}});
yl.value=MINY;yh.value=MAXY;if(yl2){yl2.value=MINY;yh2.value=MAXY;}
function updYr(){const a=activeYears(),t=a[0]+'–'+a[1];
 ['yrlab','yrlab2'].forEach(id=>{const e=document.getElementById(id);if(e)e.textContent=t;});}
function yearChanged(fromTable){
 if(fromTable){yl.value=yl2.value;yh.value=yh2.value;}
 else if(yl2){yl2.value=yl.value;yh2.value=yh.value;}
 updYr();
 build(+thr.value);                 // build() recomputes the ranking and redraws the table
}
[yl,yh].forEach(el=>el.addEventListener('input',()=>yearChanged(false)));
[yl2,yh2].forEach(el=>{if(el)el.addEventListener('input',()=>yearChanged(true));});
updYr();   // the ranking follows from build(), which runs once at the end of this script
const chemGenes={};DATA.nodes.forEach(n=>(n.chems||[]).forEach(c=>{(chemGenes[c]=chemGenes[c]||[]).push(n.label);}));const csel=document.getElementById('chemfilter');Object.keys(chemGenes).sort().forEach(c=>{const g=chemGenes[c].slice().sort();const o=document.createElement('option');o.value=c;o.textContent=c+' → '+g.join(', ');csel.appendChild(o);});csel.addEventListener('change',()=>build(+thr.value));
// disease list for the co-mention picker; the gene-only graph has none, so the row hides itself
(function(){const sel=document.getElementById('comention');const ds=DATA.nodes.filter(n=>(n.kind||'gene')==='disease').sort((a,b)=>a.label.localeCompare(b.label));
 const trow=document.getElementById('tissuerow');
 if(!ds.length){const r=document.getElementById('cmrow');if(r)r.style.display='none';if(trow)trow.style.display='none';return;}
 // say up front how much of the disease list the name-based tissue reading actually covers
 const tg={};ds.forEach(n=>{const t=tissueOf(n.label);if(t)(tg[t]=tg[t]||[]).push(n.label);});
 const multi=Object.keys(tg).filter(t=>tg[t].length>1);
 const hit=multi.reduce((s,t)=>s+tg[t].length,0);
 const el=document.getElementById('tissuen');
 if(el)el.textContent=hit+' of '+ds.length+' disease nodes fall into '+multi.length+' tissues; the rest are left where the layout puts them';
 // also a position control: replay from BASEPOS instead of rebuilding the whole network
 document.getElementById('tissuestack').addEventListener('change',applyLayoutShape);
 ds.forEach(n=>{const o=document.createElement('option');o.value=n.id;o.textContent=n.label;sel.appendChild(o);});
 sel.addEventListener('change',()=>build(+thr.value));})();
const drugBox=document.getElementById('drugsearch');function findDrug(q){q=(q||'').trim().toLowerCase();if(!q)return;const info=document.getElementById('info');const opts=[...csel.options].filter(o=>o.value);const m=opts.find(o=>o.value.toLowerCase()===q)||opts.find(o=>o.value.toLowerCase().indexOf(q)===0)||opts.find(o=>o.value.toLowerCase().indexOf(q)>=0);if(m){csel.value=m.value;build(+thr.value);info.innerHTML='Drug filter: <b>'+esc(m.value)+'</b>';}else{info.innerHTML='No drug matching "'+esc(q)+'"';}}drugBox.addEventListener('keydown',ev=>{if(ev.key==='Enter')findDrug(drugBox.value);});drugBox.addEventListener('change',()=>findDrug(drugBox.value));
// The left panel is exactly as wide as its zoom row: measured at runtime rather than guessed in
// CSS, because the buttons' width depends on the font that actually resolved. Fit's right edge
// then lands on the content edge, so the gap to the border is the padding -- the same gap the
// "genes" select has on the left. Skipped on narrow screens, where the media query takes over.
// Measured UNCONSTRAINED, which is the whole trick. Measuring the row inside the width this
// function set last time is a feedback loop: let the buttons wrap onto a second line and Fit's
// right edge is measured on that second line, so the width collapses to roughly one button --
// which guarantees the wrap next time. It ratchets narrower on every call and never recovers,
// which is exactly what a large graph did, since its redraws fire enough resize events to keep
// calling this. Clearing the width first lets the row lay out at its natural size; the same-line
// and range checks refuse a nonsense measurement rather than acting on it.
const FIT_MIN=160, FIT_MAX=520;
function fitLeftPanel(){
 const lp=document.getElementById('lpanel'), z=document.getElementById('zoom');
 if(!lp||!z)return;
 if(window.innerWidth<=700||lp.style.display==='none'){lp.style.width='';lp.style.maxWidth='';return;}
 const b=z.getElementsByTagName?z.getElementsByTagName('button'):[];
 if(!b.length||!b[0].getBoundingClientRect)return;   // nothing measurable: keep the CSS width
 const prevW=lp.style.width, prevM=lp.style.maxWidth;
 lp.style.maxWidth='none';lp.style.width='auto';     // shrink-to-fit: the row cannot wrap here
 const first=b[0].getBoundingClientRect(), last=b[b.length-1].getBoundingClientRect();
 const w=Math.ceil(last.right-first.left);
 if(Math.abs(last.top-first.top)<2&&w>=FIT_MIN&&w<=FIT_MAX){lp.style.width=w+'px';}
 else {lp.style.width=prevW;lp.style.maxWidth=prevM;}
}
fitLeftPanel();
window.addEventListener('resize',()=>{fitLeftPanel();layoutTable();if(network)network.redraw();});
document.querySelectorAll('.kindf').forEach(c=>c.addEventListener('change',()=>build(+thr.value)));
// gene-only runs hide every row the left column holds; an empty white box is worse than none
(function(){const lp=document.getElementById('lpanel');
 const kids=lp?Array.prototype.slice.call(lp.children||[]):[];   // HTMLCollection, not an array
 if(kids.length&&!kids.some(c=>!c.style||c.style.display!=='none'))lp.style.display='none';})();
buildCatFilters();buildSrcButtons();build(1);
</script></body></html>"""


KIND_ROW = (' <div class="row">Node type <button class="ihelp" data-help="kind" '
            'aria-label="About node types" aria-expanded="false">i</button></div>\n'
            ' <div class="row legend" id="kindfilters">'
            '<label><input type=checkbox class=kindf value="gene" checked> '
            '<span class="sw" style="background:#cfe3ff;border:1px solid #2b6cb0;'
            'border-radius:50%"></span> gene</label> '
            '<label><input type=checkbox class=kindf value="disease" checked> '
            '<span class="sw" style="background:#ffe0e0;border:1px solid #b3243b;'
            'transform:rotate(45deg)"></span> disease</label> '
            '<label><input type=checkbox class=kindf value="chemical" checked> '
            '<span class="sw" style="background:#ece0f8;border:1px solid #6b3fa0"></span> chemical</label>'
            '</div>\n'
            ' <div class="row mut help" data-help="kind">Shapes: gene &#9679; &middot; disease &#9670; &middot; chemical &#9632;. '
            'Only <b>gene</b> names are drawn on the canvas &mdash; disease and chemical names are long '
            'and would bury the symbols, so <b>hover</b> (or click) those nodes to read them. '
            'An edge is shown only when BOTH its endpoint types are ticked &mdash; and where the '
            'ticked types have no relations at all between them (nothing in the corpus relates a '
            '<b>disease</b> to a disease), the nodes are drawn on their own instead of leaving a '
            'blank canvas: every one that still has evidence through a partner you unticked, '
            'sized by that evidence and pooled by tissue.</div>')


def render_graph(payload, lib, miny, maxy, nxml=None, pubmed_query="", background=None):
    # nxml is kept for the caller's signature; the panel no longer prints a corpus-coverage line
    if lib:
        libtag = "<script>\n" + lib.replace("</script>", "<\\/script>") + "\n</script>"
    else:
        libtag = f'<script src="{VIS_URL}"></script>'
    # heads the left column: the PubMed query this corpus came from (empty -> nothing shown).
    # The query itself can run to several lines of boolean, so it sits behind the section's "i"
    # like every other explanation in this column -- the heading alone says where the data is from.
    qrow = (f'<div class="row" id="pubmedq">PubMed query '
            f'<button class="ihelp" aria-label="Show the PubMed query" aria-expanded="false">i</button>'
            f'<div class="mut help" style="word-break:break-word">{html.escape(pubmed_query)}</div></div>'
            if pubmed_query else "")
    title = (html.escape(pubmed_query) if pubmed_query
             else "High-confidence gene / disease / chemical relations")
    # these edges sit mostly in 0.5-0.8 (BioRED gene-disease relations), so opening the
    # slider at the old gene-only default of >=0.99 would show an almost empty canvas
    return (GRAPH_TEMPLATE.replace("__LIBTAG__", libtag)
            .replace("__TITLE__", title)
            .replace("__KINDROW__", KIND_ROW)
            .replace("__CONFDEF__", "6")
            .replace("__PUBMED_QUERY__", qrow)
            .replace("__PAYLOAD__", json.dumps(payload, ensure_ascii=False))
            .replace("__BACKGROUND__", json.dumps(background or {"n": 0, "ne": 0, "doc": {}, "par": {}}))
            .replace("__CCOLOR__", json.dumps(RCOLOR))
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
    ap.add_argument("--merge-scale", choices=MERGE_SCALES, default="rank",
                    help="which model's score represents a pair BOTH models kept: rank "
                         "(default) compares each score against its own model's distribution, "
                         "so a saturating calibrator cannot win by construction; calibrated "
                         "takes a plain max() of the calibrated probabilities (pre-2026-08 "
                         "behaviour). The kept value is the winner's calibrated probability "
                         "either way -- only the choice of winner changes")
    ap.add_argument("--typed-model", default=None,
                    help="checkpoint name supplying the typed/signed labels (default: the one "
                         "whose name contains 'biored')")
    ap.add_argument("--gate-model", default=None,
                    help="checkpoint name used as the existence gate (default: the other one)")
    args = ap.parse_args()

    set_data_root(args.data_root)
    flags = load_type_flags()
    def keep_fn(t, T):
        return qualifies_multi(t, T, flags)
    if not RE_FILE.exists():
        raise SystemExit(f"ERROR: {RE_FILE} not found under data root {DATA_ROOT} "
                         f"(run the gpu_bundle pipeline first, or pass --data-root).")

    raw = json.loads(RE_FILE.read_text(encoding="utf-8"))
    d, merge_note = merge_models(raw, args.merge, args.typed_model, args.gate_model, args.merge_scale)
    print(merge_note)
    # Ceilings are read from the PRE-merge triples: after the merge a pair keeps the typed
    # model's name but may carry the other model's score, so post-merge attribution is wrong.
    # Every absolute cutoff this run applies, including the in-browser slider's stops.
    for line in ceiling_report(raw, [args.score] + [float(x) for x in args.thresholds.split(",")]
                               + ([GRAPH_BASE, 0.95, 0.99] if not args.no_graph else [])):
        print(line)
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
        # one denominator for every lung_* analysis: this run's slice, merged with the
        # siblings' cached slices (see shared_background)
        bg = shared_background(DATA_ROOT, corpus_contrib(d, flags), RE_FILE)
        print(f"  corpus for enrichment: {bg['n']:,} publications, {bg['ne']:,} entities "
              f"from {', '.join(bg['src'])}"
              + (f"; not yet run: {', '.join(bg['missing'])}" if bg['missing'] else ""))
        universe = kept if args.score <= GRAPH_BASE else [t for t in d if keep_fn(t, GRAPH_BASE)]
        payload = graph_payload_multi(universe, flags)
        yrs = [s["yr"] for e in payload["edges"] for s in e["sents"] if s.get("yr")]
        miny, maxy = (min(yrs), max(yrs)) if yrs else (2000, 2026)
        lib = get_vis_lib()
        # corpus size = # source documents; experimental_ner/ is empty in the bundle
        # (it was a runtime symlink), so fall back to the per-document sentences/ files.
        nxml = len(list(XML_DIR.glob("*.xml"))) or len(list(SENT_DIR.glob("*.json")))
        GRAPH_OUT.parent.mkdir(parents=True, exist_ok=True)
        pubmed_query = read_pubmed_query()
        GRAPH_OUT.write_text(render_graph(payload, lib, miny, maxy, nxml, pubmed_query, bg),
                             encoding="utf-8")
        n99 = sum(1 for e in payload["edges"] if any(s["sc"] >= 0.99 for s in e["sents"]))
        npubs = len({s["pmid"] for e in payload["edges"] for s in e["sents"]})
        kinds = collections.Counter(n.get("kind", "gene") for n in payload["nodes"])
        what = ", ".join(f"{c:,} {k}" for k, c in kinds.most_common())
        print(f"wrote graph -> {GRAPH_OUT}  ({what}; {len(payload['edges']):,} edges; "
              f"{n99:,} edges have a >=0.99 sentence; {npubs:,} of {nxml:,} input XMLs produced triples)")
        cats = collections.Counter(e["cat"] for e in payload["edges"])
        signed = sum(c for k, c in cats.items()
                     if (k[4:] if k.startswith("not ") else k) in ("activates", "inhibits"))
        print(f"  edge relations: {', '.join(f'{k} {c:,}' for k, c in cats.most_common())}"
              f"  ({signed:,} signed)")
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

        # also drop a copy in the root dir, named "<current directory>_YYYY_MM_DD_M.html"
        # (illegal filename characters underscored so the name is always valid).
        dirname = re.sub(r'[\\/:*?"<>|\s]+', "_", ROOT.name)
        today = datetime.date.today().strftime("%Y_%m_%d")
        dest = ROOT / f"{dirname}_{today}_M.html"
        shutil.copy2(GRAPH_OUT, dest)
        print(f"copied graph -> {dest}")

    for r in rows:
        print(f"  score>={r['T']}: {r['triples']:,} triples / {r['sentences']:,} sent; "
              f"G {r['filt']:,} triples / {r['filt_sentences']:,} sent"
              f"; M {sum(1 for t in d if keep_fn(t, r['T'])):,} triples")


if __name__ == "__main__":
    main()
