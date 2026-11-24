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
import csv
import datetime
import html
import json
import re
import shutil
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "gpu_bundle"))
import direction as D  # noqa: E402  (shared with relation_extraction.py)

ROOT = Path(__file__).resolve().parent
# The pipeline outputs (gpu.py's writable run dir) live under kaggle_working/ by default;
# every input below is resolved beneath this data root. Override with --data-root.
DATA_ROOT = ROOT / "kaggle_working"

# module-level paths; (re)bound to DATA_ROOT by set_data_root() so --data-root can retarget them
OUT_DIR = XML_DIR = SENT_DIR = None
RE_FILE = PMC_YEARS = TARGET_FILE = CANCER_GENES = None
DISEASE_LIB = CHEM_LIB = None
JSON_OUT = GRAPH_OUT = None


def set_data_root(data_root):
    """Point every input/output path at `data_root` (the pipeline's output tree)."""
    global DATA_ROOT, OUT_DIR, XML_DIR, SENT_DIR, RE_FILE, PMC_YEARS, CANCER_GENES
    global TARGET_FILE, DISEASE_LIB, CHEM_LIB, JSON_OUT, GRAPH_OUT
    DATA_ROOT = Path(data_root).resolve()
    OUT_DIR = DATA_ROOT / "TRIPLES"
    XML_DIR = DATA_ROOT / "experimental_ner"   # input XML corpus (may be empty in the bundle)
    SENT_DIR = DATA_ROOT / "sentences"         # one JSON per source document (corpus-size fallback)
    RE_FILE = OUT_DIR / "triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json"
    PMC_YEARS = DATA_ROOT / "databases" / "pmc_years.json"
    # OncoKB Cancer Gene List (https://www.oncokb.org/api/v1/utils/cancerGeneList.txt), the
    # gene -> oncogene / tumor suppressor roles behind the Oncogene and Tumor suppressor buttons
    CANCER_GENES = DATA_ROOT / "databases" / "oncokb_cancerGeneList.tsv"
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
#   sd        a SURFACE: "SD" is stable disease (a RECIST response category, next to PR/CR),
#             which MONDO's exact synonym maps to "Salla disease", e.g. "a MEK2 mutation ...
#             in one patient with SD" (PMC3296083). Matching the surface leaves a genuine
#             "Salla disease" mention, should one ever appear, untouched.
DISEASE_IGNORE = {"neoplasm", "cancer", "os", "sd"}

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
#   mm   disease.py maps it to "Miyoshi myopathy" (413 mentions); in this cancer corpus MM is
#        multiple myeloma. DISEASE_MERGE below folds MONDO's "plasma cell myeloma" into the same
#        node, so "MM", "myeloma" and "multiple myeloma" all land on one.
#
# The rest close two holes that hid the commonest cancers (added 2026-11-24): abbreviations the
# normalizer has no id for ("PDAC" x734, "LUAD" x906), and plain names it resolves to SEVERAL
# MONDO terms at once -- "pancreatic cancer" to exocrine pancreatic carcinoma + malignant
# pancreatic neoplasm, "CRC" to colorectal cancer + colorectal carcinoma -- which single() then
# throws away. Before this, 5% of pancreatic mentions reached the graph and KRAS-PDAC (KRAS G12
# being the disease's defining driver) was all but absent.
# Only surfaces that mean one thing in a cancer corpus are listed. Left out on purpose: ALL
# (the word "all"), BC (breast / bladder), GC (gastric / germinal centre), PCa (also principal
# component analysis), PD, EC, OC, UC, um (micrometre), fl (floxed, fl/fl), MF, ACC, RMS, PDA
# (patent ductus arteriosus).
DISEASE_SURFACE = {
    "lcc": "lung large cell carcinoma", "mm": "multiple myeloma",
    # pancreas
    "pdac": "pancreatic ductal adenocarcinoma", "pdacs": "pancreatic ductal adenocarcinoma",
    "pdac tumors": "pancreatic ductal adenocarcinoma", "mpdac": "pancreatic ductal adenocarcinoma",
    "pancreatic cancer": "pancreatic cancer", "pancreatic cancers": "pancreatic cancer",
    "panin": "pancreatic intraepithelial neoplasia",
    # colon / rectum
    "crc": "colorectal cancer", "crcs": "colorectal cancer", "mcrc": "colorectal cancer",
    "colon cancer": "colon carcinoma", "colon cancers": "colon carcinoma",
    # lung
    "luad": "lung adenocarcinoma", "lusc": "squamous cell lung carcinoma",
    # breast
    "breast tumor": "breast cancer", "breast tumors": "breast cancer", "mbc": "breast cancer",
    "tnbcs": "triple-negative breast carcinoma",
    # others
    "bladder cancer": "urinary bladder carcinoma",
    "rcc": "renal cell carcinoma", "prcc": "papillary renal cell carcinoma",
    "hnsc": "head and neck squamous cell carcinoma", "oscc": "oral squamous cell carcinoma",
    "scc": "squamous cell carcinoma", "sccs": "squamous cell carcinoma",
    "crpc": "castration-resistant prostate carcinoma", "hgsoc": "high grade serous ovarian carcinoma",
    "mpm": "malignant pleural mesothelioma", "mtc": "medullary thyroid gland carcinoma",
    "ptc": "papillary thyroid carcinoma", "npc": "nasopharyngeal carcinoma",
    "cca": "cholangiocarcinoma", "gbc": "gallbladder cancer",
    "t - all": "T-cell acute lymphoblastic leukemia",
}

# Gene surfaces, the same two holes as for diseases (added 2026-11-24). The normalizer often
# resolves a surface to SEVERAL HGNC symbols -- "AR" to AR + FDXR + AKR1B1 + AREG, because the
# approved symbol is also an old alias of three other genes -- and single() then dropped it: all
# 366 plain "AR" triples were lost, leaving the AR node only glued forms (ARQ62L) and an alias
# mistake. 17,937 gene mentions went this way. Two fixes, in gene_symbol():
#  1. EXACT SYMBOL: when the surface IS one of the candidate symbols (case aside), take it --
#     "AR", "HGF", "RAC1", "BAP1", "ATR", "CDH1" (~4,900 mentions, ~200 genes). Not for the
#     surfaces in GENE_EXACT_SKIP, whose text is mostly a word or another thing: TF (transcription
#     factor), cAMP, si(RNA), me2/me3 (histone marks), mOS (median OS), SFN (sulforaphane), SAG
#     (Smoothened agonist), C6 (glioma line), TTF1 (= NKX2-1, below), Pgp (= ABCB1), ...
#  2. GENE_SURFACE: names that mean one gene in a cancer corpus. Families are deliberately left
#     out -- ERK (MAPK3 / MAPK1), NF-kB (five subunits), FGFR (four) -- rather than pick a member.
# GENE_IGNORE drops surfaces whose single mapping is wrong here: "AIS" is an AR alias (androgen
# insensitivity syndrome) but in this corpus adenocarcinoma in situ.
GENE_EXACT_SKIP = {"TF", "CAMP", "SI", "ME2", "ME3", "MOS", "GAL", "PC", "SAG", "C6", "C7", "CPD",
                   "ARC", "SMS", "HP", "GPT", "LAT", "PGP", "TPO", "SFN", "TTF1", "SRM"}
GENE_SURFACE = {
    "ar-v7": "AR", "ar-fl": "AR", "arv7": "AR",
    "p21": "CDKN1A", "p21waf1": "CDKN1A", "p21cip1": "CDKN1A", "p21waf1/cip1": "CDKN1A",
    "p16": "CDKN2A", "p16ink4a": "CDKN2A",
    "p27": "CDKN1B", "p27kip1": "CDKN1B",
    "nrf2": "NFE2L2",
    "p38": "MAPK14", "p-p38": "MAPK14", "p38 mapk": "MAPK14", "p38mapk": "MAPK14", "phospho-p38": "MAPK14",
    "p62": "SQSTM1",
    "er": "ESR1",
    "taz": "WWTR1",
    "s6": "RPS6", "ps6": "RPS6", "p-s6": "RPS6", "phospho-s6": "RPS6",
    "ttf1": "NKX2-1", "ttf-1": "NKX2-1",
}
GENE_IGNORE = {"ais"}


def gene_symbol(e, txt):
    """HGNC symbol for a GENETIC endpoint (txt = its lowercased surface), or None."""
    if txt in GENE_IGNORE:
        return None
    if txt in GENE_SURFACE:
        return GENE_SURFACE[txt]
    v = e.get("hgnc_symbol")
    if isinstance(v, list):
        up = txt.upper()
        return up if up in v and up not in GENE_EXACT_SKIP else None
    return single(v)


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
    "plasma cell myeloma": "multiple myeloma",   # MONDO's label; MM is overridden to the common name
    "pancreatic neoplasm": "pancreatic cancer",     # "pancreatic tumor(s)", as lung neoplasm above
    "colorectal carcinoma": "colorectal cancer",
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
        # the gate (PPI) model's OWN reading of a pair both models kept. The typed label above
        # is BioRED's; PPI only ever says "interacts", so "kept by both" means both found a
        # relation here, not that both said "binds" -- the graph shows each model's reading.
        if t["corroborated"] and g is not None and ty is not None:
            t["gate_cat"] = rel_cat(g)
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


# '-i' inhibitor shorthand (BRAFi, MEKi) that triples.py recast as a drug class: one node per
# class, case-insensitively ("eIF4F inhibitor" / "eiF4F inhibitor"), labelled as first seen
INHIB_LABEL = {}


def node_of(e, flags):
    """(kind, normalized_id) for one endpoint, or None when it cannot be a node.

    Dropped: endpoints with no normalized id (an un-normalized surface is not an entity we
    can merge across documents), lab controls / MKI67 (genes), process-or-phenotype
    surfaces (diseases), gene-or-process surfaces mislabelled CHEMICAL, and over-generic
    disease labels."""
    t, txt = e.get("type"), (e.get("text") or "").strip().lower()
    ic = e.get("inhibitor_class")
    if ic:          # "BRAFi" is a BRAF INHIBITOR, a drug class -- not the BRAF gene
        return ("chemical", INHIB_LABEL.setdefault(ic.lower(), ic))
    if t == "GENETIC":
        s = gene_symbol(e, txt)
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
    merged, inhib = collections.Counter(), collections.Counter()
    for t in triples:
        for side in ("subject", "object"):
            e = t[side]
            typ, txt = e.get("type"), (e.get("text") or "").strip()
            key = txt.lower()
            if e.get("inhibitor_class"):
                inhib[e["inhibitor_class"]] += 1
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
                drops["gene: no HGNC symbol" if not gene_symbol(e, key)
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
    if inhib:
        out.append(f"  '-i' inhibitor shorthand -> drug-class (chemical) nodes: {sum(inhib.values()):,} "
                   f"mention(s), " + ", ".join(f"{k} x{c}" for k, c in inhib.most_common(top)))
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
    genes = [e for e in (t["subject"], t["object"])
             if e.get("type") == "GENETIC" and not e.get("inhibitor_class")]
    if genes and not any(e.get("control") == "no" for e in genes):
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


# ----- false drug-target pairs (manual curation) --------------------------------------
# DGIdb ingests clinical-context sources -- MyCancerGenome(ClinicalTrial), TALC,
# ClearityFoundationClinicalTrial, CancerCommons -- whose claim is "this drug is used where
# this gene is altered", and records them with interaction_type "inhibitor". The type gate in
# load_dgidb therefore cannot see them, and they arrive here indistinguishable from a binding
# assertion: rituximab (anti-CD20) "targeting" EGFR, gefitinib targeting ERBB2, regorafenib
# and selumetinib targeting KRAS, panobinostat (a class I/II HDAC inhibitor) targeting the
# class III sirtuins. On the reference run 78 of 271 mappings behind coloured nodes rest only
# on such sources, and seven nodes were GREEN -- claiming an approved anti-neoplastic -- on
# that basis alone.
#
# Only the ones a drug's known pharmacology CONTRADICTS are listed. Judgement calls stay in:
# JAK2 <- ruxolitinib and MTOR <- everolimus/sirolimus are clinical-sourced and true, and
# attributions to a partner or subunit rather than the direct target (EML4 <- crizotinib, BCR <- dasatinib, PIK3R1 <-
# alpelisib, RICTOR <- dactolisib, NOTCH1 <- RO4929097) are defensible readings of "target",
# so they are deliberately NOT dropped. The identifier-namespace artifacts are a different
# problem with a different fix: target_pharm.py drops those rows at the source.
# Matched case-insensitively against the ChEBI label carried in chemical_to_target.json.
FALSE_TARGETS = {
    "EGFR": {"rituximab", "pertuzumab", "tucatinib"},      # CD20 and HER2 drugs
    "ERBB2": {"gefitinib", "rituximab"},                   # EGFR-selective / anti-CD20
    "LAG3": {"nivolumab"},                                 # nivolumab is anti-PD-1; relatlimab is the LAG-3 arm
    "KRAS": {"regorafenib", "selumetinib"},                # KRAS is the biomarker; MEK is the target
    "BRAF": {"dasatinib (anhydrous)", "dasatinib monohydrate", "dasatinib(1+)"},
    "FLT3": {"bortezomib", "clofarabine"},                 # proteasome / nucleoside analogue
    "ALK": {"ganetespib"},                                 # HSP90
    "AKT1": {"everolimus", "azd4547"},                     # mTOR / FGFR
    "SIRT1": {"panobinostat"}, "SIRT3": {"panobinostat"},  # class I/II HDAC drug vs
    "SIRT6": {"panobinostat"}, "SIRT5": {"panobinostat"},  # class III NAD+ sirtuins
    "SIRT7": {"panobinostat"},
    "PDGFB": {"axitinib", "sunitinib"},                    # they inhibit the RECEPTOR, not the ligand
    "PTCH1": {"vismodegib"},                               # vismodegib binds SMO
    "CFLAR": {"dovitinib"},                                # FGFR/VEGFR TKI
    "MTOR": {"bgj-398"}, "PIK3CA": {"bgj-398"},            # infigratinib is FGFR-selective
    "PIK3CB": {"bgj-398"}, "PIK3R1": {"bgj-398"},
}
N_FALSE_DROPPED = 0


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
        block = FALSE_TARGETS.get(g)
        if block:
            keep = [c for c in chems if (c.get("chebi_label") or "").casefold() not in block]
            globals()["N_FALSE_DROPPED"] += len(chems) - len(keep)
            chems = keep
            if not chems:                       # every claim on this gene was a false one
                continue
        tgt[g] = len(chems)                     # recount: n_chemicals is the pre-filter figure
        chems_by_gene[g] = sorted({c.get("chebi_label") for c in chems if c.get("chebi_label")})
        vrs = [str(r) for c in chems for r in (c.get("via_roles") or [])]
        has_db = any(r.startswith("DGIdb") for r in vrs)
        has_ch = any(not r.startswith("DGIdb") for r in vrs)
        tsrc[g] = "dgidb" if (has_db and not has_ch) else ("chebi+dgidb" if has_db else "chebi")
        # Colour class, read off the DGIdb drug class carried in via_roles (the ChEBI roles
        # carry no approval status at all): green = approved & anti-neoplastic, amber = approved
        # but not anti-neoplastic, pink = everything else -- a real target with NO approved drug
        # behind it, i.e. DGIdb-investigational and/or a ChEBI `has role` target role. Pink is
        # therefore not a source label: a pink gene may well be in DGIdb (check tsource), it just
        # has nothing approved pointing at it.
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
    gene_role = _cancer_gene_roles()     # not `roles`: the loop below reuses that name
    dir_sent = collections.defaultdict(set)            # (a,b,cat) -> sentences
    pair_sent = collections.defaultdict(dict)          # pair -> {sentence: [score, pmid, cat, spec, src]}
    pair_src = collections.defaultdict(set)            # pair -> model roles behind ANY of its triples
    node_sent = collections.defaultdict(dict)          # node -> {sentence: maxscore}
    # gene node -> {surface as written: sentences}: what the text actually said for the symbol
    # ("p21" for CDKN1A, "androgen receptor" for AR), shown on the node and marked in sentences
    node_surf = collections.defaultdict(lambda: collections.defaultdict(set))
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
        for nd, side in ((a, "subject"), (b, "object")):
            if nd[0] == "gene":
                node_surf[nd[1]][(t[side].get("text") or "").strip()].add(sent)
        dir_sent[(ka, kb, cat)].add(sent)
        pair_src[frozenset((ka, kb))] |= roles
        cur = pair_sent[frozenset((ka, kb))].get(sent)
        if cur is None or sc > cur[0]:
            pair_sent[frozenset((ka, kb))][sent] = [sc, pm, cat, spec, src_tag(roles), t.get("gate_cat")]
        for nd in (ka, kb):
            if node_sent[nd].get(sent, -1) < sc:
                node_sent[nd][sent] = sc
    # group the directed (a, b, cat) records by their unordered pair ONCE: scanning all of
    # dir_sent for every pair was quadratic (~125k pairs x ~200k records = hours on the
    # full corpus). Each pair gets the same candidates as before, and the sort key
    # (count, from, to, cat) is unique, so the chosen edge is identical.
    by_pair = collections.defaultdict(list)
    for (f, to, cat), ss in dir_sent.items():
        by_pair[frozenset((f, to))].append((len(ss), f, to, cat))
    edges = []
    for pr, sd in pair_sent.items():
        cands = by_pair[pr]
        cands.sort(reverse=True)
        _, ff, ft, fcat = cands[0]
        sents = [{"pmid": pm, "text": sent[:300], "sc": round(sc, 4), "yr": years.get(pm),
                  "rc": cat, "sp": sp, "sr": sr, **({"pc": pc} if sr == "both" and pc else {})}
                 for sent, (sc, pm, cat, sp, sr, pc) in sd.items()]
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
                      "tcat": tcat.get(nd, "other"),
                      **({"role": gene_role[nd]} if k == "gene" and nd in gene_role else {}),
                      # surfaces by sentence count, most frequent first (at most 12)
                      **({"surf": [[w, len(ss)] for w, ss in sorted(node_surf[nd].items(),
                                    key=lambda z: (-len(z[1]), z[0]))[:12] if w]}
                         if k == "gene" and nd in node_surf else {})})
    nr = collections.Counter(n.get("role") for n in nodes if n.get("role"))
    print(f"  OncoKB cancer-gene roles on graph genes: oncogene {nr['onc']}, tumor suppressor "
          f"{nr['tsg']}, both {nr['both']}")
    return {"nodes": nodes, "edges": edges}


def _cancer_gene_roles():
    """{HGNC symbol: 'onc' | 'tsg' | 'both'} from the OncoKB Cancer Gene List's "Gene Type"
    column (ONCOGENE / TSG / ONCOGENE_AND_TSG; NEITHER, INSUFFICIENT_EVIDENCE and blank carry no
    role). Missing file -> {} with a warning: the two buttons then keep nothing."""
    role = {"ONCOGENE": "onc", "TSG": "tsg", "ONCOGENE_AND_TSG": "both"}
    try:
        with open(CANCER_GENES, encoding="utf-8", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="	"))
    except OSError:
        print(f"  WARNING: {CANCER_GENES} not found -- Oncogene / Tumor suppressor buttons will be empty")
        return {}
    out = {}
    for r in rows:
        g, t = (r.get("Hugo Symbol") or "").strip(), (r.get("Gene Type") or "").strip().upper()
        if g and t in role:
            out[g] = role[t]
    return out


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
 #lpanel{left:12px;scrollbar-gutter:stable}
 /* an 'on' Enrich button turns bold: reserve the bold width up front so the grid never widens */
 .egrid button::after{content:attr(data-t);display:block;height:0;overflow:hidden;visibility:hidden;font-weight:600}
 /* each column can be shut to free the canvas: the x sits in the panel's top corner, and a
    closed panel leaves a small tab at the same corner that brings it back */
 .pclose{float:right;margin:-8px -10px 4px 8px;width:22px;height:22px;padding:0;border:1px solid #cdd5e0;border-radius:50%;background:#eef2f7;color:#1c2330;font:600 15px/1 Segoe UI,Arial,sans-serif;cursor:pointer}
 .pclose:hover,.popen:hover{background:#dde4ee}
 .popen{position:absolute;top:12px;z-index:6;display:none;width:34px;height:34px;padding:0;border:1px solid #cdd5e0;border-radius:8px;background:rgba(255,255,255,.97);color:#1c2330;font:18px/1 Segoe UI,Arial,sans-serif;cursor:pointer;box-shadow:0 2px 12px rgba(0,0,0,.18)}
 #lopen{left:12px} #ropen{right:12px}
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
 /* the Enrich buttons, two per row: each .row dissolves into the grid (display:contents) so its
    button and "i" are two cells, and its help note spans the full width under the pair when
    opened ("dense" lets the next button back-fill the row the note would otherwise break) */
 .egrid{display:grid;grid-template-columns:max-content max-content max-content max-content;column-gap:0;row-gap:4px;align-items:center;grid-auto-flow:row dense;margin:2px 0 8px;border:1px solid #cdd5e0;border-radius:6px;padding:4px 6px}   /* boxed like #catfilters */
 .egrid>.row{display:contents}
 .egrid>.row>button:first-child{margin:0}
 .egrid .ihelp{margin:0 10px 0 3px}
 .egrid .help{grid-column:1/-1;margin:0 0 4px;width:0;min-width:100%}   /* adds no width of its own: wraps to the grid */
 .help{display:none;margin-top:4px} .help.open{display:block}
 #conf{width:190px;cursor:pointer}
 select,#search,#genefilter,#drugsearch,#textfilter{max-width:100%;background:#fff;border:1px solid #cdd5e0;color:#1c2330;border-radius:5px;padding:3px 6px;font-size:13px}
 #search,#genefilter,#drugsearch,#textfilter{width:200px}
 mark{background:#ffe680;color:inherit;border-radius:2px;padding:0 1px}
 .gm{font-weight:700;color:#1d4e89;border-bottom:1px solid #8fb3df}
 #catfilters label{display:block;cursor:pointer;white-space:nowrap;font-size:12px;margin:1px 0}
 /* no inner scroller: the relation list is a dozen rows at most, and a box that scrolls inside
    a panel that also scrolls hides ticked types from anyone who does not think to scroll it */
 #catfilters{border:1px solid #cdd5e0;border-radius:6px;padding:4px 6px}
 #catfilters .cnt{color:#5b6677;font-size:11px}
 #zoom button,#srcbtns button,#labelbtns button,#catall,#drvbtn,#fusbtn,#itdbtn,#oexbtn,#ampbtn,#delbtn,#trubtn,#exobtn,#cnbtn,#metbtn,#resbtn,#acqbtn,#slbtn,#bmkbtn,#prgbtn,#clnbtn,#vivbtn,#oncbtn,#tsgbtn,#immbtn,#mtsbtn{background:#eef2f7;color:#1c2330;border:1px solid #cdd5e0;border-radius:6px;padding:4px 10px;cursor:pointer;margin-right:6px;font-size:13px}
 #zoom button:hover,#srcbtns button:hover,#labelbtns button:hover,#catall:hover,#drvbtn:hover,#fusbtn:hover,#itdbtn:hover,#oexbtn:hover,#ampbtn:hover,#delbtn:hover,#trubtn:hover,#exobtn:hover,#cnbtn:hover,#metbtn:hover,#resbtn:hover,#acqbtn:hover,#slbtn:hover,#bmkbtn:hover,#prgbtn:hover,#clnbtn:hover,#vivbtn:hover,#oncbtn:hover,#tsgbtn:hover,#immbtn:hover,#mtsbtn:hover{background:#dde4ee}
 #srcbtns button,#labelbtns button{margin-bottom:4px}
 #disfilters{display:block;max-height:148px;overflow-y:auto;margin-top:0}
 .dbox{border:1px solid #cdd5e0;border-radius:6px;padding:4px 6px;margin:2px 0 8px}   /* the Enrich boxes' frame (.egrid) */
 #disfilters label{display:block;margin:1px 0}
 #disfilters .dist{font-weight:600;margin-top:3px}
 #disfilters .dism{padding-left:18px}
 #disclear:hover{color:#0969da}
 /* the zoom row sets the left panel's width (see fitLeftPanel), so the last button must not
    carry a trailing margin -- 6px of it would push Fit onto a second line */
 #zoom button:last-child{margin-right:0}
 #srcbtns button.on,#labelbtns button.on,#drvbtn.on,#fusbtn.on,#itdbtn.on,#oexbtn.on,#ampbtn.on,#delbtn.on,#trubtn.on,#exobtn.on,#cnbtn.on,#metbtn.on,#resbtn.on,#acqbtn.on,#slbtn.on,#bmkbtn.on,#prgbtn.on,#clnbtn.on,#vivbtn.on,#oncbtn.on,#tsgbtn.on,#immbtn.on,#mtsbtn.on{background:#0969da;border-color:#0969da;color:#fff;font-weight:600}
 /* the name buttons must not move when one is toggled: bold text is wider, so each button reserves
    its bold width up front (a hidden bold copy of its own label) and switching weight changes nothing */
 /* ...and 3px less side padding (no trailing margin) pays for that reserve, so the pair still
    fits on one line in the left column, which is sized to the zoom row */
 #labelbtns button{display:inline-flex;flex-direction:column;align-items:center;vertical-align:top;padding:4px 7px}
 #labelbtns button:last-child{margin-right:0}
 #labelbtns button::after{content:attr(data-label);font-weight:600;height:0;overflow:hidden;visibility:hidden;pointer-events:none}
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
  #ropen{top:auto;bottom:12px}
  .vis-tooltip{max-width:88vw!important}
 }
</style></head><body>
<button id="lopen" class="popen" title="Show the left menu" aria-label="Show the left menu">&#9776;</button>
<button id="ropen" class="popen" title="Show the right menu" aria-label="Show the right menu">&#9776;</button>
<div id="lpanel">
 <button class="pclose" data-panel="lpanel" title="Close this menu" aria-label="Close the left menu">&times;</button>
 __PUBMED_QUERY__
 <div class="row">Year: <b id="yrlab"></b><br><input id="yrlo" type="range" style="width:74px"> <input id="yrhi" type="range" style="width:74px"></div>
 <div class="row" id="lblrow">Draw names for: <button class="ihelp" aria-label="About drawn names" aria-expanded="false">i</button><br><span id="labelbtns"><button class="lblb" data-kind="disease" data-label="Disease names">Disease names</button><button class="lblb" data-kind="chemical" data-label="Drug names">Drug names</button></span>
  <div class="mut help">Gene symbols are always drawn. Disease and drug names are long, repeat across many edges and bury the symbols by sheer wordiness, so they start off &mdash; switch them on to read a neighbourhood, off to see its shape. When on they are set in a <b>condensed</b> face at <b>80%</b> of the size the same weight of gene would take, and drug names are <i>italic</i> as well: the extra cues let a long name sit beside a symbol without out-shouting it, and let you tell the three kinds apart in one glance without reading them. Relative size still tracks evidence within each kind, and the font slider scales everything together. Applies in place; the layout is not recomputed.</div></div>
 <div class="row">Shrink periphery: <b id="shv">90%</b> <button class="ihelp" aria-label="About shrink periphery" aria-expanded="false">i</button><br><input id="shrink" type="range" min="0" max="90" step="5" value="90" aria-label="Shrink periphery">
  <div class="mut help">Pulls everything past the middle distance inward, so the small clusters physics flung to the edges come back where you can read them without zooming out. The centre is left alone.</div></div>
 <div class="row">Expand center: <b id="exv">100%</b> <button class="ihelp" aria-label="About expand center" aria-expanded="false">i</button><br><input id="expand" type="range" min="0" max="100" step="5" value="100" aria-label="Expand center">
  <div class="mut help">Blows the crowded core outward and carries the rest along, thinning the hairball without moving anything past its neighbours. Both reshape the drawn layout only &mdash; radially, from the stabilized positions, so nothing changes order and returning a slider to 0% restores the layout exactly.</div></div>
 <div class="row" id="tissuerow"><label><input type=checkbox id="tissuestack" checked> Stack same-tissue diseases</label> <button class="ihelp" aria-label="About tissue stacking" aria-expanded="false">i</button>
  <div class="mut help">Drops the disease nodes naming one tissue onto a single spot, overlapping, so <em>lung cancer</em>, <em>lung adenocarcinoma</em> and <em>non-small cell lung carcinoma</em> read as one place on the canvas instead of three. They stay separate nodes with their own edges and tooltips &mdash; only their positions are pooled, after the layout settles. Tissue is read from the name (<span id="tissuen"></span>).</div></div>
 <div class="row" id="disrow"><b style="color:#1c2330">Keep diseases:</b> <button class="ihelp" aria-label="About the disease keep-list" aria-expanded="false">i</button>
  <span id="disclear" class="mut" style="cursor:pointer;text-decoration:underline;float:right">clear</span>
  <div class="dbox"><div class="legend" id="disfilters"></div>
  <label style="display:block;margin-top:4px"><input type=checkbox id="disconn"> <b style="color:#1c2330">Only their direct neighbours</b></label></div>
  <div class="mut help">Ticking a disease <b>keeps it</b>: the canvas then shows those diseases only, and drops the rest along with the edges that led to them. Edges between genes and drugs keep only the sentences from <b>papers about a kept disease</b> (a paper&rsquo;s diseases are those its own sentences relate to anything), so a gene pair under <i>brain</i> shows glioma papers, not melanoma ones; an edge left with too few sentences for <i>Min unique sentences / publications</i> is dropped. Edges to a kept disease are untouched. Diseases of the same tissue (read from the name, as for <em>Stack same-tissue diseases</em>) are listed together under a bold tissue heading whose box ticks or unticks them all &mdash; half-ticked when only some are kept; the number beside it is how many diseases it holds and their edges in all. <b>Tick nothing and nothing is filtered</b>, which is the default; untick the last one (or press <em>clear</em>) to go back. The list is rebuilt on every redraw from the diseases actually in front of you, with the number of drawn edges each one carries &mdash; so it follows the score, year, text and every other control. It is built from the view <em>before</em> this filter is applied: otherwise ticking one disease would empty the menu you are ticking from. A disease you have ticked stays listed even when the current filters leave it nothing, shown <span style="color:#b3243b">(0)</span>, so you can always untick it. <b>Only their direct neighbours</b> switches on by itself when you tick the first disease (untick it to see the whole graph with only those diseases), and goes further: it keeps just the ticked diseases, the genes and drugs with a <em>drawn</em> edge straight to one of them, and the edges among those nodes &mdash; everything else drops out. Every other setting counts: a gene whose own link to the disease is hidden (its relation unticked, an Enrich button, min sentences, min connections, cluster size) is not a neighbour, however it connects to the ones that are. The <b>Match text</b> box is the exception: it searches the neighbourhood without redrawing it &mdash; the neighbours are chosen as if the box were empty, and the text then keeps the edges among them that match. Genes are linked to diseases almost only by <b>associated</b> (drugs also by <b>inhibits</b>), so with that unticked the neighbourhood shrinks to a few drugs, or to the diseases themselves. A ticked disease left with no link is drawn on its own, and the line under the graph says why. It does nothing until a disease is ticked.</div></div>
 <div class="row" id="enrichhead" style="margin-bottom:2px"><b style="color:#1c2330">Enrich for alteration:</b></div>
 <div id="enrichgrid" class="egrid">
 <div class="row"><button id="drvbtn" aria-pressed="false">Candidate driver</button> <button class="ihelp" aria-label="About enrich for candidate drivers" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence naming a <b>single amino acid change</b>, one-letter or three-letter, case-sensitive (<i>V600E</i>, <i>p.Gly12Cys</i>), with any punctuation around it (<i>G12C,</i> <i>V600E-mutant</i> <i>KRAS-G12V</i> <i>L858R/T790M</i> <i>p.L747P</i>) or glued to a gene symbol (<i>KRASG12C</i>, <i>p110&alpha;H1047R</i>), and excluding the figure, table, gene-name and cell-line look-alikes.</div></div>
 <div class="row"><button id="fusbtn" aria-pressed="false">Fusions</button> <button class="ihelp" aria-label="About enrich for fusions" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about a <b>gene fusion</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. A sentence counts if it says <i>fusion</i> (<i>EML4-ALK fusion</i>, <i>RET fusion-positive</i>; not <i>infusion</i>), or uses <i>translocation</i> in the chromosomal sense: after a fusion pair or a gene that drives fusions (<i>EML4-ALK translocation</i>, <i>ROS1 translocation</i>, <i>ALK-translocated</i>), with a karyotype (<i>t(11;22) translocation</i>), after <i>chromosomal</i>, <i>balanced</i>, <i>reciprocal</i> or <i>gene</i>, or as <i>translocation involving MLL</i>. Karyotype notation counts on its own, with no other word: <i>t(4;14)</i>, <i>t(11;14)</i>, <i>t(9;22)(q34;q11)</i>, <i>inv(16)</i>, <i>der(22)</i> &mdash; how myeloma and leukemia papers write their translocations. So does a fusion named by its genes in HGNC notation, <i>NUP98::NSD1</i>, <i>AGK::BRAF</i> (both sides genes of this graph); the older single colon, <i>TMPRSS2:ERG</i>, only with a fusion word in the sentence, since <i>EGFR:GRB2</i> or <i>MAPK8:JUN</i> are not fusions. A protein moving inside the cell &mdash; <i>nuclear translocation of p65</i>, <i>translocation of YAP to the nucleus</i> &mdash; does not count.</div></div>
 <div class="row"><button id="oexbtn" aria-pressed="false">Overexpression</button> <button class="ihelp" aria-label="About enrich for overexpression" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>overexpression</b> &mdash; any word starting <i>overexpres</i> or <i>over-expres</i>, in any case (<i>overexpression</i>, <i>overexpressed</i>, <i>over-expressing</i>) &mdash; and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. Note it catches both senses: a gene found overexpressed in tumours, and a gene overexpressed in an experiment.</div></div>
 <div class="row"><button id="ampbtn" aria-pressed="false">Amplification</button> <button class="ihelp" aria-label="About enrich for amplification" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>amplification</b> &mdash; any word starting <i>amplif</i>, in any case (<i>amplification</i>, <i>amplified</i>, <i>ERBB2-amplified</i>) &mdash; and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. The laboratory sense is not counted &mdash; a method or tool attached to the word (<i>PCR-amplified</i>, <i>amplified by PCR</i>, <i>amplified using primers</i>, <i>whole-genome amplification</i>), an assay noun after it (<i>amplification curve</i>, <i>amplification primers</i>), or a lab object amplified (<i>amplification of the cDNA</i>, <i>amplified fragments</i>) &mdash; while a gene amplification merely <em>measured</em> by PCR (<i>MET amplification detected by ddPCR</i>) still is.</div></div>
 <div class="row"><button id="delbtn" aria-pressed="false">Deletion / loss</button> <button class="ihelp" aria-label="About enrich for deletion / loss" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about a gene <b>deleted or lost</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>deletion</i> / <i>deleted</i> / <i>co-deletion</i> (<i>PTEN homozygous deletion</i>, <i>SMAD4-deleted</i>, <i>1p19q co-deletion</i>), a gene's <i>del</i> (<i>CDKN2A del</i>, <i>del(17p)</i>), <i>loss of heterozygosity</i> / <i>LOH</i>, <i>loss of function</i>, and <i>loss</i> tied to a gene or chromosome (<i>PTEN loss</i>, <i>loss of SMAD4</i>, <i>copy-number loss</i>, <i>9p loss</i>). Not counted: activating in-frame deletions (<i>EGFR exon 19 deletion</i>), engineered constructs (<i>domain deletion mutant</i>, <i>deletion constructs</i>), generic loss (<i>weight loss</i>, <i>cell loss</i>, <i>loss of viability</i>), and a lost signal rather than a lost gene (<i>loss of ERK phosphorylation</i>, <i>loss of AKT activity</i>).</div></div>
 <div class="row"><button id="trubtn" aria-pressed="false">Truncation</button> <button class="ihelp" aria-label="About enrich for truncation" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about a <b>truncating variant</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>truncating</i>, <i>truncation</i> as a mutation (<i>ASXL3 truncations</i>, <i>truncation mutations</i>, <i>premature truncation</i>), <i>truncated</i> with a mutation word (<i>truncated mutations</i>), <i>frameshift</i> / <i>frame-shift</i>, <i>nonsense</i>, <i>stop-gain</i>, <i>premature stop / termination codon</i>, splice-site variants (<i>NF2 splice-site variants</i>, <i>splicing mutation</i>), and stop or frameshift notation (<i>Q55*</i>, <i>R213Ter</i>, <i>P191fs</i>). Not counted: engineered truncated proteins (<i>truncated form of</i>, <i>truncation mutants</i>), splice isoforms (<i>AR splice variant AR-V7</i>), splicing biology and splicing-factor mutations (<i>alternative splicing</i>, <i>splicing factor mutations</i>), and <i>X</i> notation (<i>G719X</i>, where X is any amino acid).</div></div>
 <div class="row"><button id="exobtn" aria-pressed="false">Exon variants</button> <button class="ihelp" aria-label="About enrich for exon variants" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence naming a variant <b>by its exon</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>exon N</i> with a variant word close after it (<i>EGFR exon 19 deletion</i>, <i>HER2 exon 20 insertion</i>, <i>MET exon 14 skipping</i>, <i>PIK3CA exon 20 mutations</i>, <i>exon 9 PIK3CA mutations</i>, <i>exon 21 L858R</i>, <i>PDGFRA exon 18 (D842V)</i>), and the compact forms (<i>METex14</i>, <i>MET&Delta;ex14</i>, <i>ERBB2&Delta;Ex16</i>, <i>ex19del</i>, <i>EGFRex20ins</i>, <i>Ex20Ins</i>, <i>EGFR ex20</i>). Not counted: an exon named without a variant (<i>in exons 1 and 2</i>, <i>IRAK4 exon 4 by PRMT5</i>), fusion breakpoints (<i>FGFR3-TACC3 Ex14-16</i>), and look-alikes (the drug <i>EX527</i>, the gene <i>Deltex1</i>, guide RNAs such as <i>sgPpm1dexon6</i>).</div></div>
 <div class="row"><button id="cnbtn" aria-pressed="false">Copy number</button> <button class="ihelp" aria-label="About enrich for copy number" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>copy-number change</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>copy number</i> in any form (<i>copy number gain / loss / alteration / variation</i>, <i>copy-number</i>), the abbreviations <i>CNV</i>, <i>CNA</i>, <i>SCNA</i> and <i>CN gain / loss</i>, chromosome-arm and whole-chromosome changes (<i>1q gain</i>, <i>gain of 1q</i>, <i>8p loss</i>, <i>17q12 amplification</i>, <i>11q deletion</i>, <i>chromosome 7 gain</i>, <i>loss of chromosome 10</i>), and <i>copy-neutral LOH</i>. Not counted: the copy number of something that is not the tumour genome (<i>mtDNA / mitochondrial DNA / plasmid / viral / transgene copy number</i>).</div></div>
 <div class="row"><button id="metbtn" aria-pressed="false">Methylation</button> <button class="ihelp" aria-label="About enrich for methylation" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>DNA methylation or epigenetic silencing</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>hyper-</i> / <i>hypo-</i> / <i>de-methylation</i> (<i>BCAT1 promoter hypermethylation</i>, <i>hypomethylated CpGs</i>, <i>demethylating agent</i>), <i>DNA</i> / <i>promoter</i> / <i>CpG</i> / <i>cytosine</i> / <i>gene-body</i> methylation, <i>unmethylated</i> and <i>differentially methylated</i>, <i>methylome</i>, <i>DNA methyltransferase</i>, <i>CpG island</i>, and <i>epigenetic(ally) silenced</i>; a bare <i>methylation</i> (<i>MGMT methylation</i>, <i>methylation in SFRP2</i>) only when the sentence also speaks of DNA &mdash; a promoter, CpG, DNA, MGMT, bisulfite, 5-aza or decitabine. Not counted: RNA methylation (<i>m6A methylation</i>), histone and protein methylation (<i>H3K27 methylation</i>, <i>arginine methylation</i>, <i>PRMT5-mediated</i>, <i>p53 methylation</i>), and look-alikes such as <i>methylcellulose</i>.</div></div>
 <div class="row"><button id="itdbtn" aria-pressed="false">ITD</button> <button class="ihelp" aria-label="About enrich for ITD" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about an <b>internal tandem duplication</b> &mdash; <i>ITD</i> in capitals as a word of its own (<i>FLT3-ITD</i>, <i>FLT3/ITD</i>, <i>ITD-positive</i>) or <i>internal tandem duplication</i> spelled out &mdash; and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted.</div></div>
 </div>
 <div class="row" id="therapyhead" style="margin-bottom:2px"><b style="color:#1c2330">Enrich for therapy:</b> <button class="ihelp" aria-label="About enrich for therapy" aria-expanded="false">i</button>
  <div class="mut help">These buttons combine with each other in <b>OR</b>, and in <b>AND</b> with the other Enrich boxes: with <i>Candidate driver</i> and <i>Resistance</i> both on, a sentence must name a point mutation <em>and</em> speak of resistance (<i>T790M confers resistance to gefitinib</i>).</div></div>
 <div id="therapygrid" class="egrid">
 <div class="row"><button id="acqbtn" aria-pressed="false">Acquired resistance</button> <button class="ihelp" aria-label="About enrich for acquired resistance" aria-expanded="false">i</button>
  <div class="mut help">A narrower view than <i>Resistance</i>: keeps only edges with a sentence about <b>how resistance arises</b>. It counts <i>acquired</i>, <i>secondary</i>, <i>intrinsic</i>, <i>primary</i>, <i>adaptive</i>, <i>innate</i>, <i>de novo</i> and <i>emerging</i> resistance (<i>acquired erlotinib resistance</i>, <i>secondary resistant mutation</i>), <i>acquisition</i> or <i>development of resistance</i>, <i>mechanisms of resistance</i>, <i>resistance mechanisms / mutations</i>, <i>gatekeeper mutations</i>, and <i>bypass</i> signalling in a sentence that speaks of resistance or an inhibitor (<i>MET amplification bypasses EGFR inhibition</i>). With <i>Resistance</i> also on, this adds nothing.</div></div>
 <div class="row"><button id="resbtn" aria-pressed="false">Resistance</button> <button class="ihelp" aria-label="About enrich for resistance" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>resistance to a treatment</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>resistance</i> / <i>resistant</i> in any form (<i>gefitinib-resistant</i>, <i>acquired resistance</i>, <i>multidrug-resistant</i>), <i>chemoresistance</i>, <i>radioresistance</i> and <i>refractory</i>. Not counted: a cell surviving something that is not a treatment &mdash; <i>apoptosis</i>, <i>anoikis</i>, <i>ferroptosis</i> or <i>cell-death resistance</i>, <i>resistance to oxidative stress</i>, <i>insulin resistance</i>, antibiotic selection (<i>puromycin-resistant</i>) &mdash; unless the sentence also names a treatment (<i>resistance to apoptosis induced by EGFR-TKIs</i> counts).</div></div>
 <div class="row"><button id="slbtn" aria-pressed="false">Synthetic lethality</button> <button class="ihelp" aria-label="About enrich for synthetic lethality" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>synthetic lethality</b>: <i>synthetic lethal</i>, <i>synthetically lethal</i>, <i>synthetic lethality</i>, <i>synthetic sickness</i>, <i>collateral lethality</i>, and <i>SL</i> as a pairing (<i>SL interaction</i>, <i>SL gene pairs</i>, <i>MYC-SL</i>; case-sensitive, so the compound <i>SL-176</i> does not count).</div></div>
 <div class="row"><button id="bmkbtn" aria-pressed="false">Biomarker</button> <button class="ihelp" aria-label="About enrich for biomarker" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about a <b>biomarker</b>: <i>biomarker</i>, <i>predictive marker</i>, <i>companion diagnostic</i>, and a feature said to predict response, sensitivity, benefit, resistance or efficacy (<i>PBRM1 inactivation may predict benefit from anti-PD-1</i>, <i>predictive of everolimus response</i>). Prognostic biomarkers count too.</div></div>
 </div>
 <div class="row" id="clinhead" style="margin-bottom:2px"><b style="color:#1c2330">Enrich for clinical relevance:</b> <button class="ihelp" aria-label="About enrich for clinical relevance" aria-expanded="false">i</button>
  <div class="mut help">These buttons combine with each other in <b>OR</b>, and in <b>AND</b> with the other Enrich boxes: with <i>Candidate driver</i> and <i>Clinical evidence</i> both on, a sentence must name a point mutation <em>and</em> report patients or a trial.</div></div>
 <div id="clingrid" class="egrid">
 <div class="row"><button id="clnbtn" aria-pressed="false">Clinical evidence</button> <button class="ihelp" aria-label="About enrich for clinical evidence" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence reporting <b>patients in a study or on a treatment</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>clinical trial / study / response / benefit</i>, a trial phase (<i>phase II trial</i>, <i>phase I/II study</i>, <i>Phase 3</i>; only beside trial or study wording, so <i>G2/M phase</i> does not count), <i>randomised</i>, <i>open-label</i>, <i>multicentre</i>, <i>prospective</i> / <i>retrospective</i> studies, NCT numbers, response rates (<i>ORR</i>, <i>objective response rate</i>, <i>complete / partial response</i>), and <i>patients treated with / receiving</i> a drug. Plain <i>patients</i> is not enough: most such sentences are database cohorts.</div></div>
 <div class="row"><button id="prgbtn" aria-pressed="false">Prognosis</button> <button class="ihelp" aria-label="About enrich for prognosis" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>patient prognosis or survival</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>prognosis</i> / <i>prognostic</i>, named endpoints (<i>overall</i>, <i>progression-free</i>, <i>disease-free</i>, <i>relapse-free</i>, <i>event-free survival</i>; <i>PFS</i>, <i>DFS</i>, <i>RFS</i>, <i>EFS</i>), survival or outcome as a patient result (<i>poor survival</i>, <i>shorter overall survival</i>, <i>median survival</i>, <i>adverse outcome</i>, <i>survival rate</i>), hazard ratios and Kaplan-Meier curves. <i>OS</i> counts only beside a survival cue, since it is also osteosarcoma, and <i>HR</i> only as a number (<i>HR = 0.64</i>), since it is also homologous recombination or hormone receptor. Not counted: survival of cells (<i>promotes cell survival</i>, <i>increased survival of tumour cells</i>).</div></div>
 <div class="row"><button id="vivbtn" aria-pressed="false">In vivo</button> <button class="ihelp" aria-label="About enrich for in vivo" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>animal experiments</b>, the level of evidence above cell culture, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>in vivo</i>, xenografts and allografts (<i>PDX</i>, <i>patient-derived xenograft</i>, <i>CDX</i>), <i>GEMM</i>, orthotopic / syngeneic / subcutaneous models, <i>tumour-bearing</i>, <i>mice</i> and <i>rats</i>, <i>mouse / murine / animal / zebrafish models</i>, and <i>knock-in</i> / <i>transgenic</i> mice. Not counted: <i>mouse</i> as a source of cells or reagents (<i>mouse embryonic fibroblasts</i>, <i>mouse anti-</i>), and patient-derived organoids or cells, which are in vitro.</div></div>
 </div>
 <div class="row" id="biohead" style="margin-bottom:2px"><b style="color:#1c2330">Enrich for gene role and biology:</b> <button class="ihelp" aria-label="About enrich for gene role and biology" aria-expanded="false">i</button>
  <div class="mut help">These buttons combine with each other in <b>OR</b>, and in <b>AND</b> with the other Enrich boxes: with <i>Candidate driver</i> and <i>Oncogene</i> both on, a sentence must name a point mutation <em>and</em> speak of an oncogene.</div></div>
 <div id="biogrid" class="egrid">
 <div class="row"><button id="tsgbtn" aria-pressed="false">Tumor suppressor</button> <button class="ihelp" aria-label="About enrich for tumor suppressor" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with at least one end that the <b>OncoKB Cancer Gene List</b> classes as a <b>tumor suppressor</b> (<i>TSG</i>, or both oncogene and TSG, such as <i>NOTCH1</i> or <i>TP63</i>), and the nodes those edges connect. The gene decides, not the wording: every sentence of such an edge counts, and nothing is highlighted. OncoKB list of 2026-11-24, 408 tumor suppressors.</div></div>
 <div class="row"><button id="oncbtn" aria-pressed="false">Oncogene</button> <button class="ihelp" aria-label="About enrich for oncogene" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with at least one end that the <b>OncoKB Cancer Gene List</b> classes as an <b>oncogene</b> (<i>ONCOGENE</i>, or both oncogene and TSG), and the nodes those edges connect. The gene decides, not the wording: every sentence of such an edge counts, and nothing is highlighted. OncoKB list of 2026-11-24, 518 oncogenes.</div></div>
 <div class="row"><button id="mtsbtn" aria-pressed="false">Metastasis / EMT</button> <button class="ihelp" aria-label="About enrich for metastasis / EMT" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>metastasis or the epithelial-mesenchymal transition</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>metastasis</i>, <i>metastases</i>, <i>metastatic</i>, <i>metastasize</i>, <i>EMT</i> (capitals) and <i>epithelial(-to)-mesenchymal transition</i>. Note <i>metastatic</i> also marks a disease stage (<i>metastatic NSCLC patients</i>); those sentences count too.</div></div>
 <div class="row"><button id="immbtn" aria-pressed="false">Immune</button> <button class="ihelp" aria-label="About enrich for immune" aria-expanded="false">i</button>
  <div class="mut help">Keeps only edges with a sentence about <b>tumor immunology</b>, and the nodes those edges connect; the matching sentences are the ones shown, the hits highlighted. It counts <i>immune</i>, <i>immunity</i>, <i>immunotherapy</i>, <i>immunosuppressive</i>, <i>immunogenic(ity)</i>, <i>immunomodulatory</i>, <i>immunosurveillance</i>; checkpoint inhibitors or blockade; <i>PD-1</i>, <i>PD-L1</i>, <i>CTLA-4</i>, <i>LAG-3</i>, <i>TIM-3</i>, <i>TIGIT</i>, <i>ICI</i> / <i>ICB</i>, <i>CAR-T</i>; T and NK cells, <i>CD8+</i> / <i>CD4+</i>, TILs, macrophages and TAMs, dendritic cells, MDSCs, Tregs, neoantigens and interferon gamma. Not counted: the lab methods named <i>immuno-</i> (<i>immunoblot</i>, <i>immunohistochemistry</i>, <i>immunofluorescence</i>, <i>immunoprecipitation</i>, <i>immunostaining</i>), cell-cycle and DNA-damage checkpoints (<i>G2/M checkpoint</i>), T- and NK-cell lymphomas, and a singular <i>TAM</i>, which is also tamoxifen.</div></div>
 </div>
 <div class="row" id="zoom"><button id="zin">+ Zoom in</button><button id="zout">&minus; Zoom out</button><button id="zfit">Fit</button></div>
 <div class="row mut" id="info">Click a node or edge for details.</div>
</div>
<div id="panel">
 <button class="pclose" data-panel="panel" title="Close this menu" aria-label="Close the right menu">&times;</button>
 <div class="row legend"><b style="background:#cfe3ff;border-color:#2b6cb0"></b>gene <b style="background:#1b7837;border-color:#145a28"></b>approved anti-neoplastic <b style="background:#e08600;border-color:#9a6700"></b>approved (other) <b style="background:#c2185b;border-color:#7a0f3a"></b>no approved drug <button class="ihelp" data-help="col" aria-label="About node colours" aria-expanded="false">i</button></div>
 <div class="row mut help" data-help="col">Drug-target genes (corpus chemicals): deeper colour = more chemicals. <b style="color:#1b7837">Green</b> = DGIdb approved anti-neoplastic, <b style="color:#e08600">amber</b> = DGIdb approved (non-anti-neoplastic), <b style="color:#c2185b">pink</b> = no approved drug: DGIdb investigational and/or a ChEBI target role only (the gene may still be in DGIdb).</div>
 <div class="row">Min unique sentences/edge: <b id="thv">1</b><br><input id="thr" type="range" min="1" max="10" value="1"></div>
 <div class="row">Min unique publications: <b id="mpv">1</b> <button class="ihelp" aria-label="About min unique publications" aria-expanded="false">i</button><br><input id="minpub" type="range" min="1" max="10" value="1">
  <div class="mut help">Distinct PMIDs behind an edge; raise it to drop relations that rest on one paper repeating itself.</div></div>
 <div class="row">Font size: <b id="fsv">100%</b><br><input id="fscale" type="range" min="10" max="100" step="5" value="100" aria-label="Label font size"></div>
 <div class="row">Min cluster size: <select id="mincluster"><option value="2" selected>2 (off)</option><option>3</option><option>4</option><option>5</option><option>6</option><option>7</option><option>8</option><option>9</option><option>10</option><option>11</option><option>12</option></select> <button class="ihelp" aria-label="About min cluster size" aria-expanded="false">i</button><div class="mut help">Drops connected components smaller than this. <b>2 is off, not a floor</b>: the components are built from the edges that survive every other filter, so the smallest one possible already has two nodes and a pair joined by a single edge always passes. Raise it to 3 or more to peel the small islands a loose score leaves behind &mdash; and note it runs AFTER min-connections, so a node that drops below that bar can take its whole component with it.</div></div>
 <div class="row">Min connections: <select id="mindeg"><option>1</option><option selected>2</option><option>3</option><option>4</option><option>5</option><option>6</option></select> <button class="ihelp" aria-label="About min connections" aria-expanded="false">i</button><div class="mut help">Hides genes linked to fewer than this many others; thins the hairball's single-link fringe. Applied until it holds: removing a node can push its neighbours below the bar, so they go too, and every node left has at least this many links <em>on screen</em>. <b>1</b> turns it off and keeps that fringe &mdash; the only way to see a node whose whole evidence is one edge. Expect a denser canvas: it is the setting that admits every hapax in the corpus.</div></div>
 <div class="row">Search gene: <input id="search" placeholder="e.g. EGFR" autocomplete="off"></div>
 <div class="row">Filter to gene:<br><input id="genefilter" placeholder="e.g. EGFR (+neighbors)" autocomplete="off"> <select id="hops"><option value="1">1 hop</option><option value="2">2 hops</option></select></div>
 <div class="row">Search drug: <input id="drugsearch" placeholder="e.g. nivolumab" autocomplete="off"></div>
 <div class="row">Filter to drug:<br><select id="chemfilter"><option value="">(all drugs)</option></select></div>
 <div class="row">Match text in sentence: <button class="ihelp" aria-label="About the text filter" aria-expanded="false">i</button><br><input id="textfilter" placeholder="e.g. phosphorylat or /inhibit(s|ed)?/" autocomplete="off">
  <div class="mut help">Case-insensitive substring; wrap in / / for a regex. Keeps only edges with a matching sentence, and shows just those sentences. Regular expressions are the way to find edges that describe <b>variants</b>: for example, <code>/\s+[A-Z]\d+[A-Z]\s+/</code> finds edges that describe missense and nonsense variants. It is a search: an edge that passes <b>Min unique sentences</b> and <b>Min unique publications</b> on all its sentences is kept if any one of them matches (the Enrich buttons, by contrast, count only their own matching sentences against those minimums). <b>Min connections</b> and <b>Min cluster size</b> are then applied to what the query leaves.</div></div>
__KINDROW__
 <div class="row">Relation type <button class="ihelp" data-help="rel" aria-label="About relation types" aria-expanded="false">i</button>
  <div class="mut help" data-help="rel">As predicted by the RE model; &ldquo;not X&rdquo; = negated statement, drawn dashed. Unticking one hides <em>sentences</em> with that label, and any edge left without support.</div></div><div id="catfilters"></div><div class="row" style="margin-top:4px"><button id="catall" type="button">Check all</button></div>
 <div class="row mut help" data-help="rel">Edge colour = the relation the model predicted. <b>activates</b>/<b>inhibits</b> are signed and come from the BioRED checkpoint; <b>interacts</b> is the unsigned PPI verdict. Only <b>activates</b>/<b>inhibits</b> edges have an arrowhead (agent &rarr; target); <b>interacts</b>, <b>associated</b> and <b>binds</b> are symmetric and drawn without one. An edge takes its best-supported direction, and is drawn as the relation most of its sentences <em>in view</em> carry &mdash; so narrowing the filters can recolour an edge. Hover for the per-sentence labels. Thickness and arrowhead size follow the number of <b>independent publications</b> behind the edge, not its sentence count &mdash; one paper repeating itself never thickens a line.</div>
 <div class="row mut help" data-help="rel">Counts read <em>total &middot; in view</em>: the total is every edge in the file carrying at least one sentence of that type (an edge with mixed readings counts under each, so the totals exceed the edge count), &ldquo;in view&rdquo; is how many survive the current score, year, text, min-publications, min-connections and min-cluster settings. <span style="color:#b3243b">A red 0</span> means the type is ticked but everything of it is pruned &mdash; usually its edges sit in components smaller than <b>Min cluster size</b>, so lower that (or the score) to see them.</div>
 <div class="row">Training set behind the edge <button class="ihelp" data-help="src" aria-label="About training sets" aria-expanded="false">i</button></div>
 <div class="row" id="srcbtns"></div>
 <div class="row mut help" data-help="src" id="srchint">Which model kept each <em>sentence</em>. <b>PPI model</b> = every sentence the BioInfer/PPI model kept, shown with <em>its</em> reading &mdash; it has one label, <b>interacts</b>, so every other relation type reads 0 there. <b>BioRED model</b> = every sentence the BioRED model kept, with its typed (often signed) label. A sentence both models kept appears under each button with that model&rsquo;s reading (the tooltip gives both); there is no &ldquo;both&rdquo; view, because the two models do not share relation types. Under <b>All</b> each such sentence carries both readings and matches either relation box, so ticking just <b>interacts</b> shows exactly the PPI-model view, and a typed box exactly the BioRED-model view. Each count is the number of <b>edges holding at least one such sentence</b>, over the whole graph rather than the current filters; an edge counts under every model behind it, so the counts sum past All.</div>
 <div class="row">Relationship score: <b id="scval">&ge;0.99</b><br><input id="conf" type="range" min="0" max="13" step="1" value="__CONFDEF__" aria-label="Minimum relationship score"></div>
 <div class="row mut" id="stats"></div>
</div>
<div id="net"></div>
<script>
const DATA=__PAYLOAD__;
const CCOLOR=__CCOLOR__;
const DIRECTED=new Set(['activates','inhibits']);   // the only relations with a direction
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
// No "kept by both" button: the two models have different relation vocabularies (PPI says only
// "interacts", BioRED types the relation), so the overlap has no relation type of its own. A
// sentence both kept shows under each model's button with that model's reading.
const SRC_BTN=[['all','All'],['ppi','PPI model'],['biored','BioRED model']];
function buildSrcButtons(){
 // counted like the relation types: edges holding at least one sentence that model kept, so an
 // edge counts under every model behind it (the parts sum past the edge total)
 const n={all:DATA.edges.length,ppi:0,biored:0,both:0};
 DATA.edges.forEach(e=>{const r={};e.sents.forEach(s=>['ppi','biored'].forEach(k=>{if(srcPass(e,s,k))r[k]=1;}));for(const k in r)n[k]=(n[k]||0)+1;});
 document.getElementById('srcbtns').innerHTML=SRC_BTN.map(([k,lab])=>
   '<button class="srcb'+(k===SRC_MODE?' on':'')+'" data-src="'+k+'"'+(n[k]?'':' disabled')+'>'
   +lab+' ('+(n[k]||0)+' edges)</button>').join('');
 document.querySelectorAll('.srcb').forEach(b=>b.addEventListener('click',()=>{
   SRC_MODE=b.getAttribute('data-src');
   document.querySelectorAll('.srcb').forEach(x=>x.classList.toggle('on',x===b));
   build(+thr.value);}));
 // one checkpoint in the run -> nothing to separate; keep the row but say so
 if(!n.biored)document.getElementById('srchint').innerHTML=
   'Only one RE checkpoint stands behind this graph, so every sentence comes from the <b>PPI model</b>; '
   +'run step&nbsp;2 with both models (<code>--route-mode additive</code>) to split them.';
}
function activeKinds(){const b=[...document.querySelectorAll('.kindf')];return b.length?new Set(b.filter(c=>c.checked).map(c=>c.value)):null;}
const net=document.getElementById('net'); let network=null, NODEDS=null;
// the layout exactly as physics left it, plus the tissue groups of what is on screen; every
// position control replays from these, so they compose and none of them accumulates drift
let BASEPOS=null, TGROUPS=[];
// Every node's last stabilized position, kept across rebuilds (and across the times a filter
// hides it), so changing a filter does not re-roll the whole layout: nodes seen before go back
// where they were and are held there while physics places only the new ones, next to their
// neighbours. A position, once learned, is only replaced when the node is laid out afresh.
const LAYOUT={};
// Fit, but never closer than FIT_MAX_SCALE (the size nodes and labels are designed for): a view of
// three nodes used to be blown up to fill the whole screen. Frames the drawn positions, with margin.
const FIT_MAX_SCALE=1.0;
function safeFit(anim){
 const P=Object.values(network.getPositions());if(!P.length)return;
 let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity;
 P.forEach(p=>{x0=Math.min(x0,p.x);x1=Math.max(x1,p.x);y0=Math.min(y0,p.y);y1=Math.max(y1,p.y);});
 const W=net.clientWidth||800,H=net.clientHeight||600,pad=80;
 const sc=Math.min(FIT_MAX_SCALE,0.92*Math.min(W/(x1-x0+2*pad),H/(y1-y0+2*pad)));
 network.moveTo({position:{x:(x0+x1)/2,y:(y0+y1)/2},scale:sc,animation:anim?{duration:400}:false});
}
// every drawn node inside the current viewport?
function allInView(){
 const W=net.clientWidth||800,H=net.clientHeight||600;
 const a=network.DOMtoCanvas({x:0,y:0}),b=network.DOMtoCanvas({x:W,y:H});
 return Object.values(network.getPositions()).every(p=>p.x>=a.x&&p.x<=b.x&&p.y>=a.y&&p.y<=b.y);
}
function idJitter(id,r){let h=0;for(const c of String(id))h=(h*31+c.charCodeAt(0))|0;
 const a=(h%360)*Math.PI/180;return [Math.cos(a)*r,Math.sin(a)*r];}
function applyLayoutShape(){
 if(!network||!BASEPOS)return;
 reshapeRadial(network,BASEPOS,activeShrink(),activeExpand());
 if(activeTissueStack()&&TGROUPS.length)stackTissues(network,TGROUPS);
 network.redraw();
}
// Labels once faded in as you zoom: small graphs always showed every symbol; dense views revealed labels as the
// zoom scale climbs from LABEL_LO to LABEL_HI. Opacity is driven through the shared node-font colour, so
// one setOptions call recolours all labels (per-node font carries only size, inheriting this colour).
const LABEL_SMALL=60, LABEL_LO=0.05, LABEL_HI=0.25;
let LABEL_N=0, LABEL_A=-1;
// ...no longer: a label is never faded out. Zoomed out on a dense view the fade took labels to
// opacity 0 and left bare dots, so every label is now drawn opaque at every zoom.
function labelOpacity(sc){ return 1; }
// Labels are sized in SCREEN pixels, not canvas units: a node's canvas font is its on-screen
// size divided by the zoom, re-applied whenever the zoom changes. They used to live in canvas
// units, so what you saw was size x zoom -- a large view zoomed to 0.05 drew everything at the
// floor, a small one at zoom 1 drew busy genes at 96px, and the same 50% slider looked nothing
// alike from one filter to the next. ZSCALE is the zoom the fonts were last sized for.
const LABEL_MIN_PX=8;
let ZSCALE=1, zoomTimer=null;
// on-screen px for a base size (fontSize(), 11-20 for genes): the slider stretches the part
// above the 8px floor, so 50% shows the base sizes, 10% runs ~8.6-10.4px and 100% 14-32px --
// every step changes every label and none is ever under 8px
function screenPx(base){return LABEL_MIN_PX+FSCALE*Math.max(0,base-LABEL_MIN_PX);}
function refloorLabels(){ if(!network||!NODEDS)return; const sc=network.getScale(); if(!(sc>0))return;
 if(Math.abs(sc/ZSCALE-1)<0.01)return; ZSCALE=sc;
 NODEDS.update(NODEDS.get().map(n=>({id:n.id,font:nodeFont(KIND[n.id]||'gene',n._fs||11)})));}
function updateLabels(){ if(!network)return;
 if(LABEL_A!==1){LABEL_A=1;network.setOptions({nodes:{font:{color:'rgba(26,26,26,1)'}}});}
 // re-size on every wheel tick for a small view (labels stay put on screen), after the wheel
 // settles for a big one (one DataSet update per gesture, not one per tick)
 clearTimeout(zoomTimer);
 if(NODEDS&&NODEDS.length<=1500)refloorLabels();else zoomTimer=setTimeout(refloorLabels,120);}
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
 const re=new RegExp(src,flags),g=new RegExp(src,flags+'g');
 return {test:t=>re.test(t||''),hits:t=>reHits(g,t)};
}
// every [start,end) span a global regex matches in t
function reHits(g,t){const iv=[];g.lastIndex=0;let m;
 while((m=g.exec(t))!==null){
  if(!m[0].length){g.lastIndex++;continue;}   // zero-length match (e.g. /x*/): skip, never loop
  iv.push([m.index,m.index+m[0].length]);
 }
 return iv;}
// "Candidate driver" (Enrich for alteration): a fixed, case-SENSITIVE point-mutation lens (V600E, G12C as a word of
// their own, or in parentheses: "KRAS (G12C)"), ANDed with the text box -- textMatcher would
// add the i flag and let "a12b" match.
// The change itself is the token -- one-letter (G12C) or three-letter (p.Gly12Cys) -- and what
// may stand on either side is checked separately (driverLeftOK). It once had to sit between
// whitespace or inside parentheses, which lost "KRAS G12C, ...", "V600E-mutant", "KRAS-G12V",
// "L858R/T790M", "p.L747P", "KRASG12V" and "p110αH1047R" -- about a third of all mentions.
// Right side: no letter or digit may follow, nor "." plus one ("H2A.X" is a histone).
const AA3='Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val';
const DRIVER_RE=new RegExp('(?<![A-Za-z0-9])(?:p\\.)?(?:'+AA3+')\\d{1,4}(?:'+AA3+'|Ter)(?![A-Za-z0-9])|[A-Z]\\d{1,4}[A-Z](?![A-Za-z0-9])(?!\\.[A-Za-z0-9])','g');
// Left side: start, space, bracket or any punctuation ("-", "/", ",", "+", ";", ":"); a "." only
// in HGVS "p." / "c." ("Fig.S2F" is a panel); glued to letters only when they are a whole gene
// symbol in this graph (KRASG12C, TP53R273H, KrasG12D -- not CDKN2A, INPP4B, HIF1A, KMT2A, whose
// front is no gene) or end in a Greek letter (p110αH1047R, ERαY537S). "NCI-H295R" is a cell line.
let GENESYM=null;
function driverLeftOK(t,a){
 const c=t[a-1];
 if(c===undefined)return true;
 if(c==='.')return /(?:^|[^A-Za-z0-9])[pc]$/.test(t.slice(Math.max(0,a-3),a-1));
 if(/[A-Za-z0-9Ͱ-Ͽ]/.test(c)){
  if(/[Ͱ-Ͽ]/.test(c))return true;
  const run=/[A-Za-z0-9Ͱ-Ͽ]+$/.exec(t.slice(Math.max(0,a-20),a))[0];
  if(!GENESYM)GENESYM=new Set(DATA.nodes.filter(n=>n.kind==='gene').map(n=>String(n.id).toUpperCase()));
  // ...and the glued whole must not itself be a gene: PPARG + C1A is PPARGC1A
  const word=(run+(/^[A-Za-z0-9]+/.exec(t.slice(a))||[''])[0]).toUpperCase();
  return run.length>=3&&GENESYM.has(run.toUpperCase())&&!GENESYM.has(word);
 }
 if(/NCI-$/i.test(t.slice(Math.max(0,a-4),a)))return false;
 return true;
}
// Figure/table panels have the same shape -- "Figure S6D and E", "Supplementary Fig. S7G and
// S7H", "Figs. 6A and B, S3A", "Table S2A" -- so a hit is dropped when the text before it
// ends in a figure or table reference, optionally followed by a list of other panels (or, PLOS
// style, when "Fig" follows it). Only the
// last 80 characters are looked at: a reference is never far from its panels.
const FIG_CTX=/\b(?:fig(?:ure)?s?|tables?)\.?\s*(?:(?:[A-Z]?\d+[A-Za-z]{0,2}|[A-Za-z])\s*(?:[-–,;&]|\band\b|\bor\b)?\s*)*$/i;
// "fig. X1Y" with nothing in between: any letter-digits-letter label straight after a figure /
// table word is a panel, not a mutation (case-insensitive, "fig" / "fig." / "figure(s)" / "table")
const FIG_DIRECT=/\b(?:fig(?:ure)?s?|tables?)\.?\s*$/i;
// PLOS puts the word AFTER the panel: "S8B Fig", "(S2A and S2B Fig)", "S5B–S5D Figs"
const FIG_POST=/^(?:\s*(?:[-–,;&]|and|or)?\s*S\d{1,2}[A-Za-z]{0,2})*\s*Figs?\b/;
// Words of the same letter-digits-letter shape that are not amino-acid changes; a sentence
// matched only by these does not count (one naming a real mutation as well still does).
//   S6K  ribosomal S6 kinase (RPS6KB1)    E2F  transcription-factor family
//   T47D breast-cancer cell line          T98G, M059J  glioblastoma cell lines
//   G2M  the G2/M cell-cycle phase        S1P  sphingosine-1-phosphate, a lipid
//   E2A  TCF3 (E2A-PBX1, E2A-HLF)      G9A  EHMT2               H2A, H2B  histones
//   S6R  phospho-S6 ribosomal protein   S100P  a gene            B2M  beta-2-microglobulin
//   H3B  the H3B-6545 drug code
const DRIVER_STOP=new Set(['S6K','E2F','T47D','T98G','M059J','G2M','S1P','E2A','G9A','H2A','H2B','S6R','S100P','B2M','H3B']);
let DRIVERS=false;
function driverHits(t){t=t||'';
 return reHits(DRIVER_RE,t).filter(([a,b])=>{
  if(DRIVER_STOP.has(t.slice(a,b).replace(/[\s()]/g,'')))return false;
  if(!driverLeftOK(t,a))return false;
  while(/[\s(]/.test(t[a]))a++;            // the hit's own leading space / parenthesis
  // a label RIGHT after the word is a panel whatever its letter -- "Fig. S6K", "Fig. E3B"
  // (extended data), "Figure A1C" (appendix), "Table R2A": no mutation is written that way
  if(FIG_DIRECT.test(t.slice(Math.max(0,a-40),a)))return false;
  // further along a panel list, only a supplementary panel can be a figure label -- S, one or two digits, a letter (S6K,
  // S12B; a main-text panel is "3B", which the pattern never matches) -- so "Figure 3 and V600E
  // cells" keeps its V600E and "Fig. 2B and S768I mutants" its S768I
  const m=/^S\d{1,2}[A-Z]\b/.exec(t.slice(a));
  if(!m)return true;
  return !FIG_CTX.test(t.slice(Math.max(0,a-80),a))&&!FIG_POST.test(t.slice(a+m[0].length,a+m[0].length+80));
 });}
function driverMatcher(){return {test:t=>driverHits(t).length>0,hits:driverHits};}
// "Enrich for fusions": a sentence about a gene fusion. "fusion(s)" as a word of its own counts
// as it stands (so "infusion" / "transfusion" do not); "translocation" only in the CHROMOSOMAL
// sense -- in this corpus it is mostly a protein moving ("nuclear translocation of p65"), and a
// hit counts only with a fusion cue next to it and no compartment word before it.
const FUSION_WORD=/\bfusions?\b/gi;
const TRANSLOC_WORD=/\btransloca(?:tions?|ted)\b/gi;
// genes whose translocations are the fusion drivers of this literature (partners included)
const FUSION_GENES='ALK|ROS-?1|RET|NTRK[123]|NRG1|FGFR[123]|MYC|MLL|KMT2A|BCR|ABL1?|BCL2|BCL6|CCND1|IGH|ETV6|RUNX1|EWSR1|FLI1|PAX[37]|FOXO1|TFE3|TFEB|NUTM1|BRD4|PML|RARA|CBFB|MYH11|SS18|SSX[12]?|CIC|DUX4|NUP98|MECOM|NAB2|STAT6|MAML2|CRTC1|TMPRSS2|ERG|ETV[145]|EML4|KIF5B|CD74|CCDC6|NCOA4|CAMTA1|WWTR1|TAZ|YAP1|MA9|AF4|AF9|MLLT3|AFF1|PDGFR[AB]|JAK2|CRLF2|ZNF384|PBX1|TCF3|RSPO[23]';
// a compartment / direction before it: the protein-trafficking sense
const TL_NEG_PRE=/\b(?:nuclear|nucleus|cytoplasm(?:ic)?|cytosolic|mitochondrial|membrane|plasma|lysosomal|subcellular|golgi)\s*-?\s*$/i;
// the chromosomal cues: a karyotype or chromosomal word just before ("balanced", "truncating",
// "fusion/"), a fusion pair "X-Y" / "X::Y" / "X/Y" or a fusion gene just before ("ALK
// translocation", "RSPO3-translocated"), or a fusion gene with a sibling alteration ("MYC
// amplification or translocation")
// (String.raw keeps every backslash as written: in a plain '...' string "\b" is a backspace)
const TL_POS_PRE=new RegExp(String.raw`(?:\bt\(\d+;\d+\)|\b(?:chromosom\w*|balanced|reciprocal|genomic|gene|recurrent|truncating|fusion\s*/)|\b(?:`+FUSION_GENES+String.raw`)\s+(?:\w+\s+){0,2}(?:or|and)|\b[A-Z][A-Za-z0-9]{1,9}\s*(?:-|–|::|/)\s*[A-Z][A-Za-z0-9]{1,9}|\b(?:`+FUSION_GENES+String.raw`))[\s-]*$`);
// ... or a karyotype / fusion gene just after: "translocation t(9;22)", "translocations of MYC",
// "a translocation involving MLL"
const TL_POS_POST=new RegExp(String.raw`^\s*(?:\(?t\(\d+;\d+\)|(?:of|involving|between|affecting)\s+(?:the\s+)?(?:`+FUSION_GENES+String.raw`)\b)`);
// Karyotype notation counts on its own, as the chromosomal event it names: t(4;14), t(11;14),
// t(9;22)(q34;q11), t(X;18), ins(3;3), and the one-chromosome inv(16) / der(22) (a deletion,
// del(17p), is a loss -- the Deletion button's business, not this one).
// Myeloma's translocations are written almost only this way -- "myeloma cells carrying
// t(4;14)" -- and none of those sentences said "fusion". A semicolon is required for t / ins,
// so a t-statistic "t(28) = 3.1" never matches; "non-t(4;14)" still does (it names the event).
const KARYO_RE=/\b(?:t|ins)\(\s*[0-9XY]{1,2}\s*;\s*[0-9XY]{1,2}\s*\)(?:\s*\(\s*[pq][\d.]*\s*;\s*[pq][\d.]*\s*\))?|\b(?:inv|der)\(\s*[0-9XY]{1,2}\s*\)(?:\s*\(\s*[pq][\d.;pq\s]*\))?/g;
// A fusion named by its two genes. HGNC's double colon (NUP98::NSD1, AGK::BRAF, EWSR1::WT1) means
// nothing else, so two gene symbols of this graph around "::" count on their own -- 62 sentences
// wrote one without ever saying "fusion" ("AGK::BRAF induces robust ... activation"). The older
// single colon (TMPRSS2:ERG, KIAA1549:BRAF) is also proximity-assay pairs (EGFR:GRB2 PLA),
// kinase:substrate (MAPK8:JUN), statistical interactions (KRAS:CDKN2A) and transgenes (rag2:Myc),
// so it counts only with no space around it AND a fusion cue somewhere in the sentence.
const FUS_PAIR2=/\b([A-Za-z][A-Za-z0-9-]{1,14})\s?::\s?([A-Za-z][A-Za-z0-9-]{1,14})\b/g;
const FUS_PAIR1=/\b([A-Za-z][A-Za-z0-9-]{1,14})(?<!:):(?!:)([A-Za-z][A-Za-z0-9-]{1,14})\b/g;
const FUS_CUE=/\b(?:fusions?|fused|transloca\w*|rearrang\w*|chimer\w*)\b/i;
function fusPairHits(t){
 if(!GENESYM)GENESYM=new Set(DATA.nodes.filter(n=>n.kind==='gene').map(n=>String(n.id).toUpperCase()));
 const both=(a,b)=>GENESYM.has(a.toUpperCase())&&GENESYM.has(b.toUpperCase());
 const iv=[];let m;
 FUS_PAIR2.lastIndex=0;while((m=FUS_PAIR2.exec(t))!==null)if(both(m[1],m[2]))iv.push([m.index,m.index+m[0].length]);
 if(FUS_CUE.test(t)){FUS_PAIR1.lastIndex=0;while((m=FUS_PAIR1.exec(t))!==null)if(both(m[1],m[2]))iv.push([m.index,m.index+m[0].length]);}
 return iv;}
function fusionHits(t){t=t||'';
 const iv=reHits(FUSION_WORD,t).concat(reHits(KARYO_RE,t),fusPairHits(t));
 reHits(TRANSLOC_WORD,t).forEach(([a,b])=>{
   const pre=t.slice(Math.max(0,a-40),a),post=t.slice(b,b+40);
   if(TL_NEG_PRE.test(pre))return;
   if(TL_POS_PRE.test(pre)||TL_POS_POST.test(post))iv.push([a,b]);
 });
 return iv;}
let FUSIONS=false;
function fusionMatcher(){return {test:t=>fusionHits(t).length>0,hits:fusionHits};}
// "Enrich for ITD": internal tandem duplications, nearly all FLT3-ITD in this corpus. "ITD" is
// matched case-SENSITIVELY as a word ("FLT3-ITD", "FLT3/ITD", "ITD-positive"), the long form
// in any case.
const ITD_RE=/\bITDs?\b|\binternal tandem duplications?\b/gi;
let ITDS=false;
function itdHits(t){t=t||'';return reHits(ITD_RE,t).filter(([a,b])=>{const w=t.slice(a,b);return w.length>4||/^ITDs?$/.test(w);});}
function itdMatcher(){return {test:t=>itdHits(t).length>0,hits:itdHits};}
// "Enrich for overexpression": "overexpres" / "over-expres" (any hyphen glyph), any case, as the
// start of a word -- the rest of the word is highlighted with it
const OEX_RE=/\bover[-‐‑–]?expres\w*/gi;
let OEXS=false;
function oexHits(t){return reHits(OEX_RE,t||'');}
function oexMatcher(){return {test:t=>oexHits(t).length>0,hits:oexHits};}
// "Enrich for amplification": any word starting "amplif" (also after a hyphen: "ERBB2-amplified")
const AMP_RE=/\bamplif\w*/gi;
let AMPS=false;
// The laboratory sense -- copying DNA, not a gene gained in the genome -- is recognised from what
// is attached to the word itself, never from lab words elsewhere in the sentence: "MET
// amplification detected by FISH and plasma ddPCR" is a gene amplification MEASURED by PCR and
// must count. A hit is the lab sense when
//   just before it: a method ("PCR-amplified", "RT-PCR amplification", "whole-genome
//     amplification", "isothermal / rolling-circle / tyramide signal amplification"; a bare
//     "signal amplification" is not one -- under FISH it IS the gene amplification, and
//     "death signal amplification" is biology), or
//   just after it: the method or the tool ("amplified by PCR", "amplified using primers"), a
//     run-of-the-assay noun ("amplification curve / efficiency / cycle / primers / kit /
//     reaction / products"), or a lab object amplified ("amplification of the cDNA / DNA
//     fragment / PCR product / exon 19", "amplified fragments").
const AMP_LAB_PRE=/\b(?:(?:RT-?|q|dd|nested\s|real-time\s|multiplex\s)?PCR|whole[- ]genome|isothermal|rolling[- ]circle|tyramide(?:\s+signal)?|loop-mediated|library)[\s-]*$/i;
const AMP_LAB_POST=/^\s*(?:(?:was|were|is|are)\s+)?(?:(?:by|via|through)\s+(?:RT-?|q|dd|nested\s|real-time\s)?PCR|(?:using|with)\s+(?:the\s+|specific\s+|gene-specific\s+)?primers?|(?:curves?|efficienc\w*|plots?|cycles?|primers?|kits?|reactions?|steps?|products?|fragments?|amplicons?)\b|of\s+(?:the\s+)?(?:cDNA|DNA\s+fragments?|fragments?|PCR\s+products?|amplicons?|templates?|exons?\s*\d+)\b)/i;
function ampHits(t){t=t||'';return reHits(AMP_RE,t).filter(([a,b])=>!AMP_LAB_PRE.test(t.slice(Math.max(0,a-40),a))&&!AMP_LAB_POST.test(t.slice(b,b+60)));}
function ampMatcher(){return {test:t=>ampHits(t).length>0,hits:ampHits};}
// "Enrich for deletion / loss": a gene, region or allele deleted or lost.
//   deletion words: deletion(s) / deleted / deleting / co-deletion -- but not "deleterious"; not
//     an activating in-frame deletion ("exon 19 deletion" is an EGFR driver, not a loss); not an
//     engineered construct ("domain deletion mutant", "deletion constructs").
//   "del": only as a gene's call ("CDKN2A del", "PTEN del") or a karyotype ("del(17p)").
//   loss: "loss of heterozygosity" / LOH / "loss of function" always; otherwise only tied to a
//     gene or chromosome -- a gene-like symbol before it ("PTEN loss", "Smad4 loss", "p53 loss")
//     or after "loss of" ("loss of PTEN", "loss of the SMAD4 gene"), or a genomic word before it
//     ("copy-number / allelic / biallelic / genomic / chromosomal / 9p loss"). "weight loss",
//     "cell loss", "loss of viability" therefore do not count.
const DEL_WORD=/\b(?:co-?)?delet(?:ions?|ed|ing)\b/gi;
const DEL_EXON_PRE=/\bexons?\s*\d+\s*(?:in-frame\s+)?$/i;
const DEL_CONSTRUCT_PRE=/\b(?:domain|region|motif|site|binding|terminal|truncation)\s*$/i;
const DEL_CONSTRUCT_POST=/^\s*(?:mutants?|constructs?|variants?\s+of\s+the\s+protein)\b/i;
const DEL_CALL=/\b[A-Z][A-Z0-9]{1,9}\s+del\b|\bdel\s*\(\s*\d+[pq]/g;
// a gene-like symbol: capitals with a digit or 2+ capitals (PTEN, SMAD4, RB1), Smad4-style, or p53-style
const GENE_TOKEN=String.raw`(?:[A-Z][A-Z0-9-]*[0-9][A-Za-z0-9-]*|[A-Z]{2,}[A-Z0-9-]*|[A-Z][a-z]+[0-9]+[a-z]?|p[0-9]{2,3}(?:ARF|INK4a)?)`;
const LOSS_WORD=/\bloss(?:es)?\b|\bLOH\b/g;
const LOSS_ALWAYS_POST=/^\s*of\s+(?:heterozygosity|function)\b/i;
const LOSS_PRE_GENE=new RegExp(GENE_TOKEN+String.raw`\s*-?\s*$`);
const LOSS_PRE_GENOMIC=/\b(?:copy[- ]number|allelic|biallelic|monoallelic|genomic|genetic|chromosom\w*|arm|\d{1,2}[pq]\d*(?:\.\d+)?)\s*-?\s*$/i;
const LOSS_POST_GENE=new RegExp(String.raw`^\s*of\s+(?:the\s+)?(?:wild-type\s+)?(?:`+GENE_TOKEN+String.raw`\s+)?`+GENE_TOKEN+String.raw`\b`);
// ...unless what is lost is the protein's SIGNAL, not the gene: "loss of ERK phosphorylation",
// "loss of S473 AKT phosphorylation", "loss of AKT activity" are pathway read-outs
const LOSS_SIGNAL_POST=new RegExp(String.raw`^\s*of\s+(?:the\s+)?(?:\S+\s+){0,3}?(?:phospho\w*|activit\w*|activation|signal\w*|binding|dimeri\w*|kinase\s+activit\w*)\b`,'i');
// words that look like symbols but are not genes, whose "loss" is generic
const LOSS_NOT_GENE=/\b(?:DNA|RNA|ATP|ROS|MMP|ECM|EMT|PFS|OS|II|III|IV)\s*-?\s*$/;
const LOSS_NOT_GENE_POST=/^\s*of\s+(?:the\s+)?(?:DNA|RNA|ATP|ROS|MMP|ECM|EMT|PFS|OS)\b/;
let DELS=false;
function delHits(t){t=t||'';const iv=[];
 reHits(DEL_WORD,t).forEach(([a,b])=>{const pre=t.slice(Math.max(0,a-30),a),post=t.slice(b,b+40);
   if(DEL_EXON_PRE.test(pre)||DEL_CONSTRUCT_PRE.test(pre)||DEL_CONSTRUCT_POST.test(post))return;iv.push([a,b]);});
 reHits(DEL_CALL,t).forEach(x=>iv.push(x));
 reHits(LOSS_WORD,t).forEach(([a,b])=>{const w=t.slice(a,b),pre=t.slice(Math.max(0,a-30),a),post=t.slice(b,b+50);
   if(w==='LOH'||LOSS_ALWAYS_POST.test(post)||(LOSS_PRE_GENE.test(pre)&&!LOSS_NOT_GENE.test(pre))||LOSS_PRE_GENOMIC.test(pre)||(LOSS_POST_GENE.test(post)&&!LOSS_SIGNAL_POST.test(post)&&!LOSS_NOT_GENE_POST.test(post)))iv.push([a,b]);});
 return iv;}
function delMatcher(){return {test:t=>delHits(t).length>0,hits:delHits};}
// "Enrich for truncation": truncating variants -- the protein cut short.
//   always: truncating, frameshift / frame-shift, nonsense, stop-gain, premature stop / termination
//   "truncation(s)": as a mutation -- a gene symbol before it, "premature" before it, or a
//     mutation word after it; "truncated" only with a mutation word after it ("truncated form
//     of", "truncated PTCH1" constructs are engineered proteins)
//   splice: only a splice-SITE variant ("splice-site", "splice mutation", "splicing mutation",
//     "NF2 splice-site variants"); alternative splicing and splicing-factor mutations (SF3B1, U2AF1) are
//     biology, not a truncation of the gene in hand
//   notation: Q55* / R213Ter / P191fs; never "X", which in G719X / C797X means ANY amino acid
const TRU_ALWAYS=/\btruncating\b|\bframe[- ]?shifts?\b|\bnonsense\b|\bstop[- ]?gain\w*|\bpremature\s+(?:stop|termination)(?:\s+codons?)?/gi;
const TRU_NOUN=/\btruncat(?:ions?|ed)\b/gi;
// (not "mutants" / "variant": "TLOC1 truncation mutants", "CAD-truncated variant" are constructs)
const TRU_MUT_POST=/^\s*(?:mutations?|alleles?|alterations?|events?)\b/i;
const TRU_PRE_GENE=new RegExp(GENE_TOKEN+String.raw`\s*-?\s*$|\bpremature\s*$`);
// ("splice variant" is left out: "AR splice variant AR-V7", "a shorter splice variant of TLOC1"
// name protein isoforms, not mutations)
const TRU_SPLICE=/\bsplice[- ]site\w*|\bsplic(?:e|ing)\s+(?:site\s+)?(?:mutations?|alterations?)/gi;
const TRU_SPLICE_FACTOR=/^\s*(?:factor|machinery|regulator)/i;
const TRU_NOTATION=/\b[A-Z]\d{1,4}(?:\*|Ter\b|[A-Z]?fs\b(?:\*\d+)?)/g;
let TRUS=false;
function truHits(t){t=t||'';const iv=reHits(TRU_ALWAYS,t);
 reHits(TRU_NOUN,t).forEach(([a,b])=>{const w=t.slice(a,b).toLowerCase(),pre=t.slice(Math.max(0,a-30),a),post=t.slice(b,b+30);
   if(/^\s*(?:mutants?|constructs?|forms?|proteins?)\b/i.test(post))return;   // engineered, not found in a tumour
   if(TRU_MUT_POST.test(post)||(w!=='truncated'&&TRU_PRE_GENE.test(pre)))iv.push([a,b]);});
 reHits(TRU_SPLICE,t).forEach(([a,b])=>{if(!TRU_SPLICE_FACTOR.test(t.slice(b,b+20)))iv.push([a,b]);});
 reHits(TRU_NOTATION,t).forEach(x=>iv.push(x));
 return iv;}
function truMatcher(){return {test:t=>truHits(t).length>0,hits:truHits};}
// "Enrich for exon variants": a variant named by its exon.
//   spelled out: "exon(s) N" (with ranges / lists: "exons 19 and 20") and a variant word within
//     two words after it -- deletion / del / insertion / ins / skipping / duplication / dup /
//     mutation / mutant / mutated / alteration / variant / substitution / hotspot, or a protein
//     change ("exon 21 L858R", "exon 18 (D842V)") -- so "exon 9 PIK3CA mutations" counts and
//     "in exons 1 and 2" does not
//   compact: a gene glued to the exon ("METex14", "METΔex14", "ERBB2ΔEx16", "EGFRex20ins",
//     "PIK3CAEx9"), an exon with a variant suffix ("ex19del", "Ex20Ins", "exon19del"), or a
//     one-/two-digit "ex20" right after a gene symbol ("EGFR ex20") -- but not a breakpoint range
//     ("FGFR3-TACC3 Ex14-16"), and never 3+ digits ("EX527" is a sirtuin inhibitor)
const EXO_SPELLED=/\bexons?\s*\d{1,2}(?:\s*(?:-|–|,|\/|and|or|to)\s*(?:exons?\s*)?\d{1,2})*/gi;
const EXO_VARIANT_POST=/^\s*(?:\S+\s+){0,2}?\(?(?:deletions?|del\b|insertions?|ins\b|skipping|skip\b|duplications?|dup\b|mutations?|mutants?|mutated|alterations?|altered|variants?|substitutions?|hotspot|[A-Z]\d+[A-Z]\b)/i;
const EXO_GLUED=/\b[A-Z][A-Z0-9]{1,9}-?(?:Δ|Delta)?[Ee]x(?:on)?\d{1,2}(?!\d)(?:del|ins|dup|skip|mut)?\w*/g;
const EXO_SUFFIX=/(?:^|[^A-Za-z])(?:Δ)?[Ee][Xx](?:on)?\d{1,2}(?:[Dd]el|[Ii]ns|[Dd]up|[Ss]kip|[Mm]ut)\w*/g;
const EXO_BARE=/\b[Ee][Xx](?:on)?\d{1,2}(?!\d)/g;
const EXO_GENE_PRE=new RegExp(GENE_TOKEN+String.raw`\s*/?\s*(?:`+GENE_TOKEN+String.raw`\s+)?$`);
let EXOS=false;
function exoHits(t){t=t||'';const iv=[];
 // spelled out: highlight through the variant word ("exon 19 deletion", not just "exon 19")
 reHits(EXO_SPELLED,t).forEach(([a,b])=>{const m=EXO_VARIANT_POST.exec(t.slice(b,b+50));if(m)iv.push([a,b+m[0].length]);});
 const range=b=>/^\s*[-–]\s*\d/.test(t.slice(b,b+4));   // "Ex14-16": a fusion breakpoint span
 // a HYPHEN-glued exon needs a variant suffix: "ERBB2-Ex20mut" counts, "FLI1-EX9" (an editing
 // construct) does not; glued without a hyphen ("METex14", "PIK3CAEx9") counts as it stands
 reHits(EXO_GLUED,t).forEach(([a,b])=>{const w=t.slice(a,b);if(range(b))return;
   if(/-(?:Δ|Delta)?[Ee]x/.test(w)&&!/\d(?:del|ins|dup|skip|mut)/i.test(w))return;iv.push([a,b]);});
 reHits(EXO_SUFFIX,t).forEach(([a,b])=>{while(a<b&&!/[ΔEe]/.test(t[a]))a++;iv.push([a,b]);});
 reHits(EXO_BARE,t).forEach(([a,b])=>{if(range(b))return;
   // the gene must stand apart ("EGFR ex20"); hyphened on ("FLI1-EX9") it is a construct name
   if((t[a-1]!=='-'&&EXO_GENE_PRE.test(t.slice(Math.max(0,a-20),a)))||EXO_VARIANT_POST.test(t.slice(b,b+50)))iv.push([a,b]);});
 return iv;}
function exoMatcher(){return {test:t=>exoHits(t).length>0,hits:exoHits};}
// "Enrich for copy number": copy-number change of the tumour genome.
//   "copy number(s)" / "copy-number" in any phrasing, highlighted with a following gain / loss /
//     alteration / variation / change / aberration / amplification / deletion; CNV / CNA / SCNA
//     (case-sensitive) and "CN gain / loss"; chromosome-arm or whole-chromosome gains and losses
//     ("1q gain", "17q12 amplification", "gain of 1q", "chromosome 7 gain", "loss of chromosome
//     10"); "copy-neutral" (LOH)
//   not: mtDNA / mitochondrial DNA / plasmid / viral / transgene copy number
const CN_WORD=/\bcopy[- ]?numbers?(?:[- ](?:gains?|loss(?:es)?|alterations?|variations?|variants?|changes?|aberrations?|amplifications?|deletions?|status|profiles?))?|\bcopy[- ]neutral\b/gi;
const CN_NOT_PRE=/\b(?:mtDNA|mitochondrial\s+DNA|mitochondrial|plasmid|viral|virus|transgene|vector)\s*$/i;
const CN_ABBR=/\b(?:S?CNAs?|CNVs?|CN[- ](?:gains?|loss(?:es)?|amplifications?))\b/g;
const CN_ARM=/\b(?:[1-9]|1\d|2[0-2]|X|Y)[pq](?:\d+(?:\.\d+)?)?\s*-?\s*(?:gains?|loss(?:es)?|amplifications?|deletions?)\b|\b(?:gains?|loss(?:es)?)\s+of\s+(?:chromosomes?\s+)?(?:[1-9]|1\d|2[0-2]|X|Y)[pq]?\b|\bchromosomes?\s+(?:[1-9]|1\d|2[0-2]|X|Y)[pq]?\s+(?:gains?|loss(?:es)?)\b|\bchromosom(?:e|al)\s+(?:gains?|loss(?:es)?)\b/gi;
let CNS=false;
function cnHits(t){t=t||'';const iv=[];
 reHits(CN_WORD,t).forEach(([a,b])=>{if(!CN_NOT_PRE.test(t.slice(Math.max(0,a-30),a)))iv.push([a,b]);});
 reHits(CN_ABBR,t).forEach(x=>iv.push(x));
 reHits(CN_ARM,t).forEach(x=>iv.push(x));
 return iv;}
function cnMatcher(){return {test:t=>cnHits(t).length>0,hits:cnHits};}
// "Enrich for methylation": DNA methylation and epigenetic silencing of the tumour genome.
//   strong (count as they stand): hyper- / hypo- / de-methylat*, "DNA / promoter / CpG /
//     cytosine / gene-body / global methylation", (un / differentially) methylated CpGs or genes,
//     methylome, DNA methyltransferase, CpG island, epigenetic(ally) silenc*
//   bare "methylation / methylated": only if the sentence also speaks of DNA -- promoter, CpG,
//     DNA, MGMT, bisulfite, 5-aza / azacitidine / decitabine
//   never near an RNA, histone or protein cue: m6A / m5C / RNA, histone / H3K27 / H4K20 / H3R,
//     arginine / lysine, PRMT / EZH2 / SETD / SMYD / KMT-mediated, "protein methylation",
//     "p53 methylation" (protein); "methylcellulose" never matches (only methylat* is read)
const MET_STRONG=/\b(?:hyper|hypo|de)-?methylat\w*|\b(?:DNA|promoter|CpG|cytosine|gene[- ]body|global|genome[- ]wide)\s+(?:hyper|hypo)?methylat\w*|\b(?:un|differentially\s+)methylated\b|\bmethylomes?\b|\bDNA\s+methyltransferases?\b|\bCpG\s+islands?\b|\bepigenetic(?:ally)?\s+silenc\w*/gi;
const MET_BARE=/\bmethylat(?:ion|ed)\b/gi;
const MET_DNA_CUE=/\bpromoters?\b|\bCpG|\bc?t?DNA\b|\bcfDNA\b|\bMGMT\b|bisulfite|\b5-aza|azacitidine|decitabine/i;
// a bare methylation also counts when it is plainly about a gene's DNA: a biomarker word after it
// ("methylation status / levels / markers / pattern"), "gene methylation", a gene symbol right
// before it ("TBX15 methylation", "Shox2 and Rassf1a methylation"), or "methylation of (the) GENE"
const MET_QUAL_POST=/^\s*(?:status|levels?|markers?|biomarkers?|patterns?|profil\w*|signatures?|analys[ie]s|-?driven)\b/i;
// (a trailing "p" is the promoter: "TERTp methylation")
const MET_GENE_PRE=new RegExp(String.raw`(?:\bgenes?|\bdifferential|`+GENE_TOKEN+String.raw`p?|\b[A-Z][a-z]+[0-9]+[a-z]?)\s*$`);
// "methylation of (the) GENE", "Methylated Shox2"
const MET_GENE_POST=new RegExp(String.raw`^\s*(?:(?:levels?\s+)?of\s+(?:the\s+)?(?:specific\s+)?)?(?:`+GENE_TOKEN+String.raw`|[A-Z][a-z]+[0-9]+[a-z]?)\b`);
// RNA / histone / protein methylation: the RNA marks are case-SENSITIVE ("m6A"; with /i, "me"
// and "my" matched too), the rest case-insensitive
// (an RNA MARK, not the word RNA: "RNA abundance and DNA methylation" is DNA methylation)
const MET_NOT_RNA=/\bm[1-9][A-Z]\b|2['’‘]?[-‐]?O[-‐]methyl|\bm?RNA\s+methylat/;
const MET_NOT_CTX=/histone|\bH[234][A-Z]?\d*[KR]\d+|arginine|lysine|\bPRMT\d*|\bEZH2\b|\bPRC2\b|\bSETD\w*|\bSET7|\bSMYD\d*|\bKMT\w*-mediated|protein\s+methyl|\bp53\s+methylation|\bSTAT3\s+methylation|methylat\w*\s+STAT3/i;
let METS=false;
function metHits(t){t=t||'';const iv=[];
 const near=(a,b)=>{const w=t.slice(Math.max(0,a-60),Math.min(t.length,b+25));return MET_NOT_CTX.test(w)||MET_NOT_RNA.test(w);};
 reHits(MET_STRONG,t).forEach(([a,b])=>{if(!near(a,b))iv.push([a,b]);});
 const dna=MET_DNA_CUE.test(t);
 reHits(MET_BARE,t).forEach(([a,b])=>{if(near(a,b))return;
   const pre=t.slice(Math.max(0,a-30),a),post=t.slice(b,b+40);
   if(dna||MET_QUAL_POST.test(post)||MET_GENE_PRE.test(pre)||MET_GENE_POST.test(post))iv.push([a,b]);});
 return iv;}
function metMatcher(){return {test:t=>metHits(t).length>0,hits:metHits};}
// ---- "Enrich for therapy": a second box. Its buttons combine in OR among themselves, and the
// box combines in AND with the other boxes: Drivers + Resistance keeps a sentence that names
// a point mutation AND speaks of resistance (T790M confers resistance to gefitinib).
// Resistance: resist* (resistance, resistant, TKI-resistant, multidrug-resistant), chemo- /
// radio-resist*, refractory. Not counted: a cell surviving something that is not a treatment
// (apoptosis / anoikis / ferroptosis / cell-death / oxidative-stress resistance, insulin
// resistance, antibiotic selection markers) -- unless the sentence also names a treatment
// ("resistance to apoptosis induced by EGFR-TKIs" stays).
// no word boundary demanded where PDF text glues words (erlotinibresistant, radioresistance48)
const RES_RE=/(?:\bchemo-?|\bradio-?)?resist(?:ance|ances|ant|ants|ed|ing|s)?(?![a-z])|(?:\bchemo-?)?refractor(?:y|iness)(?![a-z])/gi;
const RES_NOT_PRE=/(?:apopto\w*|anoikis|ferropto\w*|pyropto\w*|cell[- ]death|death|oxidative[- ]stress|stress|insulin|leptin|antibiotic|puromycin|neomycin|G418|hygromycin|blasticidin|zeocin|ampicillin|kanamycin|colonization|vascular|electrical)[- ]?$/i;
const RES_NOT_POST=/^\s+(?:to\s+)?(?:the\s+)?(?:apopto\w*|anoikis|ferropto\w*|pyropto\w*|cell[- ]death|death|oxidative|stress|insulin|leptin|puromycin|neomycin|G418|hygromycin|blasticidin|zeocin|ampicillin|kanamycin|antibiotics?|training|exercise)\b/i;
const THER_CUE=/\b(?:inhibit(?:or|ors|ion)|TKIs?|therap\w*|chemo\w*|drugs?|treat\w*|blockade|agents?|radiation|radiotherapy|irradiation|\w+(?:tinib|rafenib|metinib|mab|platin|taxel|rubicin|ciclib|lisib|parib|zomib|mustine|fen|rolimus|tecan|citabine|uracil))\b/i;
function resNotTherapy(t,a,b){return !THER_CUE.test(t)&&(RES_NOT_PRE.test(t.slice(Math.max(0,a-25),a))||RES_NOT_POST.test(t.slice(b,b+40)));}
let RESS=false;
function resHits(t){t=t||'';return reHits(RES_RE,t).filter(([a,b])=>!resNotTherapy(t,a,b));}
function resMatcher(){return {test:t=>resHits(t).length>0,hits:resHits};}
// Acquired resistance: how resistance arises, a narrower view than Resistance --
// acquired / secondary / intrinsic / primary / adaptive / innate / de novo / emerging resistance
// (one word may sit between: "acquired gefitinib-resistant cells"), acquisition or development of
// resistance, mechanisms of resistance, resistance mechanisms / mutations, gatekeeper mutations,
// and bypass signalling in a sentence that speaks of resistance or an inhibitor.
const ACQ_RE=new RegExp([
 String.raw`\b(?:acquired|secondary|intrinsic|primary|adaptive|innate|de[- ]novo|emerg\w*|on-target|off-target)\s+(?:[\w/-]+[- ])?(?:chemo-?|radio-?)?resist\w*`,
 String.raw`\bacquisition\s+of\s+(?:[\w/-]+[- ]){0,2}resist\w*`,
 String.raw`\bdevelop\w*\s+(?:[\w/-]+[- ]){0,3}resist\w*`,
 String.raw`\bmechanisms?\s+(?:of|underlying|for|driving|behind|conferring)\s+(?:[\w/-]+[- ]){0,3}resist\w*`,
 String.raw`\bresist(?:ance|ant)[- ](?:mechanisms?|mutations?|alleles?|variants?|clones?)\b`,
 String.raw`\bgatekeeper\s+(?:mutations?|residues?|variants?)\b`].join('|'),'gi');
const BYPASS_RE=/\bbypass(?:es|ed|ing)?\b/gi, BYPASS_CTX=/resist|\bTKIs?\b|inhibit/i;
let ACQS=false;
function acqHits(t){t=t||'';
 const h=reHits(ACQ_RE,t).filter(([a,b])=>!resNotTherapy(t,a,b));
 return BYPASS_CTX.test(t)?h.concat(reHits(BYPASS_RE,t)):h;}
function acqMatcher(){return {test:t=>acqHits(t).length>0,hits:acqHits};}
// Synthetic lethality: synthetic(ally) lethal*, synthetic sick*, collateral lethality, and the
// abbreviation SL as a pairing (SL interaction / pair / partner / gene / effect, MYC-SL) --
// case-sensitive, so the compound SL-176 does not count.
const SL_RE=/\bsynthetic(?:al(?:ly)?)?[- ]?(?:lethal\w*|sick\w*)|\bcollateral[- ]lethal\w*|\bselective(?:ly)?[- ]lethal\w*/gi;
const SL_ABBR=/\bSL(?=[- ](?:interactions?|pairs?|partners?|genes?|effects?|vulnerabilit\w*|relationships?|screens?|targets?)\b)|\b[A-Z][A-Za-z0-9]*-SL\b(?!-)/g;
let SLS=false;
function slHits(t){t=t||'';return reHits(SL_RE,t).concat(reHits(SL_ABBR,t));}
function slMatcher(){return {test:t=>slHits(t).length>0,hits:slHits};}
// Biomarker: biomarker(s), predictive marker, companion diagnostic, and a feature said to predict
// response / sensitivity / benefit / resistance / efficacy ("PBRM1 inactivation may predict
// benefit from anti-PD-1", "predictive of everolimus response"). Prognostic biomarkers count too.
const BMK_RE=/\bbio-?markers?\b|\bpredictive\s+(?:bio-?)?markers?\b|\bcompanion\s+diagnostics?\b/gi;
const BMK_PRED=/\bpredict(?:s|ed|ive|ing|ors?)?\s+(?:of\s+)?(?:[\w/-]+\s+){0,2}(?:response|responsiveness|sensitivity|benefit|resistance|efficacy)\b/gi;
let BMKS=false;
function bmkHits(t){t=t||'';return reHits(BMK_RE,t).concat(reHits(BMK_PRED,t));}
function bmkMatcher(){return {test:t=>bmkHits(t).length>0,hits:bmkHits};}
// ---- "Enrich for clinical relevance": a third box, OR inside, AND with the other two boxes.
// Prognosis: prognos* (prognosis, prognostic), named survival endpoints (overall / progression-
// free / disease-free / relapse-free / event-free / metastasis-free survival), survival or outcome
// qualified as a patient result (poor / shorter / median ... survival, adverse outcome, survival
// rate / time / benefit), hazard ratios and Kaplan-Meier; the abbreviations PFS / DFS / RFS / EFS /
// MFS case-sensitive, and OS only beside a survival cue (OS is also osteosarcoma). Cell survival
// ("promotes cell survival", "survival of tumour cells") is biology, not prognosis, and does not count.
const PRG_RE=new RegExp([
 String.raw`\bprognos\w*`,
 String.raw`\b(?:overall|progression[- ]free|disease[- ]free|relapse[- ]free|recurrence[- ]free|event[- ]free|metastasis[- ]free|distant[- ]metastasis[- ]free|cancer[- ]specific|disease[- ]specific|failure[- ]free|leukemia[- ]free|leukaemia[- ]free)\s+survival\b`,
 String.raw`\b(?:poor(?:er|est)?|worse|shorter|longer|reduced|improved|better|inferior|superior|prolonged|median|decreased|increased|lower|higher|unfavou?rable|favou?rable|patients?'?|5-year|five-year|3-year|10-year)\s+(?:(?:overall|clinical|patient|disease)\s+)?(?:survival|outcomes?)\b`,
 String.raw`\b(?:adverse|dismal|clinical|patient|treatment)\s+outcomes?\b`,
 String.raw`\bsurvival\s+(?:rates?|times?|benefit|outcomes?|analys[ie]s|curves?|probability|advantage)\b`,
 String.raw`\bhazard\s+ratios?\b|\bKaplan[- ]Meier\b`].join('|'),'gi');
const PRG_ABBR=/\b(?:PFS|DFS|RFS|EFS|MFS|DMFS|CSS)\b|\bHR\s*[=:,(]?\s*\d+\.\d/g;
const PRG_OS=/\bOS\b/g, PRG_OS_PRE=/\b(?:poor(?:er|est)?|worse|shorter|longer|inferior|superior|improved|better|prolonged|median|reduced|decreased|increased|favou?rable|unfavou?rable)\s+$/i, PRG_OS_CUE=/surviv|\bPFS\b|\bDFS\b|hazard|\bHR\s*[=:,(]?\s*\d|months|Kaplan|prognos|median/i;
const CELL_SURV=/\b(?:cell|cellular|clonogenic|tumou?r[- ]cell|neuronal)\s+$/i;
let PRGS=false;
function prgHits(t){t=t||'';
 const h=reHits(PRG_RE,t).filter(([a,b])=>{const w=t.slice(a,b);
   // "increased survival" / "reduced survival rate" of cells: drop when "cell(s)" sits right before
   return !(/surviv/i.test(w)&&!/free|overall|specific|patient|median|year/i.test(w)&&(CELL_SURV.test(t.slice(Math.max(0,a-20),a))||/^\s+of\s+(?:\w+\s+){0,2}cells?\b/i.test(t.slice(b,b+40))));});
 const cue=PRG_OS_CUE.test(t);
 return h.concat(reHits(PRG_ABBR,t),reHits(PRG_OS,t).filter(([a])=>cue||PRG_OS_PRE.test(t.slice(Math.max(0,a-15),a))));}
function prgMatcher(){return {test:t=>prgHits(t).length>0,hits:prgHits};}
// Clinical evidence: patients in a study or on a treatment -- clinical trial(s), a trial phase
// (phase I / II / III, 1-3, Ib, I/II; only beside trial / study / clinical / randomised wording,
// so a cell-cycle "G2/M phase" or "S phase" does not count), randomised / prospective /
// retrospective studies, NCT numbers, response rates (ORR, objective / overall response rate,
// complete / partial response), and patients treated with / receiving / given a drug.
const CLN_RE=new RegExp([
 String.raw`\bclinical\s+(?:trials?|studies|study|testing|evaluation|responses?|benefit|efficacy|activity|data|setting|practice|use)\b`,
 String.raw`\b(?:randomi[sz]ed|open[- ]label|pivotal|multicent(?:er|re)|single[- ]arm|first[- ]in[- ]human|dose[- ]escalation)\b`,
 String.raw`\b(?:prospective|retrospective)\s+(?:\w+\s+)?(?:study|studies|analys[ie]s|cohorts?|trials?|series|review)\b`,
 String.raw`\bphase\s*(?:I{1,3}|[1-3])(?:[ab]|\s*\/\s*(?:I{1,3}|[1-3])[ab]?)?\s+(?:(?:\w+[- ])?\w+\s+)?(?:trials?|stud(?:y|ies)|clinical|testing|cohorts?|portion|part|expansion|dose)`,
 String.raw`\b(?:objective|overall)\s+response\s+rates?\b|\b(?:complete|partial)\s+(?:responses?|remissions?)\b|\bORR\b|\bNCT\d{6,}`,
 String.raw`\bpatients?\s+(?:\w+\s+){0,2}(?:treated|receiving|received|given|enrolled|randomi[sz]ed)\b`].join('|'),'gi');
let CLNS=false;
function clnHits(t){return reHits(CLN_RE,t||'');}
function clnMatcher(){return {test:t=>clnHits(t).length>0,hits:clnHits};}
// In vivo: in vivo, xenografts (PDX, patient-derived / cell-line-derived xenografts), allografts,
// orthotopic / syngeneic / transgenic / knock-in / nude / tumour-bearing models, mice and rats
// (plural: the animals in an experiment), mouse / murine / animal models, zebrafish models.
// Not counted: "mouse" as a source of cells or reagents (mouse embryonic fibroblasts, mouse
// anti-), and patient-derived organoids or cells, which are in vitro.
const VIV_RE=new RegExp([
 String.raw`\bin[- ]vivo(?![a-z])`,   // also "in vivo32", a citation glued on
 String.raw`\b(?:mouse|murine)\s+survival\b`,
 String.raw`\b(?:xeno|allo)[- ]?grafts?\w*|\bPDXs?\d*\b|\bCDXs?\d*\b|\bGEMMs?\b`,
 String.raw`\b(?:orthotopic|syngeneic|subcutaneous|intracranial)\s+(?:\w+\s+)?(?:models?|tumou?rs?|implantation|injection|xenografts?|mice)\b`,
 String.raw`\b(?:mice|rats|animals)\b`,
 String.raw`\b(?:mouse|murine|animal|rat|zebrafish|rodent)\s+(?:[\w-]+\s+){0,2}(?:models?|xenografts?|experiments?|studies)\b`,
 // tissue or tumours taken from the animal -- not "mouse cells" / "mouse cell lines", which are in vitro
 String.raw`\b(?:mouse|murine)\s+(?:[\w-]+\s+)?(?:tumou?rs?|lungs?|livers?|pancreas|tissues?|lymphomas|melanomas|carcinomas|HCCs?|PCa)\b(?!\s+cells?)`,
 String.raw`\btumou?r[- ]bearing\b|\bknock[- ]?in\s+mice\b|\btransgenic\s+(?:mice|mouse|animals?|models?)\b`].join('|'),'gi');
let VIVS=false;
function vivHits(t){return reHits(VIV_RE,t||'');}
function vivMatcher(){return {test:t=>vivHits(t).length>0,hits:vivHits};}
// ---- "Enrich for gene role and biology": a fourth box, OR inside, AND with the other boxes.
// Oncogene / Tumor suppressor: not sentence wording but the GENE -- an edge counts when either
// end is a gene the OncoKB Cancer Gene List classes as an oncogene (ONCOGENE or
// ONCOGENE_AND_TSG) / a tumor suppressor (TSG or ONCOGENE_AND_TSG). The role rides on each gene
// node ("role": onc / tsg / both), baked in by _cancer_gene_roles(). Every sentence of such an
// edge passes, so the matcher takes the edge as a second argument and highlights nothing.
// (They once matched the words "oncogene" / "tumor suppressor", which few sentences use: the
// edge's genes say more than the sentence's vocabulary.)
const ROLE={};DATA.nodes.forEach(n=>{if(n.role)ROLE[n.id]=n.role;});
function hasRole(id,r){const x=ROLE[id];return x===r||x==='both';}
let ONCS=false, TSGS=false;
function roleMatcher(r){return {test:(t,e)=>!!e&&(hasRole(e.from,r)||hasRole(e.to,r)),hits:()=>[]};}
function oncMatcher(){return roleMatcher('onc');}
function tsgMatcher(){return roleMatcher('tsg');}
// Immune: immune (immune checkpoint, immune evasion, immune cells), immunity, immunotherapy,
// immunosuppressive, immunogenic(ity), immunomodulatory, immunocompetent, immunosurveillance;
// checkpoint inhibitors / blockade; macrophages, dendritic cells, MDSCs, neoantigens, interferon
// gamma; and, case-sensitive, PD-1 / PD-L1 / PD-L2 / CTLA-4 / LAG-3 / TIM-3 / TIGIT, ICI / ICB,
// CAR-T, TILs, TAMs, Tregs, T / NK cells, CD8+ / CD4+, IFN-g. Not counted: the lab methods named
// immuno- (immunoblot, immunohistochemistry, immunofluorescence, immunoprecipitation,
// immunostaining, immunoreactivity, immunoglobulin) -- the word list simply leaves them out --
// cell-cycle and DNA-damage checkpoints (G2/M checkpoint, checkpoint kinase), T- / NK-cell
// lymphomas and leukaemias (diseases), and a singular "TAM" (also tamoxifen).
const IMM_WORD=/\bimmun(?:e|ity|ities|otherap\w*|osuppress\w*|ogenic\w*|omodulat\w*|ocompetent|odeficient|ological|oinhibitory|oediting|oevasion|o-?oncology|ochemotherapy|osurveillance)(?![a-z])|\bcheckpoint\s+(?:inhibit\w*|blockade|therap\w*|immunotherap\w*)|\bmacrophages?\b|\bneoantigens?\b|\bdendritic[- ]cells?\b|\bmyeloid[- ]derived\s+suppressor\s+cells?\b|\binterferon[- ]gamma\b/gi;
const IMM_ABBR=/\bPD-?L?[12](?![0-9])|\bCTLA-?4\b|\bLAG-?3\b|\bTIM-?3\b|\bTIGIT\b|\bIC[IB]s?\b|\bCAR[- ]?T\b|\bTILs?\b|\bTAMs\b|\bTregs?\b|\b(?:T|NK)[- ]cells?\b(?![- ](?:lymphoma|leuka?emia|neoplasm))|\bCD[48]\+|\bIFN-?(?:γ|g\b|gamma\b)/g;
let IMMS=false;
function immHits(t){t=t||'';return reHits(IMM_WORD,t).concat(reHits(IMM_ABBR,t));}
function immMatcher(){return {test:t=>immHits(t).length>0,hits:immHits};}
// Metastasis / EMT: metasta* (metastasis, metastases, metastatic, metastasize), EMT and
// epithelial(-to)-mesenchymal transition. Note "metastatic" also marks a disease stage
// (metastatic NSCLC patients); those sentences count too.
// EMT markers count too (mesenchymal markers / phenotype / properties, epithelial markers); the
// mesenchymal GBM subtype does not
const MTS_RE=/\bmetasta\w*|\bepithelial[- ](?:to[- ])?mesenchymal\w*|\bmesenchymal[- ]transition\b|\b(?:mesenchymal|epithelial)\s+(?:markers?|phenotypes?|propert(?:y|ies)|profiles?|features|states?|traits?|characteristics)\b/gi, EMT_RE=/\bEMT\b/g;
let MTSS=false;
function mtsHits(t){t=t||'';return reHits(MTS_RE,t).concat(reHits(EMT_RE,t));}
function mtsMatcher(){return {test:t=>mtsHits(t).length>0,hits:mtsHits};}
// the stats line: each box's switched-on buttons joined by 'or', the boxes by 'and'
function enrichLabel(){const alt=[DRIVERS&&'candidate drivers',FUSIONS&&'fusions',ITDS&&'ITD',OEXS&&'overexpression',AMPS&&'amplification',DELS&&'deletion / loss',TRUS&&'truncation',EXOS&&'exon variants',CNS&&'copy number',METS&&'methylation'].filter(Boolean).join(' or '),
  ther=[RESS&&'resistance',ACQS&&'acquired resistance',SLS&&'synthetic lethality',BMKS&&'biomarker'].filter(Boolean).join(' or '),
  clin=[PRGS&&'prognosis',CLNS&&'clinical evidence',VIVS&&'in vivo'].filter(Boolean).join(' or '),
  bio=[ONCS&&'oncogene',TSGS&&'tumor suppressor',IMMS&&'immune',MTSS&&'metastasis / EMT'].filter(Boolean).join(' or ');
 const g=[alt,ther,clin,bio].filter(Boolean);
 return g.length>1?g.map(x=>'('+x+')').join(' and '):g.join('');}
function andMatchers(a,b){if(!a||!b)return a||b;return {test:(t,e)=>a.test(t,e)&&b.test(t,e),hits:t=>[...a.hits(t),...b.hits(t)]};}
// the two Enrich buttons widen each other: a sentence naming a mutation OR a fusion counts
function orMatchers(a,b){if(!a||!b)return a||b;return {test:(t,e)=>a.test(t,e)||b.test(t,e),hits:t=>[...a.hits(t),...b.hits(t)]};}
let TM=null;   // matcher in force for the current view; used to highlight hits in sentence text
// mark hits on the RAW text (so a query containing <, > or & still highlights), escaping each piece as we go
// The words the text used for a gene node (payload "surf": [[surface, sentences], ...]), e.g.
// CDKN1A <- "p21", "CDKN1A", "p21Cip1"; AR <- "AR", "androgen receptor", "ARQ62L".
const SURF={};DATA.nodes.forEach(n=>{if(n.surf)SURF[n.id]=n.surf;});
function surfLine(id){const v=SURF[id];if(!v||!v.length)return '';
 return 'found in the text as: '+v.map(([w,c])=>w+' ('+c+')').join(', ');}
// spans in t where an endpoint gene is written, by its own surfaces (case-sensitive, not inside
// a longer word -- but a glued surface such as "ARQ62L" is itself one of the surfaces)
function geneIv(t,e){const iv=[];if(!e)return iv;
 [e.from,e.to].forEach(id=>{(SURF[id]||[]).forEach(([w])=>{if(!w||w.length<2)return;let i=t.indexOf(w);
   while(i>=0){const p=t[i-1],q=t[i+w.length];
     if(!(p&&/[A-Za-z0-9]/.test(p))&&!(q&&/[A-Za-z0-9]/.test(q)))iv.push([i,i+w.length]);
     i=t.indexOf(w,i+1);}});});
 return iv;}
// mark hits on the RAW text (so a query containing <, > or & still highlights), escaping each piece
// as we go: <mark> for the text box / Enrich hits, .gm for the gene as the text wrote it (both
// may apply to one stretch)
function hl(t,e){
 t=t||'';
 const iv=TM?TM.hits(t):[], gv=geneIv(t,e);
 if(!iv.length&&!gv.length)return esc(t);
 const cut=new Set([0,t.length]);iv.concat(gv).forEach(([a,b])=>{cut.add(a);cut.add(b);});
 const pts=[...cut].sort((x,y)=>x-y);let out='';
 for(let k=0;k<pts.length-1;k++){const a=pts[k],b=pts[k+1];if(a>=b)continue;
   const m=iv.some(([x,y])=>x<=a&&b<=y), g=gv.some(([x,y])=>x<=a&&b<=y);let piece=esc(t.slice(a,b));
   if(g)piece='<span class=gm>'+piece+'</span>';if(m)piece='<mark>'+piece+'</mark>';out+=piece;}
 return out;
}
// The relation-type boxes filter SENTENCES, not whole edges. An edge's payload "cat" is the
// label that won over the whole file, so once a score/year/text filter narrows the support the
// sentences still in view can all carry a different label -- SOCS1-MALAT1 is 'interacts' over
// its nine sentences, but the three that mention "lung" all say 'activates'. Filtering,
// counting and colouring on the per-sentence label is what keeps the panel honest: "activates
// 0 in view" then means no visible sentence says activates, which is what the reader sees.
// A sentence BOTH models kept carries BioRED's typed label (rc) and PPI's own reading (pc,
// always "interacts" / "not interacts"). They agree that a relation exists, not on its type,
// so under All such a sentence holds BOTH readings and matches a relation box if either is
// ticked: "interacts" then gives the same sentences (hence edges, nodes, PMIDs) under All as
// under PPI model, and "binds" the same under All as under BioRED model. A single-model view
// shows only that model's reading; the tooltip gives both.
function ppiCat(s){return s.pc||((s.rc||'').indexOf('not ')===0?'not interacts':'interacts');}
function sentCats(e,s){const b=s.rc||e.cat;if(sentSrc(e,s)!=='both')return [b];
 if(SRC_MODE==='ppi')return [ppiCat(s)];if(SRC_MODE==='biored')return [b];
 const p=ppiCat(s);return p===b?[b]:[b,p];}
// the reading a sentence is SHOWN as: its first one among the ticked boxes (BioRED's typed
// label first), so an edge drawn under just "interacts" is coloured interacts, not binds
let VIEW_CATS=null;
function sentCat(e,s){const c=sentCats(e,s);if(VIEW_CATS)for(const x of c)if(VIEW_CATS.has(x))return x;return c[0];}
// Same story for the training set: e.src is the union of model roles over ALL the pair's
// triples, so a pair where PPI said "interacts" in one sentence and BioRED said "activates" in
// another is tagged 'both' -- and an edge-level filter then hid its PPI-only sentences from the
// PPI view. Filtering per sentence makes the buttons agree with the relation boxes.
function sentSrc(e,s){return s.sr||e.src||'ppi';}
// a model's button shows EVERY sentence that model kept, corroborated ones included
function srcPass(e,s,mode){if(!mode||mode==='all')return true;const r=sentSrc(e,s);return r===mode||r==='both';}
function visSents(e,conf,lo,hi,tm,cats,mode){return e.sents.filter(s=>s.sc>=conf&&passYear(s.yr,lo,hi)&&(!tm||tm.test(s.text,e))&&(!cats||sentCats(e,s).some(c=>cats.has(c)))&&srcPass(e,s,mode));}
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
function relTag(s){if(!s.rc)return '';const c=sentCat({},s);
 const src=s.sr==='both'?('BioRED: '+esc(s.rc)+' · PPI: '+esc(ppiCat(s))):(s.sr?esc(SRCLAB[s.sr]||s.sr):'');
 return ' <span class=mut style="color:'+(CCOLOR[c]||'#888')+'">['+esc(c)+(s.sp?' ?':'')+']</span>'+(src?' <span class=mut>'+src+'</span>':'');}
function edgeTip(e,vis,cat){const d=document.createElement('div');let h=edgeHead(e,vis,cat);const lim=20;vis.slice(0,lim).forEach(s=>{h+='<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span>'+relTag(s)+' '+hl(s.text,e)+'</div>';});if(vis.length>lim)h+='<div class=more>+'+(vis.length-lim)+' more</div>';d.innerHTML=h;return d;}
// --- the disease keep-list ------------------------------------------------------------------
// A WHITELIST over disease nodes: empty means no constraint (the default), and any tick means
// "these diseases only". It touches nothing else -- a gene-gene edge has no disease endpoint and
// cannot be cut by it -- so this narrows which diseases the picture is allowed to contain,
// rather than narrowing the picture.
const DIS_KEEP=new Set();
// pmid -> Set of the disease nodes that paper's sentences relate to (any score, built once)
let PAPER_DIS=null;
function paperDiseases(){
 if(PAPER_DIS)return PAPER_DIS;
 PAPER_DIS={};
 DATA.edges.forEach(e=>{const dz=[e.from,e.to].filter(x=>KIND[x]==='disease');if(!dz.length)return;
   e.sents.forEach(s=>{const z=PAPER_DIS[s.pmid]||(PAPER_DIS[s.pmid]=new Set());dz.forEach(d=>z.add(d));});});
 return PAPER_DIS;}
function disPass(id){return !DIS_KEEP.size||KIND[id]!=='disease'||DIS_KEEP.has(id);}
// Rebuilt on every redraw from the edges the rest of the controls left, so the menu is always
// the diseases actually in front of you. Two things it must not do: it must not be built from
// its OWN leavings, or ticking one disease would empty the list you are ticking from (so the
// caller passes the pre-filter edges); and it must not drop a ticked disease that the current
// filters have starved, or that tick could never be undone --
// those stay listed at (0).
// Keeping a disease is almost always a question about ITS neighbourhood, so the first tick (the
// keep-list going from empty to non-empty) switches "Only their direct neighbours" on. Unticking
// that box afterwards is respected until the keep-list is emptied again.
function disConnDefault(was){if(was===0&&DIS_KEEP.size>0){const d=document.getElementById('disconn');if(d)d.checked=true;}}
function buildDiseaseFilter(edges){
 const box=document.getElementById('disfilters'); if(!box)return;
 const seen={};
 edges.forEach(o=>[o.e.from,o.e.to].forEach(nd=>{if(KIND[nd]==='disease')seen[nd]=(seen[nd]||0)+1;}));
 const ids=Object.keys(seen);
 DIS_KEEP.forEach(id=>{if(seen[id]===undefined)ids.push(id);});
 const lab=id=>String(labelById[id]||id);
 const byCount=(a,b)=>(seen[b]||0)-(seen[a]||0)||lab(a).localeCompare(lab(b));
 // same-tissue diseases stacked under one heading (the tissue read off the name, as for the
 // canvas stacking); a tissue with a single disease, or none, stays a plain row. Blocks are
 // ordered by their edge total, members inside by their own count.
 const by={};ids.forEach(id=>{const t=tissueOf(lab(id));(by[t||'~'+id]=by[t||'~'+id]||[]).push(id);});
 const blocks=Object.keys(by).map(k=>({tissue:by[k].length>1?k:'',ids:by[k].sort(byCount),n:by[k].reduce((a,id)=>a+(seen[id]||0),0)}));
 blocks.sort((a,b)=>b.n-a.n||lab(a.ids[0]).localeCompare(lab(b.ids[0])));
 const row=(id,ind)=>{const n=seen[id]||0;
   return '<label'+(ind?' class=dism':'')+'><input type=checkbox class=disf value="'+esc(id)+'"'+(DIS_KEEP.has(id)?' checked':'')
     +'> '+esc(lab(id))+' <span class=cnt'+(n?'':' style="color:#b3243b"')+'>('+n+')</span></label>';};
 box.innerHTML=ids.length
   ? blocks.map(bk=>bk.tissue
       ? '<label class=dist><input type=checkbox class=disg data-ids="'+esc(bk.ids.join('|'))+'"> '+esc(bk.tissue)
         +' <span class=cnt>('+bk.ids.length+' diseases, '+bk.n+')</span></label>'+bk.ids.map(id=>row(id,true)).join('')
       : row(bk.ids[0],false)).join('')
   : '<span class="mut">no diseases in view</span>';
 // a heading is ticked when all its diseases are kept, half-ticked when some are
 box.querySelectorAll('.disg').forEach(g=>{const m=g.dataset.ids.split('|'),k=m.filter(id=>DIS_KEEP.has(id)).length;
   g.checked=k===m.length;g.indeterminate=k>0&&k<m.length;
   g.addEventListener('change',()=>{const was=DIS_KEEP.size;m.forEach(id=>{if(g.checked)DIS_KEEP.add(id);else DIS_KEEP.delete(id);});
     disConnDefault(was);build(+document.getElementById('thr').value);});});
 box.querySelectorAll('.disf').forEach(c=>c.addEventListener('change',()=>{
   const was=DIS_KEEP.size;
   if(c.checked)DIS_KEEP.add(c.value);else DIS_KEEP.delete(c.value);
   disConnDefault(was);build(+document.getElementById('thr').value);
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
// Spring rest length from the edge's sentence count in view: the better-supported a link, the
// closer its two ends sit. Log scale, so 1 -> 200, 2 -> 165, 4 -> 130, 8 -> 95, 16 -> 60 and
// a floor of 45 from about 22 sentences on (a hub's hundreds would otherwise stack its partners).
function edgeLength(w){return Math.max(45,Math.min(200,200-35*Math.log2(Math.max(1,w))));}
function arrowScale(np){return Math.min(0.5+0.18*Math.sqrt(Math.max(0,np-1)),1.6);}
// heads the details box whenever it lists sentences -- above the "A -> B" line, not in the
// hover tooltip, which is transient and already framed by the edge you are pointing at
const INFO_HEAD='<div class="ihead">Sentences of interest</div>';
function scaleNode(s){return 6+Math.sqrt(s)*3.4;}
// on-screen px at the 50% slider, by unique sentences in view: log scale, 11px for 1 sentence to
// 20px from 64 on (2 -> 12.5, 8 -> 15.5, 16 -> 17). It once jumped 13 -> 26 between 4 and 5
// sentences and ran to 96, so one ticked box could double half the labels.
function fontSize(c){return 11+9*Math.min(1,Math.log2(Math.max(1,c))/6);}
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
 const f={size:screenPx(base)/Math.max(ZSCALE,1e-3)};   // screen px -> canvas units at this zoom
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
 const conf=activeConf(), cats=activeCats(); VIEW_CATS=cats; const [ylo,yhi]=activeYears(); const mc=activeMinCluster(); const md=activeMinDegree(); const mp=activeMinPub(); FSCALE=activeFontScale();
 // Enrich buttons combine in OR within a box and in AND across boxes: Drivers + Oncogene keeps a
 // sentence that names a point mutation AND speaks of an oncogene
 const txt=activeText(), tx=textMatcher(txt), lens=andMatchers(andMatchers(andMatchers(orMatchers(orMatchers(orMatchers(orMatchers(orMatchers(orMatchers(orMatchers(orMatchers(orMatchers(DRIVERS?driverMatcher():null,FUSIONS?fusionMatcher():null),ITDS?itdMatcher():null),OEXS?oexMatcher():null),AMPS?ampMatcher():null),DELS?delMatcher():null),TRUS?truMatcher():null),EXOS?exoMatcher():null),CNS?cnMatcher():null),METS?metMatcher():null),
   orMatchers(orMatchers(orMatchers(RESS?resMatcher():null,ACQS?acqMatcher():null),SLS?slMatcher():null),BMKS?bmkMatcher():null)),
   orMatchers(orMatchers(PRGS?prgMatcher():null,CLNS?clnMatcher():null),VIVS?vivMatcher():null)),
   orMatchers(orMatchers(orMatchers(ONCS?oncMatcher():null,TSGS?tsgMatcher():null),IMMS?immMatcher():null),MTSS?mtsMatcher():null)), tm=andMatchers(tx,lens); TM=tm;
 const kinds=activeKinds();
 let edges=[];
 // An edge's support is what survives the score, year, relation and source settings, and the
 // thresholds below judge that first. The Enrich lenses come further down and the support
 // thresholds are judged AGAIN on the sentences they leave, so "min 2 sentences" holds for the
 // lens evidence drawn: with the lens threshold-blind, an edge with 4 sentences of which 1 named
 // a mutation was drawn showing that one sentence under a floor of 2. The text box is exempt --
 // it is a search for edges that mention something, judged on the edge's own support (see the
 // lens step below); the min-connections / min-cluster pruning after it judges what is left.
 // Some node types have no relations OF THEIR OWN: the RE routing pairs a disease with a gene
 // or a chemical, never with another disease, so ticking disease alone can never draw an edge.
 // Such a view is left empty and the stats line says why. (It once fell back to drawing those
 // nodes unconnected, but that fallback skipped min-connections, min-cluster and the disease
 // keep-list, so unticking gene and chemical brought back diseases the full view had filtered
 // out. Unticking a node type now only ever removes nodes.) kindOK: any edge fits the ticks.
 let kindOK=false;
 const disConnOn=DIS_KEEP.size>0&&document.getElementById('disconn').checked;
 DATA.edges.forEach(e=>{
   const kf=!kinds||(kinds.has(KIND[e.from])&&kinds.has(KIND[e.to]));
   if(kf)kindOK=true;
   const sup=visSents(e,conf,ylo,yhi,null,cats,SRC_MODE); if(sup.length<thr)return;
   const np=new Set(sup.map(s=>s.pmid)).size;
   if(mp>1&&np<mp)return;
   if(kf)edges.push({e:e,sup:sup,vis:sup,w:sup.length,np:np,cat:viewCat(e,sup)});
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
 // the lens: keep the edges with at least one sentence matching the text box and the Enrich
 // buttons, and show those sentences only. The support thresholds are judged again only for the
 // Enrich buttons -- at least `thr` lens sentences from at least `mp` publications. The text box
 // is a search, not a lens on the evidence: an edge that passes the thresholds on its own support
 // and has one sentence with the word is found (CTNNB1-SNAI1, five sentences, one naming OXCT1,
 // vanished under "OXCT1" at min 2 sentences when the text match was counted against the floor).
 // The late steps -- lens / text, focus re-cut, structure, disease keep-list -- as one function,
 // so they can run twice: "Only their direct neighbours" needs the picture WITHOUT the text box
 // to decide who the neighbours are (see below). `m` is the matcher in force; `menu` rebuilds the
 // Keep diseases list (only the real pass does).
 function shape(es,m,menu){
  if(m)es=es.filter(o=>{const v=o.sup.filter(s=>m.test(s.text,o.e));if(!v.length)return false;
    if(lens){const L=o.sup.filter(s=>lens.test(s.text,o.e));if(L.length<thr)return false;if(mp>1&&new Set(L.map(s=>s.pmid)).size<mp)return false;}
    o.vis=v;o.w=v.length;o.np=new Set(v.map(s=>s.pmid)).size;o.cat=viewCat(o.e,v);return true;});
  // and re-cut the neighbourhood on what the lens left, so a focus view never shows the seed's
  // neighbours to each other with the seed itself missing
  if(m&&focusActive)es=focusKeep(es,focusSeeds,focusHops);
  // Structure comes LAST, on the edges that are actually drawn: min-connections peels the
  // single-link fringe, then min-cluster drops the small components -- so "min cluster size 6"
  // is a statement about the picture, not about some graph behind it. (It once ran before the
  // lens, and typing a word then shattered a 6-node floor into pairs and triples.)
  // Skipped under gene/chem focus, where you asked for a neighborhood and want all of it.
  if(!focusActive){
    // degree floor, peeled until it holds (the md-core): drop the edges of every node with fewer
    // than md links, re-count, repeat. A single pass let a removed node still vote -- a gene with
    // one gene partner and one chemical whose only link it was counted 2, kept its gene edge, and
    // ticking "chemical" changed the gene count with not one chemical drawn.
    if(md>1){
      for(;;){
        const deg={};es.forEach(o=>{deg[o.e.from]=(deg[o.e.from]||0)+1;deg[o.e.to]=(deg[o.e.to]||0)+1;});
        const next=es.filter(o=>deg[o.e.from]>=md&&deg[o.e.to]>=md);
        if(next.length===es.length)break;
        es=next;
      }
    }
    const adj={};
    es.forEach(o=>{(adj[o.e.from]=adj[o.e.from]||[]).push(o.e.to);(adj[o.e.to]=adj[o.e.to]||[]).push(o.e.from);});
    const comp={};let cid=0;
    for(const n in adj){if(comp[n]!==undefined)continue;const sk=[n];comp[n]=cid;while(sk.length){const x=sk.pop();(adj[x]||[]).forEach(y=>{if(comp[y]===undefined){comp[y]=cid;sk.push(y);}});}cid++;}
    const csz={};for(const n in comp)csz[comp[n]]=(csz[comp[n]]||0)+1;
    es=es.filter(o=>csz[comp[o.e.from]]>=mc);
  }
  // the disease keep-list, built from what the controls above left and applied after it, so the
  // menu never filters itself out of existence (see buildDiseaseFilter)
  if(menu)buildDiseaseFilter(es);
  if(DIS_KEEP.size)es=es.filter(o=>disPass(o.e.from)&&disPass(o.e.to));
  return es;
 }
 // Direct neighbours, strictly: the ticked diseases, every node with a DRAWN edge to one of
 // them, and the drawn edges among that set -- under all settings: relation boxes, Enrich
 // buttons, min sentences / publications / connections / cluster size. (It was once chosen
 // relation- and lens-blind, so a gene whose own disease link was hidden stayed on through its
 // links to other neighbours, two steps out; with "associated" unticked that was most of them.)
 // One step, not "connected to": the gene graph is one giant component through hubs like TP53,
 // so walking further reaches nearly everything whatever is ticked. The price: a disease whose
 // links are all hidden (say, only relations it never carries are ticked) stands alone.
 // The text box is the one exception: it SEARCHES the neighbourhood, it does not redraw it. The
 // neighbours are picked from the picture without the text, and the text then keeps the edges
 // among them that match. (Judged after the text, typing OXCT1 under "skin" emptied the canvas:
 // CTNNB1-SNAI1 says OXCT1, but CTNNB1's and SNAI1's own links to melanoma do not, so neither
 // was a neighbour any more and their matching edge went with them.)
 const disConn=disConnOn;
 let near=null;
 if(disConn){
   const base=(tx?shape(edges.map(o=>Object.assign({},o)),lens,false):null);
   edges=shape(edges,tm,true);
   near=new Set(DIS_KEEP);
   (base||edges).forEach(o=>{if(DIS_KEEP.has(o.e.from))near.add(o.e.to);if(DIS_KEEP.has(o.e.to))near.add(o.e.from);});
   edges=edges.filter(o=>near.has(o.e.from)&&near.has(o.e.to));
 } else edges=shape(edges,tm,true);
 // Disease context for the edges that do NOT touch a kept disease (gene-gene, gene-drug): only
 // the sentences from papers about a kept disease. The endpoints being neighbours of glioma says
 // nothing about where each of their sentences comes from, and CTNNB1-SNAI1 under "brain" showed
 // melanoma papers. A paper's diseases are the ones its own sentences relate to anything
 // (paperDiseases); a paper with none -- about a third -- is placed nowhere and drops out of
 // these edges. Edges to a kept disease stay as they are. The edge is then re-judged on what is
 // left, against min sentences / publications like the Enrich buttons.
 if(DIS_KEEP.size){
   const pd=paperDiseases(), inScope=s=>{const z=pd[s.pmid];if(!z)return false;for(const d of z)if(DIS_KEEP.has(d))return true;return false;};
   edges=edges.filter(o=>{
     if(DIS_KEEP.has(o.e.from)||DIS_KEEP.has(o.e.to))return true;
     const v=o.vis.filter(inScope);
     if(!v.length||v.length<thr)return false;
     const np=new Set(v.map(s=>s.pmid)).size; if(mp>1&&np<mp)return false;
     if(v.length!==o.vis.length){o.vis=v;o.w=v.length;o.np=np;o.cat=viewCat(o.e,v);}
     return true;});
 }
 // A kept disease the other settings left without a single link -- typically because only
 // relation types it never carries are ticked (diseases are linked by "associated", almost never
 // by "binds") -- used to vanish, and under "Only their direct neighbours" it took the whole
 // canvas with it, unexplained. It is drawn on its own instead, and the stats line says why.
 const loneDis=[];
 if(DIS_KEEP.size&&(!kinds||kinds.has('disease'))){
   const linked=new Set();edges.forEach(o=>{linked.add(o.e.from);linked.add(o.e.to);});
   DIS_KEEP.forEach(id=>{if(!linked.has(id))loneDis.push(id);});
 }
 updateCatCounts(edges,cats);   // edges is final here (category, score, year, degree, cluster, text, disease)
 const keep=new Set();edges.forEach(o=>{keep.add(o.e.from);keep.add(o.e.to);});
 loneDis.forEach(id=>keep.add(id));
 const nss={};edges.forEach(o=>{o.vis.forEach(s=>{(nss[o.e.from]=nss[o.e.from]||new Set()).add(s.text);(nss[o.e.to]=nss[o.e.to]||new Set()).add(s.text);});});
 const nsz=id=>(nss[id]?nss[id].size:0);
 const nodes=DATA.nodes.filter(n=>keep.has(n.id)).map(n=>{const k=n.kind||'gene';
  const fs=fontSize(nsz(n.id))*(KIND_FS[k]||1);
  return {id:n.id,label:nodeLabel(n),value:nsz(n.id),size:scaleNode(nsz(n.id)),shape:n.shape||'dot',title:n.label+((n.kind&&n.kind!=='gene')?'  ['+n.kind+']':'')+' — '+nsz(n.id)+' unique sentences (in view)'+(n.target?' · drug target: '+n.target+' chemicals'+(n.tcat==='green'?' (approved anti-neoplastic)':(n.tcat==='amber'?' (approved, non-anti-neoplastic)':' (no approved drug: investigational / ChEBI role)')):'')+(SURF[n.id]?String.fromCharCode(10)+surfLine(n.id):''),color:nodeColor(n),_fs:fs,font:nodeFont(k,fs),// the one visual channel still free: fill is drug-target shading, shape is node kind
   borderWidth:1,borderWidthSelected:2};});
 // no `value`: vis would then scale the width itself and ignore edgeWidth()
 // arrowheads only on directed (signed) relations; interacts/associated/binds are symmetric,
 // so an arrow there would only show which entity the sentence happened to mention first
 const eds=edges.map((o,i)=>({id:i,from:o.e.from,to:o.e.to,width:edgeWidth(o.np),length:edgeLength(o.w),color:{color:CCOLOR[o.cat]||o.e.color,opacity:0.6},dashes:o.cat.indexOf('not ')===0,arrows:{to:{enabled:DIRECTED.has(o.cat.replace(/^not /,'')),scaleFactor:arrowScale(o.np)}},title:edgeTip(o.e,o.vis,o.cat)}));
 const vpub=new Set();edges.forEach(o=>o.vis.forEach(s=>vpub.add(s.pmid)));
 const nkinds=new Set(nodes.map(n=>KIND[n.id]));
 // name what is actually on screen: "378 diseases" beats "378 genes" in a disease-only view
 const KLAB={gene:'genes',disease:'diseases',chemical:'chemicals'};
 document.getElementById('stats').innerHTML='Showing <b>'+nodes.length+'</b> '+(nkinds.size===1?(KLAB[[...nkinds][0]]||'nodes'):'nodes').replace(/s$/,nodes.length===1?'':'s')+', <b>'+edges.length+'</b> edges, <b>'+vpub.size+'</b> publications (&ge;'+conf+')'+(txt?' &middot; text: <b>'+esc(txt)+'</b>':'')+(enrichLabel()?' &middot; <b>enriched for '+enrichLabel()+'</b>':'')+(focusActive?' &middot; focus: <b>'+esc(focusLabel)+'</b>':'')+(kinds&&!kindOK?' &middot; <b>nothing relates '+(kinds.size===1?'a '+[...kinds][0]+' to another '+[...kinds][0]:'these node types to each other')+'</b> &mdash; '+(kinds.has('disease')&&kinds.size===1?'diseases link only to genes and chemicals, so tick one of them':'tick another node type'):'')
   +(cats.size?'':' &middot; <b>no relation type is ticked</b>')
   +(loneDis.length?' &middot; <b>'+loneDis.length+'</b> kept disease'+(loneDis.length>1?'s have':' has')+' no link under these settings, drawn on '+(loneDis.length>1?'their':'its')+' own'
     +(cats.has('associated')?'':' &mdash; diseases are linked by <b>associated</b>, so tick it to see '+(loneDis.length>1?'their':'its')+' genes and drugs'):'');
 // seed positions: known nodes where they were (held fixed while physics settles the rest), new
 // nodes at the mean of their already-placed neighbours, nudged apart so they do not coincide
 const known=nodes.filter(n=>LAYOUT[n.id]);
 known.forEach(n=>{n.x=LAYOUT[n.id].x;n.y=LAYOUT[n.id].y;n.fixed={x:true,y:true};});
 const nbr={};eds.forEach(e=>{(nbr[e.from]=nbr[e.from]||[]).push(e.to);(nbr[e.to]=nbr[e.to]||[]).push(e.from);});
 nodes.forEach(n=>{if(LAYOUT[n.id])return;const p=(nbr[n.id]||[]).map(m=>LAYOUT[m]).filter(Boolean);
   if(p.length){const [jx,jy]=idJitter(n.id,40);n.x=p.reduce((a,q)=>a+q.x,0)/p.length+jx;n.y=p.reduce((a,q)=>a+q.y,0)/p.length+jy;}});
 const nNew=nodes.length-known.length;
 // the camera before this rebuild: kept if the new picture still fits inside it (see settled)
 const view=network?{position:network.getViewPosition(),scale:network.getScale()}:null;
 const data={nodes:new vis.DataSet(nodes),edges:new vis.DataSet(eds)};
 NODEDS=data.nodes;
 // springConstant 0.2, not 0.02: each edge's length is its sentence count (edgeLength), and at
 // 0.02 the repulsion swamped the springs -- every link stretched to ~300 whatever its support.
 const options={layout:{improvedLayout:false},physics:{enabled:nNew>0,stabilization:{iterations:known.length?120:200},barnesHut:{gravitationalConstant:-14000,springLength:130,springConstant:0.2,avoidOverlap:0.3}},interaction:{hover:true,tooltipDelay:120},nodes:{shape:'dot',scaling:{min:6,max:60,label:{enabled:false,drawThreshold:1}},font:{color:'rgba(26,26,26,0)'}},edges:{smooth:false,arrowStrikethrough:false,hoverWidth:0,selectionWidth:0,arrows:{to:{enabled:true,scaleFactor:0.6}}}};
 if(network)network.destroy();
 network=new vis.Network(net,data,options);
 LABEL_N=nodes.length; LABEL_A=-1;
 // reshape once physics is off, so nothing pulls the moved nodes back, and fit afterwards so the
 // new positions are inside the view. BASEPOS is captured first: it is the layout to replay from.
 TGROUPS=tissueGroups(nodes);
 const settled=()=>{network.setOptions({physics:false});
   if(known.length)data.nodes.update(known.map(n=>({id:n.id,fixed:false})));   // draggable again
   BASEPOS=network.getPositions();
   Object.assign(LAYOUT,BASEPOS);
   applyLayoutShape();
   // Camera: stay put when every node is still inside the old view; otherwise frame them with
   // safeFit. A plain fit zoomed straight in on a handful of nodes until they filled the screen,
   // and re-fitting on every tick made the view lurch in and out as relation boxes were ticked.
   if(view){network.moveTo({position:view.position,scale:view.scale,animation:false});
     if(!allInView())safeFit(false);}
   else safeFit(false);
   LABEL_A=-1;updateLabels();refloorLabels();network.redraw();};   // fonts for the final camera, at once
 // nothing new to place: no physics runs, so no stabilization event will come
 if(nNew>0)network.on('stabilizationIterationsDone',settled);else settled();
 network.on('zoom',updateLabels);
 network.on('animationFinished',updateLabels);
 fitLeftPanel();    // self-correcting: a panel narrowed by anything recovers on the next redraw
 const _e=edges;
 network.on('click',p=>{const info=document.getElementById('info');
   if(p.nodes.length){const n=DATA.nodes.find(x=>x.id===p.nodes[0]);info.innerHTML='<b>'+n.label+'</b>: '+nsz(n.id)+' unique sentences (in view)'+(SURF[n.id]?'<div class=mut style="margin-top:3px">'+esc(surfLine(n.id))+'</div>':'');}
   else if(p.edges.length){const o=_e[p.edges[0]];info.innerHTML=INFO_HEAD+edgeHead(o.e,o.vis,o.cat)+o.vis.map(s=>'<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span> '+hl(s.text,o.e)+'</div>').join('');}});
}
const thr=document.getElementById('thr');
let CATTOT={};
// One edge counts under EVERY relation its sentences carry (a nine-sentence edge holding both
// 'interacts' and 'activates' readings counts once under each), so the totals add up to more
// than the edge count -- they answer "how many edges can show me this label", which is the
// question the tick boxes and the "in view" half answer too.
function catsOf(e){const c={};e.sents.forEach(s=>sentCats(e,s).forEach(x=>c[x]=1));return Object.keys(c);}
function buildCatFilters(){CATTOT={};DATA.edges.forEach(e=>catsOf(e).forEach(c=>CATTOT[c]=(CATTOT[c]||0)+1));const cats=Object.keys(CATTOT).sort((a,b)=>CATTOT[b]-CATTOT[a]);document.getElementById('catfilters').innerHTML=cats.map(c=>'<label><input type=checkbox class=catf value="'+esc(c)+'"'+(c==='binds'?' checked':'')+'> <span class=sw style="background:'+(CCOLOR[c]||'#888')+'"></span> '+esc(c)+' <span class=cnt data-cat="'+esc(c)+'">('+CATTOT[c]+')</span></label>').join('');document.querySelectorAll('.catf').forEach(c=>c.addEventListener('change',()=>build(+thr.value)));}
// Counts are LIVE: "(total · N in view)" is recomputed from the sentences actually drawn, so
// "N in view" is exactly the number of drawn edges that can show you a sentence tagged with
// that relation -- 0 means no visible sentence carries the label, never "0 but there it is in
// the tooltip". The total is the whole payload (everything qualifying at >=0.5), which is NOT
// what you see: the score slider, year range, text filter, min-publications, min-connections
// and min-cluster size all prune afterwards. A category pruned to nothing reads 0 (in red)
// instead of looking available: ticking it alone would leave the canvas blank, typically
// because its edges form components smaller than "Min cluster size".
function updateCatCounts(edges,cats){const seen={};edges.forEach(o=>{const c={};o.vis.forEach(s=>sentCats(o.e,s).forEach(x=>c[x]=1));for(const k in c)seen[k]=(seen[k]||0)+1;});
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
 if(NODEDS)NODEDS.update(NODEDS.get().map(n=>({id:n.id,font:nodeFont(KIND[n.id]||'gene',n._fs||11)})));});
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
(function(){const c=document.getElementById('disclear');
 if(c)c.addEventListener('click',()=>{if(!DIS_KEEP.size)return;DIS_KEEP.clear();build(+thr.value);});
 const dc=document.getElementById('disconn');
 if(dc)dc.addEventListener('change',()=>{if(DIS_KEEP.size)build(+thr.value);});})();
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
document.getElementById('metbtn').addEventListener('click',ev=>{METS=!METS;ev.currentTarget.classList.toggle('on',METS);ev.currentTarget.setAttribute('aria-pressed',METS);build(+thr.value);});
document.getElementById('resbtn').addEventListener('click',ev=>{RESS=!RESS;ev.currentTarget.classList.toggle('on',RESS);ev.currentTarget.setAttribute('aria-pressed',RESS);build(+thr.value);});
document.getElementById('acqbtn').addEventListener('click',ev=>{ACQS=!ACQS;ev.currentTarget.classList.toggle('on',ACQS);ev.currentTarget.setAttribute('aria-pressed',ACQS);build(+thr.value);});
document.getElementById('slbtn').addEventListener('click',ev=>{SLS=!SLS;ev.currentTarget.classList.toggle('on',SLS);ev.currentTarget.setAttribute('aria-pressed',SLS);build(+thr.value);});
document.getElementById('bmkbtn').addEventListener('click',ev=>{BMKS=!BMKS;ev.currentTarget.classList.toggle('on',BMKS);ev.currentTarget.setAttribute('aria-pressed',BMKS);build(+thr.value);});
document.getElementById('prgbtn').addEventListener('click',ev=>{PRGS=!PRGS;ev.currentTarget.classList.toggle('on',PRGS);ev.currentTarget.setAttribute('aria-pressed',PRGS);build(+thr.value);});
document.getElementById('clnbtn').addEventListener('click',ev=>{CLNS=!CLNS;ev.currentTarget.classList.toggle('on',CLNS);ev.currentTarget.setAttribute('aria-pressed',CLNS);build(+thr.value);});
document.getElementById('vivbtn').addEventListener('click',ev=>{VIVS=!VIVS;ev.currentTarget.classList.toggle('on',VIVS);ev.currentTarget.setAttribute('aria-pressed',VIVS);build(+thr.value);});
document.getElementById('oncbtn').addEventListener('click',ev=>{ONCS=!ONCS;ev.currentTarget.classList.toggle('on',ONCS);ev.currentTarget.setAttribute('aria-pressed',ONCS);build(+thr.value);});
document.getElementById('tsgbtn').addEventListener('click',ev=>{TSGS=!TSGS;ev.currentTarget.classList.toggle('on',TSGS);ev.currentTarget.setAttribute('aria-pressed',TSGS);build(+thr.value);});
document.getElementById('immbtn').addEventListener('click',ev=>{IMMS=!IMMS;ev.currentTarget.classList.toggle('on',IMMS);ev.currentTarget.setAttribute('aria-pressed',IMMS);build(+thr.value);});
document.getElementById('mtsbtn').addEventListener('click',ev=>{MTSS=!MTSS;ev.currentTarget.classList.toggle('on',MTSS);ev.currentTarget.setAttribute('aria-pressed',MTSS);build(+thr.value);});
document.getElementById('cnbtn').addEventListener('click',ev=>{CNS=!CNS;ev.currentTarget.classList.toggle('on',CNS);ev.currentTarget.setAttribute('aria-pressed',CNS);build(+thr.value);});
document.getElementById('exobtn').addEventListener('click',ev=>{EXOS=!EXOS;ev.currentTarget.classList.toggle('on',EXOS);ev.currentTarget.setAttribute('aria-pressed',EXOS);build(+thr.value);});
document.getElementById('trubtn').addEventListener('click',ev=>{TRUS=!TRUS;ev.currentTarget.classList.toggle('on',TRUS);ev.currentTarget.setAttribute('aria-pressed',TRUS);build(+thr.value);});
document.getElementById('delbtn').addEventListener('click',ev=>{DELS=!DELS;ev.currentTarget.classList.toggle('on',DELS);ev.currentTarget.setAttribute('aria-pressed',DELS);build(+thr.value);});
document.getElementById('ampbtn').addEventListener('click',ev=>{AMPS=!AMPS;ev.currentTarget.classList.toggle('on',AMPS);ev.currentTarget.setAttribute('aria-pressed',AMPS);build(+thr.value);});
document.getElementById('oexbtn').addEventListener('click',ev=>{OEXS=!OEXS;ev.currentTarget.classList.toggle('on',OEXS);ev.currentTarget.setAttribute('aria-pressed',OEXS);build(+thr.value);});
document.getElementById('itdbtn').addEventListener('click',ev=>{ITDS=!ITDS;ev.currentTarget.classList.toggle('on',ITDS);ev.currentTarget.setAttribute('aria-pressed',ITDS);build(+thr.value);});
// tick every relation type at once (one rebuild, not one per box)
document.getElementById('catall').addEventListener('click',()=>{document.querySelectorAll('.catf').forEach(c=>c.checked=true);build(+thr.value);});
document.getElementById('fusbtn').addEventListener('click',ev=>{FUSIONS=!FUSIONS;ev.currentTarget.classList.toggle('on',FUSIONS);ev.currentTarget.setAttribute('aria-pressed',FUSIONS);build(+thr.value);});
document.getElementById('drvbtn').addEventListener('click',ev=>{DRIVERS=!DRIVERS;ev.currentTarget.classList.toggle('on',DRIVERS);ev.currentTarget.setAttribute('aria-pressed',DRIVERS);build(+thr.value);});
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
document.getElementById('zfit').addEventListener('click',()=>{if(network)safeFit(true);});
// the year range: dragging either handle relabels it and rebuilds the graph
const yl=document.getElementById('yrlo'),yh=document.getElementById('yrhi');
[yl,yh].forEach(el=>{el.min=MINY;el.max=MAXY;});
yl.value=MINY;yh.value=MAXY;
function updYr(){const a=activeYears();document.getElementById('yrlab').textContent=a[0]+'–'+a[1];}
[yl,yh].forEach(el=>el.addEventListener('input',()=>{updYr();build(+thr.value);}));
updYr();
const chemGenes={};DATA.nodes.forEach(n=>(n.chems||[]).forEach(c=>{(chemGenes[c]=chemGenes[c]||[]).push(n.label);}));const csel=document.getElementById('chemfilter');Object.keys(chemGenes).sort().forEach(c=>{const g=chemGenes[c].slice().sort();const o=document.createElement('option');o.value=c;o.textContent=c+' → '+g.join(', ');csel.appendChild(o);});csel.addEventListener('change',()=>build(+thr.value));
// disease list for the tissue controls; the gene-only graph has none, so the row hides itself
(function(){const ds=DATA.nodes.filter(n=>(n.kind||'gene')==='disease').sort((a,b)=>a.label.localeCompare(b.label));
 const trow=document.getElementById('tissuerow');
 if(!ds.length){if(trow)trow.style.display='none';return;}
 // say up front how much of the disease list the name-based tissue reading actually covers
 const tg={};ds.forEach(n=>{const t=tissueOf(n.label);if(t)(tg[t]=tg[t]||[]).push(n.label);});
 const multi=Object.keys(tg).filter(t=>tg[t].length>1);
 const hit=multi.reduce((s,t)=>s+tg[t].length,0);
 const el=document.getElementById('tissuen');
 if(el)el.textContent=hit+' of '+ds.length+' disease nodes fall into '+multi.length+' tissues; the rest are left where the layout puts them';
 // also a position control: replay from BASEPOS instead of rebuilding the whole network
 document.getElementById('tissuestack').addEventListener('change',applyLayoutShape);
})();
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
// measured once, then locked: the panel keeps one static width whatever the filters do
let LP_W=0;
function fitLeftPanel(){
 const lp=document.getElementById('lpanel'), z=document.getElementById('zoom');
 if(!lp||!z)return;
 if(window.innerWidth<=700||lp.style.display==='none'){lp.style.width='';lp.style.maxWidth='';return;}
 if(LP_W){lp.style.maxWidth='none';lp.style.width=LP_W+'px';return;}
 const b=z.getElementsByTagName?z.getElementsByTagName('button'):[];
 if(!b.length||!b[0].getBoundingClientRect)return;   // nothing measurable: keep the CSS width
 const prevW=lp.style.width, prevM=lp.style.maxWidth;
 lp.style.maxWidth='none';lp.style.width='auto';     // shrink-to-fit: the row cannot wrap here
 const first=b[0].getBoundingClientRect(), last=b[b.length-1].getBoundingClientRect();
 // the panel is as wide as the zoom row, or as the two-column Enrich grid if that is wider
 const eg=document.getElementById('enrichgrid');
 // ...and wide enough for an indented 'non-small cell lung carcinoma (9999)' row in Keep diseases
 // on one line, plus the list's scrollbar: a fixed probe, so the width does not follow the data
 const df=document.getElementById('disfilters');let dw=0;
 if(df){const p=document.createElement('label');p.className='dism';
   p.style.cssText='position:absolute;visibility:hidden;white-space:nowrap;display:inline-block;box-sizing:content-box';
   p.innerHTML='<input type=checkbox> non-small cell lung carcinoma <span class=cnt>(9999)</span>';
   df.appendChild(p);dw=Math.ceil(p.getBoundingClientRect().width)+18;df.removeChild(p);
   const bx=df.closest('.dbox');if(bx){const c=getComputedStyle(bx);dw+=Math.ceil(parseFloat(c.paddingLeft)+parseFloat(c.paddingRight)+parseFloat(c.borderLeftWidth)+parseFloat(c.borderRightWidth));}}
 const gw=[...document.querySelectorAll('.egrid')].reduce((m,g)=>Math.max(m,Math.ceil(g.offsetWidth)),0);   // every Enrich box
 const w=Math.max(Math.ceil(last.right-first.left),gw,dw);   // offsetWidth: the box's padding and border included
 if(Math.abs(last.top-first.top)<2&&w>=FIT_MIN&&w<=FIT_MAX){lp.style.width=w+'px';
   // the reserved scrollbar gutter comes out of the width set here, so add it back: without it
   // the content area is narrower than measured and the Enrich grids' notes spill out of their box
   const cs=getComputedStyle(lp), inner=lp.clientWidth-parseFloat(cs.paddingLeft)-parseFloat(cs.paddingRight);
   LP_W=w+3+Math.max(0,Math.ceil(w-inner));lp.style.width=LP_W+'px';}   // +3px of breathing room
 else {lp.style.width=prevW;lp.style.maxWidth=prevM;}
}
document.querySelectorAll('.egrid button').forEach(b=>{b.dataset.t=b.textContent;});
fitLeftPanel();
window.addEventListener('resize',()=>{fitLeftPanel();if(network)network.redraw();});
document.querySelectorAll('.kindf').forEach(c=>c.addEventListener('change',()=>build(+thr.value)));
// gene-only runs hide every row the left column holds; an empty white box is worse than none
(function(){const lp=document.getElementById('lpanel');
 const kids=lp?Array.prototype.slice.call(lp.children||[]).filter(c=>(c.className||'').indexOf('pclose')<0):[];   // HTMLCollection, not an array; the close button is not content
 if(kids.length&&!kids.some(c=>!c.style||c.style.display!=='none'))lp.style.display='none';})();
// close / reopen each column; the reopen tab only ever appears for a panel the user shut, so a
// gene-only run whose left column hid itself above gets no tab for an empty box
function setPanel(id,open){const p=document.getElementById(id),t=document.getElementById(id==='lpanel'?'lopen':'ropen');
 if(!p||!t)return;p.style.display=open?'':'none';t.style.display=open?'none':'block';
 if(open&&id==='lpanel')fitLeftPanel();}
document.querySelectorAll('.pclose').forEach(b=>b.addEventListener('click',()=>setPanel(b.getAttribute('data-panel'),false)));
document.getElementById('lopen').addEventListener('click',()=>setPanel('lpanel',true));
document.getElementById('ropen').addEventListener('click',()=>setPanel('panel',true));
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


def render_graph(payload, lib, miny, maxy, nxml=None, pubmed_query=""):
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
    # the score slider opens at its lowest stop (SCORE_STEPS[0] = 0.5): these edges sit mostly
    # in 0.5-0.8 (BioRED gene-disease relations), so a high default shows an almost empty canvas
    return (GRAPH_TEMPLATE.replace("__LIBTAG__", libtag)
            .replace("__TITLE__", title)
            .replace("__KINDROW__", KIND_ROW)
            .replace("__CONFDEF__", "0")
            .replace("__PUBMED_QUERY__", qrow)
            .replace("__PAYLOAD__", json.dumps(payload, ensure_ascii=False))
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
    # subject = the agent, not merely the first-mentioned entity ("lactate ... reduced by
    # AZD3965"); a no-op on triples relation_extraction.py already oriented.
    n_flip = sum(D.orient(t) for t in raw)
    print(f"direction: {n_flip:,} signed triple(s) reversed to put the agent first "
          f"(passive/agent cue -- gpu_bundle/direction.py)")
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
        universe = kept if args.score <= GRAPH_BASE else [t for t in d if keep_fn(t, GRAPH_BASE)]
        payload = graph_payload_multi(universe, flags)
        if N_FALSE_DROPPED:
            print(f"  drug targets: dropped {N_FALSE_DROPPED} curated false gene-chemical "
                  f"pair(s) (FALSE_TARGETS)")
        yrs = [s["yr"] for e in payload["edges"] for s in e["sents"] if s.get("yr")]
        miny, maxy = (min(yrs), max(yrs)) if yrs else (2000, 2026)
        lib = get_vis_lib()
        # corpus size = # source documents; experimental_ner/ is empty in the bundle
        # (it was a runtime symlink), so fall back to the per-document sentences/ files.
        nxml = len(list(XML_DIR.glob("*.xml"))) or len(list(SENT_DIR.glob("*.json")))
        GRAPH_OUT.parent.mkdir(parents=True, exist_ok=True)
        pubmed_query = read_pubmed_query()
        GRAPH_OUT.write_text(render_graph(payload, lib, miny, maxy, nxml, pubmed_query),
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
