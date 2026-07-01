#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""relation_extraction.py -- learned RE step that puts a *score* on each triple.

Where triples.py keeps the raw text between two consecutive entities as a
"predicate" iff a regex cue matches (no model, no score), this module runs
BioBERT sequence-classification RE model(s) over entity pairs and emits, for each
pair it keeps, a predicted relation *label* and a softmax *confidence* -- the
per-triple score BioBERT NER itself cannot give you.

METHOD (standard BioBERT RE, entity-marker / "entity blinding" scheme)
----------------------------------------------------------------------
For a sentence with typed entity spans (the same sentences/*.json BioBERT NER
output that triples.py consumes), and for each ordered entity pair (a before b,
by character offset) within a char window:

  1. Blind the two target entities in the sentence with type markers, e.g.
       "TP53 inhibited MDM2 in ..."  ->  "@GENE$ inhibited @GENE$ in ..."
     so the model attends to the pair, not the surface names.
  2. Feed the marked sentence to AutoModelForSequenceClassification (a BioBERT
     checkpoint FINE-TUNED for RE -- ChemProt / GAD / your own; see train_re.py).
  3. softmax over the label set -> (argmax label, raw prob). If a calibrator
     (calibration.json, fit by train_re.py) sits next to the checkpoint, the raw
     prob is mapped to a calibrated `p_rel` (RE_CALIBRATE=0 disables); the raw
     value is kept as predicate.score, the calibrated as predicate.score_calibrated.
  4. COMPOSITE SCORING (triples_strategy.html section 7): the kept `score` fuses
     `p_rel` with independent evidence the relation softmax doesn't see --
        score = p_rel x min(subj,obj NER) x cue_factor x margin_factor
                      x section_factor x self_pair_factor
     NER is near-saturated in this corpus, so it gates rather than discriminates;
     `p_rel` carries the weight. cue_factor rewards a relational cue (sentence
     `relation` cue or a RELCUE stem in the connecting text); margin_factor uses the
     sentence `result_margin`; section_factor weights results/abstract > discussion
     > intro; self_pair_factor penalizes a gene paired with itself unless the
     relation is reflexive-plausible (RE_REFLEXIVE_OK). Each factor is stored in
     `score_components` for ablation.

Pairs whose predicted label is a "no-relation" class (auto-detected, or set via
RE_NEG_LABELS) are dropped, as are pairs whose COMPOSITE score < RE_MIN_SCORE.

POLARITY / MODALITY (triples_strategy.html section 6.1)
------------------------------------------------------
A detector tags each triple with polarity in {positive, negated} and modality in
{asserted, speculated} from grammatical-negation / hedge cues near the pair.
Negation has two backends (RE_NEG_BACKEND): "cuelist" (default, zero-dependency
window) or "negspacy" (NegEx/ConText sentence-scope; `pip install spacy negspacy`,
falls back to cue-list if unavailable). Speculation/modality always uses the
cue-list window (NegEx covers negation only). Negated or speculated relations are
KEPT and flagged ("X does not bind Y" is a fact), not dropped or down-scored; the
tags (+ matched assertion_cues, "negex" when negspaCy fires) support filtering.

TYPE-PAIR ROUTING
-----------------
No single biomedical RE model covers every entity-type pair, so each pair is
routed to the model trained for it (markers are produced per endpoint type by
TYPE_MARKER, which is correct for both ChemProt and GAD):

    {GENETIC, CHEMICAL}  ->  $RE_MODEL_CHEMPROT   (CPR:3/4/5/6/9 + false)  @CHEMICAL$/@GENE$
    {GENETIC, DISEASE}   ->  $RE_MODEL_GAD        (associated + false)     @GENE$/@DISEASE$
    {GENETIC, GENETIC}   ->  $RE_MODEL_PPI        (interacts + false)      @GENE$/@GENE$
    {CHEMICAL, CHEMICAL} ->  $RE_MODEL_DDI        (mechanism/effect/...)   @DRUG$/@DRUG$
    any pair             ->  $RE_MODEL            (global fallback, if set)

Each route blinds with the markers its checkpoint was TRAINED on (ROUTE_MARKERS):
DDI uses @DRUG$ for both CHEMICAL entities, while the ChemProt route uses @CHEMICAL$
for the same NER type -- so blinding is deferred until after routing.

A pair with no applicable model is skipped and counted. Set at least one of the
env vars to a fine-tuned checkpoint -- a base LM like dmis-lab/biobert-base-cased-v1.1
has a RANDOM head and is useless here. Train them with train_re.py
(--task chemprot / --task gad / --task ppi).

DEPENDENCY-PATH OVERLAY (triples_strategy.html section 5.2 step 4)
-----------------------------------------------------------------
All-pairs generation explodes in high-entity sentences (up to 27 entities = 351
pairs), most of them spurious. Set RE_DEP_OVERLAY=1 to prune, in sentences with
>= RE_DEP_MIN_ENTITIES entities (default 5), any pair whose two entity head tokens
are more than RE_DEP_MAX_PATH (default 4) edges apart on the spaCy dependency tree.
Needs a parser model (`python -m spacy download en_core_web_sm`); fails OPEN (keeps
the pair) if the model is missing or a span can't be aligned. The run prints how
many pairs were pruned.

    $env:RE_MODEL_CHEMPROT = (Resolve-Path .\\chemprot-biobert-re)
    $env:RE_MODEL_GAD      = (Resolve-Path .\\gad-biobert-re)
    $env:RE_MODEL_PPI      = (Resolve-Path .\\ppi-biobert-re)

Input  : sentences/*.json                         (BioBERT NER spans)
Output : TRIPLES/triples_re.json                  scored triples; predicate carries
           {text:<label>, type:"relation", score:<p_rel>, model:<name>} and the triple
           carries the COMPOSITE score, score_2nd, score_components, self_pair, cues,
           result_margin (+ subject/object NER score)
         TRIPLES/relation_extraction.html          summary (label & score dists,
           per-model counts, samples)

Optionally also writes the GENETIC/DISEASE/CHEMICAL-normalized variant when
--normalize is passed, by reusing triples.py's normalization chain.

Run::  python relation_extraction.py            (set the RE_MODEL_* env vars first)
       python relation_extraction.py --normalize --limit 200
"""

import argparse
import glob
import html
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

# reuse triples.py's helpers so surfaces / normalization stay identical
import triples as T
import calibration as Cal

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent
SENT_DIR = ROOT / "sentences"
OUT_DIR = ROOT / "TRIPLES"
JSON_OUT = OUT_DIR / "triples_re.json"
NORM_OUT = OUT_DIR / "triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json"
HTML_OUT = OUT_DIR / "relation_extraction.html"

# --- config (override via env) ------------------------------------------------
BATCH_SIZE = int(os.environ.get("RE_BATCH", "32"))
MAX_LEN = int(os.environ.get("RE_MAX_LEN", "256"))
# keep a pair only if its endpoints are within this many characters of each other
# (0 = no limit). Caps the all-pairs combinatorics on long sentences.
MAX_PAIR_CHARS = int(os.environ.get("RE_MAX_PAIR_CHARS", "300"))
MIN_SCORE = float(os.environ.get("RE_MIN_SCORE", "0.5"))
# apply a per-checkpoint probability calibrator (calibration.json next to the
# model, fit by train_re.py) to p_rel before composite scoring. RE_CALIBRATE=0 disables.
CALIBRATE = os.environ.get("RE_CALIBRATE", "1").strip().lower() in ("1", "true", "yes", "on")
# on long full-corpus runs, flush partial results to disk every N batches so a crash
# or machine sleep doesn't lose hours of scoring (0 = write only at the end)
FLUSH_BATCHES = int(os.environ.get("RE_FLUSH_BATCHES", "200"))

# --- dependency-path overlay (triples_strategy.html section 5.2 step 4) --------
# In high-entity sentences (up to 27 entities -> 351 all-pairs), most pairs are
# spurious. When enabled, a sentence with >= RE_DEP_MIN_ENTITIES entities keeps a
# pair only if the two entity head tokens are within RE_DEP_MAX_PATH edges on the
# dependency tree. Opt-in (needs a spaCy parser model); fails OPEN (keeps the pair)
# if the model is missing or a span can't be aligned -- never silently loses recall.
DEP_OVERLAY = os.environ.get("RE_DEP_OVERLAY", "0").strip().lower() in ("1", "true", "yes", "on")
DEP_MIN_ENTITIES = int(os.environ.get("RE_DEP_MIN_ENTITIES", "5"))
DEP_MAX_PATH = int(os.environ.get("RE_DEP_MAX_PATH", "4"))
DEP_SPACY_MODEL = os.environ.get("RE_DEP_SPACY_MODEL", "en_core_web_sm")
# NER label -> entity-marker type in the blinded sentence. MUST match the markers
# the routed checkpoint was fine-tuned on (ChemProt: GENE/CHEMICAL; GAD: GENE/DISEASE).
TYPE_MARKER = {"GENETIC": "GENE", "DISEASE": "DISEASE", "CHEMICAL": "CHEMICAL"}
# explicit negative-class label names (comma-sep env), else auto-detected per model
NEG_LABELS = {s.strip().lower() for s in os.environ.get("RE_NEG_LABELS", "").split(",") if s.strip()}

# --- composite scoring (triples_strategy.html section 7) ----------------------
# The kept `score` fuses the relation-model softmax with independent evidence:
#   score = p_rel x min(ner) x cue_factor x margin_factor x section_factor x self_pair_factor
# NER is near-saturated in this corpus, so it acts as a gate; p_rel carries the
# discriminative weight. Each factor is recorded in score_components for ablation.
CUE_PRESENT, CUE_ABSENT = 1.0, 0.85          # relational cue in sentence / connecting text?
# self-pair (both endpoints same gene): penalize unless the relation is reflexive-plausible
SELF_PAIR_PENALTY = float(os.environ.get("RE_SELF_PAIR_PENALTY", "0.5"))
REFLEXIVE_OK = {s.strip().lower() for s in os.environ.get(
    "RE_REFLEXIVE_OK", "interacts,binds,phosphorylates,dimerizes,oligomerizes").split(",") if s.strip()}

# --- polarity / modality (triples_strategy.html section 6.1) ------------------
# A lightweight cue-list-in-a-window detector (no spaCy dependency): scan a window
# around the entity pair for grammatical-negation and speculation/hedge cues.
# A negated or speculated relation is still emitted (it IS a fact: "X does not bind
# Y") -- polarity/modality are annotations, not score factors or drop conditions.
# backend for negation: "cuelist" (default, zero-dependency window) or "negspacy"
# (NegEx/ConText scope resolution; needs `pip install spacy negspacy`). negspaCy
# covers NEGATION only -- speculation/modality always uses the cue-list window.
NEG_BACKEND = os.environ.get("RE_NEG_BACKEND", "cuelist").strip().lower()
NEG_WINDOW_L = int(os.environ.get("RE_NEG_WINDOW_LEFT", "40"))    # chars before subject
NEG_WINDOW_R = int(os.environ.get("RE_NEG_WINDOW_RIGHT", "30"))   # chars after object
# grammatical negation only (NOT biological "loss of"/"deficient", which don't
# negate the stated relation)
import re as _re
NEG_RE = _re.compile(
    r"\b(?:no|not|n't|cannot|can't|without|never|neither|nor|none|fails?|failed|"
    r"unable|absen(?:ce|t)|lacks?|lacking|unaffected|unchanged|unaltered|negative)\b", _re.I)
SPEC_RE = _re.compile(
    r"\b(?:may|might|could|would|possibl[ey]|potential(?:ly)?|putative|suggest\w*|"
    r"likely|presumabl\w*|appears?|appeared|seems?|seemed|hypothes\w*|speculat\w*|"
    r"propos\w*|candidate|predict\w*|probabl[ey]|presum\w*|conceivabl\w*|assum\w*)\b", _re.I)


_NEGEX_NLP = None
_NEGEX_WARNED = False


def _get_negex_nlp():
    """Lazy blank-English spaCy pipeline (tokenizer + sentencizer + negex). Blank,
    not a full model -- NegEx needs only tokens, sentence bounds, and the target
    entities (which we set by hand), not a parser/NER."""
    global _NEGEX_NLP
    if _NEGEX_NLP is None:
        import spacy
        import negspacy  # noqa: F401  (registers the "negex" factory)
        try:
            from negspacy.negation import Negex  # noqa: F401
        except Exception:
            pass
        nlp = spacy.blank("en")
        nlp.add_pipe("sentencizer")
        try:
            nlp.add_pipe("negex")
        except Exception:                       # older negspacy needs an explicit termset
            from negspacy.termsets import termset
            nlp.add_pipe("negex", config={"neg_termset": termset("en_clinical").get_patterns()})
        _NEGEX_NLP = nlp
    return _NEGEX_NLP


def _negspacy_negated(text, a, b):
    """True if NegEx marks either target entity as negated (sentence-scoped)."""
    nlp = _get_negex_nlp()
    doc = nlp.get_pipe("sentencizer")(nlp.make_doc(text))
    spans = []
    for e in (a, b):
        sp = doc.char_span(e["start"], e["end"], label=e["label"], alignment_mode="expand")
        if sp is not None:
            spans.append(sp)
    if not spans:
        return False
    doc.set_ents(spans)
    doc = nlp.get_pipe("negex")(doc)
    return any(ent._.negex for ent in doc.ents)


def negation_modality(text, a, b):
    """Polarity / modality of the relation. Returns (polarity, modality, sorted_cues).
    Modality (speculation) always uses the cue-list window; negation uses the
    selected backend (cue-list window or negspaCy NegEx scope)."""
    region = text[max(0, a["start"] - NEG_WINDOW_L):min(len(text), b["end"] + NEG_WINDOW_R)]
    spec = [m.lower() for m in SPEC_RE.findall(region)]
    if NEG_BACKEND == "negspacy":
        global _NEGEX_WARNED
        try:
            neg = ["negex"] if _negspacy_negated(text, a, b) else []
        except Exception as e:
            if not _NEGEX_WARNED:
                print(f"  [warn] RE_NEG_BACKEND=negspacy unavailable ({type(e).__name__}: {e}); "
                      f"falling back to cue-list negation. pip install spacy negspacy", flush=True)
                _NEGEX_WARNED = True
            neg = [m.lower() for m in NEG_RE.findall(region)]
    else:
        neg = [m.lower() for m in NEG_RE.findall(region)]
    return ("negated" if neg else "positive",
            "speculated" if spec else "asserted",
            sorted(set(neg + spec)))

# type-pair -> env var naming its fine-tuned checkpoint. RE_MODEL is a global
# fallback applied to any pair without a specific route.
ROUTE_ENV = {
    frozenset({"GENETIC", "CHEMICAL"}): "RE_MODEL_CHEMPROT",
    frozenset({"GENETIC", "DISEASE"}): "RE_MODEL_GAD",
    # same-type pairs collapse to a single-element frozenset
    frozenset({"GENETIC"}): "RE_MODEL_PPI",
    frozenset({"CHEMICAL"}): "RE_MODEL_DDI",       # drug-drug interaction
}
FALLBACK_ENV = "RE_MODEL"

# per-route marker overrides (merged onto TYPE_MARKER). A route's checkpoint must
# be blinded with the markers it was TRAINED on -- DDI uses @DRUG$ for both drugs,
# so on the CHEMICAL-CHEMICAL route a CHEMICAL entity is blinded as @DRUG$ (not the
# @CHEMICAL$ that the ChemProt route uses for the same NER type).
ROUTE_MARKERS = {
    frozenset({"CHEMICAL"}): {**TYPE_MARKER, "CHEMICAL": "DRUG"},
}


def route_markers(key):
    return ROUTE_MARKERS.get(key, TYPE_MARKER)

_MODEL_CACHE = {}   # model_name -> (tok, model, device, id2label, neg_ids)
_CALIBRATORS = {}   # model_name -> calibration spec (or None)


def blind(text, a, b, markers=None):
    """Sentence with entities a (earlier) and b (later) replaced by type markers.
    `markers` is a {NER label -> marker name} map (default the global TYPE_MARKER);
    a route may override it (e.g. DDI blinds CHEMICAL as @DRUG$)."""
    markers = markers or TYPE_MARKER

    def mk(label):
        return f"@{markers.get(label, label)}$"
    return (text[:a["start"]] + mk(a["label"]) + text[a["end"]:b["start"]]
            + mk(b["label"]) + text[b["end"]:])


_PARSER_NLP = None        # None=untried, False=tried&failed, else the spaCy nlp
_PARSER_WARNED = False
_DEP_STATS = {"sentences": 0, "examined": 0, "pruned": 0}


def _get_parser():
    """Lazy spaCy pipeline WITH a dependency parser (NER disabled). Returns the nlp,
    or None (warning once) if the model is unavailable -> overlay fails open."""
    global _PARSER_NLP, _PARSER_WARNED
    if _PARSER_NLP is None:
        try:
            import spacy
            nlp = spacy.load(DEP_SPACY_MODEL, disable=["ner", "lemmatizer", "textcat"])
            if not nlp.has_pipe("parser"):
                raise RuntimeError(f"{DEP_SPACY_MODEL} has no parser")
            _PARSER_NLP = nlp
        except Exception as e:
            _PARSER_NLP = False
            if not _PARSER_WARNED:
                print(f"  [warn] dependency overlay unavailable ({type(e).__name__}: {e}); "
                      f"high-entity pairs NOT pruned. python -m spacy download {DEP_SPACY_MODEL}",
                      flush=True)
                _PARSER_WARNED = True
    return _PARSER_NLP or None


def _dep_path_len(tok_a, tok_b):
    """Edge count of the shortest path between two tokens on the (undirected)
    dependency tree, via BFS over head+children links."""
    if tok_a.i == tok_b.i:
        return 0
    from collections import deque
    seen, q = {tok_a.i}, deque([(tok_a, 0)])
    while q:
        tok, d = q.popleft()
        neighbors = list(tok.children)
        if tok.head.i != tok.i:
            neighbors.append(tok.head)
        for nb in neighbors:
            if nb.i == tok_b.i:
                return d + 1
            if nb.i not in seen:
                seen.add(nb.i)
                q.append((nb, d + 1))
    return 10 ** 6            # disconnected (e.g. spans in different sentences)


def _dep_path_ok(doc, a, b):
    """True if the two entities' head tokens are within DEP_MAX_PATH dep edges.
    Fails OPEN (True) if either span can't be aligned to tokens."""
    sa = doc.char_span(a["start"], a["end"], alignment_mode="expand")
    sb = doc.char_span(b["start"], b["end"], alignment_mode="expand")
    if sa is None or sb is None:
        return True
    return _dep_path_len(sa.root, sb.root) <= DEP_MAX_PATH


def candidate_pairs(sentences, pmid):
    """Yield `meta` for every in-window entity pair in a doc. Blinding is deferred
    to run() (per-route markers), so meta carries the entity offsets (_a/_b) and the
    sentence text. Pairs are ordered by reading position (a before b). In high-entity
    sentences the dependency-path overlay (if enabled) prunes distant pairs."""
    for s in sentences:
        text = s.get("text", "")
        section = s.get("section_type", "")
        ents = [e for e in s.get("entities", []) if e.get("start") is not None]
        if len(ents) < 2:
            continue
        ents.sort(key=lambda e: e["start"])
        # parse once per sentence only when the overlay applies to it
        dep_doc = None
        if DEP_OVERLAY and len(ents) >= DEP_MIN_ENTITIES:
            nlp = _get_parser()
            if nlp is not None:
                dep_doc = nlp(text)
                _DEP_STATS["sentences"] += 1
        for i, a in enumerate(ents):
            for b in ents[i + 1:]:
                if b["start"] < a["end"]:        # overlapping spans -> skip
                    continue
                if MAX_PAIR_CHARS and (b["start"] - a["end"]) > MAX_PAIR_CHARS:
                    break                        # later b's are only farther
                if dep_doc is not None:          # dependency-path overlay
                    _DEP_STATS["examined"] += 1
                    if not _dep_path_ok(dep_doc, a, b):
                        _DEP_STATS["pruned"] += 1
                        continue
                connecting = " ".join(text[a["end"]:b["start"]].split())
                polarity, modality, asn_cues = negation_modality(text, a, b)
                yield {
                    "subject": {"text": T.entity_text(a), "type": a["label"], "score": a.get("score")},
                    "object": {"text": T.entity_text(b), "type": b["label"], "score": b.get("score")},
                    "connecting": connecting, "cues": s.get("cues") or [],
                    "result_margin": s.get("result_margin"),
                    "polarity": polarity, "modality": modality, "assertion_cues": asn_cues,
                    "pmid": pmid, "section": section, "sentence": text,
                    "_a": {"start": a["start"], "end": a["end"], "label": a["label"]},
                    "_b": {"start": b["start"], "end": b["end"], "label": b["label"]},
                }


def resolve_routes():
    """type-pair frozenset -> checkpoint name, plus a global fallback. Exits if
    nothing is configured."""
    routes = {}
    for key, var in ROUTE_ENV.items():
        v = os.environ.get(var, "").strip()
        if v:
            routes[key] = v
    fallback = os.environ.get(FALLBACK_ENV, "").strip() or None
    if not routes and not fallback:
        sys.exit("No RE checkpoint configured. Set RE_MODEL_CHEMPROT and/or RE_MODEL_GAD "
                 "(or RE_MODEL as a global fallback) to a FINE-TUNED checkpoint -- train "
                 "with train_re.py. See the module docstring.")
    return routes, fallback


def get_model(name):
    """Lazy-load + cache (tok, model, device, id2label, neg_ids) for a checkpoint."""
    if name in _MODEL_CACHE:
        return _MODEL_CACHE[name]
    try:
        import torch  # noqa: F401
        from transformers import AutoTokenizer, AutoModelForSequenceClassification
    except Exception as e:
        sys.exit(f"need torch + transformers ({type(e).__name__}: {e}).  pip install torch transformers")
    try:
        tok = AutoTokenizer.from_pretrained(name)
    except (ValueError, OSError):
        from transformers import BertTokenizer    # BioBERT ships only vocab.txt (slow WordPiece)
        tok = BertTokenizer.from_pretrained(name)
    model = AutoModelForSequenceClassification.from_pretrained(name)
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(device).eval()
    id2label = {int(i): str(l) for i, l in model.config.id2label.items()}
    if len(id2label) < 2:
        sys.exit(f"{name} has {len(id2label)} label(s); not an RE classifier.")
    neg_ids = set()
    for i, l in id2label.items():
        ll = l.lower()
        if (ll in NEG_LABELS if NEG_LABELS
                else ll in {"0", "false", "none", "no_relation", "negative", "not_related"}
                or ll.startswith("no_") or "no-relation" in ll):
            neg_ids.add(i)
    _CALIBRATORS[name] = Cal.load(Path(name) / "calibration.json") if CALIBRATE else None
    _MODEL_CACHE[name] = (tok, model, device, id2label, neg_ids)
    return _MODEL_CACHE[name]


def score_batch(tok, model, device, marked):
    """softmax probabilities for a batch of marked sentences. (B, num_labels)."""
    import torch
    enc = tok(marked, padding=True, truncation=True, max_length=MAX_LEN,
              return_tensors="pt").to(device)
    with torch.no_grad():
        return torch.softmax(model(**enc).logits, dim=-1).cpu().tolist()


def section_factor(section):
    s = (section or "").lower()
    if "result" in s or "abstract" in s or "figure" in s:
        return 1.0
    if "discussion" in s:
        return 0.97
    if "intro" in s:
        return 0.90
    return 0.95


def composite_score(p_rel, meta, label, self_pair):
    """Fuse the relation softmax with NER / cue / margin / section / self-pair evidence.
    Returns (score, components_dict). See triples_strategy.html section 7."""
    s_subj = meta["subject"].get("score") or 1.0
    s_obj = meta["object"].get("score") or 1.0
    ner = min(s_subj, s_obj)
    relational = ("relation" in (meta.get("cues") or [])) or bool(T.RELCUE.search(meta.get("connecting", "")))
    cue_f = CUE_PRESENT if relational else CUE_ABSENT
    m = meta.get("result_margin")
    m = m if isinstance(m, (int, float)) else 0.0
    margin_f = 1.0 + max(0.0, min(0.1, (m - 0.01) / 0.05))
    sec_f = section_factor(meta.get("section"))
    self_f = SELF_PAIR_PENALTY if (self_pair and label.lower() not in REFLEXIVE_OK) else 1.0
    comps = {"p_rel": round(p_rel, 4), "ner": round(ner, 4), "cue": cue_f,
             "margin": round(margin_f, 4), "section": sec_f, "self_pair": self_f}
    # margin_factor can boost slightly above 1.0; clamp the product to [0,1]
    return round(min(1.0, p_rel * ner * cue_f * margin_f * sec_f * self_f), 4), comps


def run(limit=None):
    routes, fallback = resolve_routes()
    # collect candidate pairs, bucketed by the checkpoint that should score them
    buckets = defaultdict(list)       # model_name -> [(marked, meta), ...]
    n_files = n_sent = n_cand = skipped = 0
    skip_pairs = Counter()
    _DEP_STATS.update(sentences=0, examined=0, pruned=0)
    for fp in sorted(glob.glob(str(SENT_DIR / "*.json"))):
        try:
            rec = json.load(open(fp, encoding="utf-8"))
        except Exception:
            continue
        n_files += 1
        pmid = Path(rec.get("source_file", Path(fp).name)).stem
        sents = rec.get("sentences", [])
        n_sent += len(sents)
        for meta in candidate_pairs(sents, pmid):
            n_cand += 1
            key = frozenset({meta["subject"]["type"], meta["object"]["type"]})
            name = routes.get(key, fallback)
            if not name:
                skipped += 1
                skip_pairs[tuple(sorted((meta["subject"]["type"], meta["object"]["type"])))] += 1
                continue
            # blind NOW, with the markers the routed checkpoint was trained on
            marked = blind(meta["sentence"], meta["_a"], meta["_b"], route_markers(key))
            buckets[name].append((marked, meta))
            if limit and (n_cand - skipped) >= limit:
                break
        if limit and (n_cand - skipped) >= limit:
            break

    print(f"routes: " + ", ".join(f"{'/'.join(sorted(k))}->{v}" for k, v in routes.items())
          + (f"  fallback={fallback}" if fallback else ""))
    print(f"candidate pairs: {n_cand:,} from {n_sent:,} sentences / {n_files:,} files; "
          f"routed to {len(buckets)} model(s), skipped {skipped:,} "
          f"({', '.join(f'{a}-{b}:{c}' for (a, b), c in skip_pairs.most_common(5)) or 'none'})")
    if DEP_OVERLAY and _DEP_STATS["examined"]:
        d = _DEP_STATS
        print(f"dependency overlay: pruned {d['pruned']:,}/{d['examined']:,} pairs "
              f"({100 * d['pruned'] / d['examined']:.1f}%) across {d['sentences']:,} "
              f"high-entity (>={DEP_MIN_ENTITIES}) sentences (max path {DEP_MAX_PATH})")

    def _flush(ts):
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        JSON_OUT.write_text(json.dumps(ts, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    triples = []
    total_to_score = sum(len(v) for v in buckets.values())
    scored_n = batch_n = 0
    for name, items in buckets.items():
        tok, model, device, id2label, neg_ids = get_model(name)
        cal = _CALIBRATORS.get(name)
        print(f"  [{name}] labels={list(id2label.values())}"
              f"{' calibrated(' + cal['method'] + ')' if cal else ''} -> scoring {len(items):,} pairs",
              flush=True)
        kept = 0
        for i in range(0, len(items), BATCH_SIZE):
            chunk = items[i:i + BATCH_SIZE]
            probs = score_batch(tok, model, device, [m for m, _ in chunk])
            for (m, meta), p in zip(chunk, probs):
                order = sorted(range(len(p)), key=lambda k: p[k], reverse=True)
                top, second = order[0], (order[1] if len(p) > 1 else order[0])
                if top in neg_ids:
                    continue
                label = id2label[top]
                # surface self-pair (same gene mention twice); refined post-normalization
                self_pair = (meta["subject"]["type"] == meta["object"]["type"]
                             and meta["subject"]["text"] == meta["object"]["text"])
                p_raw = p[top]
                p_rel = Cal.apply(cal, p_raw) if cal else p_raw   # calibrated relation prob
                score, comps = composite_score(p_rel, meta, label, self_pair)
                if score < MIN_SCORE:                 # threshold on the COMPOSITE score
                    continue
                kept += 1
                pred = {"text": label, "type": "relation", "score": round(p_raw, 4), "model": name}
                if cal:
                    pred["score_calibrated"] = round(p_rel, 4)
                triples.append({
                    "subject": meta["subject"],
                    "predicate": pred,
                    "object": meta["object"],
                    "score": score, "score_2nd": round(p[second], 4),
                    "score_components": comps, "self_pair": self_pair,
                    "polarity": meta["polarity"], "modality": meta["modality"],
                    "assertion_cues": meta["assertion_cues"],
                    "cues": meta.get("cues"), "result_margin": meta.get("result_margin"),
                    "pmid": meta["pmid"], "section": meta["section"], "sentence": meta["sentence"],
                })
            scored_n += len(chunk)
            batch_n += 1
            if FLUSH_BATCHES and batch_n % FLUSH_BATCHES == 0:
                _flush(triples)
                print(f"  ... scored {scored_n:,}/{total_to_score:,}  kept {len(triples):,}  (checkpointed)",
                      flush=True)
        print(f"  [{name}] kept {kept:,}", flush=True)

    _flush(triples)
    print(f"kept {len(triples):,} scored triples (>= {MIN_SCORE}) -> {JSON_OUT}")
    return triples


def normalize(triples):
    """Run triples.py's GENETIC->HGNC, DISEASE->MONDO, CHEMICAL->ChEBI chain so the
    scored triples carry the same id-normalization as triples_*_normalized.json."""
    rm, gk, ge = T.build_genetic_maps()
    dm, cm = T.build_disease_map(), T.build_chemical_map()

    def _sym(el):
        hs = el.get("hgnc_symbol")
        return hs if isinstance(hs, str) else ("|".join(hs) if isinstance(hs, list) else None)

    out, refined = [], 0
    for t in triples:
        t = T.normalize_genetic(t, rm, gk, ge)
        t = T.normalize_disease(t, dm)
        t = T.normalize_chemical(t, cm)
        # refine self-pair: same HGNC gene under different surfaces wasn't caught at
        # scoring time -- re-apply the penalty (unless reflexive-plausible) via components
        su, ob = _sym(t["subject"]), _sym(t["object"])
        if su is not None and su == ob and not t.get("self_pair"):
            t["self_pair"] = True
            if t["predicate"]["text"].lower() not in REFLEXIVE_OK:
                old = t["score_components"]["self_pair"]
                t["score"] = round(t["score"] / old * SELF_PAIR_PENALTY, 4)
                t["score_components"]["self_pair"] = SELF_PAIR_PENALTY
                refined += 1
        out.append(t)
    if refined:
        print(f"  self-pair refinement: penalized {refined:,} same-HGNC pairs found post-normalization")
    NORM_OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"normalized -> {NORM_OUT}")
    return out


def render_html(triples):
    esc = html.escape
    labels = Counter(t["predicate"]["text"] for t in triples)
    pairs = Counter((t["subject"]["type"], t["object"]["type"]) for t in triples)
    models = Counter(t["predicate"]["model"] for t in triples)
    buckets = Counter()
    for t in triples:
        lo = int(t["score"] * 10) * 10
        buckets[f"{lo}-{lo + 10}%"] += 1
    samp = sorted(triples, key=lambda t: t["score"], reverse=True)[:60]

    def tbl(rows):
        return "".join(f'<tr><td>{esc(str(a))}</td><td class="num">{b:,}</td></tr>' for a, b in rows)

    lab_rows = tbl(labels.most_common())
    pair_rows = tbl([(f"{a}->{b}", c) for (a, b), c in pairs.most_common()])
    mdl_rows = tbl(models.most_common())
    bkt_rows = tbl(sorted(buckets.items()))
    pm_rows = tbl([(f"{pol} / {mod}", c) for (pol, mod), c in
                   Counter((t["polarity"], t["modality"]) for t in triples).most_common()])
    def pm_tag(t):
        bits = []
        if t["polarity"] != "positive":
            bits.append(t["polarity"])
        if t["modality"] != "asserted":
            bits.append(t["modality"])
        return esc(", ".join(bits))
    samp_rows = "".join(
        f'<tr><td>{esc(t["subject"]["text"])}</td><td><code>{esc(t["predicate"]["text"])}</code></td>'
        f'<td>{esc(t["object"]["text"])}</td><td class="num">{t["score"]:.3f}</td>'
        f'<td>{pm_tag(t)}</td><td><code>{esc(t["pmid"])}</code></td></tr>' for t in samp)
    style = (" body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;"
             "margin:2rem auto;padding:0 1rem;color:#1a1a1a;} h1{font-size:1.45rem;}"
             " h2{font-size:1.15rem;margin-top:1.6rem;border-bottom:1px solid #ddd;padding-bottom:.3rem;}"
             " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;}"
             " th,td{border:1px solid #ccc;padding:.3rem .6rem;text-align:left;}"
             " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
             " .big{font-size:2rem;font-weight:700;} code{background:#f3f3f3;padding:1px 5px;border-radius:3px;}"
             " .headline{background:#eef4fb;border:1px solid #cdddf0;border-radius:8px;padding:.8rem 1rem;margin:1rem 0;}")
    HTML_OUT.write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>relation_extraction &mdash; scored RE triples</title><style>{style}</style></head><body>
<h1>Scored relation triples (BioBERT RE, type-pair routed)</h1>
<p>Entity pairs from <code>sentences/*.json</code> routed by entity-type pair to a fine-tuned
BioBERT checkpoint and classified under the entity-marker scheme; each kept triple carries a
softmax <code>score</code>. Data: <code>TRIPLES/triples_re.json</code>. Produced by
<code>relation_extraction.py</code>.</p>
<div class="headline"><span class="big">{len(triples):,}</span> scored triples
(score &ge; {MIN_SCORE}) across {len(models)} model(s)</div>
<h2>By model</h2><table><tr><th>checkpoint</th><th class="num">triples</th></tr>{mdl_rows}</table>
<h2>Predicted relation labels</h2><table><tr><th>label</th><th class="num">triples</th></tr>{lab_rows}</table>
<h2>By entity-type pair</h2><table><tr><th>subj-&gt;obj</th><th class="num">triples</th></tr>{pair_rows}</table>
<h2>Score distribution</h2><table><tr><th>bucket</th><th class="num">triples</th></tr>{bkt_rows}</table>
<h2>Polarity / modality</h2>
<p>Negation &amp; speculation detected by a cue-list window around each pair. Negated/speculated
relations are <em>kept</em> &mdash; "X does not bind Y" is a fact &mdash; and flagged for downstream use.</p>
<table><tr><th>polarity / modality</th><th class="num">triples</th></tr>{pm_rows}</table>
<h2>Top 60 by score</h2><table><tr><th>subject</th><th>relation</th><th>object</th>
<th class="num">score</th><th>polarity/modality</th><th>pmid</th></tr>{samp_rows}</table>
</body></html>""", encoding="utf-8")
    print(f"Wrote {HTML_OUT}")


def main():
    ap = argparse.ArgumentParser(description="learned, type-pair-routed RE step with per-triple scores")
    ap.add_argument("--normalize", action="store_true",
                    help="also write the HGNC/MONDO/ChEBI-normalized variant")
    ap.add_argument("--limit", type=int, default=None, help="cap candidate pairs (smoke test)")
    args = ap.parse_args()
    triples = run(limit=args.limit)
    if args.normalize:
        normalize(triples)
    render_html(triples)


if __name__ == "__main__":
    main()
