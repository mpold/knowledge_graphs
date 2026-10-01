#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bigbio_to_re.py -- convert a BigBIO KB dataset into train_re.py marked-TSV format.

BigBIO (https://huggingface.co/bigbio) distributes biomedical corpora under a
unified KB schema (passages / entities / relations). train_re.py instead wants
entity-blinded sentences with a label, e.g.

    @GENE$ directly binds @GENE$ .    1

This script bridges them: per sentence, it forms entity pairs of the task's type,
blinds the two targets with @GENE$/@CHEMICAL$/@DISEASE$ markers, and labels each
pair positive (a relation is annotated between them) or negative (co-occurring,
no relation). It writes train.tsv / dev.tsv / test.tsv ready for

    python train_re.py --task ppi --data <out>

TASKS (marker scheme + which pairs are candidates + label policy)
    ppi      : both endpoints -> @GENE$ ;   GENE-GENE pairs ; binary 1/0
    chemprot : @CHEMICAL$ / @GENE$       ;   CHEMICAL-GENE   ; CPR:3/4/5/6/9 + false
    gad      : @GENE$ / @DISEASE$         ;   GENE-DISEASE    ; binary 1/0
    ddi      : both endpoints -> @DRUG$  ;   DRUG-DRUG pairs ; mechanism/effect/advise/int + false
    biored   : @GENE$/@DISEASE$/@CHEMICAL$/@VARIANT$ ; the pair types BioRED
               annotates ; typed+signed labels (Positive_Correlation /
               Negative_Correlation / Bind / Association) + false

BioRED (Luo et al. 2022) annotates relations at the DOCUMENT level, so a related
pair yields a positive row in every sentence where both endpoints co-occur -- the
standard distant-supervision artefact. Two knobs address it:
    --require-cue          demote a document-level positive to negative unless the
                           connecting text carries a relational stem (RELCUE)
    --sample-positives N   print N random positive rows (blinded, with their
                           connecting text) for a hand-read; writes no files
Its eight relation types are also steeply skewed, so the four rare chemical-chemical
ones (Cotreatment / Comparison / Drug_Interaction / Conversion) are folded into
Association by default; --biored-all-types keeps all eight.

EXTRA BIND SUPERVISION (biored only)
BioRED's `Bind` class is ~1.4% of its rows, so the trained model almost never predicts
it. --extra-bind appends Bind rows from BioNLP-ST event corpora, which annotate binding
as a mention-level EVENT (trigger + Theme/Theme2/... arguments) rather than a relation:
    bigbio/bionlp_st_2013_ge   GENIA Event 2013 -- full-text papers, Results included
                               (~270 Bind rows, gene-gene only)
    bigbio/bionlp_st_2013_pc   Pathway Curation 2013 -- abstracts (~510 Bind rows,
                               ~45 of them chemical-gene)
(bigbio/bionlp_st_2011_ge does not load: its BigBIO script declares an empty
validation split, which `datasets` rejects.)
Every pair of themes of one Binding event is a `Bind` row. A co-occurring pair linked
by no event at all is `false`. A pair linked only by some OTHER event (regulation,
phosphorylation, ...) is skipped: BioRED would call it a correlation, so it must not
be taught as false. All annotated splits of these corpora go into train.tsv only --
BioRED's dev/test stay untouched, so dev F1 and calibration stay comparable across runs.

INPUT
    --dataset bigbio/bioinfer [--config bioinfer_bigbio_kb]   (needs `datasets`)
  or
    --input-json docs.json   (a list of BigBIO-KB documents, or {split: [docs]})
                             -- offline; also how this script is unit-tested.

Negative pairs explode in entity-dense sentences; --neg-ratio caps negatives to
N x positives (random, seeded). Datasets without a dev/test split can be carved
with --val-frac / --test-frac -- by DOCUMENT, so no sentence is in both train and dev.

Sanity-check the entity/relation type mapping for a corpus BEFORE converting:
    python bigbio_to_re.py --task ppi --dataset bigbio/bioinfer --print-types

Run::  python bigbio_to_re.py --task ppi --dataset bigbio/bioinfer --out ppi_data
       python bigbio_to_re.py --task ppi --input-json docs.json --out ppi_data --neg-ratio 3
       python bigbio_to_re.py --task biored --dataset bigbio/biored --out biored_data --val-frac 0
       python bigbio_to_re.py --task biored --dataset bigbio/biored --out biored_data --val-frac 0 \
           --extra-bind bigbio/bionlp_st_2013_ge,bigbio/bionlp_st_2013_pc
"""

import argparse
import json
import random
import re
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

CHEMPROT_EVAL = {"CPR:3", "CPR:4", "CPR:5", "CPR:6", "CPR:9"}   # the evaluated CPR groups
DDI_TYPES = {"mechanism", "effect", "advise", "int"}            # DDIExtraction-2013 positive types

# BioRED relation types (Luo et al. 2022). Unknown/other annotated relations fall
# back to the generic unsigned bucket rather than being dropped.
BIORED_TYPES = {"positive_correlation": "Positive_Correlation",
                "negative_correlation": "Negative_Correlation",
                "association": "Association", "bind": "Bind",
                "cotreatment": "Cotreatment", "comparison": "Comparison",
                "drug_interaction": "Drug_Interaction", "conversion": "Conversion"}

# The four chemical-chemical bookkeeping types are rare (plausibly tens of instances
# each across ~6,500 relations) and none of them adds sign or direction to the graph,
# so by default they collapse into Association: per-class support rises and the
# metric stops being dominated by classes we never write to the graph. Keep all
# eight with --biored-all-types.
BIORED_COLLAPSE = {"Cotreatment", "Comparison", "Drug_Interaction", "Conversion"}

# Entity-type pairs BioRED actually annotates. Same-type pairs collapse to a
# 1-element frozenset on both sides, so GENE-GENE is written as one member.
BIORED_PAIRS = {frozenset(("GENE",)), frozenset(("CHEMICAL",)), frozenset(("VARIANT",)),
                frozenset(("GENE", "DISEASE")), frozenset(("CHEMICAL", "GENE")),
                frozenset(("CHEMICAL", "DISEASE")), frozenset(("VARIANT", "DISEASE")),
                frozenset(("CHEMICAL", "VARIANT"))}

# Relational stems, mirroring triples.py's RELCUE. Used by --require-cue to demote a
# document-level BioRED positive whose connecting text asserts nothing.
RELCUE = re.compile(
    r"induc|inhibit|activat|regulat|express|\bbind|bound|target|associat|interact|"
    r"phosphorylat|suppress|promot|mediat|overexpress|knock\s?down|knock\s?out|"
    r"silenc|deplet|mutat|mutant|treat|encod|\bcaus|block|modulat|stimulat|increas|"
    r"decreas|reduc|elevat|attenuat|amelior|enhanc|repress|sensiti|resist|cleav|"
    r"degrad|secret|recruit|antagon|agoni|deficien|depend|correlat|involv|signal|"
    r"\brole\b|\beffect|abolish|abrogat|disrupt|restor|rescu|prevent|trigger|driv|"
    r"confer|\bloss\b|deletion|amplif|fusion|translocat|methylat|acetylat|ubiquitin|"
    r"glycosylat|inhibitor|activator|agonist|antagonist", re.I)

_SENT_NLP = None


def _sentencizer():
    """Cached blank-English spaCy sentencizer; None -> regex fallback."""
    global _SENT_NLP
    if _SENT_NLP is None:
        try:
            import spacy
            nlp = spacy.blank("en")
            nlp.add_pipe("sentencizer")
            _SENT_NLP = nlp
        except Exception:
            _SENT_NLP = False
    return _SENT_NLP or None


def sent_spans(text):
    """[(sentence_text, start_char, end_char), ...] over `text`."""
    nlp = _sentencizer()
    if nlp is not None:
        return [(s.text, s.start_char, s.end_char) for s in nlp(text).sents]
    spans, pos = [], 0
    for seg in re.findall(r"[^.!?]*[.!?]+(?:\s+|$)|[^.!?]+$", text):
        if seg.strip():
            i = text.find(seg, pos)
            spans.append((seg.rstrip(), i, i + len(seg.rstrip())))
            pos = i + len(seg)
    return spans or [(text, 0, len(text))]


def marker_for(task, etype):
    t = (etype or "").lower()
    if task == "ddi":
        return "DRUG"                               # ddi: every entity is a drug
    if task == "chemprot":
        return "CHEMICAL" if "chem" in t or "drug" in t else "GENE"
    if task == "gad":
        return "DISEASE" if "dis" in t else "GENE"
    if task == "biored":
        # substring matching (the existing idiom) so minor naming variants in the
        # loader do not break it. Order matters: `variant` is tested first.
        if "variant" in t:                       return "VARIANT"
        if "chem" in t or "drug" in t:           return "CHEMICAL"
        if "dis" in t or "phenotyp" in t:        return "DISEASE"
        if "cell" in t:                          return "CELLLINE"
        if "taxon" in t or "organism" in t or "species" in t: return "SPECIES"
        return "GENE"
    return "GENE"                                   # ppi: every entity is a protein


def valid_pair(task, ma, mb):
    s = {ma, mb}
    if task == "ppi":
        return ma == "GENE" and mb == "GENE"
    if task == "ddi":
        return ma == "DRUG" and mb == "DRUG"
    if task == "chemprot":
        return s == {"CHEMICAL", "GENE"}
    if task == "gad":
        return s == {"GENE", "DISEASE"}
    if task == "biored":
        # SPECIES and CELLLINE are deliberately absent from BIORED_PAIRS: BioRED
        # annotates them as entities but they are not relation endpoints we want as
        # graph edges. They still get blinded correctly if they appear as a target.
        return frozenset((ma, mb)) in BIORED_PAIRS
    return True


def label_for(task, rel_type, chemprot_eval_only, biored_collapse=True):
    if rel_type is None:                            # no annotated relation -> negative
        return "false" if task in ("chemprot", "ddi", "biored") else "0"
    if task == "chemprot":
        return rel_type if (not chemprot_eval_only or rel_type in CHEMPROT_EVAL) else "false"
    if task == "ddi":                               # normalize "DDI-effect"/"effect" -> the 4 types
        r = rel_type.lower().replace("ddi-", "").strip()
        return r if r in DDI_TYPES else "int"       # an annotated DDI with no/other subtype -> generic int
    if task == "biored":
        lab = BIORED_TYPES.get(rel_type.strip().lower(), "Association")
        return "Association" if (biored_collapse and lab in BIORED_COLLAPSE) else lab
    return "1"


_WORD = re.compile(r"\w")
# entity mentions that could not be placed cleanly in their passage text
_SPAN_STATS = {"realigned": 0, "dropped": 0}


def _boundary_ok(ptext, i, j, txt):
    """False if ptext[i:j] sits INSIDE a longer word, e.g. `dex` within `dexamethasone`.
    Only word characters on the entity's own edges are checked, so a mention that
    legitimately starts or ends on punctuation (`-catenin`, `(RELN)`) is unaffected."""
    if not txt:
        return False
    if _WORD.match(txt[0]) and i > 0 and _WORD.match(ptext[i - 1]):
        return False
    if _WORD.match(txt[-1]) and j < len(ptext) and _WORD.match(ptext[j]):
        return False
    return True


def _find_aligned(ptext, txt):
    """First occurrence of txt in ptext that is not embedded inside a longer word, else -1."""
    i = ptext.find(txt)
    while i != -1:
        if _boundary_ok(ptext, i, i + len(txt), txt):
            return i
        i = ptext.find(txt, i + 1)
    return -1


def _local_span(ent, ptext, poff):
    """(start,end) of an entity inside one passage's text, or None.

    Uses global offsets shifted by the passage offset; falls back to a text search.
    Both routes must land on WORD BOUNDARIES: a span that cuts into a longer word
    produces a corrupt blinded sentence (`@CHEMICAL$amethasone-induced` from an entity
    matched at `dex` inside `dexamethasone`), which is worse than no row at all. A
    misaligned annotated span retries as a text search; if nothing clean is found the
    mention is dropped and counted in _SPAN_STATS."""
    offs = ent.get("offsets") or []
    txt = (ent.get("text") or [""])[0] if isinstance(ent.get("text"), list) else (ent.get("text") or "")
    belongs_here = not offs                    # no offsets -> only the text search can place it
    if offs:
        ls = min(o[0] for o in offs) - poff
        le = max(o[1] for o in offs) - poff
        if 0 <= ls < le <= len(ptext) and (not txt or ptext[ls:le] == txt or txt in ptext[ls:le]):
            belongs_here = True                # the annotation puts this mention in THIS passage
            if _boundary_ok(ptext, ls, le, ptext[ls:le]):
                return ls, le
            _SPAN_STATS["realigned"] += 1      # stale/shifted offsets -> try the text search
    if txt:
        i = _find_aligned(ptext, txt)
        if i != -1:
            return i, i + len(txt)
    if belongs_here:                           # every entity is tried against every passage, so
        _SPAN_STATS["dropped"] += 1            # only count the ones that should have been placeable
    return None


def _sentence_pairs(doc, entities, mark, valid):
    """Yield (ea, eb, marked_sentence, connecting) for every pair of `entities` that
    co-occur in one sentence of `doc`: ordered by reading position, blinded with the
    markers mark(type) gives, and kept only if valid(marker_a, marker_b)."""
    for psg in doc.get("passages", []):
        ptext = (psg.get("text") or [""])[0] if isinstance(psg.get("text"), list) else (psg.get("text") or "")
        offs = psg.get("offsets") or [[0, len(ptext)]]
        poff = offs[0][0]
        located = []
        for e in entities:
            sp = _local_span(e, ptext, poff)
            if sp:
                located.append((e, sp))
        for s_text, s0, s1 in sent_spans(ptext):
            here = [(e, (ls - s0, le - s0)) for e, (ls, le) in located if ls >= s0 and le <= s1]
            for i in range(len(here)):
                for j in range(i + 1, len(here)):
                    ea, (as0, ae0) = here[i]
                    eb, (bs0, be0) = here[j]
                    ma, mb = mark(ea.get("type")), mark(eb.get("type"))
                    if not valid(ma, mb):
                        continue
                    if as0 > bs0:                   # order by reading position
                        ea, eb = eb, ea
                        as0, ae0, bs0, be0 = bs0, be0, as0, ae0
                        ma, mb = mb, ma
                    if bs0 < ae0:                   # overlapping spans -> skip
                        continue
                    connecting = s_text[ae0:bs0]
                    marked = (s_text[:as0] + f"@{ma}$" + connecting
                              + f"@{mb}$" + s_text[be0:])
                    yield ea, eb, " ".join(marked.split()), connecting


def iter_instances(doc, task, chemprot_eval_only, biored_collapse=True, require_cue=False,
                   with_connecting=False):
    """Yield (marked_sentence, label) for every candidate entity pair in a doc.

    `with_connecting` additionally yields the raw connecting text (for --sample-positives).
    `require_cue` demotes a positive whose connecting text carries no relational stem --
    the mitigation for document-level annotation (BioRED_task_mapping.md section 5)."""
    rels = {}
    for r in doc.get("relations", []):
        a1, a2 = r.get("arg1_id"), r.get("arg2_id")
        if a1 and a2:
            rels[frozenset((a1, a2))] = r.get("type", "") or ""
    for ea, eb, marked, connecting in _sentence_pairs(
            doc, doc.get("entities", []), lambda t: marker_for(task, t),
            lambda ma, mb: valid_pair(task, ma, mb)):
        rel = rels.get(frozenset((ea.get("id"), eb.get("id"))))
        label = label_for(task, rel, chemprot_eval_only, biored_collapse)
        # BioRED annotates at the DOCUMENT level, so a related pair emits a
        # positive row in EVERY sentence where both endpoints co-occur --
        # including sentences that merely mention them. --require-cue keeps
        # only those whose connecting text asserts something.
        if require_cue and rel is not None and not RELCUE.search(connecting):
            label = label_for(task, None, chemprot_eval_only, biored_collapse)
        yield (marked, label, " ".join(connecting.split())) if with_connecting \
            else (marked, label)


def event_marker(etype):
    """BioNLP-ST entity type -> BioRED marker, or None if it is not a relation endpoint.
    GE has Protein; PC has Gene_or_gene_product / Simple_chemical. Binding sites and
    promoters (GE `Entity`), Complex and Cellular_component are not endpoints BioRED
    would relate, so they stay as plain text in the sentence."""
    t = (etype or "").lower()
    if "protein" in t or "gene" in t:
        return "GENE"
    if "chem" in t:
        return "CHEMICAL"
    return None


def _event_entities(ev, events, seen=None):
    """Entity ids an event reaches through its arguments, following nested events."""
    seen = set() if seen is None else seen
    if ev.get("id") in seen:
        return set()
    seen.add(ev.get("id"))
    out = set()
    for a in ev.get("arguments", []):
        ref = a.get("ref_id")
        if ref in events:
            out |= _event_entities(events[ref], events, seen)
        elif ref:
            out.add(ref)
    return out


def iter_bind_instances(doc):
    """Yield (marked_sentence, label) Bind/false rows from a BigBIO-KB EVENT document.

    Bind  : the two mentions are both Themes of one Binding event.
    false : the pair is linked by no event at all (directly or through nesting).
    A pair linked only by a non-Binding event yields nothing (see EXTRA BIND above).
    Annotation is per mention, so no document-level artefact and no --require-cue."""
    ents = [e for e in doc.get("entities", []) if event_marker(e.get("type"))]
    ids = {e.get("id") for e in ents}
    events = {ev.get("id"): ev for ev in doc.get("events", [])}
    bind, linked = set(), set()
    for ev in events.values():
        reach = _event_entities(ev, events) & ids
        linked |= {frozenset(p) for p in combinations(reach, 2)}
        if (ev.get("type") or "").lower() == "binding":
            themes = {a.get("ref_id") for a in ev.get("arguments", [])
                      if (a.get("role") or "").lower().startswith("theme") and a.get("ref_id") in ids}
            bind |= {frozenset(p) for p in combinations(themes, 2)}
    for ea, eb, marked, _ in _sentence_pairs(
            doc, ents, event_marker, lambda ma, mb: frozenset((ma, mb)) in BIORED_PAIRS):
        key = frozenset((ea.get("id"), eb.get("id")))
        if key in bind:
            yield marked, "Bind"
        elif key not in linked:
            yield marked, "false"


NEG_LABELS = {"0", "false"}


def convert_split(docs, task, chemprot_eval_only, neg_ratio, rng,
                  biored_collapse=True, require_cue=False):
    pos, neg = [], []
    neg_labels = NEG_LABELS
    for doc in docs:
        for sent, label in iter_instances(doc, task, chemprot_eval_only,
                                          biored_collapse, require_cue):
            (neg if label in neg_labels else pos).append((sent, label))
    if neg_ratio is not None and pos:
        cap = int(neg_ratio * len(pos))
        if len(neg) > cap:
            neg = rng.sample(neg, cap)
    rows = pos + neg
    rng.shuffle(rows)
    return rows


def print_types(splits, task, chemprot_eval_only, biored_collapse=True):
    """Dump distinct entity/relation types + how this task maps them. No files written."""
    ent_types, rel_types, n_docs = Counter(), Counter(), 0
    for docs in splits.values():
        for doc in docs:
            n_docs += 1
            for e in doc.get("entities", []):
                ent_types[e.get("type")] += 1
            for r in doc.get("relations", []):
                rel_types[r.get("type")] += 1
    print(f"docs: {n_docs:,} across splits {list(splits)}  (task={task})")
    print(f"\nentity types ({len(ent_types)}) -> marker:")
    for t, c in ent_types.most_common():
        print(f"  {str(t):28.28s} {c:>9,}  -> @{marker_for(task, t)}$")
    print(f"\nrelation types ({len(rel_types)}) -> label:")
    for t, c in rel_types.most_common():
        lab = label_for(task, t if t is not None else "", chemprot_eval_only, biored_collapse)
        kind = "filtered->false" if lab in ("false", "0") else "POSITIVE"
        print(f"  {str(t):28.28s} {c:>9,}  -> {lab:8s} ({kind})")
    ev_types = Counter(ev.get("type") for docs in splits.values() for doc in docs
                       for ev in doc.get("events", []))
    if ev_types:
        print(f"\nevent types ({len(ev_types)}) -- only Binding is used, via --extra-bind:")
        for t, c in ev_types.most_common():
            print(f"  {str(t):28.28s} {c:>9,}")
        print("event-corpus endpoints (--extra-bind): "
              + ", ".join(f"{t} -> @{event_marker(t)}$" for t in ent_types if event_marker(t)))
    markers = sorted({marker_for(task, t) for t in ent_types})
    pairs = [f"{a}-{b}" for i, a in enumerate(markers) for b in markers[i:] if valid_pair(task, a, b)]
    print(f"\nvalid candidate type-pairs for task={task}: "
          f"{', '.join(pairs) if pairs else '(NONE -- entity types do not map to this task!)'}")
    print("\nIf a marker mapping looks wrong, the source uses different type names than "
          "marker_for() expects -- adjust marker_for or pick the right --task before converting.")


def sample_positives(splits, task, chemprot_eval_only, n, rng,
                     biored_collapse=True, require_cue=False):
    """Print n random POSITIVE rows (blinded sentence + connecting text) for a hand-read.

    BioRED_task_mapping.md section 5.1: with document-level annotation the only way to
    know the label quality is to read positives and count how many actually assert the
    relation. That number is the ceiling on what the model can learn. No files written."""
    rows = []
    for split, docs in splits.items():
        for doc in docs:
            for sent, label, conn in iter_instances(doc, task, chemprot_eval_only,
                                                    biored_collapse, require_cue,
                                                    with_connecting=True):
                if label not in NEG_LABELS:
                    rows.append((split, label, sent, conn))
    if not rows:
        print("no positive rows produced -- check the entity/relation types (--print-types).")
        return
    picked = rng.sample(rows, min(n, len(rows)))
    print(f"{len(rows):,} positive rows total; showing {len(picked)} at random "
          f"(task={task}, require_cue={require_cue})\n")
    for i, (split, label, sent, conn) in enumerate(picked, 1):
        print(f"[{i:>3}] {label}   ({split})")
        print(f"      {sent}")
        print(f"      connecting: {conn!r}\n")
    print("Count how many of these sentences actually ASSERT the labelled relation.\n"
          "A low fraction means the document-level artefact is biting: re-run with\n"
          "--require-cue (and re-sample) before training.")


def write_tsv(path, rows):
    lines = ["index\tsentence\tlabel"]
    lines += [f"{i}\t{s}\t{l}" for i, (s, l) in enumerate(rows)]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


SPLIT_FILE = {"train": "train.tsv", "validation": "dev.tsv", "valid": "dev.tsv",
              "dev": "dev.tsv", "test": "test.tsv"}


def load_json_docs(path):
    """{split_name: [docs]} from a local BigBIO-KB JSON (a list, or {split: [docs]})."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else {"train": data}


def load_docs(args):
    """Return {split_name: [docs]}."""
    if args.input_json:
        return load_json_docs(args.input_json)
    return load_bigbio(args.dataset, args.config)


def load_bigbio(dataset, config=None):
    """{split_name: [docs]} from a BigBIO HF dataset in its *_bigbio_kb schema."""
    try:
        import datasets as hfds
        from datasets import load_dataset
    except Exception as e:
        sys.exit(f"need the `datasets` library for --dataset ({type(e).__name__}). pip install 'datasets<4'")

    # BigBIO corpora are distributed as dataset *loading scripts*. `datasets` >= 4.0
    # removed loading-script support and no longer honors `trust_remote_code` (it warns
    # "not supported anymore" and ignores it), so the *_bigbio_kb configs cannot be
    # fetched there. Only pass the flag on versions that still run scripts; otherwise
    # the load below fails and we point at the version pin.
    major = int(hfds.__version__.split(".")[0])
    kw = {"trust_remote_code": True} if major < 4 else {}
    cfg = config or f"{dataset.split('/')[-1]}_bigbio_kb"
    try:
        ds = load_dataset(dataset, name=cfg, **kw)
    except Exception as e:
        hint = ""
        if major >= 4:
            hint = (f"\n`datasets` {hfds.__version__} no longer runs dataset loading scripts, which "
                    f"BigBIO ({dataset}) relies on. Install a compatible version first:\n"
                    f"    pip install 'datasets<4'")
        # A BigBIO loading script pulls in the parser for the corpus's native format --
        # e.g. biored is BioC XML and imports `bioc`, which is not on the Kaggle image.
        # That failure surfaces here as a plain ModuleNotFoundError about a package the
        # user never named, so say what to install.
        if isinstance(e, ModuleNotFoundError) and e.name:
            hint += (f"\nThe loading script for {dataset} needs the `{e.name}` package "
                     f"(it parses the corpus's native format). Install it and re-run:\n"
                     f"    pip install {e.name}")
        sys.exit(f"could not load {dataset} (config {cfg}): {type(e).__name__}: {e}{hint}\n"
                 f"Check the config name (try --config) and that the dataset has a *_bigbio_kb schema.")
    return {split: list(ds[split]) for split in ds.keys()}


def extra_bind_rows(sources, neg_ratio, rng, max_docs=None):
    """Bind/false rows from every annotated split of each event corpus in `sources`.

    A split with no events at all is skipped: BioNLP-ST withholds test annotations,
    and every pair there would otherwise be written as a false negative."""
    rows = []
    for src in sources:
        splits = load_json_docs(src) if src.lower().endswith(".json") else load_bigbio(src)
        pos, neg = [], []
        for split, docs in splits.items():
            docs = docs[:max_docs] if max_docs else docs
            if not any(doc.get("events") for doc in docs):
                print(f"  [extra-bind] {src} {split}: no events annotated -- skipped")
                continue
            for doc in docs:
                for sent, label in iter_bind_instances(doc):
                    (neg if label in NEG_LABELS else pos).append((sent, label))
        if neg_ratio is not None and len(neg) > int(neg_ratio * len(pos)):
            neg = rng.sample(neg, int(neg_ratio * len(pos)))
        print(f"  [extra-bind] {src}: {len(pos):,} Bind + {len(neg):,} false -> train.tsv")
        rows += pos + neg
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", required=True, choices=["ppi", "chemprot", "gad", "ddi", "biored"])
    ap.add_argument("--dataset", help="BigBIO HF dataset, e.g. bigbio/bioinfer")
    ap.add_argument("--config", help="BigBIO config (default <name>_bigbio_kb)")
    ap.add_argument("--input-json", help="local BigBIO-KB docs (list, or {split: [docs]})")
    ap.add_argument("--out", help="output dir (default <task>_re_data)")
    ap.add_argument("--neg-ratio", type=float, default=None,
                    help="cap negatives to N x positives per split (default: keep all)")
    ap.add_argument("--val-frac", type=float, default=0.0, help="carve a dev split from train if none exists")
    ap.add_argument("--test-frac", type=float, default=0.0, help="carve a test split from train if none exists")
    ap.add_argument("--chemprot-all-cpr", action="store_true",
                    help="keep all CPR relation types as positives (default: only CPR:3/4/5/6/9)")
    ap.add_argument("--biored-all-types", action="store_true",
                    help="biored: keep all 8 relation types (default folds the 4 rare "
                         "chemical-chemical types into Association)")
    ap.add_argument("--require-cue", action="store_true",
                    help="biored: demote an annotated positive to negative unless its connecting "
                         "text carries a relational stem (document-level annotation mitigation)")
    ap.add_argument("--sample-positives", type=int, default=None, metavar="N",
                    help="print N random positive rows for a hand-read, then exit (no files written)")
    ap.add_argument("--print-types", action="store_true",
                    help="dump distinct entity/relation types + their mapping, then exit (no files written)")
    ap.add_argument("--extra-bind", type=lambda s: [x.strip() for x in s.split(",") if x.strip()],
                    default=None, metavar="SRC[,SRC...]",
                    help="biored: append Bind/false rows from BioNLP-ST event corpora to train.tsv "
                         "(BigBIO names, e.g. bigbio/bionlp_st_2013_ge,bigbio/bionlp_st_2013_pc, "
                         "or local BigBIO-KB .json files)")
    ap.add_argument("--extra-bind-neg-ratio", type=float, default=1.0,
                    help="cap the extra corpora's false rows to N x their Bind rows (default 1.0; "
                         "BioRED already supplies plenty of negatives)")
    ap.add_argument("--max-docs", type=int, default=None, help="limit docs per split (quick test)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if not args.dataset and not args.input_json:
        ap.error("provide --dataset or --input-json")
    if args.extra_bind and args.task != "biored":
        ap.error("--extra-bind adds BioRED `Bind` rows; it needs --task biored")

    rng = random.Random(args.seed)
    splits = load_docs(args)
    if args.max_docs:
        splits = {k: v[:args.max_docs] for k, v in splits.items()}

    collapse = not args.biored_all_types
    if args.print_types:
        print_types(splits, args.task, not args.chemprot_all_cpr, collapse)
        return
    if args.sample_positives:
        sample_positives(splits, args.task, not args.chemprot_all_cpr, args.sample_positives,
                         rng, collapse, args.require_cue)
        return

    out = Path(args.out or f"{args.task}_re_data")
    out.mkdir(parents=True, exist_ok=True)

    # group docs by output file
    by_file = {}
    for split, docs in splits.items():
        by_file.setdefault(SPLIT_FILE.get(split.lower(), "train.tsv"), []).extend(docs)

    # carve dev/test from train BY DOCUMENT, before any rows exist. One sentence yields a row
    # per entity pair, so carving rows put near-copies of training sentences in dev (95% of
    # BioInfer's dev rows) -- inflating dev F1 and the threshold/calibrator fit on it.
    if (args.val_frac or args.test_frac) and by_file.get("train.tsv"):
        pool = list(by_file["train.tsv"])
        rng.shuffle(pool)
        n = len(pool)
        n_test = int(args.test_frac * n)
        n_val = int(args.val_frac * n)
        by_file.setdefault("test.tsv", []).extend(pool[:n_test])
        by_file.setdefault("dev.tsv", []).extend(pool[n_test:n_test + n_val])
        by_file["train.tsv"] = pool[n_test + n_val:]
        print(f"carved by document: {n_val:,} dev + {n_test:,} test docs out of {n:,} train docs")

    converted = {fname: convert_split(docs, args.task, not args.chemprot_all_cpr, args.neg_ratio,
                                      rng, collapse, args.require_cue)
                 for fname, docs in by_file.items()}

    # after the carve, so event-corpus rows never leak into a carved dev/test split
    if args.extra_bind:
        extra = extra_bind_rows(args.extra_bind, args.extra_bind_neg_ratio, rng, args.max_docs)
        pool = converted.setdefault("train.tsv", [])
        pool.extend(extra)
        rng.shuffle(pool)

    total = 0
    for fname, rows in converted.items():
        write_tsv(out / fname, rows)
        dist = Counter(l for _, l in rows)
        total += len(rows)
        print(f"{fname:10s} {len(rows):>7,} rows  labels={dict(dist)}")
    print(f"-> {total:,} instances in {out}/  (task={args.task})")
    if _SPAN_STATS["realigned"] or _SPAN_STATS["dropped"]:
        print(f"  [spans] {_SPAN_STATS['realigned']:,} mention(s) had an annotated span cutting into a "
              f"word and were re-located by text search; {_SPAN_STATS['dropped']:,} could not be placed "
              f"cleanly and were skipped (they would have blinded a partial word).")
    if "train.tsv" not in converted:
        print("  [warn] no train.tsv produced -- check entity types / relations in the source.")


if __name__ == "__main__":
    main()
