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

INPUT
    --dataset bigbio/bioinfer [--config bioinfer_bigbio_kb]   (needs `datasets`)
  or
    --input-json docs.json   (a list of BigBIO-KB documents, or {split: [docs]})
                             -- offline; also how this script is unit-tested.

Negative pairs explode in entity-dense sentences; --neg-ratio caps negatives to
N x positives (random, seeded). Datasets without a dev/test split can be carved
with --val-frac / --test-frac.

Sanity-check the entity/relation type mapping for a corpus BEFORE converting:
    python bigbio_to_re.py --task ppi --dataset bigbio/bioinfer --print-types

Run::  python bigbio_to_re.py --task ppi --dataset bigbio/bioinfer --out ppi_data
       python bigbio_to_re.py --task ppi --input-json docs.json --out ppi_data --neg-ratio 3
"""

import argparse
import json
import random
import re
import sys
from collections import Counter
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

CHEMPROT_EVAL = {"CPR:3", "CPR:4", "CPR:5", "CPR:6", "CPR:9"}   # the evaluated CPR groups
DDI_TYPES = {"mechanism", "effect", "advise", "int"}            # DDIExtraction-2013 positive types
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
    return True


def label_for(task, rel_type, chemprot_eval_only):
    if rel_type is None:                            # no annotated relation -> negative
        return "false" if task in ("chemprot", "ddi") else "0"
    if task == "chemprot":
        return rel_type if (not chemprot_eval_only or rel_type in CHEMPROT_EVAL) else "false"
    if task == "ddi":                               # normalize "DDI-effect"/"effect" -> the 4 types
        r = rel_type.lower().replace("ddi-", "").strip()
        return r if r in DDI_TYPES else "int"       # an annotated DDI with no/other subtype -> generic int
    return "1"


def _local_span(ent, ptext, poff):
    """(start,end) of an entity inside one passage's text, or None.
    Uses global offsets shifted by the passage offset; falls back to a text search."""
    offs = ent.get("offsets") or []
    txt = (ent.get("text") or [""])[0] if isinstance(ent.get("text"), list) else (ent.get("text") or "")
    if offs:
        ls = min(o[0] for o in offs) - poff
        le = max(o[1] for o in offs) - poff
        if 0 <= ls < le <= len(ptext) and (not txt or ptext[ls:le] == txt or txt in ptext[ls:le]):
            return ls, le
    if txt:
        i = ptext.find(txt)
        if i != -1:
            return i, i + len(txt)
    return None


def iter_instances(doc, task, chemprot_eval_only):
    """Yield (marked_sentence, label) for every candidate entity pair in a doc."""
    rels = {}
    for r in doc.get("relations", []):
        a1, a2 = r.get("arg1_id"), r.get("arg2_id")
        if a1 and a2:
            rels[frozenset((a1, a2))] = r.get("type", "") or ""
    for psg in doc.get("passages", []):
        ptext = (psg.get("text") or [""])[0] if isinstance(psg.get("text"), list) else (psg.get("text") or "")
        offs = psg.get("offsets") or [[0, len(ptext)]]
        poff = offs[0][0]
        located = []
        for e in doc.get("entities", []):
            sp = _local_span(e, ptext, poff)
            if sp:
                located.append((e, sp))
        for s_text, s0, s1 in sent_spans(ptext):
            here = [(e, (ls - s0, le - s0)) for e, (ls, le) in located if ls >= s0 and le <= s1]
            for i in range(len(here)):
                for j in range(i + 1, len(here)):
                    ea, (as0, ae0) = here[i]
                    eb, (bs0, be0) = here[j]
                    ma, mb = marker_for(task, ea.get("type")), marker_for(task, eb.get("type"))
                    if not valid_pair(task, ma, mb):
                        continue
                    if as0 > bs0:                   # order by reading position
                        ea, eb = eb, ea
                        as0, ae0, bs0, be0 = bs0, be0, as0, ae0
                        ma, mb = mb, ma
                    if bs0 < ae0:                   # overlapping spans -> skip
                        continue
                    marked = (s_text[:as0] + f"@{ma}$" + s_text[ae0:bs0]
                              + f"@{mb}$" + s_text[be0:])
                    marked = " ".join(marked.split())
                    rel = rels.get(frozenset((ea.get("id"), eb.get("id"))))
                    yield marked, label_for(task, rel, chemprot_eval_only)


def convert_split(docs, task, chemprot_eval_only, neg_ratio, rng):
    pos, neg = [], []
    neg_labels = {"0", "false"}
    for doc in docs:
        for sent, label in iter_instances(doc, task, chemprot_eval_only):
            (neg if label in neg_labels else pos).append((sent, label))
    if neg_ratio is not None and pos:
        cap = int(neg_ratio * len(pos))
        if len(neg) > cap:
            neg = rng.sample(neg, cap)
    rows = pos + neg
    rng.shuffle(rows)
    return rows


def print_types(splits, task, chemprot_eval_only):
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
        lab = label_for(task, t if t is not None else "", chemprot_eval_only)
        kind = "filtered->false" if lab in ("false", "0") else "POSITIVE"
        print(f"  {str(t):28.28s} {c:>9,}  -> {lab:8s} ({kind})")
    markers = sorted({marker_for(task, t) for t in ent_types})
    pairs = [f"{a}-{b}" for i, a in enumerate(markers) for b in markers[i:] if valid_pair(task, a, b)]
    print(f"\nvalid candidate type-pairs for task={task}: "
          f"{', '.join(pairs) if pairs else '(NONE -- entity types do not map to this task!)'}")
    print("\nIf a marker mapping looks wrong, the source uses different type names than "
          "marker_for() expects -- adjust marker_for or pick the right --task before converting.")


def write_tsv(path, rows):
    lines = ["index\tsentence\tlabel"]
    lines += [f"{i}\t{s}\t{l}" for i, (s, l) in enumerate(rows)]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


SPLIT_FILE = {"train": "train.tsv", "validation": "dev.tsv", "valid": "dev.tsv",
              "dev": "dev.tsv", "test": "test.tsv"}


def load_docs(args):
    """Return {split_name: [docs]}."""
    if args.input_json:
        data = json.loads(Path(args.input_json).read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {"train": data}
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
    cfg = args.config or f"{args.dataset.split('/')[-1]}_bigbio_kb"
    try:
        ds = load_dataset(args.dataset, name=cfg, **kw)
    except Exception as e:
        hint = ""
        if major >= 4:
            hint = (f"\n`datasets` {hfds.__version__} no longer runs dataset loading scripts, which "
                    f"BigBIO ({args.dataset}) relies on. Install a compatible version first:\n"
                    f"    pip install 'datasets<4'")
        sys.exit(f"could not load {args.dataset} (config {cfg}): {type(e).__name__}: {e}{hint}\n"
                 f"Check the config name (try --config) and that the dataset has a *_bigbio_kb schema.")
    return {split: list(ds[split]) for split in ds.keys()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", required=True, choices=["ppi", "chemprot", "gad", "ddi"])
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
    ap.add_argument("--print-types", action="store_true",
                    help="dump distinct entity/relation types + their mapping, then exit (no files written)")
    ap.add_argument("--max-docs", type=int, default=None, help="limit docs per split (quick test)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if not args.dataset and not args.input_json:
        ap.error("provide --dataset or --input-json")

    rng = random.Random(args.seed)
    splits = load_docs(args)
    if args.max_docs:
        splits = {k: v[:args.max_docs] for k, v in splits.items()}

    if args.print_types:
        print_types(splits, args.task, not args.chemprot_all_cpr)
        return

    out = Path(args.out or f"{args.task}_re_data")
    out.mkdir(parents=True, exist_ok=True)

    # build per-output-split rows
    converted = {}
    for split, docs in splits.items():
        rows = convert_split(docs, args.task, not args.chemprot_all_cpr, args.neg_ratio, rng)
        fname = SPLIT_FILE.get(split.lower(), "train.tsv")
        converted.setdefault(fname, []).extend(rows)

    # carve dev/test from train if requested and missing
    if (args.val_frac or args.test_frac) and "train.tsv" in converted:
        pool = converted["train.tsv"]
        rng.shuffle(pool)
        n = len(pool)
        n_test = int(args.test_frac * n)
        n_val = int(args.val_frac * n)
        converted["test.tsv"] = converted.get("test.tsv", []) + pool[:n_test]
        converted["dev.tsv"] = converted.get("dev.tsv", []) + pool[n_test:n_test + n_val]
        converted["train.tsv"] = pool[n_test + n_val:]

    total = 0
    for fname, rows in converted.items():
        write_tsv(out / fname, rows)
        dist = Counter(l for _, l in rows)
        total += len(rows)
        print(f"{fname:10s} {len(rows):>7,} rows  labels={dict(dist)}")
    print(f"-> {total:,} instances in {out}/  (task={args.task})")
    if "train.tsv" not in converted:
        print("  [warn] no train.tsv produced -- check entity types / relations in the source.")


if __name__ == "__main__":
    main()
