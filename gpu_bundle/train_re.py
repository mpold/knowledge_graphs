#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""train_re.py -- fine-tune BioBERT for biomedical relation extraction.

One trainer, several tasks (a registry -- add more by extending TASKS). Each
produces a sequence-classification checkpoint that relation_extraction.py loads
as-is: trained on the same entity-marker scheme and keeping a "false" negative
class, both of which that script already understands (and routes by entity-type
pair -- see RE_MODEL_CHEMPROT / RE_MODEL_GAD there).

TASKS
-----
  chemprot : BioCreative VI CHEMPROT. CHEMICAL<->GENE relations, markers
             @CHEMICAL$ / @GENE$. 6 classes: CPR:3/4/5/6/9 + false.
  gad      : Genetic Association Database. GENE<->DISEASE association, markers
             @GENE$ / @DISEASE$. Binary: associated (1) + false (0).
  ppi      : gene-gene / protein-protein interaction (AIMed / BioInfer / ...).
             Both endpoints blinded to @GENE$. Binary: interacts (1) + false (0).
  ddi      : drug-drug interaction (DDIExtraction 2013). Both endpoints -> @DRUG$.
             4 classes: mechanism / effect / advise / int (+ false).
The official metric for all is micro-F1 over the POSITIVE classes (false
excluded); that is what --data dev reports and what selects the best checkpoint.

DATA
----
--data is a directory of BioBERT-preprocessed TSVs (train.tsv / dev.tsv /
test.tsv). Each line is the entity-blinded sentence and its label, tab-separated;
column order and an optional index/header are auto-detected (the field carrying
an @MARKER$ is the sentence; the label is the last known-label field). These TSVs
ship with the BioBERT RE release (dmis-lab/biobert, datasets/RE/{ChemProt,GAD}).

No data handy?  `python train_re.py --task gad --smoke` writes a tiny synthetic
set to ./gad_smoke and trains 1 epoch on CPU just to prove the pipeline (NOT a
real model).

USAGE
-----
    pip install torch transformers accelerate
    python train_re.py --task chemprot --data path/to/ChemProt --out chemprot-biobert-re
    python train_re.py --task gad      --data path/to/GAD      --out gad-biobert-re
    # then wire both into relation_extraction.py:
    $env:RE_MODEL_CHEMPROT = (Resolve-Path .\\chemprot-biobert-re)
    $env:RE_MODEL_GAD      = (Resolve-Path .\\gad-biobert-re)
    python relation_extraction.py --normalize

Run `python train_re.py -h` for all flags.
"""

import argparse
import re
import sys
from pathlib import Path

NEG = "false"                                  # canonical negative class (after mapping)
MARKER_RE = re.compile(r"@[A-Z]+\$")           # entity marker in a blinded sentence

# task registry: label_names maps RAW corpus labels -> readable predicate (the
# negative raw label MUST map to NEG); known_raw is the set of raw label strings
# used to locate the label column; markers are added as special tokens with
# --add-marker-tokens.
TASKS = {
    "chemprot": {
        "label_names": {"CPR:3": "upregulator/activator", "CPR:4": "downregulator/inhibitor",
                        "CPR:5": "agonist", "CPR:6": "antagonist",
                        "CPR:9": "substrate/product-of", "false": NEG},
        "known_raw": {"CPR:3", "CPR:4", "CPR:5", "CPR:6", "CPR:9", "false", "true"},
        "markers": ["@CHEMICAL$", "@GENE$"],
    },
    "gad": {
        "label_names": {"1": "associated", "0": NEG, "true": "associated", "false": NEG},
        "known_raw": {"0", "1", "true", "false"},
        "markers": ["@GENE$", "@DISEASE$"],
    },
    # gene-gene / protein-protein interaction. BOTH endpoints are blinded to the
    # same @GENE$ marker (the two targets; other mentions stay as text), matching
    # how PPI corpora (AIMed / BioInfer / HPRD50 / IEPA / LLL) are reformatted for
    # BioBERT. Binary: interacts (1) + false (0).
    "ppi": {
        "label_names": {"1": "interacts", "0": NEG, "true": "interacts", "false": NEG},
        "known_raw": {"0", "1", "true", "false"},
        "markers": ["@GENE$"],
    },
    # drug-drug interaction (DDIExtraction 2013). BOTH endpoints -> @DRUG$. 4 positive
    # types + false. Handles both the BioBERT-release labels (DDI-mechanism / DDI-false)
    # and the bare forms (mechanism / false) emitted by bigbio_to_re.py.
    "ddi": {
        "label_names": {"DDI-mechanism": "mechanism", "DDI-effect": "effect",
                        "DDI-advise": "advise", "DDI-int": "int", "DDI-false": NEG,
                        "mechanism": "mechanism", "effect": "effect", "advise": "advise",
                        "int": "int", "false": NEG, "true": "int"},
        "known_raw": {"DDI-mechanism", "DDI-effect", "DDI-advise", "DDI-int", "DDI-false",
                      "mechanism", "effect", "advise", "int", "false", "true"},
        "markers": ["@DRUG$"],
    },
}


def read_examples(path, cfg):
    """Parse one BioBERT-style RE TSV -> [(sentence, mapped_label), ...].

    Robust to column order, a leading index column, and a header row: the field
    containing an @MARKER$ is the sentence; the label is the last field that is a
    known raw label (last-column-first avoids matching an index column like "0"/"1"),
    else any other known-label field.
    """
    known, names = cfg["known_raw"], cfg["label_names"]
    rows = []
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        fields = ln.split("\t")
        sent = next((f for f in fields if MARKER_RE.search(f)), None)
        if sent is None:                       # header / malformed -> skip
            continue
        last = fields[-1].strip()
        if last in known and fields[-1] is not sent:
            raw = last
        else:
            raw = next((f.strip() for f in fields if f is not sent and f.strip() in known), last)
        rows.append((sent, names.get(raw, raw)))
    if not rows:
        sys.exit(f"no usable rows in {path} (need an @MARKER$ sentence + a label column)")
    return rows


def write_smoke(task, d):
    """Tiny synthetic set in the task's format so the pipeline can be tested
    without downloading the corpus. Not a meaningful model."""
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    if task == "chemprot":
        pos = [("@CHEMICAL$ is a potent inhibitor of @GENE$ .", "CPR:4"),
               ("@CHEMICAL$ strongly activates @GENE$ in cells .", "CPR:3"),
               ("@CHEMICAL$ acts as an agonist of @GENE$ .", "CPR:5"),
               ("@CHEMICAL$ is an antagonist at the @GENE$ receptor .", "CPR:6"),
               ("@GENE$ metabolizes @CHEMICAL$ to its active form .", "CPR:9")]
        neg = [("@CHEMICAL$ was measured alongside @GENE$ in plasma .", "false"),
               ("Levels of @CHEMICAL$ and @GENE$ were reported .", "false")]
    elif task == "gad":
        pos = [("@GENE$ mutations are strongly associated with @DISEASE$ .", "1"),
               ("Variants in @GENE$ confer increased risk of @DISEASE$ .", "1"),
               ("@GENE$ polymorphism is linked to @DISEASE$ susceptibility .", "1")]
        neg = [("@GENE$ and @DISEASE$ were both described in the cohort .", "0"),
               ("No relationship between @GENE$ and @DISEASE$ was found .", "0")]
    elif task == "ppi":
        pos = [("@GENE$ directly interacts with @GENE$ in the complex .", "1"),
               ("@GENE$ binds @GENE$ to form a heterodimer .", "1"),
               ("@GENE$ phosphorylates @GENE$ during signaling .", "1")]
        neg = [("@GENE$ and @GENE$ were both quantified by qPCR .", "0"),
               ("Expression of @GENE$ was normalized to @GENE$ .", "0")]
    elif task == "ddi":
        pos = [("@DRUG$ increases the plasma concentration of @DRUG$ .", "mechanism"),
               ("@DRUG$ enhances the anticoagulant effect of @DRUG$ .", "effect"),
               ("Patients taking @DRUG$ should avoid @DRUG$ .", "advise"),
               ("A pharmacokinetic interaction between @DRUG$ and @DRUG$ occurs .", "int")]
        neg = [("@DRUG$ and @DRUG$ were administered in the trial .", "false"),
               ("No interaction between @DRUG$ and @DRUG$ was found .", "false")]
    else:
        sys.exit(f"no smoke template for task {task}")
    body = (pos * 8) + (neg * 8)
    header = "index\tsentence\tlabel\n"
    for name in ("train.tsv", "dev.tsv", "test.tsv"):
        lines = [f"{i}\t{s}\t{l}" for i, (s, l) in enumerate(body)]
        (d / name).write_text(header + "\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote synthetic {task} set -> {d}")
    return str(d)


def build_training_args(out, epochs, bs, lr, seed=42):
    """TrainingArguments across transformers versions (eval_strategy vs
    evaluation_strategy)."""
    from transformers import TrainingArguments
    common = dict(output_dir=out, num_train_epochs=epochs,
                  per_device_train_batch_size=bs, per_device_eval_batch_size=bs * 2,
                  learning_rate=lr, weight_decay=0.01, warmup_ratio=0.1,
                  logging_steps=50, save_total_limit=1, seed=seed, data_seed=seed,
                  load_best_model_at_end=True, metric_for_best_model="f1",
                  greater_is_better=True, report_to=[])
    for key in ("eval_strategy", "evaluation_strategy"):
        try:
            return TrainingArguments(**{key: "epoch"}, save_strategy="epoch", **common)
        except TypeError:
            continue
    return TrainingArguments(**common)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", required=True, choices=sorted(TASKS), help="RE task / dataset")
    ap.add_argument("--data", help="dir with train.tsv / dev.tsv / test.tsv (BioBERT format)")
    ap.add_argument("--smoke", action="store_true",
                    help="generate a tiny synthetic set and train 1 epoch (pipeline test only)")
    ap.add_argument("--model", default="dmis-lab/biobert-base-cased-v1.1",
                    help="base LM to fine-tune (default BioBERT)")
    ap.add_argument("--out", help="output checkpoint dir (default <task>-biobert-re)")
    ap.add_argument("--epochs", type=float, default=3.0)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--lr", type=float, default=2e-5)
    ap.add_argument("--max-len", type=int, default=128)
    ap.add_argument("--add-marker-tokens", action="store_true",
                    help="add the task's @MARKER$ tokens as special tokens (default off, matches BioBERT)")
    ap.add_argument("--calibration", choices=["isotonic", "platt", "none"], default="isotonic",
                    help="fit a probability calibrator on the dev split (default isotonic)")
    ap.add_argument("--seed", type=int, default=42,
                    help="random seed for reproducibility (data shuffling + classifier-head init)")
    args = ap.parse_args()
    cfg = TASKS[args.task]
    out = args.out or f"{args.task}-biobert-re"

    if args.smoke:
        args.data = write_smoke(args.task, f"{args.task}_smoke")
        args.epochs = min(args.epochs, 1.0)
    if not args.data:
        ap.error("--data is required (or use --smoke)")

    data_dir = Path(args.data)
    train = read_examples(data_dir / "train.tsv", cfg)
    dev_p, test_p = data_dir / "dev.tsv", data_dir / "test.tsv"
    dev = read_examples(dev_p, cfg) if dev_p.exists() else None
    test = read_examples(test_p, cfg) if test_p.exists() else None

    # label space: positives sorted, negative last; consistent id2label/label2id
    labels = sorted({l for _, l in train if l != NEG})
    if any(l == NEG for _, l in train):
        labels.append(NEG)
    label2id = {l: i for i, l in enumerate(labels)}
    id2label = {i: l for l, i in label2id.items()}
    neg_id = label2id.get(NEG)
    print(f"task={args.task}  labels ({len(labels)}): {labels}")
    print(f"train={len(train):,}  dev={len(dev) if dev else 0:,}  test={len(test) if test else 0:,}")

    # heavy imports here so --help / --smoke parsing stay fast and dependency-light
    try:
        import numpy as np
        import torch
        from torch.utils.data import Dataset
        from transformers import (AutoTokenizer, AutoModelForSequenceClassification,
                                  DataCollatorWithPadding, Trainer, set_seed)
    except Exception as e:
        sys.exit(f"need torch + transformers + accelerate ({type(e).__name__}: {e}).  "
                 f"pip install torch transformers accelerate")

    # Seed BEFORE building the model: the classifier head is reinitialized at
    # from_pretrained() time, so seeding here (not just inside Trainer) is what
    # makes the fresh head -- the main source of run-to-run variation -- reproducible.
    set_seed(args.seed)
    print(f"seed={args.seed}")

    # Some BioBERT checkpoints (e.g. dmis-lab/biobert-base-cased-v1.1) ship a
    # minimal config.json with no "model_type", so the Auto* loaders can't infer
    # the architecture. Fall back to the explicit BERT classes in that case.
    try:
        tok = AutoTokenizer.from_pretrained(args.model)
    except (ValueError, KeyError):
        from transformers import BertTokenizerFast
        tok = BertTokenizerFast.from_pretrained(args.model)
    if args.add_marker_tokens:
        tok.add_special_tokens({"additional_special_tokens": cfg["markers"]})
    try:
        model = AutoModelForSequenceClassification.from_pretrained(
            args.model, num_labels=len(labels), id2label=id2label, label2id=label2id,
            ignore_mismatched_sizes=True)  # reinit a fresh head when fine-tuning from one
    except (ValueError, KeyError):
        from transformers import BertForSequenceClassification
        model = BertForSequenceClassification.from_pretrained(
            args.model, num_labels=len(labels), id2label=id2label, label2id=label2id,
            ignore_mismatched_sizes=True)
    if args.add_marker_tokens:
        model.resize_token_embeddings(len(tok))

    class DS(Dataset):
        def __init__(self, rows):
            self.enc = tok([s for s, _ in rows], truncation=True, max_length=args.max_len)
            self.labels = [label2id[l] for _, l in rows]

        def __len__(self):
            return len(self.labels)

        def __getitem__(self, i):
            item = {k: torch.tensor(v[i]) for k, v in self.enc.items()}
            item["labels"] = torch.tensor(self.labels[i])
            return item

    def compute_metrics(p):
        preds = np.argmax(p.predictions, axis=-1)
        gold = p.label_ids
        # micro P/R/F1 over POSITIVE classes only (official: exclude the negative class)
        tp = fp = fn = 0
        for pr, gd in zip(preds, gold):
            if pr == neg_id and gd == neg_id:
                continue
            if pr == gd and pr != neg_id:
                tp += 1
            else:
                if pr != neg_id:
                    fp += 1
                if gd != neg_id:
                    fn += 1
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        return {"precision": prec, "recall": rec, "f1": f1, "accuracy": float((preds == gold).mean())}

    targs = build_training_args(out, args.epochs, args.batch_size, args.lr, args.seed)
    tkw = dict(model=model, args=targs, train_dataset=DS(train),
               eval_dataset=DS(dev) if dev else None, compute_metrics=compute_metrics,
               data_collator=DataCollatorWithPadding(tok))
    try:                                  # transformers>=4.46 renamed tokenizer-> processing_class
        trainer = Trainer(processing_class=tok, **tkw)
    except TypeError:
        trainer = Trainer(tokenizer=tok, **tkw)

    trainer.train()
    trainer.save_model(out)
    tok.save_pretrained(out)
    print(f"saved checkpoint -> {out}")

    if test:
        m = trainer.evaluate(DS(test))
        print("TEST  " + "  ".join(f"{k.replace('eval_', '')}={v:.4f}"
                                   for k, v in m.items() if k.startswith("eval_")))

    # ---- probability calibration (triples_strategy.html section 7) ----
    # Fit on dev: for each POSITIVE-predicted example, (softmax of predicted label,
    # 1 iff prediction correct). Save calibration.json next to the checkpoint;
    # relation_extraction.py applies it to p_rel before composite scoring.
    def _pos_pairs(rows):
        pr = trainer.predict(DS(rows))
        z = pr.predictions - pr.predictions.max(axis=1, keepdims=True)
        probs = np.exp(z) / np.exp(z).sum(axis=1, keepdims=True)
        xs, ys = [], []
        for p, g in zip(probs, pr.label_ids):
            k = int(p.argmax())
            if k == neg_id:                  # only positive predictions become triples
                continue
            xs.append(float(p[k])); ys.append(1.0 if k == g else 0.0)
        return xs, ys

    if args.calibration != "none" and dev:
        import calibration as Cal
        xs, ys = _pos_pairs(dev)
        spec = Cal.fit(xs, ys, method=args.calibration)
        if spec is None:
            print(f"calibration: skipped ({len(xs)} positive-predicted dev examples / single class)")
        else:
            Cal.save(f"{out}/calibration.json", spec)
            exs, eys = _pos_pairs(test) if test else (xs, ys)
            where = "test" if test else "dev(in-sample)"
            e_raw, e_cal = Cal.ece(exs, eys), Cal.ece([Cal.apply(spec, x) for x in exs], eys)
            print(f"calibration: {args.calibration} on {len(xs)} dev positives -> {out}/calibration.json")
            print(f"  ECE on {where} ({len(exs)} positives): raw {e_raw:.4f} -> calibrated {e_cal:.4f}")

    var = {"chemprot": "RE_MODEL_CHEMPROT", "gad": "RE_MODEL_GAD",
           "ppi": "RE_MODEL_PPI", "ddi": "RE_MODEL_DDI"}.get(args.task, "RE_MODEL")
    print(f"\nUse it:  $env:{var} = (Resolve-Path .\\{out})")
    print("then:    python relation_extraction.py --normalize")


if __name__ == "__main__":
    main()
