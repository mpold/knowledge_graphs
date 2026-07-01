#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_re_pipeline.py -- one command for the whole biomedical RE pipeline.

Unifies the three steps we otherwise run by hand:

    1. CONVERT    bigbio_to_re.py   -- BigBIO corpus  -> marked TSVs (train/dev/test)
    2. TRAIN      train_re.py       -- fine-tune BioBERT into a checkpoint
    3. CALIBRATE  calibration.py     -- evaluate on test (metrics + confusion matrix),
                                       fit a probability calibrator on dev, write
                                       <model>/calibration.json AND summaries/calibration.html

Steps 1 and 2 are run as subprocesses (so they stay byte-for-byte the same scripts,
with all their existing flags and fixes). Training is run with --calibration none so
calibration happens exactly once, here, as an explicit and re-runnable step.

DEFAULTS reproduce the run we did:
    python run_re_pipeline.py
is equivalent to
    python bigbio_to_re.py --task ppi --dataset bigbio/bioinfer --out ppi_data \
        --val-frac 0.1 --seed 0
    python train_re.py --task ppi --data ppi_data --out ppi-biobert-re \
        --seed 42 --calibration none
    # + calibration on the saved checkpoint -> ppi-biobert-re/calibration.json

Note: --val-frac defaults to 0.1 because train_re.py requires a dev split (the bare
conversion command produces only train/test, which makes training stop with an error).

USAGE
    python run_re_pipeline.py                         # ppi / bioinfer, all defaults
    python run_re_pipeline.py --task gad --dataset bigbio/gad
    python run_re_pipeline.py --skip-convert          # reuse existing TSVs
    python run_re_pipeline.py --skip-convert --skip-train   # (re)calibrate only

GPU / KAGGLE
    The training and calibration steps use HuggingFace Trainer, which moves the
    model to CUDA automatically when a GPU is visible -- no code change needed.
    On Kaggle: enable the accelerator (Settings -> Accelerator -> GPU), then just
    run the script. By default it pins to ONE GPU (--gpus 0), which is best for
    Kaggle's 2x T4 (Trainer would otherwise DataParallel across both and add
    overhead for this small model). The script prints a [device] line at startup
    so you can confirm the GPU is active.
      python run_re_pipeline.py                 # single GPU (default)
      python run_re_pipeline.py --gpus 0,1      # use both Kaggle T4s
      python run_re_pipeline.py --gpus cpu      # force CPU
    Kaggle's torch already ships with CUDA, so no torch reinstall is required.

    NOTE: step 1 (bigbio_to_re.py) fetches the BigBIO corpus, which is distributed
    as a dataset loading script, so install `datasets` once per session before running:
        !pip install -q datasets
    (Skip this when --skip-convert reuses TSVs that are already on disk.)
"""

import argparse
import datetime
import html
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SUMMARY_PATH = ROOT / "summaries" / "calibration.html"


def _fmt_bytes(n):
    step = 1024.0
    for unit in ("B", "KB", "MB", "GB"):
        if n < step or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= step
    return f"{n:.1f} GB"


def write_summary_html(ctx):
    """Render the run summary (metrics + confusion matrix + outputs) to summaries/calibration.html."""
    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)   # mkdir('summaries') if needed
    e = html.escape
    labels = ctx["labels"]
    cm = ctx["confusion"]
    totals = [sum(row) for row in cm]
    col_tot = [sum(cm[r][c] for r in range(len(labels))) for c in range(len(labels))]

    # confusion-matrix table (rows = gold/true, cols = predicted)
    head_cells = "".join(f"<th>pred: {e(l)}</th>" for l in labels)
    cm_rows = []
    for r, lbl in enumerate(labels):
        cells = []
        for c in range(len(labels)):
            diag = ' class="diag"' if r == c else ""
            cells.append(f"<td{diag}>{cm[r][c]:,}</td>")
        cm_rows.append(f"<tr><th>true: {e(lbl)}</th>{''.join(cells)}<td class='tot'>{totals[r]:,}</td></tr>")
    cm_foot = "".join(f"<td class='tot'>{c:,}</td>" for c in col_tot)
    cm_table = (f"<table class='cm'><thead><tr><th></th>{head_cells}<th class='tot'>total</th></tr></thead>"
                f"<tbody>{''.join(cm_rows)}<tr><th class='tot'>total</th>{cm_foot}"
                f"<td class='tot'>{sum(totals):,}</td></tr></tbody></table>")

    m = ctx["metrics"]
    thr = m.get("threshold")
    thr_rows = ""
    if thr is not None:
        thr_rows = (f"<tr><td>Positive-decision threshold (tuned on dev)</td><td>{thr:.4f}</td></tr>"
                    f"<tr><td>F1 at plain argmax (baseline)</td><td>{m.get('f1_argmax', m['f1']):.4f}</td></tr>")
    metrics_table = (
        "<table><thead><tr><th>Metric</th><th>Value</th></tr></thead><tbody>"
        f"<tr><td>Precision (micro, positives)</td><td>{m['precision']:.4f}</td></tr>"
        f"<tr><td>Recall (micro, positives)</td><td>{m['recall']:.4f}</td></tr>"
        f"<tr><td><strong>F1 (micro, positives)</strong></td><td><strong>{m['f1']:.4f}</strong></td></tr>"
        f"<tr><td>Accuracy (all classes)</td><td>{m['accuracy']:.4f}</td></tr>"
        f"{thr_rows}"
        "</tbody></table>")

    cal = ctx["calibration"]
    if cal["done"]:
        ecal = f"{cal['ece_cal']:.4f}" if cal["ece_cal"] is not None else "&mdash;"
        cal_html = (
            f"<p>Method: <code>{e(cal['method'])}</code>, fit on {cal['n_dev_pos']:,} dev positives "
            f"&rarr; <code>{e(cal['json'])}</code></p>"
            "<table><thead><tr><th></th><th>Raw</th><th>Calibrated</th></tr></thead><tbody>"
            f"<tr><td>ECE on {e(cal['eval_name'])} ({cal['n_eval_pos']:,} positives)</td>"
            f"<td>{cal['ece_raw']:.4f}</td><td><strong>{ecal}</strong></td></tr></tbody></table>")
    else:
        cal_html = f"<p class='note'>Calibration not produced: {e(cal['reason'])}</p>"

    out_rows = "".join(
        f"<tr><td><code>{e(name)}</code></td><td>{e(size)}</td><td>{e(desc)}</td></tr>"
        for name, size, desc in ctx["outputs"])
    outputs_table = ("<table><thead><tr><th>File</th><th>Size</th><th>What</th></tr></thead>"
                     f"<tbody>{out_rows}</tbody></table>")

    splits = ctx["splits"]
    splits_html = ", ".join(f"{k} = {v:,}" for k, v in splits.items() if v is not None)

    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RE pipeline summary &mdash; {e(ctx['task'])}</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
         line-height: 1.55; max-width: 900px; margin: 2rem auto; padding: 0 1.2rem; color: #1b1b1b; }}
  h1 {{ font-size: 1.6rem; border-bottom: 2px solid #ddd; padding-bottom: .4rem; }}
  h2 {{ font-size: 1.2rem; margin-top: 1.8rem; border-bottom: 1px solid #eee; padding-bottom: .3rem; }}
  code {{ background: #f3f3f3; padding: .1rem .35rem; border-radius: 4px; font-size: .9em; }}
  table {{ border-collapse: collapse; margin: .8rem 0; }}
  th, td {{ border: 1px solid #ddd; padding: .4rem .7rem; text-align: right; }}
  th {{ background: #f7f7f7; text-align: center; }}
  table.cm td.diag {{ background: #e6f4ea; font-weight: 600; }}
  table.cm td.tot, table.cm th.tot {{ background: #fafafa; color: #555; }}
  .note {{ background: #fff8e1; border-left: 4px solid #f0c040; padding: .6rem 1rem; border-radius: 4px; }}
  .meta {{ color: #666; font-size: .9em; }}
</style>
</head>
<body>
<h1>Relation-extraction pipeline summary</h1>
<p class="meta">task = <code>{e(ctx['task'])}</code> &middot; dataset = <code>{e(ctx['dataset'])}</code>
 &middot; model = <code>{e(ctx['model'])}</code> &middot; seed = {ctx['seed']}
 &middot; device = {e(ctx['device'])} &middot; generated {e(ctx['timestamp'])}</p>
<p class="meta">data splits: {e(splits_html)} &middot; labels: {", ".join(f"<code>{e(l)}</code>" for l in labels)}
 &middot; evaluated on: <strong>{e(ctx['metrics']['eval_name'])}</strong> ({ctx['metrics']['n']:,} examples)</p>

<h2>Test metrics</h2>
{metrics_table}
<p class="meta">F1 is micro-averaged over the positive classes only (the negative class
<code>{e(ctx['neg'])}</code> is excluded) &mdash; the official RE metric.{
" A positive is predicted when 1&minus;P(" + e(ctx['neg']) + ") reaches the threshold tuned on dev to maximize F1; the metrics and confusion matrix below reflect that operating point. (ECE is reported over the argmax positive predictions, independent of this threshold.)" if thr is not None else ""}</p>

<h2>Confusion matrix</h2>
<p class="meta">Rows = true label, columns = predicted label. Diagonal (correct) highlighted.</p>
{cm_table}

<h2>Calibration</h2>
{cal_html}

<h2>Outputs</h2>
{outputs_table}
</body>
</html>
"""
    SUMMARY_PATH.write_text(doc, encoding="utf-8")
    print(f"summary written -> {SUMMARY_PATH}")


def setup_device(gpus):
    """Select/report the compute device and export it to child processes.

    `gpus` becomes CUDA_VISIBLE_DEVICES so both this process and the
    train/convert subprocesses see the same devices:
      "0"      -> use the first GPU only (default; best for Kaggle's 2x T4, since
                  HF Trainer otherwise auto-enables DataParallel across both, which
                  adds overhead and can break for a small BERT model)
      "0,1"    -> expose both GPUs (Trainer will DataParallel across them)
      "" / cpu -> force CPU
      "all"    -> leave CUDA_VISIBLE_DEVICES untouched (use whatever is visible)
    """
    sel = (gpus or "").strip().lower()
    if sel in ("", "cpu", "none"):
        os.environ["CUDA_VISIBLE_DEVICES"] = ""        # hide all GPUs -> CPU
    elif sel != "all":
        os.environ["CUDA_VISIBLE_DEVICES"] = gpus      # pin to the requested device(s)
    # else "all": leave inherited CUDA_VISIBLE_DEVICES as-is

    try:
        import torch
    except Exception:
        print("[device] torch not importable yet; GPU check skipped.")
        return
    if torch.cuda.is_available():
        names = ", ".join(torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count()))
        print(f"[device] GPU ENABLED -- {torch.cuda.device_count()} visible: {names}  "
              f"(CUDA {torch.version.cuda}, torch {torch.__version__})")
    else:
        msg = "[device] running on CPU"
        if sel in ("", "cpu", "none"):
            msg += " (forced via --gpus cpu)"
        elif "+cpu" in torch.__version__:
            msg += f" -- torch {torch.__version__} is a CPU-only build (no CUDA support compiled in)"
        else:
            msg += " -- no CUDA device visible. On Kaggle, enable a GPU accelerator in "
            msg += "Notebook settings (Settings -> Accelerator -> GPU)."
        print(msg)


def run_step(name, cmd):
    """Run a subprocess step, streaming its output; abort the pipeline on failure.
    Child inherits os.environ (incl. CUDA_VISIBLE_DEVICES set by setup_device)."""
    print(f"\n{'=' * 70}\n[{name}]  {' '.join(cmd)}\n{'=' * 70}", flush=True)
    r = subprocess.run(cmd, cwd=str(ROOT), env=os.environ.copy())
    if r.returncode != 0:
        sys.exit(f"[{name}] failed with exit code {r.returncode} -- pipeline aborted.")


def convert(args):
    cmd = [sys.executable, "bigbio_to_re.py", "--task", args.task,
           "--dataset", args.dataset, "--out", args.data, "--seed", str(args.convert_seed)]
    if args.config:
        cmd += ["--config", args.config]
    if args.val_frac:
        cmd += ["--val-frac", str(args.val_frac)]
    if args.test_frac:
        cmd += ["--test-frac", str(args.test_frac)]
    if args.neg_ratio is not None:
        cmd += ["--neg-ratio", str(args.neg_ratio)]
    run_step("1/3 CONVERT", cmd)


def train(args):
    # --calibration none: calibration is done explicitly in step 3 below.
    # --max-len matches the calibration step so train/eval tokenize identically.
    cmd = [sys.executable, "train_re.py", "--task", args.task, "--data", args.data,
           "--out", args.model, "--seed", str(args.seed), "--epochs", str(args.epochs),
           "--max-len", str(args.max_len), "--calibration", "none"]
    if args.base_model:
        cmd += ["--model", args.base_model]
    if args.add_marker_tokens:
        cmd += ["--add-marker-tokens"]
    if args.lr is not None:
        cmd += ["--lr", str(args.lr)]
    if args.batch_size is not None:
        cmd += ["--batch-size", str(args.batch_size)]
    run_step("2/3 TRAIN", cmd)


def calibrate(args):
    """Evaluate the checkpoint, fit calibration.py on dev, and write the HTML summary.

    Runs predictions once per split, then: (a) optionally tunes the positive-decision
    threshold on dev to maximize micro-F1 and computes test metrics + a confusion
    matrix at that operating point (plain argmax with --no-tune-threshold), (b) fits
    a calibrator on dev positives (unless --calibration none) and reports ECE on test,
    and (c) renders everything to summaries/calibration.html."""
    print(f"\n{'=' * 70}\n[3/3 CALIBRATE + SUMMARY]  {args.calibration} on {args.model}\n{'=' * 70}", flush=True)

    sys.path.insert(0, str(ROOT))                       # train_re + calibration importable
    import numpy as np
    import torch
    from torch.utils.data import Dataset
    from transformers import (AutoTokenizer, AutoModelForSequenceClassification,
                              Trainer, DataCollatorWithPadding)
    import calibration as Cal
    from train_re import TASKS, read_examples, NEG

    cfg = TASKS[args.task]
    model_dir = Path(args.model)
    tok = AutoTokenizer.from_pretrained(str(model_dir))
    model = AutoModelForSequenceClassification.from_pretrained(str(model_dir))
    label2id = model.config.label2id
    n_lbl = len(label2id)
    id2label = {i: l for l, i in label2id.items()}
    labels = [id2label[i] for i in range(n_lbl)]
    neg_id = label2id.get(NEG)

    data = Path(args.data)
    dev = read_examples(data / "dev.tsv", cfg)
    test = read_examples(data / "test.tsv", cfg) if (data / "test.tsv").exists() else None
    train_n = len(read_examples(data / "train.tsv", cfg)) if (data / "train.tsv").exists() else None

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

    trainer = Trainer(model=model, data_collator=DataCollatorWithPadding(tok))

    def predict(rows):
        pr = trainer.predict(DS(rows))
        z = pr.predictions - pr.predictions.max(axis=1, keepdims=True)
        probs = np.exp(z) / np.exp(z).sum(axis=1, keepdims=True)
        return probs.argmax(axis=1), np.asarray(pr.label_ids), probs

    def pos_pairs(preds, gold, probs):
        xs, ys = [], []
        for k, g, p in zip(preds, gold, probs):
            if int(k) == neg_id:                         # only positive predictions become triples
                continue
            xs.append(float(p[int(k)]))
            ys.append(1.0 if int(k) == int(g) else 0.0)
        return xs, ys

    # predictions: dev (for calibration fit) + an evaluation split (test if present)
    dev_preds, dev_gold, dev_probs = predict(dev)
    if test:
        ev_preds, ev_gold, ev_probs, ev_name = (*predict(test), "test")
    else:
        ev_preds, ev_gold, ev_probs, ev_name = dev_preds, dev_gold, dev_probs, "dev(in-sample)"

    # ---- micro P/R/F1 over positive classes (official RE metric) ----
    def micro_prf(preds, gold):
        tp = fp = fn = 0
        for p, g in zip(preds, gold):
            if p == neg_id and g == neg_id:
                continue
            if p == g and p != neg_id:
                tp += 1
            else:
                if p != neg_id:
                    fp += 1
                if g != neg_id:
                    fn += 1
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        return prec, rec, f1

    # ---- optional decision threshold tuned on dev, applied to the eval split ----
    # Predict a positive iff its score 1-P(neg) clears the threshold; the positive
    # class is the argmax over non-negative classes. Sweeping the dev pos-scores
    # finds the F1-maximizing operating point, then we report eval at that point.
    def best_pos(p):
        bk, bv = None, -1.0
        for k in range(n_lbl):
            if k != neg_id and p[k] > bv:
                bk, bv = k, p[k]
        return bk if bk is not None else int(p.argmax())

    def predict_at(probs, t):
        return np.array([best_pos(p) if (1.0 - p[neg_id]) >= t else neg_id for p in probs])

    def tune_threshold(probs, gold):
        best_t, best_f1 = 0.0, -1.0
        for t in [0.0] + sorted({float(1.0 - p[neg_id]) for p in probs}):
            f1 = micro_prf(predict_at(probs, t), gold)[2]
            if f1 > best_f1:
                best_f1, best_t = f1, t
        return best_t

    f1_argmax = micro_prf(ev_preds, ev_gold)[2]
    threshold, ev_used = None, ev_preds
    if args.tune_threshold and neg_id is not None:
        threshold = tune_threshold(dev_probs, dev_gold)
        ev_used = predict_at(ev_probs, threshold)
        print(f"threshold tuned on dev: {threshold:.4f}  "
              f"(argmax F1={f1_argmax:.4f} -> tuned F1={micro_prf(ev_used, ev_gold)[2]:.4f} on {ev_name})")

    # ---- confusion matrix + metrics at the chosen operating point ----
    cm = [[0] * n_lbl for _ in range(n_lbl)]
    for g, p in zip(ev_gold, ev_used):
        cm[int(g)][int(p)] += 1
    prec, rec, f1 = micro_prf(ev_used, ev_gold)
    acc = float((ev_used == ev_gold).mean())
    metrics = {"precision": prec, "recall": rec, "f1": f1, "accuracy": acc,
               "eval_name": ev_name, "n": len(ev_gold),
               "threshold": threshold, "f1_argmax": f1_argmax}
    print(f"metrics on {ev_name}: P={prec:.4f} R={rec:.4f} F1={f1:.4f} acc={acc:.4f} (n={len(ev_gold):,})")

    # ---- calibration ----
    cal = {"done": False, "reason": "", "method": args.calibration}
    if args.calibration == "none":
        cal["reason"] = "--calibration none"
        print("calibration: skipped (--calibration none)")
    else:
        xs, ys = pos_pairs(dev_preds, dev_gold, dev_probs)
        spec = Cal.fit(xs, ys, method=args.calibration)
        if spec is None:
            cal["reason"] = f"too few/one-class dev positives ({len(xs)})"
            print(f"calibration: skipped ({len(xs)} positive-predicted dev examples / single class)")
        else:
            out_json = model_dir / "calibration.json"
            Cal.save(str(out_json), spec)
            exs, eys = pos_pairs(ev_preds, ev_gold, ev_probs)
            e_raw = Cal.ece(exs, eys)
            e_cal = Cal.ece([Cal.apply(spec, x) for x in exs], eys)
            cal.update(done=True, n_dev_pos=len(xs), n_eval_pos=len(exs), eval_name=ev_name,
                       ece_raw=e_raw, ece_cal=e_cal, json=f"{model_dir.name}/calibration.json")
            print(f"calibration: {args.calibration} on {len(xs)} dev positives -> {out_json}")
            print(f"  ECE on {ev_name} ({len(exs)} positives): raw {e_raw:.4f} -> calibrated {e_cal:.4f}")

    # ---- outputs listing (model dir files + data TSVs + this summary) ----
    desc = {"model.safetensors": "trained model weights", "config.json": "model config + label map",
            "tokenizer.json": "fast tokenizer", "tokenizer_config.json": "tokenizer config",
            "vocab.txt": "WordPiece vocab", "special_tokens_map.json": "special tokens",
            "training_args.bin": "training hyperparameters", "calibration.json": "probability calibrator"}
    outputs = []
    for f in sorted(p for p in model_dir.glob("*") if p.is_file()):
        outputs.append((f"{model_dir.name}/{f.name}", _fmt_bytes(f.stat().st_size),
                        desc.get(f.name, "")))
    for split in ("train.tsv", "dev.tsv", "test.tsv"):
        fp_ = data / split
        if fp_.exists():
            outputs.append((f"{data.name}/{split}", _fmt_bytes(fp_.stat().st_size),
                            "entity-blinded TSV (sentence + label)"))
    outputs.append((str(SUMMARY_PATH.relative_to(ROOT)).replace("\\", "/"), "-", "this summary"))

    write_summary_html({
        "task": args.task, "dataset": args.dataset, "model": str(model_dir),
        "seed": args.seed,
        "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU",
        "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "labels": labels, "neg": NEG, "confusion": cm, "metrics": metrics, "calibration": cal,
        "outputs": outputs,
        "splits": {"train": train_n, "dev": len(dev), "test": len(test) if test else None},
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", default="ppi", choices=["ppi", "chemprot", "gad", "ddi"])
    ap.add_argument("--dataset", default="bigbio/bioinfer", help="BigBIO HF dataset")
    ap.add_argument("--config", help="BigBIO config (default <name>_bigbio_kb)")
    ap.add_argument("--data", default="ppi_data", help="TSV dir (convert output / train input)")
    ap.add_argument("--model", default="ppi-biobert-re", help="checkpoint output dir")
    ap.add_argument("--val-frac", type=float, default=0.1,
                    help="dev fraction carved from train (REQUIRED by training; default 0.1)")
    ap.add_argument("--test-frac", type=float, default=0.0, help="carve a test split if none exists")
    ap.add_argument("--neg-ratio", type=float, default=None, help="cap negatives to N x positives")
    ap.add_argument("--epochs", type=float, default=3.0)
    ap.add_argument("--max-len", type=int, default=128, help="tokenizer max length (training + calibration)")
    ap.add_argument("--base-model", default=None,
                    help="base LM to fine-tune (train_re --model; default BioBERT, e.g. "
                         "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext for PubMedBERT)")
    ap.add_argument("--add-marker-tokens", action="store_true",
                    help="add the task's @MARKER$ entity tokens as dedicated special tokens (RE gain)")
    ap.add_argument("--lr", type=float, default=None, help="learning rate (train_re default 2e-5)")
    ap.add_argument("--batch-size", type=int, default=None, help="train batch size (train_re default 16)")
    ap.add_argument("--seed", type=int, default=42, help="training seed (reproducibility)")
    ap.add_argument("--convert-seed", type=int, default=0, help="conversion seed (sampling/shuffle)")
    ap.add_argument("--calibration", choices=["isotonic", "platt", "none"],
                    default="isotonic", help="dev calibrator fit on dev positives")
    ap.add_argument("--tune-threshold", action=argparse.BooleanOptionalAction, default=True,
                    help="tune the positive-decision threshold on dev to maximize micro-F1, then "
                         "report test metrics at that operating point (default on; --no-tune-threshold "
                         "to report plain argmax)")
    ap.add_argument("--skip-convert", action="store_true", help="reuse existing TSVs in --data")
    ap.add_argument("--skip-train", action="store_true", help="reuse existing checkpoint in --model")
    ap.add_argument("--gpus", default="0",
                    help="CUDA_VISIBLE_DEVICES for training/calibration: '0' (default, single GPU "
                         "-- recommended on Kaggle's 2x T4), '0,1' (both), 'all' (inherit), "
                         "'cpu'/'' (force CPU)")
    args = ap.parse_args()

    setup_device(args.gpus)

    if not args.skip_convert:
        convert(args)
    else:
        print("[1/3 CONVERT] skipped (--skip-convert)")
    if not args.skip_train:
        train(args)
    else:
        print("[2/3 TRAIN] skipped (--skip-train)")
    calibrate(args)

    print(f"\n{'=' * 70}\npipeline done -> model in {args.model}/  "
          f"(task={args.task})\n{'=' * 70}")
    print(f"Use it:  $env:RE_MODEL_{args.task.upper()} = (Resolve-Path .\\{args.model})")
    print("then:    python relation_extraction.py --normalize")


if __name__ == "__main__":
    main()
