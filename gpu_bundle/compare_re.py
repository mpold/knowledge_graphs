#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compare_re.py -- put two RE checkpoints side by side on the SAME entity pairs.

BioRED_task_mapping.md section 8, step 5: "Compare against the current PPI checkpoint
on the same documents -- the question is not raw F1 but whether the extra edge types
are *right*." Dev-set F1 cannot answer that (the two models are scored on different
corpora with different label sets), so the comparison has to happen on OUR sentences.

INPUT
    TRIPLES/triples_re.json written by
        python relation_extraction.py --normalize --route-mode additive
    In additive mode every applicable checkpoint scores each pair, and each triple
    carries `pair_id` (stable per pmid+sentence+offsets) and `predicate.model` -- so
    the two verdicts on one pair join on pair_id.

WHAT IT REPORTS
    * coverage      : pairs kept by both models / by only one (a pair missing for a
                      model means that model called it negative OR its composite score
                      fell below RE_MIN_SCORE -- the two are not distinguishable here)
    * agreement     : label cross-tab on the pairs both kept, i.e. what the typed model
                      calls the edges the binary model calls `interacts`
    * sign          : how many edges the BioRED-style model gives a SIGNED label
                      (upregulator/activator, downregulator/inhibitor, binds) -- the
                      whole point of the corpus swap
    * score deltas  : mean composite score per model on the shared pairs
    * samples       : disagreements and BioRED-only edges, for a hand-read

OUTPUT
    summaries/compare_re.html  (+ a text summary on stdout)

Run::  python compare_re.py
       python compare_re.py --baseline ppi-biobert-re --candidate biored-biobert-re
       python compare_re.py --input TRIPLES/triples_re.json --samples 40
"""

import argparse
import html
import json
import sys
from collections import Counter
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent
DEFAULT_IN = ROOT / "TRIPLES" / "triples_re.json"
HTML_OUT = ROOT / "summaries" / "compare_re.html"

# labels that actually add sign/direction to the graph (the BioRED upgrade); anything
# else is unsigned coverage. Compared lower-cased against predicate.text.
SIGNED = {"upregulator/activator", "downregulator/inhibitor", "binds",
          "positive_correlation", "negative_correlation", "bind"}


def load(path):
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except FileNotFoundError:
        sys.exit(f"no {path} -- run relation_extraction.py first.")
    except Exception as e:
        sys.exit(f"could not read {path}: {type(e).__name__}: {e}")
    if not isinstance(data, list):
        sys.exit(f"{path} is not a list of triples.")
    return data


def pick_models(triples, baseline, candidate):
    """(baseline, candidate) checkpoint names. Defaults: the BioRED-ish model is the
    candidate, the other one the baseline."""
    models = [m for m, _ in Counter(t["predicate"].get("model") for t in triples).most_common()
              if m]
    if baseline and candidate:
        return baseline, candidate
    if len(models) < 2:
        return (models[0] if models else None), None
    cand = candidate or next((m for m in models if "biored" in m.lower()), None)
    base = baseline or next((m for m in models if m != cand), None)
    if cand is None:                      # no biored-looking name: just take the top two
        base, cand = models[0], models[1]
    return base, cand


def index_by_pair(triples, model):
    """pair_id -> triple, for one model. Keeps the highest-scoring triple if a pair
    somehow repeats (it should not)."""
    out = {}
    for t in triples:
        if t["predicate"].get("model") != model:
            continue
        pid = t.get("pair_id")
        if not pid:
            continue
        if pid not in out or t["score"] > out[pid]["score"]:
            out[pid] = t
    return out


def describe(t):
    return (f'{t["subject"]["text"]} -[{t["predicate"]["text"]}]-> {t["object"]["text"]}',
            f'{t["subject"]["type"]}-{t["object"]["type"]}')


def compare(triples, base, cand, n_samples):
    A, B = index_by_pair(triples, base), index_by_pair(triples, cand)
    shared = sorted(set(A) & set(B))
    only_a = sorted(set(A) - set(B))
    only_b = sorted(set(B) - set(A))

    cross = Counter()
    agree_pairs, disagree = [], []
    for pid in shared:
        la, lb = A[pid]["predicate"]["text"], B[pid]["predicate"]["text"]
        cross[(la, lb)] += 1
        (agree_pairs if la.lower() == lb.lower() else disagree).append(pid)

    def mean_score(idx, pids):
        vals = [idx[p]["score"] for p in pids]
        return sum(vals) / len(vals) if vals else 0.0

    signed_b = [p for p in B if B[p]["predicate"]["text"].lower() in SIGNED]
    stats = {
        "n_base": len(A), "n_cand": len(B), "shared": len(shared),
        "only_base": len(only_a), "only_cand": len(only_b),
        "mean_base_shared": mean_score(A, shared), "mean_cand_shared": mean_score(B, shared),
        "signed_cand": len(signed_b),
        "signed_shared": len([p for p in signed_b if p in A]),
        "pairs_base": Counter(f'{A[p]["subject"]["type"]}-{A[p]["object"]["type"]}' for p in A),
        "pairs_cand": Counter(f'{B[p]["subject"]["type"]}-{B[p]["object"]["type"]}' for p in B),
        "labels_base": Counter(A[p]["predicate"]["text"] for p in A),
        "labels_cand": Counter(B[p]["predicate"]["text"] for p in B),
    }
    samples = {
        "disagree": [(A[p], B[p]) for p in disagree[:n_samples]],
        "only_cand": [B[p] for p in sorted(only_b, key=lambda p: -B[p]["score"])[:n_samples]],
        "only_base": [A[p] for p in sorted(only_a, key=lambda p: -A[p]["score"])[:n_samples]],
    }
    return stats, cross, samples


def render_html(base, cand, stats, cross, samples, src):
    esc = html.escape

    def tbl(rows, h1, h2):
        body = "".join(f'<tr><td>{esc(str(a))}</td><td class="num">{b:,}</td></tr>' for a, b in rows)
        return f'<table><tr><th>{h1}</th><th class="num">{h2}</th></tr>{body}</table>'

    cross_rows = "".join(
        f'<tr><td><code>{esc(a)}</code></td><td><code>{esc(b)}</code></td>'
        f'<td class="num">{c:,}</td><td>{"same" if a.lower() == b.lower() else ""}</td></tr>'
        for (a, b), c in cross.most_common())
    dis_rows = "".join(
        f'<tr><td>{esc(a["subject"]["text"])}</td><td>{esc(a["object"]["text"])}</td>'
        f'<td><code>{esc(a["predicate"]["text"])}</code> <span class="num">{a["score"]:.2f}</span></td>'
        f'<td><code>{esc(b["predicate"]["text"])}</code> <span class="num">{b["score"]:.2f}</span></td>'
        f'<td class="sent">{esc(a["sentence"][:240])}</td></tr>' for a, b in samples["disagree"])

    def only_rows(items):
        return "".join(
            f'<tr><td>{esc(t["subject"]["text"])}</td><td><code>{esc(t["predicate"]["text"])}</code></td>'
            f'<td>{esc(t["object"]["text"])}</td><td class="num">{t["score"]:.2f}</td>'
            f'<td class="sent">{esc(t["sentence"][:240])}</td></tr>' for t in items)

    style = (" body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1100px;"
             "margin:2rem auto;padding:0 1rem;color:#1a1a1a;} h1{font-size:1.45rem;}"
             " h2{font-size:1.15rem;margin-top:1.6rem;border-bottom:1px solid #ddd;padding-bottom:.3rem;}"
             " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;width:100%;}"
             " th,td{border:1px solid #ccc;padding:.3rem .6rem;text-align:left;vertical-align:top;}"
             " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
             " .big{font-size:2rem;font-weight:700;} code{background:#f3f3f3;padding:1px 5px;border-radius:3px;}"
             " .sent{color:#444;font-size:.9em;} .headline{background:#eef4fb;border:1px solid #cdddf0;"
             "border-radius:8px;padding:.8rem 1rem;margin:1rem 0;}"
             " .note{background:#fff8e1;border-left:4px solid #f0c040;padding:.6rem 1rem;border-radius:4px;}")
    HTML_OUT.parent.mkdir(parents=True, exist_ok=True)
    HTML_OUT.write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>compare_re &mdash; {esc(base)} vs {esc(cand)}</title><style>{style}</style></head><body>
<h1>Two RE checkpoints on the same entity pairs</h1>
<p>Baseline <code>{esc(base)}</code> vs candidate <code>{esc(cand)}</code>, joined on
<code>pair_id</code> from <code>{esc(src)}</code> (written by
<code>relation_extraction.py --route-mode additive</code>). The question this answers is not
which model has the higher corpus F1 &mdash; they are evaluated on different corpora with
different label sets &mdash; but whether the candidate's extra edge types are <em>right</em>
on our sentences.</p>
<div class="headline"><span class="big">{stats['shared']:,}</span> pairs scored by both
&middot; {stats['only_base']:,} baseline-only &middot; {stats['only_cand']:,} candidate-only
&middot; <strong>{stats['signed_cand']:,}</strong> candidate edges carry a SIGNED label</div>
<p class="note">A pair missing for one model means that model predicted its negative class
<em>or</em> its composite score fell below <code>RE_MIN_SCORE</code> &mdash; those two cases are
not distinguishable from this file.</p>

<h2>Coverage</h2>
{tbl([("pairs kept by " + base, stats['n_base']), ("pairs kept by " + cand, stats['n_cand']),
      ("kept by both", stats['shared']), ("only " + base, stats['only_base']),
      ("only " + cand, stats['only_cand']),
      ("signed candidate edges", stats['signed_cand']),
      ("signed candidate edges the baseline also kept", stats['signed_shared'])],
     "measure", "pairs")}
<p>Mean composite score on the shared pairs: <code>{esc(base)}</code>
{stats['mean_base_shared']:.3f} &middot; <code>{esc(cand)}</code> {stats['mean_cand_shared']:.3f}</p>

<h2>Label agreement on the shared pairs</h2>
<table><tr><th>{esc(base)}</th><th>{esc(cand)}</th><th class="num">pairs</th><th></th></tr>{cross_rows}</table>

<h2>Entity-type pairs covered</h2>
<div style="display:flex;gap:1rem;flex-wrap:wrap">
<div style="flex:1;min-width:280px"><h3>{esc(base)}</h3>{tbl(stats['pairs_base'].most_common(), "pair", "triples")}</div>
<div style="flex:1;min-width:280px"><h3>{esc(cand)}</h3>{tbl(stats['pairs_cand'].most_common(), "pair", "triples")}</div>
</div>

<h2>Labels</h2>
<div style="display:flex;gap:1rem;flex-wrap:wrap">
<div style="flex:1;min-width:280px"><h3>{esc(base)}</h3>{tbl(stats['labels_base'].most_common(), "label", "triples")}</div>
<div style="flex:1;min-width:280px"><h3>{esc(cand)}</h3>{tbl(stats['labels_cand'].most_common(), "label", "triples")}</div>
</div>

<h2>Disagreements (hand-read these)</h2>
<table><tr><th>subject</th><th>object</th><th>{esc(base)}</th><th>{esc(cand)}</th><th>sentence</th></tr>{dis_rows}</table>

<h2>Kept only by {esc(cand)}</h2>
<table><tr><th>subject</th><th>relation</th><th>object</th><th class="num">score</th><th>sentence</th></tr>
{only_rows(samples['only_cand'])}</table>

<h2>Kept only by {esc(base)}</h2>
<table><tr><th>subject</th><th>relation</th><th>object</th><th class="num">score</th><th>sentence</th></tr>
{only_rows(samples['only_base'])}</table>
</body></html>""", encoding="utf-8")
    print(f"Wrote {HTML_OUT}")


def main():
    ap = argparse.ArgumentParser(description="compare two RE checkpoints on the same entity pairs")
    ap.add_argument("--input", default=str(DEFAULT_IN), help="scored triples (default TRIPLES/triples_re.json)")
    ap.add_argument("--baseline", default=None, help="checkpoint name (default: the non-BioRED one)")
    ap.add_argument("--candidate", default=None, help="checkpoint name (default: the BioRED one)")
    ap.add_argument("--samples", type=int, default=25, help="rows per sample table (default 25)")
    args = ap.parse_args()

    triples = load(args.input)
    base, cand = pick_models(triples, args.baseline, args.candidate)
    if not cand:
        print(f"only one checkpoint in {args.input} ({base or 'none'}) -- nothing to compare.\n"
              f"Re-run: python relation_extraction.py --normalize --route-mode additive "
              f"with both RE_MODEL_PPI and RE_MODEL_BIORED set.")
        return
    if not any(t.get("pair_id") for t in triples):
        sys.exit(f"{args.input} has no pair_id fields -- it was written by an older "
                 f"relation_extraction.py. Re-run the RE step.")

    stats, cross, samples = compare(triples, base, cand, args.samples)
    print(f"baseline={base}  candidate={cand}  (source {args.input})")
    print(f"  kept: {stats['n_base']:,} vs {stats['n_cand']:,};  both {stats['shared']:,}, "
          f"only-baseline {stats['only_base']:,}, only-candidate {stats['only_cand']:,}")
    print(f"  mean score on shared pairs: {stats['mean_base_shared']:.3f} vs {stats['mean_cand_shared']:.3f}")
    print(f"  candidate signed edges: {stats['signed_cand']:,} "
          f"({stats['signed_shared']:,} of them on pairs the baseline also kept)")
    same = sum(c for (a, b), c in cross.items() if a.lower() == b.lower())
    if stats["shared"]:
        print(f"  identical label on {same:,}/{stats['shared']:,} shared pairs "
              f"({100 * same / stats['shared']:.1f}%) -- the rest is where the typed model "
              f"says something the binary one cannot")
    for (a, b), c in cross.most_common(8):
        print(f"    {a:>28.28s} -> {b:<28.28s} {c:>7,}")
    render_html(base, cand, stats, cross, samples, args.input)


if __name__ == "__main__":
    main()
