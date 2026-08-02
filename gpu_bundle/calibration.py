#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""calibration.py -- map a relation model's raw softmax to a calibrated probability.

triples_strategy.html section 7: "Raw softmax is over-confident; fit Platt /
isotonic calibration of p_rel on the RE dev split so the score is a true
probability." Used by train_re.py (fit + save per checkpoint) and
relation_extraction.py (load + apply to p_rel before composite scoring).

Dependency-free:
  * isotonic (default): pool-adjacent-violators (PAV) -> monotonic step map,
    flexible, no parametric assumption.
  * platt: logistic calibration sigmoid(a*logit(p)+b), fit by gradient descent.

A calibrator is fit on (p, y) pairs where p is the softmax of the PREDICTED label
and y=1 iff that prediction was correct -- restricted to POSITIVE predictions
(the ones that become triples). The spec is a small JSON dict saved next to the
checkpoint as calibration.json.

NO CALIBRATOR MAY ASSERT CERTAINTY. Every output is clipped into
[CLIP, 1-CLIP]. Isotonic is a step function fit by PAV, so a pure top bin
(every dev example correct) produces a hard 1.0 and a pure bottom bin a hard
0.0 -- on the reference run ppi-biobert-re's 197 breakpoints collapsed onto 8
levels whose extremes were exactly 0.0 and 1.0, i.e. "this pair is certain",
asserted from a 780-row dev split. Clipping is applied BOTH when fitting and
when applying, so specs written before this change are corrected on load
without refitting.
"""

import bisect
import json
import math
from pathlib import Path

_EPS = 1e-6
# no calibrated probability may be 0 or 1: a finite dev split cannot justify
# certainty, and a hard 1.0 lets one model monopolize a cross-model max().
CLIP = 1e-3


def _clip(p):
    return min(1.0 - CLIP, max(CLIP, float(p)))


def _logit(p):
    p = min(1 - _EPS, max(_EPS, p))
    return math.log(p / (1 - p))


def _fit_isotonic(xs, ys):
    """PAV isotonic regression -> compact monotonic breakpoints {x:[...], y:[...]}."""
    pairs = sorted(zip(xs, ys))
    X = [float(x) for x, _ in pairs]
    Y = [float(y) for _, y in pairs]
    val, wt = [], []
    for y in Y:
        val.append(y)
        wt.append(1.0)
        while len(val) > 1 and val[-2] > val[-1]:        # pool adjacent violators
            v = (val[-2] * wt[-2] + val[-1] * wt[-1]) / (wt[-2] + wt[-1])
            w = wt[-2] + wt[-1]
            val.pop(); wt.pop()
            val[-1], wt[-1] = v, w
    fitted = []
    for v, w in zip(val, wt):
        fitted += [v] * int(round(w))
    bx, by = [], []
    for x, y in zip(X, fitted):
        if bx and x == bx[-1]:
            by[-1] = _clip(y)                            # same x -> keep monotone latest
        else:
            bx.append(x); by.append(_clip(y))
    return {"method": "isotonic", "x": bx, "y": by, "clip": CLIP}


def _fit_platt(xs, ys, iters=3000, lr=0.1):
    """Logistic calibration on the logit feature, fit by batch gradient descent."""
    f = [_logit(x) for x in xs]
    a, b, n = 1.0, 0.0, len(f)
    for _ in range(iters):
        ga = gb = 0.0
        for fi, yi in zip(f, ys):
            p = 1 / (1 + math.exp(-(a * fi + b)))
            e = p - yi
            ga += e * fi; gb += e
        a -= lr * ga / n
        b -= lr * gb / n
    return {"method": "platt", "a": a, "b": b, "clip": CLIP}


def fit(xs, ys, method="isotonic"):
    if len(xs) < 20 or len(set(ys)) < 2:
        return None                                      # too little signal to calibrate
    return _fit_platt(xs, ys) if method == "platt" else _fit_isotonic(xs, ys)


def apply(spec, p):
    """Calibrated probability for a raw softmax value p (identity if spec is None).

    The result is clipped into [CLIP, 1-CLIP] regardless of what the spec holds, so
    calibration.json files fitted before clipping existed are corrected here."""
    if not spec:
        return p
    if spec["method"] == "platt":
        return _clip(1 / (1 + math.exp(-(spec["a"] * _logit(p) + spec["b"]))))
    xs, ys = spec["x"], spec["y"]
    if p <= xs[0]:
        return _clip(ys[0])
    if p >= xs[-1]:
        return _clip(ys[-1])
    i = bisect.bisect_right(xs, p) - 1
    x0, x1, y0, y1 = xs[i], xs[i + 1], ys[i], ys[i + 1]
    return _clip(y0 if x1 == x0 else y0 + (y1 - y0) * (p - x0) / (x1 - x0))


def ceiling(spec):
    """Highest calibrated probability this calibrator can ever emit.

    A threshold above this is unreachable for the model behind `spec` -- which is how
    a 0.95 cutoff silently became gene-only once BioRED got a Platt calibrator whose
    ceiling sits below it. Used by the graph step to warn instead of quietly filtering."""
    if not spec:
        return 1.0
    return apply(spec, 1.0 - _EPS) if spec["method"] == "platt" else _clip(max(spec["y"]))


def ece(ps, ys, bins=10):
    """Expected calibration error: weighted |accuracy - confidence| over bins."""
    n = len(ps)
    if not n:
        return 0.0
    e = 0.0
    for b in range(bins):
        lo, hi = b / bins, (b + 1) / bins
        idx = [i for i, p in enumerate(ps) if (lo < p <= hi) or (b == 0 and p <= hi)]
        if not idx:
            continue
        conf = sum(ps[i] for i in idx) / len(idx)
        acc = sum(ys[i] for i in idx) / len(idx)
        e += len(idx) / n * abs(acc - conf)
    return e


def save(path, spec):
    Path(path).write_text(json.dumps(spec, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def load(path):
    p = Path(path)
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
