#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression checks for gpu_bundle/antibodies.py: [[mention]], symbol, expected treatment."""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "gpu_bundle"))
import antibodies as A
import constructs as K

CASES = [
    ("combination of [[anti-4-1BB]] and anti-PD-1", "TNFRSF9", "anti-4-1BB antibody"),
    ("deletion with [[anti–PD-1]] treatment prolonged survival", "PDCD1", "anti-PD-1 antibody"),
    ("combining AraC with [[anti-IL-4]] also extended remission", "IL4", "anti-IL-4 antibody"),
    ("cells treated with [[anti-MN1]] siRNA (Figure 7C)", "MN1", None),
    ("detected in the [[anti-menin]] immunoprecipitate (Fig. 2C)", "MEN1", None),
    ("as shown by ChIP using an [[anti-AF9]] antibody", "MLLT3", None),
    ("expression of [[α-SMA]] and collagen I", "ACTA2", None),
    ("whereas [[CD19]] CAR T cells demonstrated effects", "CD19", "CD19 CAR-T"),
    ("[[CSPG4]]-CAR T cells evinced upregulation", "CSPG4", "CSPG4 CAR-T"),
    ("bicistronic [[CD19]]/CD22 CAR T-cell therapy", "CD19", "CD19 CAR-T"),
    ("concomitant [[CD19]] and CD22 CAR-T cells", "CD19", "CD19 CAR-T"),
    ("dual targeting of [[CD19]] and FLT3, which showed", "CD19", None),
    ("successful targeting with [[MSLN]]-targeted ADCs in vitro", "MSLN", "MSLN ADC"),
    ("[[CD47]] blockade with IL4 stimulation", "CD47", None),
]


def main():
    bad = 0
    for marked, sym, want in CASES:
        m = re.search(r"\[\[(.+?)\]\]", marked)
        text = marked.replace("[[", "").replace("]]", "")
        got = A.treatment(text, m.start(), m.start() + len(m.group(1)), m.group(1), sym)
        if got != want:
            bad += 1
            print(f"FAIL want {want}, got {got}: {marked}")
    t = {"subject": {"type": "GENETIC", "text": "anti-CD33"}, "object": {"type": "CHEMICAL", "text": "lintuzumab"},
         "sentence": "anti-CD33 lintuzumab", "pair_id": "x:0:0-9:10-20"}
    if A.tag([t], lambda e: "CD33") != (1, 1):
        bad += 1
        print("FAIL anti-CD33 lintuzumab")
    s = "ST2 blockade decreased systemic interferon-γ, IL-17, and IL-23 but increased IL-10 and IL-33 plasma levels."
    for partner, want in (("IL-10", K.DOWN), ("interferon-γ", K.UP)):
        t = {"subject": {"type": "GENETIC", "text": "ST2"}, "object": {"type": "GENETIC", "text": partner},
             "predicate": {"text": K.UP}, "sentence": s, "pmid": K.ST2_PAPER}
        K.fix_names(t)
        if t["predicate"]["text"] != want:
            bad += 1
            print("FAIL", partner)
    n = len(CASES) + 3
    print(f"{n - bad}/{n} passed")
    return bad


if __name__ == "__main__":
    sys.exit(bool(main()))
