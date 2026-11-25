#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression checks for gpu_bundle/drug_direction.py: [[gene or disease]] <<drug>>."""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "gpu_bundle"))
import drug_direction as DD

# text, kind of the subject, label before, expected (drug is the subject, label)
CASES = [
    ("The selective [[BCL-2]] inhibitor, <<venetoclax>>, exhibited in vitro additivity.", "gene", "inhibits", (True, "inhibits")),
    ("Cytotoxicity of the [[BCL-2]] inhibitor venetoclax (<<ABT-199>>).", "gene", "inhibits", (True, "inhibits")),
    ("Cells are hypersensitive to the [[BRD4]] inhibitor <<JQ1>>, which lowers MYC.", "gene", "activates", (True, "inhibits")),
    ("Pharmacological inhibition of [[STAT5]] with <<pimozide>> enhanced sensitivity to imatinib.", "gene", "inhibits", (True, "inhibits")),
    ("[[CDK9]]-specific inhibitors <<BAY1143572>> and ADZ4573", "gene", "inhibits", (True, "inhibits")),
    ("[[MEF2C]] and MEIS1 were repressed genes post <<VTP50469>> treatment.", "gene", "inhibits", (True, "inhibits")),
    ("Expression of [[PU.1]] increased in response to <<ATRA>>/TCP.", "gene", "activates", (True, "activates")),
    ("We found an increase in [[BAX]] expression under <<revumenib>>.", "gene", "activates", (True, "activates")),
    ("[[RAD17]] (DNA replication and DNA damage response) were repressed upon MI-538 + <<chidamide>> treatment.", "gene", "inhibits", (True, "inhibits")),
    ("[[PD-L1]] upregulation induced by <<DS-5272>> in cells", "gene", "activates", (True, "activates")),
    ("Reduced [[leukemia]] burden with BZ1, <<A485>>, and combination.", "disease", "inhibits", (True, "inhibits")),
    ("[[MLL]]-rearranged cells were sensitive to <<pinometostat>>.", "gene", "activates", (False, "associated")),
    ("[[AML]] with higher HOXA11 had an improved response to <<Ara-C>>.", "disease", "activates", (False, "associated")),
    ("The [[leukemia]] cells were sensitive to <<chloroquine>> treatment.", "disease", "inhibits", (True, "inhibits")),
    ("High [[TK1]] expression was associated with <<cytarabine>> resistance.", "gene", "activates", (False, "associated")),
    ("Knockdown of Rassf4 in [[Gata2]]-low cells resulted in <<doxorubicin>> or nutlin-3 resistance.", "gene", "activates", (False, "associated")),
    ("[[SALL4]] blocked apoptosis that was induced by <<ATRA>>.", "gene", "activates", (False, "associated")),
    ("[[FTO]] inhibits <<ATRA>>-induced differentiation of APL cells.", "gene", "inhibits", (False, "inhibits")),
    ("[[ABCB1]] reduced <<succinate>> levels.", "gene", "inhibits", (False, "inhibits")),
]
NAME = {"activates": DD.UP, "inhibits": DD.DOWN, "associated": "associated"}


def main():
    bad = 0
    for marked, kind, before, (swap, want) in CASES:
        text, spans = re.sub(r"\[\[|\]\]|<<|>>", "", marked), {}
        for m in re.finditer(r"\[\[(.+?)\]\]|<<(.+?)>>", marked):
            at = len(re.sub(r"\[\[|\]\]|<<|>>", "", marked[:m.start()]))
            spans["g" if m.group(1) else "d"] = (at, at + len(m.group(1) or m.group(2)))
        g, d = spans["g"], spans["d"]
        t = {"subject": "g", "object": "d", "predicate": {"text": NAME[before]}, "polarity": "positive",
             "direction": "reading_order", "sentence": text,
             "pair_id": f"x:0:{g[0]}-{g[1]}:{d[0]}-{d[1]}"}
        DD.fix(t, kind, text[d[0]:d[1]])
        if (t["subject"] == "d", t["predicate"]["text"]) != (swap, NAME[want]):
            bad += 1
            print(f"FAIL want {swap, want}, got {t['subject'], t['predicate']['text']}: {marked}")
    print(f"{len(CASES) - bad}/{len(CASES)} passed")
    return bad


if __name__ == "__main__":
    sys.exit(bool(main()))
