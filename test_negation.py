"""Regression sentences for relation_extraction._neg_cues (which pairs a negation word denies).
Run: python test_negation.py      [[..]] marks the pair, {..} the sentence's other entities."""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "gpu_bundle"))
import relation_extraction as RX

CASES = [   # (negated?, sentence)
    (0, "[[AML]] patients with {JL1} expression showed higher {CD13} than those without {JL1} expression and were associated with {RUNX1}-{RUNX1T1} and {CBFB}-[[MYH11]] translocations."),
    (0, "[[AML]] patients with and without [[EZH2]] mutation were enrolled."),
    (0, "NB4 and U937 [[leukemia]] cells without an [[MLL]] translocation were less susceptible."),
    (1, "Knockdown of [[ITGB3]] reduced p-SYK without changing overall [[SYK]] levels."),
    (0, "Incubation with [[CD19]]-negative [[AML]] cell lines induced IL-2."),
    (0, "[[IFN]] signaling was enhanced in the absence of [[STAT3b]]."),
    (0, "[[PD0332991]] did not affect {CDK6} protein levels, but suppressed [[RB1]] phosphorylation."),
    (0, "Treatment with [[EPZ]] reduced [[FLT3]] expression with no observable effect on {HOXA9}."),
    (0, "Suppression of [[WDR5]] enhanced [[RPS24]] expression (data not shown)."),
    (0, "In depleted cells [[hSNF5]] could no longer activate [[INK4a]] expression."),
    (1, "Interestingly, the [[Myb]] expression was not reduced after the addition of [[SGC0946]]."),
    (1, "[[PLZF]] did not result in decreased expression of [[MYC]]."),
    (1, "[[HDAC1]] silencing exerts no significant effect on [[miR-146a]] expression."),
    (1, "Neither {AF4} nor [[AF9]] transactivated the [[EIbCAT]] reporter."),
    (1, "We did not observe a global change in [[CTCF]] levels upon [[KDM3B]] depletion."),
    (1, "[[EVI1]] high had lower rates of EFS but had no independent prognostic value for pediatric [[AML]]."),
    (1, "[[MIAT]] knockdown could reduce {MMP9} expression but has no effects on [[MMP2]]."),
    (0, "{MIAT} reduced [[MMP9]] expression but has no effects on [[MMP2]]."),
    (0, "Overexpression of [[USP7]], but not a catalytically inactive {USP7} mutant, led to increased [[MLL5]] stability."),
    (0, "[[KMT2A]] downregulation is not the only mechanism of [[MEN1]]-KO-induced toxicity."),
    (0, "Downregulation of [[BCL2]] may be in part but not entirely responsible for effects of [[MYBMIM]]."),
    (1, "[[AMPK]] deletion reduced {Glut1}, but not [[Glut4]], in AML cells."),
    (0, "[[AMPK]] deletion reduced [[Glut1]], but not {Glut4}, in AML cells."),
    (1, "[[AZD0530]] inhibited {SFKs} but not the activation of {Erk}, {Akt} and [[Stat5]]."),
]


def parse(marked):
    text, ents, pair, i = "", [], [], 0
    for m in re.finditer(r"\[\[(.+?)\]\]|\{(.+?)\}", marked):
        text += marked[i:m.start()]
        name = m.group(1) or m.group(2)
        ent = {"start": len(text), "end": len(text) + len(name), "label": "GENETIC"}
        text, i = text + name, m.end()
        ents.append(ent)
        if m.group(1):
            pair.append(ent)
    return text + marked[i:], pair, ents


if __name__ == "__main__":
    bad = 0
    for want, marked in CASES:
        text, (a, b), ents = parse(marked)
        got = RX._neg_cues(text, a, b, ents)
        if bool(got) != bool(want):
            bad += 1
            print(f"FAIL want {'negated' if want else 'positive'}, got {got}: {marked}")
    print(f"{len(CASES) - bad}/{len(CASES)} passed")
    sys.exit(bool(bad))
