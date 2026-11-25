"""Regression triples for constructs.py: the signed si/shRNA triples of the KMT2A graph, hand-checked.
Run: python test_constructs.py     (sentence, offsets of the pair, expected relation or None = unsigned)"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "gpu_bundle"))
import constructs as K
import controls as C

CASES = [
    ('A The addition of exogenous CTP rescued the shCTPS2 induced inhibition of CLL proliferation.',
     [44, 51, 74, 77], ('shCTPS2', 'activates', 'CLL')),
    ('Based on this analysis, treatments with the Prime-Fect/siFLT3 complex significantly inhibited the leukemia burden by day 13, in contrast to the Prime-Fect/C-siRNA group (Figure 8b,c).',
     [55, 61, 98, 106], ('siFLT3', 'activates', 'leukemia')),
    ('Besides, the expression of TRAP, c-Fos, and NFATc1 was suppressed by shKAT6A, whereas TRPV4 overexpression abrogated this decrease (Fig. 4D).',
     [44, 50, 69, 76], ('shKAT6A', 'activates', 'NFATc1')),
    ('Colony formation assays demonstrated a decrease in the proliferation of B-ALL cells in the MLL-AF9+shIGFBP7+DOX group compared to the MLL-AF9+shNC+DOX group.',
     [99, 107, 108, 111], None),
    ('Compared with the TNBS (sh‐MLL1) group, the TNBS (sh‐MLL1 + sh‐Gata4) group showed significantly more severe weight loss (Figure 9A) and a higher DAI score (Figure 9B).',
     [50, 57, 60, 68], None),
    ('Further, pathway analysis by gene set enrichment analysis revealed that the top gene signatures in shMBNL1 knockout are those related to FLT3 inhibitor treatment, including TP0903 and gilteritinib in AML cells (fig. S2F).',
     [99, 106, 173, 179], None),
    ('Furthermore, the expression of antioxidative genes, SOD and CAT, decreased by si-MLL1 (Figure 3D).',
     [60, 63, 78, 85], ('si-MLL1', 'activates', 'CAT')),
    ('In addition, flow cytometry analysis of the cell cycle indicated that the MLL-AF9+shIGFBP7+DOX group exhibited an increase in the G0/G1 phase and a decrease in the S and G2/M phases compared to the MLL-AF9+shNC+DOX group.',
     [82, 90, 91, 94], None),
    ('In this study, we demonstrated that si-MLL1 aggravated the apoptosis rate and increased the expression of cleaved-caspase 3 and Bax upon IPEC-J2 cell exposure to 1 μg/mL DON.',
     [36, 43, 106, 123], ('si-MLL1', 'inhibits', 'cleaved-caspase 3')),
    ('In this study, we demonstrated that si-MLL1 aggravated the apoptosis rate and increased the expression of cleaved-caspase 3 and Bax upon IPEC-J2 cell exposure to 1 μg/mL DON.',
     [36, 43, 128, 131], ('si-MLL1', 'inhibits', 'Bax')),
    ('In this study, we demonstrated that si-MLL1 aggravated the apoptosis rate and increased the expression of cleaved-caspase 3 and Bax upon IPEC-J2 cell exposure to 1 μg/mL DON.',
     [36, 43, 170, 173], None),
    ('Indeed, both MLL1 and H3K4me3 bound to the promoter of TEAD4, and the RNA-seq data revealed a 48.9% reduction in TEAD4 expression in the si-MLL1 group (Fig. 5F).',
     [113, 118, 137, 144], ('si-MLL1', 'activates', 'TEAD4')),
    ('Interestingly, knockdown of WDR5 increased AURKA protein levels, whereas siHOTTIP and siMLL1 decreased AURKA protein (Figure 4E).',
     [73, 81, 103, 116], ('siHOTTIP', 'activates', 'AURKA protein')),
    ('Interestingly, knockdown of WDR5 increased AURKA protein levels, whereas siHOTTIP and siMLL1 decreased AURKA protein (Figure 4E).',
     [86, 92, 103, 116], ('siMLL1', 'activates', 'AURKA protein')),
    ('Moreover, the sh‐HOTTIP + oe‐MLL1 group showed notable elevations in TLR4 mRNA and protein levels compared to the sh‐HOTTIP + oe‐NC group.',
     [14, 23, 69, 78], None),
    ('Results of immunofluorescence showed that OGD could upregulate MLL1 and cleaved caspase-3 levels at the same cell, while siMLL1 could reverse these effects significantly (Figure 2G).',
     [80, 89, 121, 127], None),
    ('The Cell Counting Kit-8 (CCK-8) assays demonstrated a decrease in cell viability in the MLL-AF9+shNC+DOX group compared to the MLL-AF9+shIGFBP7+DOX group.',
     [101, 104, 135, 143], None),
    ('Thus, intramuscular injection of si‐MLL1 inhibited the expansion of Pax7‐positive satellite cells through the down‐regulation of Myf5 and Cyclin D1, and subsequently blunted muscle regeneration.',
     [33, 40, 68, 72], ('si-MLL1', 'activates', 'Pax7')),
    ('Thus, intramuscular injection of si‐MLL1 inhibited the expansion of Pax7‐positive satellite cells through the down‐regulation of Myf5 and Cyclin D1, and subsequently blunted muscle regeneration.',
     [33, 40, 129, 133], ('si-MLL1', 'activates', 'Myf5')),
    ('Thus, intramuscular injection of si‐MLL1 inhibited the expansion of Pax7‐positive satellite cells through the down‐regulation of Myf5 and Cyclin D1, and subsequently blunted muscle regeneration.',
     [33, 40, 138, 147], ('si-MLL1', 'activates', 'Cyclin D1')),
    ('Transfection of Panc1 cells with siAURKA also induced Annexin V staining and PARP cleavage (Figure 4C) and inhibited Panc1 cell migration in a scratch assay (Figure 4D).',
     [33, 40, 54, 63], ('siAURKA', 'inhibits', 'Annexin V')),
    ('Transfection of Panc1 cells with siAURKA also induced Annexin V staining and PARP cleavage (Figure 4C) and inhibited Panc1 cell migration in a scratch assay (Figure 4D).',
     [33, 40, 77, 81], ('siAURKA', 'inhibits', 'PARP')),
    ('Upregulation of PRDM1 at 180 min after IR was markedly attenuated in shPTIP-transduced HL-60 cells compared with shCtrl-transduced cells (Fig. 3B).',
     [16, 21, 69, 75], ('shPTIP', 'activates', 'PRDM1')),
]
CONTROLS = [("shCtr", True), ("siCtr", True), ("sh-Ctr", True), ("shCtrl", None), ("shCtr9", False),
            ("CTR9", False), ("shCTRC", False)]     # None: left to controls.py
NAMES = [("PMC5388432.grobid.tei", "shM1", "MRLN", "HOTAIRM1"), ("PMC5388432.grobid.tei", "A6", "TWF1", None),
         ("PMC1.grobid.tei", "shM1", "MRLN", "MRLN")]
# a gene named with its loss: [[gene]] <<partner>>, partner type, label before (agent first in the
# text), expected (subject, relation) or None = untouched, "associated" = drug partner
LOSS = [
    ("[[MLL1]] knockout significantly downregulated <<LILRB4>> levels.", "GENETIC", "inhibits", ("MLL1", "activates")),
    ("Knockdown of [[HOTAIR]] induced expression of <<HOXD10>> and PCDHB5.", "GENETIC", "activates", ("HOTAIR", "inhibits")),
    ("[[Mettl3]] knockdown decreased <<Akt1>> mRNA levels.", "GENETIC", "activates", ("Mettl3", "activates")),
    ("<<Runx1>> upregulation in [[Runx2]]-deficient cells", "GENETIC", "activates", ("Runx2", "inhibits")),
    ("In turn, <<HoxA9>> expression was decreased upon loss of [[Dot1l]].", "GENETIC", "inhibits", ("Dot1l", "activates")),
    ("[[Fbxo22]] loss significantly represses development of <<AML>>.", "DISEASE", "inhibits", ("Fbxo22", "activates")),
    # the lost gene named as the cells' treatment; not when the siRNA targets a fusion
    ("Expression of <<c-Jun>>, c-Myc, Akt was downregulated in HOS cells transfected with [[MLLT3]]-siRNA.", "GENETIC", "activates", ("MLLT3", "activates")),
    ("The level of <<c-Jun>>, BCL-2 and c-Myc was significantly decreased in HOS cells with [[MLLT3]] knockdown.", "GENETIC", "activates", ("MLLT3", "activates")),
    ("<<Bax>> was increased in cells treated with [[MLLT3]] siRNA.", "GENETIC", "activates", ("MLLT3", "inhibits")),
    ("<<caspase 8>> expression was increased in cells treated with MLL–[[AFF1]] siRNA.", "GENETIC", "activates", None),
    ("<<caspase 8>> expression was increased in cells treated with MLL::[[AFF1]] siRNA.", "GENETIC", "activates", None),
    ("[[CEBPE]] depletion could reverse the activation of <<S100A9>>.", "GENETIC", "activates", None),
    ("[[ERa]] knockdown downregulates E2 and BPA induced expression of <<HOXC6>>.", "GENETIC", "activates", None),
    ("[[Mll1]] knockdown significantly increased sensitivity to <<Wnt3a>>.", "GENETIC", "activates", None),
    ("Loss of [[Necdin]] enhanced differentiation, leading to delayed development of <<leukemia>>.", "DISEASE", "activates", None),
    ("In [[BRG1]]-depleted cells hSNF5 could no longer activate <<INK4a>> expression.", "GENETIC", "activates", None),
    ("Knockdown of [[Taf1c]] decreased the expression of <<MLL>>-ENL target genes.", "GENETIC", "activates", None),
    ("Loss of [[Puma]] increased resistance to <<cytarabine>>.", "CHEMICAL", "activates", "associated"),
    ("<<doxycycline>> induced a 60% loss of [[MEN1]] expression.", "CHEMICAL", "inhibits", "associated"),
    ("The morphology induced by [[LIMK1]] depletion was markedly enhanced by <<palbociclib>>.", "CHEMICAL", "inhibits", "associated"),
    ("Cells were very sensitive to <<bortezomib>> and loss of [[Noxa]] conveyed resistance.", "CHEMICAL", "activates", "associated"),
    ("<<MRT199665>> treatment led to dose-dependent loss of [[MEF2C]] phosphorylation.", "CHEMICAL", "inhibits", None),
    ("<<VTP-50469>>-induced loss of [[IKAROS]] protein was observed.", "CHEMICAL", "inhibits", None),
    ("Knockdown of [[GPX4]] by shRNA and inhibition of GPX4 enzyme activity by <<RSL3>> induced ferroptosis.", "CHEMICAL", "inhibits", None),
]
CAT = {K.UP: "activates", K.DOWN: "inhibits"}


def run():
    bad = 0
    for sent, (a0, a1, b0, b1), want in CASES:
        ends = [{"type": "GENETIC", "text": sent[i:j]} for i, j in ((a0, a1), (b0, b1))]
        t = {"subject": ends[0], "object": ends[1], "predicate": {"text": K.UP}, "sentence": sent,
             "pair_id": f"x:0:{a0}-{a1}:{b0}-{b1}", "polarity": "positive"}
        K.fix_sign(t)
        cat = CAT.get(t["predicate"]["text"])
        flat = [c.replace("‐", "-").replace(" ", "") for c in (t["subject"]["text"], t["object"]["text"])]
        got = (flat[0], cat, flat[1]) if cat else None
        exp = want and (want[0].replace(" ", ""), want[1], want[2].replace(" ", ""))
        if got != exp:
            bad += 1
            print(f"FAIL want {exp}, got {got}: {sent[:90]}")
    for name, want in CONTROLS:
        if want is not None and (bool(K.CONTROL.search(name)) != want or (C.classify(name) == C.NTC_CLASS) != want):
            bad += 1
            print("FAIL control", name)
    for pmid, txt, before, want in NAMES:
        e = {"type": "GENETIC", "text": txt, "hgnc_symbol": before}
        K.fix_names({"pmid": pmid, "subject": e, "object": {"type": "DISEASE", "text": "x"}})
        if e["hgnc_symbol"] != want:
            bad += 1
            print("FAIL name", pmid, txt, e["hgnc_symbol"])
    for marked, kind, before, want in LOSS:
        found = {m.group(0)[0]: m for m in re.finditer(r"\[\[.+?\]\]|<<.+?>>", marked)}
        text = re.sub(r"\[\[|\]\]|<<|>>", "", marked)
        pair = []
        for k, m in sorted(found.items(), key=lambda z: z[1].start()):
            at = len(re.sub(r"\[\[|\]\]|<<|>>", "", marked[:m.start()]))
            pair.append(({"type": "GENETIC" if k == "[" else kind, "text": m.group(0)[2:-2]}, at))
        (a, i), (b, j) = pair
        t = {"subject": a, "object": b, "sentence": text, "polarity": "positive",
             "predicate": {"text": K.UP if before == "activates" else K.DOWN},
             "pair_id": f"x:0:{i}-{i + len(a['text'])}:{j}-{j + len(b['text'])}"}
        r = K.fix_loss(t)
        got = None if r is None else "associated" if r == "unsigned" else \
            (t["subject"]["text"], CAT[t["predicate"]["text"]])
        if got != want:
            bad += 1
            print(f"FAIL want {want}, got {got}: {marked}")
    n = len(CASES) + len(CONTROLS) + len(NAMES) + len(LOSS)
    print(f"{n - bad}/{n} passed")
    return bad


if __name__ == "__main__":
    sys.exit(bool(run()))
