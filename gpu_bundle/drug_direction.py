#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""drug_direction.py -- the agent of a signed relation between a drug and a gene or disease.

The models give no direction and direction.py takes the first-mentioned entity as the agent unless
the text has an agent cue, so "the BCL-2 inhibitor venetoclax" became `BCL2 -> inhibits ->
venetoclax`. For a signed, non-negated triple whose subject is a gene or disease and whose object
is a drug or chemical:

    "G inhibitor D", "inhibition of G with D"           -> D inhibits G
    drug response ("G ... resistance to D"), activates  -> associated (G does not activate D)
    drug response of a gene, inhibits                   -> left ("G inhibits D-induced apoptosis")
    natural metabolite, gene-loss or construct reading  -> left
    anything else                                       -> D is the agent, the sign stays

Two drugs (fix_pair): the models read a combination as a sign ("VTP50469 IC50 fell when combined
with MK-2206" became VTP-50469 inhibits MK-2206).

    combination, synergy, sensitivity, IC50, "A+B"      -> associated
    ... but "B hampered A-induced death" is antagonism  -> left
    "A inhibitors B", "A-induced ... blocked by B", and
    single sentences listed in PAIR_SWAP                -> B inhibits A

Used by high_confidence_g.py and audit_lists.py at graph time, after constructs.apply().
"""
import re

import direction as D

UP, DOWN = "upregulator/activator", "downregulator/inhibitor"
METABOLITES = {s.lower() for s in (
    "3',5'-cyclic AMP", "succinate", "oxalate(1-)", "calcium(2+)", "carbon dioxide", "water",
    "S-adenosyl-L-homocysteine", "S-adenosyl-L-methionine", "reactive oxygen species", "iron atom",
    "electron", "nucleoside", "hydrogen peroxide", "nitric oxide", "ATP", "ADP", "AMP", "GTP", "CTP",
    "lactate", "glucose", "cholesterol", "glutathione", "2-hydroxyglutarate", "acetyl-CoA", "dioxygen",
    "superoxide")}
_INHIB = r"(?:inhibitors?|antagonists?|degraders?|blockers?)"
# the text between gene and drug: "-specific inhibitor, ", " inhibitors such as ", " inhibitor A ("
INHIBITOR_OF = re.compile(r"[\s-]*(?:(?:specific|selective|small[- ]molecule|pathway)\s+)?" + _INHIB +
                          r"(?:\s*[,:(]\s*|\s+(?:such as|including|like|named|called)\s+|\s+)"
                          r"(?:[\w-]+\s*(?:\(|,|/|\band\b|\bor\b)\s*)?", re.I)
INHIBITION_WITH = re.compile(r"\s+(?P<n>inhibition|blockade)?\s*(?:with|by|using)\s+", re.I)
_INHIBITION_OF = re.compile(r"\b(?:inhibition|blockade) of (?:the )?$", re.I)
_RESP = (r"(?:(?:chemo)?resist\w*|(?:hyper|chemo)?-?sensiti\w*|refractor\w*|toler(?:a\w+)|suscept\w*|"
         r"responsive\w*|responds?|responded|protect\w*|rescu\w*|cytotoxicity|efficacy|survival|"
         r"(?<!in )(?<!damage )(?<!immune )respons\w*)")
RESPONSE = re.compile(r"\b" + _RESP + r"\b", re.I)
# "... doxorubicin or nutlin-3 resistance"
RESPONSE_AFTER = re.compile(r"(?:\s*(?:,|/|\b(?:and/or|and|or|plus)\b)\s*[\w-]+){0,3}[\s-]+" + _RESP
                            + r"\b", re.I)
INDUCED_BY = re.compile(r"\binduced\s+by\s+(?:[\w-]+\s+){0,3}$", re.I)
INDUCED = re.compile(r"\s?[-‐–]?\s?induced\b", re.I)

COMBINED = re.compile(r"\b(?:synerg|combin|additiv|co-?treat|sensiti[sz]|sensitivity|potentiat|IC50)"
                      r"|\b(?:abrogat|revers|overc[ao]m)\w* (?:[\w-]+ ){1,2}resistance", re.I)
PLUS = re.compile(r"\s*[+＋]\s*")                        # "Verapamil+DOX group"
ANTAGONIZE = re.compile(r"\b(?:ameliorat|attenuat|block|abolish|suppress|inhibit|reduc|disrupt|prevent|"
                        r"hamper|diminish|impair)\w*", re.I)
INHIBITORS = re.compile(r"\s+(?:inhibitors?|chelators?|antagonists?)\s+", re.I)
RESCUED_BY = re.compile(r"[-‐–\s]induced\b.*\b(?:rescued|blocked|abolished|prevented|reversed)\s+"
                        r"(?:\w+\s+)?(?:by|using|with)\s+(?:[\w-]+\s+){0,3}(?:[\w-]+,\s*)?", re.I)
# (paper, first drug, second drug): the second one blocks the first one's effect
PAIR_SWAP = {("PMC5226548", "cytarabine", "AC220"), ("PMC6726985", "Cycloheximide", "MG132")}


def fix_pair(t):
    """One signed drug-drug triple, in place. Returns None, "drug_agent" or "associated"."""
    m = re.search(r":(\d+)-(\d+):(\d+)-(\d+)$", t.get("pair_id") or "")
    if t["predicate"].get("text") not in (UP, DOWN) or t.get("polarity") == "negated" or not m \
            or t.get("direction") == "drug_agent":
        return None
    n, text = [int(x) for x in m.groups()], t.get("sentence") or ""
    a, b = sorted([(n[0], n[1]), (n[2], n[3])])
    gap = text[a[1]:b[0]]
    key = ((t.get("pmid") or "").split(".")[0], text[a[0]:a[1]], text[b[0]:b[1]])
    if INHIBITORS.fullmatch(gap) or RESCUED_BY.fullmatch(gap) or key in PAIR_SWAP:
        if t.get("direction") != "agent_cue":            # else the second drug is already first
            t["subject"], t["object"] = t["object"], t["subject"]
        t["predicate"]["text"], t["direction"] = DOWN, "drug_agent"
        return "drug_agent"
    if INDUCED.match(text, b[1]) and ANTAGONIZE.search(gap):
        return None
    if COMBINED.search(text) or PLUS.fullmatch(gap):
        t["predicate"]["text"], t["direction"] = "associated", "drug_combination"
        return "associated"
    return None


def fix(t, kind, drug):
    """Orient one triple in place: subject a gene / disease (`kind`), object the chemical node
    `drug`. Returns None (untouched), "drug_agent" or "associated"."""
    p = t["predicate"].get("text")
    m = re.search(r":(\d+)-(\d+):(\d+)-(\d+)$", t.get("pair_id") or "")
    if p not in (UP, DOWN) or t.get("polarity") == "negated" or not m or drug.lower() in METABOLITES \
            or t.get("direction") in ("loss", "construct", "drug_agent"):
        return None
    n, text = [int(x) for x in m.groups()], t.get("sentence") or ""
    a, b = sorted([(n[0], n[1]), (n[2], n[3])])
    gene_first = t.get("direction") != "agent_cue"       # direction.py swapped the mention order
    gap, c = text[a[1]:b[0]], b if gene_first else a
    w = INHIBITION_WITH.fullmatch(gap)
    if kind == "gene" and gene_first and (
            INHIBITOR_OF.fullmatch(gap) or w and (w.group("n") or _INHIBITION_OF.search(text[:a[0]]))):
        t["predicate"]["text"] = DOWN
    else:
        by = INDUCED_BY.search(gap)                      # "G blocked apoptosis induced by D"
        blocks = gene_first and (by or INDUCED.match(text, c[1])) and D.ACTIVE_VERB.search(
            gap[:by.start()] if by else gap)
        if RESPONSE.search(gap) or RESPONSE_AFTER.match(text, c[1]) or blocks:
            if p == UP:
                t["predicate"]["text"], t["direction"] = "associated", "drug_response"
                return "associated"
            if kind == "gene":
                return None
    t["subject"], t["object"], t["direction"] = t["object"], t["subject"], "drug_agent"
    return "drug_agent"
