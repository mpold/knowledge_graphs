#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""antibodies.py -- an antibody or CAR-T treatment is not its target gene.

Normalization maps "anti-4-1BB" to TNFRSF9 and the "CD19" of "CD19 CAR-T cells" to CD19, so what
the treatment did was drawn on the gene. tag() names the treatment instead, through the same
`inhibitor_class` field that makes "BRAFi" a drug class:

    "anti-PD-1 treatment"            -> anti-PD-1 antibody
    "CD19 CAR T cells", "CD19/CD22 CAR-T", "CD19 and CD22 CAR-T"  -> CD19 CAR-T, CD22 CAR-T
    "MSLN-targeted ADCs"             -> MSLN ADC
    "anti-CD33 lintuzumab"           -> one drug: the pair is dropped

Not tagged: siRNA / antisense ("anti-MN1 siRNA"), lab reagents (immunoprecipitation, ChIP, blots,
staining). "X blockade" keeps the gene. Used by high_confidence_g.py and audit_lists.py at graph
time, before constructs.apply().
"""
import re

ANTI = re.compile(r"^anti\s?[-‐–]?\s?(?=[A-Z0-9])(?P<x>.+?)(?:\s*antibod(?:y|ies))?$")
NOT_DRUG = re.compile(r"\s*[-‐–]?\s*(?:si ?RNA|shRNA|antisense|oligo|miR|immuno(?:precip|blot|stain|fluor)"
                      r"|IP\b|ChIP|western|signal|stain|beads?\b|conjugated (?:beads|resin))", re.I)
LAB = re.compile(r"(?:ChIP|immunoprecipitat\w*|immunoblot\w*|western blot\w*|immunostain\w*|stained|blotted|"
                 r"probed|detected)\b[^.;]{0,40}\b(?:with|using|by)\s+(?:an?\s+|the\s+)?$", re.I)
_MATE = r"(?:(?:\s?[/‐–-]\s?|\s+and\s+)[A-Z][A-Za-z0-9]+)?"
CAR = re.compile(_MATE + r"[\s‐–-]*(?:[Bb]ispecific\s+)?(?:CARs?\b|[Cc]himeric antigen receptor)")
ADC = re.compile(r"[-‐–\s](?:targeted|targeting|directed)\s+(?:ADCs?\b|antibody[-‐– ]drug conjugates?)")
MAB = re.compile(r"(?:mab|vedotin|ozogamicin|mafodotin)$", re.I)


def _flat(s):
    return re.sub(r"\W", "", s)


def treatment(text, start, end, surface, symbol):
    """Name of the treatment the mention at [start, end) stands for, or None."""
    m = ANTI.match(surface.strip())
    if m:
        if NOT_DRUG.match(text, end) or LAB.search(text[max(0, start - 80):start]):
            return None
        return "anti-" + re.sub(r"\s", "", m.group("x")).replace("‐", "-").replace("–", "-") + " antibody"
    if CAR.match(text, end):
        return symbol + " CAR-T"
    if ADC.match(text, end):
        return symbol + " ADC"
    return None


def tag(triples, symbol_of):
    """Tag treatment endpoints in place; drop "anti-X DRUG-mab" pairs. `symbol_of(endpoint)` gives
    the gene node's symbol or None. Returns (endpoints tagged, pairs dropped)."""
    tagged, keep = 0, []
    for t in triples:
        m = re.search(r":(\d+)-(\d+):(\d+)-(\d+)$", t.get("pair_id") or "")
        n, text, hit = [int(x) for x in m.groups()] if m else [], t.get("sentence") or "", set()
        spans = sorted([(n[0], n[1]), (n[2], n[3])]) if n else []
        ends = [t["subject"], t["object"]]
        for e in ends:
            if e.get("type") != "GENETIC" or e.get("inhibitor_class"):
                continue
            at = [k for k, (a, b) in enumerate(spans) if k not in hit
                  and _flat(text[a:b]) == _flat(e.get("text") or "")]
            sym = symbol_of(e) if at else None
            name = sym and treatment(text, *spans[at[0]], e["text"], sym)
            if name:
                hit.add(at[0])
                e["inhibitor_class"] = name
                tagged += 1
        if hit and len(spans) == 2 and re.fullmatch(r"\s*(?:antibod(?:y|ies)\s*)?[(,]?\s*",
                                                    text[spans[0][1]:spans[1][0]]) \
                and any(e.get("type") == "CHEMICAL" and MAB.search(e.get("text") or "") for e in ends):
            continue
        keep.append(t)
    dropped = len(triples) - len(keep)
    triples[:] = keep
    return tagged, dropped
