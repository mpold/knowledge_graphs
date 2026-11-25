#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""constructs.py -- signed relations of an siRNA / shRNA construct or of a gene that was lost.

Normalization maps "siMLL1" to the gene it targets (KMT2A), so "siMLL1 decreased AURKA" became
`KMT2A -> inhibits -> AURKA`, the opposite of what the knockdown shows. For a signed triple with a
construct at one end, apply() reads what the construct did to the other entity from the sentence:

    "siX decreased Y", "Y was reduced by siX", "a reduction in Y in the siX group"
        -> X activates Y  (the sign opposite to the verb's, X as the subject)
    group labels and comparisons ("shX + DOX group"), both ends constructs, no readable verb
        -> associated     (no sign can be read reliably)

The sign comes from the verb and not from flipping the model's label, because the model has
sometimes already allowed for the knockdown. Negated triples are left alone.

A gene named with its loss ("X knockdown", "loss of X", "X-deficient") gets the same reading, more
cautiously, as the model is right on most of these (fix_loss):

    gene or disease partner, verb readable       -> sign and subject from the sentence
    drug or chemical partner                     -> associated, unless the drug caused the loss ("loss of X increased resistance to
                                                    cytarabine" does not say X activates it)
    guarded sentence or no readable verb         -> left as the model had it

Also here, since the normalized files of finished runs already carry the errors: "shCtr" is a
control (but shCtr9 targets CTR9), and a few names that mean something else in one paper.
Used by high_confidence_g.py and audit_lists.py at graph time.
"""
import re

UP, DOWN = "upregulator/activator", "downregulator/inhibitor"
SIGNED = {UP, DOWN}

CONSTRUCT = re.compile(r"^(?:si|sh)[-‐– ]?(?=[A-Z0-9])|\b(?:si|sh)[- ]?RNAs?\b")
CONTROL = re.compile(r"^(?:sh|si|sg)-?ctr(?![A-Za-z0-9])", re.I)
# (paper, surface) -> HGNC symbol, None = not a gene there
PAPER_NAMES = {
    ("PMC5388432", "shM1"): "HOTAIRM1",     # an shRNA against HOTAIRM1, not MRLN
    ("PMC5388432", "A6"): None,             # HOXA6 in "HOXA4, A5, A6 and A7", not TWF1
}

_UP = (r"(?:up-?regulat|increas|enhanc|elevat|induc|activat|promot|stimulat|aggravat|augment|"
       r"rais|restor)")
_DOWN = (r"(?:down-?regulat|decreas|reduc|inhibit|suppress|attenuat|abolish|block|impair|diminish|"
         r"lower|repress|prevent|abrogat|curb|blunt)")
_END = r"(?:e[sd]?|s|es|ed|ing|t?ions?|ement)?"
EFFECT = re.compile(rf"\b(?:(?P<up>{_UP})|(?P<down>{_DOWN})){_END}\b", re.I)
_LY = r"(?:\w+ly\s+)?"
_LINK = r"(?:by|in|upon|after|following|with)\s+(?:(?:the|an?|those|cells?|of)\s+){0,2}"
HARD = re.compile(r"[;:]|\b(?:while|whereas|but|although|however|compared|than|contrast)\b", re.I)
SOFT = re.compile(r",|\b(?:and|through|via)\b", re.I)
# Y is the setting or the treatment rather than what changed
SETTING = re.compile(r"\b(?:upon|exposure|exposed|under|response|treated|treatment|presence|"
                     r"absence|following|after|during|group|related)\b", re.I)
COMBINED = re.compile(r"[+＋]")


def is_construct(e):
    return e.get("type") == "GENETIC" and bool(CONSTRUCT.search((e.get("text") or "").strip()))


def _sign(matches):
    """+1 / -1 for one effect word, or two that multiply ("induced inhibition of"); else None."""
    if not 1 <= len(matches) <= 2:
        return None
    s = 1
    for m in matches:
        s *= 1 if m.group("up") else -1
    return s


def effect(text, c, y):
    """What the construct at span c did to the entity at span y: +1, -1 or None (unreadable)."""
    if COMBINED.search(text[max(0, c[0] - 3):c[0]]) or COMBINED.search(text[c[1]:c[1] + 3]):
        return None
    if c[1] <= y[0]:                                     # construct first
        mid = text[c[1]:y[0]]
        if HARD.search(mid):
            return None
        for seg_end, seg in _segments(mid):              # nearest clause with an effect word
            found = list(EFFECT.finditer(seg))
            if found:
                if SETTING.search(mid[seg_end - len(seg) + found[-1].end():]):
                    return None
                return _sign(found)
        # "in siX cells, Y was reduced"
        m = re.match(r"[\w\s-]{0,30}?\b(?:was|were|is|are)\s+" + _LY + r"(?=" + EFFECT.pattern + ")",
                     text[y[1]:], re.I)
        return _sign([EFFECT.match(text, y[1] + m.end())]) if m else None
    mid = text[y[1]:c[0]]                                # construct second
    if HARD.search(mid):
        return None
    before = re.search(r"(" + EFFECT.pattern + r")\s+(?:of|in)\s+(?:the\s+)?$", text[:y[0]], re.I)
    noun = EFFECT.match(text, before.start(1)) if before else None
    tail = re.search(r"(?P<w>" + EFFECT.pattern + r")\s+" + _LY + _LINK + r"$", mid, re.I)
    if tail:
        verb = EFFECT.match(mid, tail.start("w"))
        return _sign([verb] + ([noun] if noun else []))
    if noun and re.fullmatch(r"\s*(?:(?:mRNA|protein|expression|levels?|activity)\s+){0,3}" + _LINK,
                             mid, re.I):
        return _sign([noun])                             # "a reduction in Y expression in the siX"
    return None


def _segments(mid):
    """(end offset, text) of the clauses of mid, the nearest to Y first."""
    cuts = [0] + [m.end() for m in SOFT.finditer(mid)]
    ends = [m.start() for m in SOFT.finditer(mid)] + [len(mid)]
    return [(e, mid[s:e]) for s, e in reversed(list(zip(cuts, ends)))]


# --- a gene named with its loss -------------------------------------------------------------
_LOSS = (r"(?:knock(?:ed)?[- ]?(?:down|out)|silenc\w+|deplet\w+|deletion|loss|ablation|"
         r"inactivation|deficien\w+|si ?RNA|shRNA)")
LOSS_AFTER = re.compile(r"[\s-]*(?:expression\s+|gene\s+)?" + _LOSS, re.I)
LOSS_BEFORE = re.compile(_LOSS + r"\s+(?:of\s+)?(?:the\s+)?$", re.I)
# the verb is not the plain outcome of the loss: rescue, conclusion, sensitivity, a second cause
GUARD = re.compile(r"\b(?:rescu\w*|revers\w*|demonstrat\w*|suggest\w*|indicat\w*|no longer|could not|"
                   r"hallmarks?|latency|genes? (?:negatively|positively)|overexpress\w*|"
                   r"restored? the growth|failed|unable|not|sensitiv\w*|suscept\w*|resist\w*|delay\w*|"
                   r"alleviat\w*|leading to|effects? on)\b"
                   r"|[-‐–](?:induced|driven|mediated)\s|\S\s+induced\s+(?:expression|levels?)\b", re.I)
# the entity is half of a fusion, or only the cause of what changed ("MLL-ENL-induced AML")
PART = re.compile(r"\s*(?:[-–/]|::)\s*[A-Z0-9]|[-‐–]\s?(?:induced|driven|mediated|rearranged|fusion)")
# the drug did it: it causes the loss ("SNDX-50469 treatment depleted Bcl-xL", "loss of X by DRUG"),
# unless it only switches a construct on
_DRUG_LOSS = re.compile(r"loss|deplet|inactivat", re.I)
_CAUSES = re.compile(r"\b(?:induc\w+|led to|leads to|result\w*|caus\w+|trigger\w+|treatment|"
                     r"associated with)\b", re.I)
_LOSS_BY = re.compile(r"\s*(?:\w+\s+)?(?:induced by|by|upon)\s+(?:\w+\s+){0,3}", re.I)
_NOT_DRUG_LOSS = re.compile(r"knock|shRNA|siRNA|CRISPR|rescu|sensitiv|resist|suscept", re.I)
SWITCH = re.compile(r"^(?:dox(?:ycycline)?|tam(?:oxifen)?|(?:4[\s-]*)?OHT)$", re.I)


def drug_effect(text, g, c, y):
    """True when the drug at span y lowers the gene at span g (c = with its loss wording)."""
    if SWITCH.match(text[y[0]:y[1]].strip()):
        return False
    gene, drug = re.escape(text[g[0]:g[1]]), re.escape(text[y[0]:y[1]])
    if re.search(rf"{gene}\s+inhibitors?\s+{drug}|inhibition of {gene}[\w\s]{{0,30}}\bby {drug}",
                          text, re.I):
        return True
    if not _DRUG_LOSS.search(text[c[0]:c[1]]):
        return False
    if y[1] <= c[0]:
        gap = text[y[1]:c[0]]
        return len(gap) <= 100 and bool(_CAUSES.search(gap)) \
            and not _NOT_DRUG_LOSS.search(text[max(0, y[0] - 40):c[0]])
    return bool(_LOSS_BY.fullmatch(text[c[1]:y[0]])) \
        and not _NOT_DRUG_LOSS.search(text[max(0, c[0] - 40):y[0]])


# "Y was decreased in HOS cells transfected with X siRNA / in cells with X knockdown"
_IN_CELLS = re.compile(r"\s+" + _LY + r"in\s+(?:[\w./+-]+\s+){0,4}?cells?\s+"
                       r"(?:(?:transfected|transduced|treated)\s+)?with\s+", re.I)
_FUSED = re.compile(r"\w\s?(?:[-‐–/]|::)\s?$")           # "MLL-AFF1 siRNA": the fusion was lost


def cells_effect(text, g, c, y):
    """Effect on Y (span y) in cells that lost the gene at span g (c = with its loss wording)."""
    if y[1] > c[0] or c[0] != g[0] or _FUSED.search(text[max(0, g[0] - 4):g[0]]):
        return None
    mid = text[y[1]:c[0]]
    found = list(EFFECT.finditer(mid))
    if HARD.search(mid) or len(found) != 1 or not _IN_CELLS.fullmatch(mid, found[0].end()):
        return None
    return _sign(found)


def loss_span(text, start, end):
    """Span of the gene at [start, end) with its loss wording, or None."""
    m = LOSS_AFTER.match(text, end)
    if m:
        return start, m.end()
    lo = max(0, start - 30)
    m = LOSS_BEFORE.search(text[lo:start])
    return (lo + m.start(), end) if m else None


def fix_loss(t):
    """As fix_sign, for a gene named with its loss. None (untouched), "signed" or "unsigned"."""
    if t["predicate"].get("text") not in SIGNED or t.get("polarity") == "negated" \
            or t.get("direction") in ("construct", "loss") \
            or any(is_construct(t[s]) for s in ("subject", "object")):
        return None
    m = re.search(r":(\d+)-(\d+):(\d+)-(\d+)$", t.get("pair_id") or "")
    if not m:
        return None
    n, text = [int(x) for x in m.groups()], t.get("sentence") or ""
    span = dict(zip(("subject", "object"), sorted([(n[0], n[1]), (n[2], n[3])])))
    if t.get("direction") == "agent_cue":                # direction.py swapped the mention order
        span = {"subject": span["object"], "object": span["subject"]}
    lost = {s: loss_span(text, *span[s]) for s in span
            if t[s].get("type") == "GENETIC" and not t[s].get("inhibitor_class")}
    lost = {s: v for s, v in lost.items() if v}
    if len(lost) != 1:
        return None
    side, c = next(iter(lost.items()))
    other = "object" if side == "subject" else "subject"
    y = span[other]
    if t[other].get("type") == "CHEMICAL" or t[other].get("inhibitor_class"):
        if drug_effect(text, span[side], c, y):
            return None
        t["direction"], t["predicate"]["text"] = "loss", "associated"
        return "unsigned"
    if GUARD.search(text[min(c[1], y[1]):max(c[0], y[0])]) or PART.match(text, y[1]) \
            or PART.match(text, c[1]):
        return None
    got = effect(text, c, y)
    if got is None:
        got = cells_effect(text, span[side], c, y)
    if got is None:
        return None
    t["direction"], t["predicate"]["text"] = "loss", DOWN if got > 0 else UP
    if side == "object":
        t["subject"], t["object"] = t["object"], t["subject"]
    return "signed"


def fix_names(t):
    """Control and per-paper names, in place. Returns the number of endpoints changed."""
    n, paper = 0, (t.get("pmid") or "").split(".")[0]
    for e in (t["subject"], t["object"]):
        txt = (e.get("text") or "").strip()
        if e.get("type") != "GENETIC":
            continue
        if CONTROL.search(txt) and e.get("control") != "yes":
            e["control"] = "yes"
            n += 1
        if (paper, txt) in PAPER_NAMES and e.get("hgnc_symbol") != PAPER_NAMES[paper, txt]:
            e["hgnc_symbol"] = PAPER_NAMES[paper, txt]
            n += 1
    return n


def fix_sign(t):
    """Re-sign one triple in place. Returns None (untouched), "signed" or "unsigned"."""
    sides = [s for s in ("subject", "object") if is_construct(t[s])]
    if not sides or t["predicate"].get("text") not in SIGNED or t.get("polarity") == "negated" \
            or t.get("direction") == "construct":
        return None
    got = None
    m = re.search(r":(\d+)-(\d+):(\d+)-(\d+)$", t.get("pair_id") or "")
    if m and len(sides) == 1:
        n = [int(x) for x in m.groups()]
        text = t.get("sentence") or ""
        spans = sorted([(n[0], n[1]), (n[2], n[3])])
        hit = [bool(CONSTRUCT.search(text[a:b])) for a, b in spans]
        if hit.count(True) == 1:
            c, y = (spans[0], spans[1]) if hit[0] else (spans[1], spans[0])
            got = effect(text, c, y)
    t["direction"] = "construct"
    if got is None:
        t["predicate"]["text"] = "associated"
        return "unsigned"
    t["predicate"]["text"] = DOWN if got > 0 else UP     # knocking X down did the opposite
    if sides[0] == "object":
        t["subject"], t["object"] = t["object"], t["subject"]
    return "signed"


def apply(triples):
    """All of the above, in place. Returns ((re-signed, unsigned) of the constructs, names changed,
    (read from the sentence, unsigned) of the lost genes)."""
    out = [fix_sign(t) for t in triples]
    names = sum(fix_names(t) for t in triples)
    lost = [fix_loss(t) for t in triples]
    return ((out.count("signed"), out.count("unsigned")), names,
            (lost.count("signed"), lost.count("unsigned")))
