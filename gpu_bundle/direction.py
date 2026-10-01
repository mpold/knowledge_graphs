#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""direction.py -- which entity of a signed relation is the agent.

The relation models predict a TYPE and SIGN only (BioRED's Positive/Negative_Correlation
are symmetric), so relation_extraction.py writes the first-mentioned entity as `subject`.
That is wrong whenever the sentence is passive: "lactate levels ... abolished by AZD3965"
became `lactate -> inhibits -> AZD3965`. This module reverses a signed triple when the text
between the two entities ends in an agent cue that introduces the SECOND entity:

    "... reduced by X"   "... treated with X"   "... in response to X"   "... after X"

and is NOT one of the look-alikes that do not introduce an agent:

    "by enhancing X" / "by increased X"           (means clause)
    "marked by X" / "accompanied by X"            (observation)
    "A reduced MCT4 expression induced by X"      (A counteracts X's effect)
    "A knockdown abolished activation ... by X"   (A is the agent of a transitive verb)

Hand-checked on 100 triples the cue fired on (two random samples of 50): 0 of 11 wrong
flips made, 82 of 85 right flips kept. On the 2026-10-01 corpus it reverses ~12% of
signed triples. Only signed relations are touched: binds / associated / interacts carry
no direction. Imported by relation_extraction.py (new runs) and high_confidence_g.py
(triples from runs made before it existed).
"""

import re

SIGNED = {"upregulator/activator", "downregulator/inhibitor"}

AGENT_CUE = re.compile(
    r"(?:\b(?:by|upon|following|after|in response to)"
    r"|\b(?:treat\w*|pretreat\w*|cotreat\w*|incubat\w*|expos\w*|stimulat\w*|supplement\w*|challeng\w*)"
    r"\s+(?:\w+ly\s+)?(?:with|to))"
    r"\s+(?:[\w-]+\s+){0,3}$", re.I)
# "by" that does NOT introduce an agent: a means clause ("by enhancing X", "by increased X")
BY_MEANS = re.compile(r"^by\s+(?:\w+ing|increased|enhanced|decreased|reduced|elevated|diminished|"
                      r"upregulated|downregulated|higher|lower)\b", re.I)
# ... or an observation ("marked by X", "accompanied by X", "as shown by X")
BY_OBSERVED = re.compile(r"\b(?:marked|accompanied|characteri[sz]ed|followed|shown|detected|measured|"
                         r"assessed|determined|confirmed|evidenced|indicated|reflected|evaluated|"
                         r"quantified|visuali[sz]ed|analy[sz]ed|replaced)\s+by\s+(?:[\w-]+\s+){0,3}$", re.I)
# the first entity counteracts an effect the second caused: "A reduced X induced by B"
CAUSED_BY = re.compile(r"\b(?:induced|caused|mediated|triggered|driven|elicited|evoked|stimulated)\s+by\s+"
                       r"(?:[\w-]+\s+){0,3}$", re.I)
VERB_STEM = (r"(?:inhibit|reduc|decreas|increas|prevent|block|abolish|attenuat|revers|suppress|"
             r"hinder|enhanc|promot|lessen|restor|abrogat|rescu|counteract|antagoni[sz]|impair|"
             r"diminish|negat|alleviat|ameliorat|eliminat|activat|induc|repress|stimulat)")
ACTIVE_VERB = re.compile(r"(?<!\bwas )(?<!\bwere )(?<!\bis )(?<!\bare )(?<!\bbe )(?<!\bbeen )"
                         r"\b" + VERB_STEM + r"(?:e|ed|es|s)?\b", re.I)
VERB_TOKEN = re.compile(r"^" + VERB_STEM + r"(?:e|ed|es|s)?$", re.I)
AUX = {"was", "were", "is", "are", "be", "been", "being"}
# a token after the verb that marks it intransitive ("Smad3 increased observably when ...")
NO_OBJECT = re.compile(r"^(?:\w+ly|when|after|in|following|upon|by|at|with|compared|during|under|"
                       r"to|from|but|and|while|as|than|,|\.|;)$", re.I)


def _subject_is_agent(c):
    """The first entity is the agent of a transitive verb right after it:
    'HIF1A reduced the induction ...', 'LKB1 knockdown abolished activation ...'."""
    toks = re.findall(r"[\w-]+|[,.;]", c)
    for k, w in enumerate(toks[:3]):
        if w.lower() in AUX:
            return False
        if VERB_TOKEN.match(w):
            nxt = toks[k + 1] if k + 1 < len(toks) else ""
            return bool(nxt) and not NO_OBJECT.match(nxt)
    return False


def agent_cue(connecting):
    """The cue (e.g. 'by', 'treated with') if the second-mentioned entity is the agent, else None.
    `connecting` is the text strictly between the first and the second entity."""
    c = " ".join(connecting.split()) + " "
    m = AGENT_CUE.search(c)
    if not m:
        return None
    cue = m.group(0).strip()
    if cue.lower().startswith("by") and (BY_MEANS.match(cue) or BY_OBSERVED.search(c)):
        return None
    cb = CAUSED_BY.search(c)
    if cb and ACTIVE_VERB.search(c[:cb.start()]):
        return None
    if _subject_is_agent(c):
        return None
    return cue


def connecting_from_pair_id(sentence, pid):
    """Text between the two entities, from the sentence-local offsets at the end of a
    relation_extraction pair_id ('pmid:digest:a0-a1:b0-b1'); None if it cannot be parsed."""
    try:
        a, b = pid.rsplit(":", 2)[-2:]
        a_end = int(a.split("-")[1])
        b_start = int(b.split("-")[0])
    except (AttributeError, ValueError, IndexError):
        return None
    return sentence[a_end:b_start] if 0 <= a_end <= b_start <= len(sentence) else None


def orient(triple, connecting=None):
    """Set triple['direction'] and, for a signed relation whose agent is the second-mentioned
    entity, swap subject and object in place. Idempotent: a triple that already carries
    'direction' is left alone. Returns True if it was swapped."""
    if "direction" in triple:
        return False
    triple["direction"] = "reading_order"
    if triple.get("predicate", {}).get("text") not in SIGNED:
        return False
    if connecting is None:
        connecting = connecting_from_pair_id(triple.get("sentence", ""), triple.get("pair_id"))
    if connecting is None:
        return False
    cue = agent_cue(connecting)
    if not cue:
        return False
    triple["subject"], triple["object"] = triple["object"], triple["subject"]
    triple["direction"] = "agent_cue"
    triple["direction_cue"] = cue
    return True
