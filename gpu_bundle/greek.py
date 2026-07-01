#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expand single-character Greek symbols to spelled-out Greek words, build a JSON
library keyed by the original clean_genetic_ne value, then split that library
into multi-protein complexes (mapped to verified HGNC symbols) vs the rest.

Inputs : GENETIC/greek_clean_genetic_ne.tsv
         (columns: clean_genetic_ne, occurrences, n_source_forms)
         databases/hgnc_complete_set_2026-05-01.json
Outputs:
  * GENETIC/full_greek.json     -- every entity, original -> {expanded, occ.}
  * GENETIC/greek_complex.json  -- the subset whose greek_expanded value names
                                     a multi-protein complex, each mapped to its
                                     verified HGNC subunit symbols (Stage A)
  * GENETIC/greek.json          -- non-complex values whose greek_expanded form
                                     is string-equal to a symbol/alias/prev/name
                                     field of exactly ONE HGNC gene (Stage B)
  * GENETIC/greek_ambiguous.json -- as greek.json but string-equal to MULTIPLE
                                     genes' fields
  * GENETIC/unmatched_greek.json -- entries matched by neither stage
  * GENETIC/greek_cosine.json   -- ADVISORY overlay: unmatched multi-word values
                                     given a suggested gene by word-order cosine
                                     vs HGNC name fields (Stage C; not a partition
                                     member -- these stay in unmatched_greek.json)

================================ MAPPING STRATEGY ============================
Goal: map each greek_expanded value to HGNC gene symbol(s) from
hgnc_complete_set, splitting into complex / single / ambiguous / unmatched.

(0) HGNC INDEXING. Build value -> {(symbol, hgnc_id)} indexes over the six fields
    symbol, alias_symbol, prev_symbol, name, alias_name, prev_name (array fields
    indexed per element), under four key normalizations:
      - literal (case-sensitive)            idx
      - casefold                            idx_ci
      - hsfold  = casefold + '-' -> ' '     idx_hs   (hyphen<->space)
      - delsep  = casefold + drop '-' & ' ' idx_del  (separators removed)
    Each gene is keyed by its approved symbol + hgnc_id.

(A) STAGE A -- multi-protein complexes. classify_complex() flags a value that
    NAMES a complex (an assembly of >=2 distinct gene products, e.g. NF-kappaB =
    NFKB1/RELA, integrin alphavbeta3) and maps it to its verified subunit genes
    -> greek_complex.json. Single subunits (IKKbeta=IKBKB) and coordinate paralog
    mentions (GSK3alpha/beta) are deliberately NOT treated as complexes here.

(B) STAGE B -- single/ambiguous match cascade (match_cascade), tried in order;
    each pass runs only on what the previous left unmatched and the winning pass
    is recorded as match_mode. Genes are unioned across the six fields; the
    highest-priority field that hits is recorded as match_field.
      1.  key case-insensitive  the ORIGINAL full_greek.json key (pre-expansion
                                surface form), casefolded, is exactly equal to a
                                symbol/alias/prev/name field (idx_ci) -- tried
                                FIRST, before any glyph expansion / transform
                                (tau->MAPT, beta-actin->ACTB). A deliberate curated
                                override still wins for its few keys (Rho->RHOA,
                                not the coincidental rhodopsin RHO).
      2.  curated synonym       hand-curated CURATED_SYNONYMS for names HGNC has
                                no field for (HIF-2alpha->EPAS1, NF-kappaB p65->
                                RELA, and whole families: PKC/PLC/Galpha/tubulins
                                /PPAR/...). str value -> single, list -> ambiguous.
      3.  case-sensitive        exact equality (idx, on the expanded value)
      4.  case-insensitive      casefold equality (idx_ci, on the expanded value)
      5.  hyphen/whitespace     hsfold equality (idx_hs): NF-kappaB == NF kappaB
      6.  greek-letter split    split a trailing Greek word + optional hyphen/
                                space: HIF1alpha / gamma-H2AX -> 'HIF1 alpha'
      7.  letter-number split   XXX-<n><greek> -> 'XXX<n> <greek>' (HIF-2alpha)
      8.  letter-number (ci)    case-insensitive-prefix variant, hyphen-folded
      9.  prefix-stripped       drop leading si/sh/p/pro, RE-CASCADE the remainder
      10. IFN expansion         IFN-gamma / IFNgamma -> 'interferon gamma'
      11. alpha-stripped        leading 'alpha' = anti-X (alphaPD-1->PD-1);
                                _STRIP_DENYLIST guards intrinsic-alpha (alphaSMA...)
      12. phospho-stripped      drop leading 'phospho-', RE-CASCADE remainder
      13. p- stripped           drop leading 'p-', RE-CASCADE remainder
      14. leading-greek strip   drop a leading Greek word (gammaH2AX->H2AX);
                                _STRIP_DENYLIST guards mis-linking remainders
      15. integrin subunit      parse a lone integrin chain -> ITGA*/ITGB*
      16. separator-deletion    delsep equality: HIF1alpha == HIF-1alpha -> HIF1A
    Then two wrappers re-run the WHOLE cascade on a cleaned form of the original:
      17. qualifier-stripped    drop leading state words (nuclear/mutant/wild-
                                type/...) and trailing molecule/assay words
                                (protein/mRNA/gene/antibody/...): HIF-1alpha
                                protein -> HIF-1alpha -> HIF1A
      18. residue-stripped      drop phospho-site / point-mutation notes
                                (GSK-3beta S9A, GSK3beta ( Ser9 ), -Tyr216) -> GSK3B
    Pass 1 keys on the raw original surface; passes 3+ key on the expanded value.
    Passes 9/12/13 recurse one level (depth guard); curated targets are verified
    against the approved-symbol set at load (a warning lists any that are not).

(C) ROUTING. distinct-gene count of the winning pass decides the bucket:
    1 gene  -> greek.json (single, "hgnc_symbol": "SYM", "hgnc_id": ...)
    >=2     -> greek_ambiguous.json ("hgnc_symbol": [..], "hgnc_ids": [..],
              "n_genes": N)
    0       -> unmatched_greek.json.
    Approved-symbol guard: before a >=2-gene hit is filed as ambiguous, if exactly
    ONE of those genes carried the winning key in its approved `symbol` field (the
    rest only via alias/prev/name), the value resolves to that gene -- the string
    IS its official symbol -- and the single entry is flagged
    "approved_symbol_guard": true. Curated-synonym and integrin-parse results are
    authoritative and bypass the guard (so deliberate ambiguous lists are kept).

(D) NON-GENE TRIAGE. High-frequency (occ>=5) unmatched values reviewed and found
    NOT to be a single HGNC gene (oncolytic viruses, PNA reagents, the metabolite
    alpha-ketoglutarate, mouse-only genes, CAR/viral constructs, bare/ambiguous
    tokens) are recorded in NON_GENE with a reason and listed in the report, so
    they read as deliberately unmapped rather than missed.

(E) STAGE C -- cosine fallback (cosine_fallback). Stage B links only by exact
    equality of a transformed key. This fuzzy fallback re-examines the UNMATCHED
    multi-word values and suggests a gene when a value's word multiset
    (near-)equals a multi-word HGNC name-field surface but the word order /
    punctuation differs (word-count cosine >= COSINE_THRESHOLD; 1.0 = a pure
    reordering). It writes greek_cosine.json keyed by the original entity with the
    suggested symbol(s), cosine score and matched surface. It is an ADVISORY
    OVERLAY: the entities stay in unmatched_greek.json and the four-way partition
    is unchanged. NON_GENE-triaged values are skipped so it never contradicts (D).
=============================================================================

full_greek.json schema:

  {
    "<original clean_genetic_ne>": {
        "greek_expanded": "<value with Greek glyphs spelled out>",
        "occurrences": <int>
    },
    ...
  }

greek_complex.json adds "complex", "hgnc_symbol", "hgnc_ids", "n_genes" and
"match_mode". The complex split is skipped (with a notice) if the HGNC file is
absent. greek_complex.json and unmatched_greek.json together partition
full_greek.json.

The glyph -> word table is built from the Unicode database by scanning the
Greek & Coptic (U+0370..U+03FF) and Greek Extended (U+1F00..U+1FFF) blocks: any
character whose Unicode name starts with "GREEK" and contains one of the 24
Greek letter names is mapped to that name. This covers capital and small
letters (Α->Alpha, α->alpha), the final sigma (ς->sigma), glyph variants
(ϕ phi symbol->phi) and accented Greek Extended letters. Case is preserved
(Δ->Delta, δ->delta). Non-letter Greek marks (the bare tonos ΄) and Coptic
letters (ϩ hori, Ϫ gangia) are NOT Greek-alphabet letters, have no Greek-word
form, and are left untouched; the script reports any such residue.
"""

import argparse
import json
import math
import re
import sys
import unicodedata
from collections import Counter, defaultdict
from html import escape as _esc

IN_PATH = "GENETIC/greek_clean_genetic_ne.tsv"
OUT_PATH = "GENETIC/full_greek.json"
HGNC_PATH = "databases/hgnc_complete_set_2026-05-01.json"
COMPLEX_OUT_PATH = "GENETIC/greek_complex.json"
UNMATCHED_OUT_PATH = "GENETIC/unmatched_greek.json"
SINGLE_OUT_PATH = "GENETIC/greek.json"
AMBIGUOUS_OUT_PATH = "GENETIC/greek_ambiguous.json"
COSINE_OUT_PATH = "GENETIC/greek_cosine.json"
HTML_OUT_PATH = "GENETIC/greek.html"

# HGNC fields tested for exact string equality, in priority order (the
# highest-priority field that matches is reported as match_field).
MATCH_FIELDS = ["symbol", "alias_symbol", "prev_symbol",
                "name", "alias_name", "prev_name"]

# Stage C (cosine fallback) compares multi-word values to multi-word HGNC *name*
# surfaces by word-set cosine; only the descriptive name fields carry multi-word
# surfaces worth this fuzzy bag-of-words match.
NAME_FIELDS = ["name", "alias_name", "prev_name"]
COSINE_THRESHOLD = 0.99

# Canonical 24 Greek letter names. Unicode spells lambda "LAMDA" in character
# names; normalise it to the conventional "lambda".
GREEK_NAMES = {
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
    "iota", "kappa", "lamda", "lambda", "mu", "nu", "xi", "omicron", "pi",
    "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega",
}


def build_translation_table():
    """Return {codepoint: word} for every Greek-alphabet glyph, derived from the
    Unicode character names."""
    table = {}
    for cp in list(range(0x0370, 0x0400)) + list(range(0x1F00, 0x2000)):
        try:
            name = unicodedata.name(chr(cp))
        except ValueError:
            continue
        if not name.startswith("GREEK"):
            continue  # skips COPTIC letters and anything non-Greek
        words = name.split()
        letter = next((w for w in words if w.lower() in GREEK_NAMES), None)
        if letter is None:
            continue  # e.g. "GREEK TONOS" -> not a letter
        word = "lambda" if letter.lower() == "lamda" else letter.lower()
        if "CAPITAL" in words:
            word = word.capitalize()
        table[cp] = word
    return table


def is_greek_block(ch):
    cp = ord(ch)
    return 0x0370 <= cp <= 0x03FF or 0x1F00 <= cp <= 0x1FFF


# --------------------------------------------------------------------------
# Multi-protein complex classification + HGNC verification
# --------------------------------------------------------------------------
# Some greek_expanded values name a MULTI-PROTEIN COMPLEX -- a physical assembly
# of two or more DISTINCT gene products acting as one functional unit (e.g.
# NF-kappaB, the NFKB1(p50)/RELA(p65) dimer) -- rather than a single gene. This
# is a curated classifier tuned to the families present in this data. It keeps
# SINGLE SUBUNITS out (IKKbeta=IKBKB, CK2alpha=CSNK2A1, integrin beta1=ITGB1,
# NF-kappaB p65=RELA) and ignores coordinate paralog mentions (GSK3alpha/beta).

GENE_SETS = {
    "NF-kappaB": ["NFKB1", "RELA"],          # canonical p50/p65 dimer
    "gamma-secretase": ["PSEN1", "PSEN2", "NCSTN", "APH1A", "APH1B", "PSENEN"],
    "IKK complex": ["CHUK", "IKBKB", "IKBKG"],
    "MutSbeta": ["MSH2", "MSH3"],
    "lymphotoxin alpha/beta": ["LTA", "LTB"],
    "beta-catenin destruction complex":
        ["CTNNB1", "APC", "AXIN1", "GSK3B", "CSNK1A1"],
    "CK2 holoenzyme": ["CSNK2A1", "CSNK2A2", "CSNK2B"],
    "AMPK complex": ["PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2",
                     "PRKAG1", "PRKAG2", "PRKAG3"],
    "beta-catenin/TCF complex": ["CTNNB1", "TCF7L2"],
    "MutLalpha": ["MLH1", "PMS2"],
}
# Hand-curated synonyms that no string normalization reaches because HGNC simply
# does not carry the spelled-out form as a field (e.g. EPAS1 has alias 'HIF2A'
# but not 'HIF-2alpha'). Each maps to a SINGLE unambiguous approved symbol;
# verified at load time. 'Abeta*' is the amyloid-beta peptide -> its parent gene
# APP. Keys are exact greek_expanded values.
CURATED_SYNONYMS = {
    "HIF-2alpha": "EPAS1", "HIF2alpha": "EPAS1",
    "HIF-3alpha": "HIF3A", "HIF3alpha": "HIF3A",
    "TGF-beta1": "TGFB1", "TGF-beta2": "TGFB2", "TGF-beta3": "TGFB3",
    "TGFbeta1": "TGFB1", "TGFbeta2": "TGFB2", "TGFbeta3": "TGFB3",
    "TGF-alpha": "TGFA", "TGFalpha": "TGFA", "TGF alpha": "TGFA",
    "Tgfalpha": "TGFA",                        # transforming growth factor alpha
    "IFN-beta": "IFNB1", "IFNbeta": "IFNB1",
    "SDF-1alpha": "CXCL12",
    "GSK3beta": "GSK3B", "GSK-3beta": "GSK3B", "GSK3-beta": "GSK3B",
    "GSK3-beta kinase": "GSK3B", "GSK-3beta-specific siRNA": "GSK3B",
    "GSK3alpha": "GSK3A",
    "GSk-3beta": "GSK3B", "GSk3beta": "GSK3B",  # lowercase-k typos
    "GSK3alpha / beta": ["GSK3A", "GSK3B"],   # coordinate paralog mention
    "PDGFRalpha": "PDGFRA", "PDGFRbeta": "PDGFRB", "hPDGFRalpha": "PDGFRA",
    "PDGFR-alpha": "PDGFRA", "PDGFR-beta": "PDGFRB", "Pdgfr-alpha": "PDGFRA",
    "CK2alpha": "CSNK2A1", "CK2beta": "CSNK2B", "CK1delta": "CSNK1D",
    "IDH3alpha": "IDH3A", "IDH3beta": "IDH3B", "IDH3gamma": "IDH3G",
    # CK2 alpha-prime (the 2nd catalytic subunit) -> CSNK2A2; the data uses a
    # prime char (U+2032) with a space.
    "CK2alpha ′": "CSNK2A2", "CK2alpha ′ subunit": "CSNK2A2",
    "CK2alpha ′-siRNA": "CSNK2A2", "CK2alpha'": "CSNK2A2",
    "CHKalpha": "CHKA", "CHKalpha I199N / F200N": "CHKA",   # choline kinase alpha
    "CHKbeta": "CHKB",                         # choline kinase beta
    "IRE1alpha": "ERN1", "LTbetaR": "LTBR", "IKKepsilon": "IKBKE",
    # PGC-1 transcriptional coactivators: alpha = PPARGC1A, beta = PPARGC1B
    "PGC-1alpha": "PPARGC1A", "PGC1alpha": "PPARGC1A", "PGC-1 alpha": "PPARGC1A",
    "Pgc1alpha": "PPARGC1A",                   # case variant
    "PGC-1beta": "PPARGC1B", "PGC1beta": "PPARGC1B", "PGC-1 beta": "PPARGC1B",
    "PGC-1beta-/-": "PPARGC1B",
    "Rho": "RHOA", "rho": "RHOA",             # canonical Rho GTPase (vs RHO/RHOD)
    # PPAR receptors: alpha=PPARA; beta and delta are the same gene PPARD.
    "PPARalpha": "PPARA", "Pparalpha": "PPARA",   # incl. lowercase-tail variant
    "Ppargamma": "PPARG",                      # case variant (PPARgamma in HGNC)
    "PPARdelta": "PPARD", "PPAR-delta": "PPARD", "PPAR delta": "PPARD",
    "PPARbeta": "PPARD", "PPARbeta / delta": "PPARD",
    "p110alpha": "PIK3CA", "p110beta": "PIK3CB",
    "AMPKalpha1": "PRKAA1", "FBW7alpha": "FBXW7",
    "alpha-SMA": "ACTA2", "beta-actin": "ACTB",
    "IL-13Ralpha2": "IL13RA2", "IL13Ralpha2": "IL13RA2", "integrin beta1": "ITGB1",
    "SA-beta-gal": "GLB1", "SA-beta-Gal": "GLB1",
    "Abeta": "APP", "Abeta42": "APP", "Abeta40": "APP", "Abeta1-42": "APP",
    "Abeta42 / 40": "APP", "Abeta 42 / 40": "APP",   # amyloid-beta peptide ratios
    "Abeta38 / 40": "APP", "Abeta42 / 38": "APP",
    # catenins by Greek letter: alpha-catenin is ambiguous across the three
    # alpha-catenin paralogs (CTNNA1/2/3 -> list value); delta-catenin = CTNND2
    # (CTNND1 is "p120 catenin"). A list value routes to greek_ambiguous.json.
    "alpha-catenin": ["CTNNA1", "CTNNA2", "CTNNA3"], "delta-catenin": "CTNND2",
    # alpha-tubulin: ambiguous across the protein-coding 'tubulin alpha N' genes
    # (pseudogenes / antisense / 'alpha-like' excluded).
    "alpha-tubulin": ["TUBA1A", "TUBA1B", "TUBA1C", "TUBA3C", "TUBA3D",
                      "TUBA3E", "TUBA4A", "TUBA4B", "TUBA8"],
    "beta-tubulin": ["TUBB", "TUBB1", "TUBB2A", "TUBB2B", "TUBB3", "TUBB4A",
                     "TUBB4B", "TUBB6", "TUBB8", "TUBB8B"],
    # gamma-tubulin: TUBG1/TUBG2 only (the TUBGCP* are gamma-TuRC complex
    # components, not gamma-tubulins). delta/epsilon-tubulin are single genes.
    "gamma-tubulin": ["TUBG1", "TUBG2"],
    "delta-tubulin": "TUBD1", "epsilon-tubulin": "TUBE1",
    # NF-kappaB named subunits -> their single genes (ambiguous as a bare token,
    # but the NF-kappaB context fixes the gene). p65=RELA, p50/p105=NFKB1,
    # p52/p100=NFKB2. ('nuclear/cytoplasmic ...' variants are reached via the
    # qualifier strip re-matching these keys.)
    "NF-kappaB p65": "RELA", "NFkappaB p65": "RELA",
    "NF-kappaB p65 subunit": "RELA",
    "p65NF-kappaB": "RELA", "NFkappaB ( p65 )": "RELA",
    "p-NF-kappaBp65": "RELA",
    "NF-kappaB p50": "NFKB1", "NFkappaB p50": "NFKB1",
    "NF-kappaB p105": "NFKB1", "NFkappaB p105": "NFKB1",
    "NF-kappaB p52": "NFKB2", "NFkappaB p52": "NFKB2",
    "NF-kappaB p100": "NFKB2", "NFkappaB p100": "NFKB2",
}

# PKC isoform family: each Greek letter -> its gene, in glued/hyphen/space
# spellings (PKCalpha / PKC-alpha / PKC alpha). iota and lambda are the same
# human gene PRKCI; mu = PRKD1 (PKD1), nu = PRKD3 (PKD3).
_PKC_ISO = {"alpha": "PRKCA", "beta": "PRKCB", "gamma": "PRKCG",
            "delta": "PRKCD", "epsilon": "PRKCE", "eta": "PRKCH",
            "theta": "PRKCQ", "zeta": "PRKCZ", "iota": "PRKCI",
            "lambda": "PRKCI", "mu": "PRKD1", "nu": "PRKD3"}
for _iso, _g in _PKC_ISO.items():
    for _f in (f"PKC{_iso}", f"PKC-{_iso}", f"PKC {_iso}"):
        CURATED_SYNONYMS.setdefault(_f, _g)
for _f in ("PKCbetaII", "PKC betaII", "PKC beta II", "PKC-beta II"):
    CURATED_SYNONYMS.setdefault(_f, "PRKCB")   # beta II isoform
CURATED_SYNONYMS.setdefault("aPKCzeta", "PRKCZ")        # atypical PKC zeta
CURATED_SYNONYMS.setdefault("rPKCepsilon", "PRKCE")     # rat PKC epsilon
CURATED_SYNONYMS.setdefault("CA-PKCalpha", "PRKCA")     # constitutively-active
# coordinate multi-isoform mentions -> ambiguous across the named isoforms
CURATED_SYNONYMS.setdefault("PKCalpha / beta", ["PRKCA", "PRKCB"])
CURATED_SYNONYMS.setdefault("PKC alpha and gamma", ["PRKCA", "PRKCG"])
CURATED_SYNONYMS.setdefault("PKC alpha, delta, and epsilon",
                            ["PRKCA", "PRKCD", "PRKCE"])

# PLC isoform family (phospholipase C): numbered isoforms -> single genes; a
# bare Greek letter -> ambiguous across that subfamily. delta has no PLCD2
# (pseudogene); epsilon/zeta are single genes.
_PLC = {"beta1": "PLCB1", "beta2": "PLCB2", "beta3": "PLCB3", "beta4": "PLCB4",
        "gamma1": "PLCG1", "gamma2": "PLCG2",
        "delta1": "PLCD1", "delta3": "PLCD3", "delta4": "PLCD4",
        "eta1": "PLCH1", "eta2": "PLCH2",
        "epsilon": "PLCE1", "epsilon1": "PLCE1", "zeta": "PLCZ1", "zeta1": "PLCZ1",
        "beta": ["PLCB1", "PLCB2", "PLCB3", "PLCB4"],
        "gamma": ["PLCG1", "PLCG2"], "delta": ["PLCD1", "PLCD3", "PLCD4"],
        "eta": ["PLCH1", "PLCH2"]}
for _iso, _g in _PLC.items():
    for _f in (f"PLC{_iso}", f"PLC-{_iso}", f"PLC {_iso}"):
        CURATED_SYNONYMS.setdefault(_f, _g)
del _iso, _g, _f

# G-protein alpha subunit family (Galpha...): subtype -> GNA* gene; a bare 'i'/'t'
# is ambiguous across that subfamily. ('0' = OCR of the 'o' subtype -> GNAO1;
# Galpha15/16 are the same gene GNA15.)
_GA = {"s": "GNAS", "i1": "GNAI1", "i2": "GNAI2", "i3": "GNAI3",
       "o": "GNAO1", "0": "GNAO1", "z": "GNAZ", "q": "GNAQ",
       "11": "GNA11", "12": "GNA12", "13": "GNA13", "14": "GNA14",
       "15": "GNA15", "16": "GNA15", "olf": "GNAL",
       "t1": "GNAT1", "t2": "GNAT2", "t3": "GNAT3",
       "i": ["GNAI1", "GNAI2", "GNAI3"], "t": ["GNAT1", "GNAT2", "GNAT3"]}
for _sub, _g in _GA.items():
    for _f in (f"Galpha{_sub}", f"Galpha-{_sub}", f"Galpha {_sub}",
               f"G-alpha{_sub}"):
        CURATED_SYNONYMS.setdefault(_f, _g)
CURATED_SYNONYMS.setdefault("Galphai1 / 3", ["GNAI1", "GNAI3"])   # coordinate
CURATED_SYNONYMS.setdefault("ITGalphaV", "ITGAV")      # integrin alpha V
CURATED_SYNONYMS.setdefault("shITGalphaV _ 150", "ITGAV")   # shRNA construct

# High-frequency (occ>=10) unmatched forms mapped to their HGNC genes.
_BATCH = {
    # alpha-synuclein
    "alpha-syn": "SNCA", "alphaSyn": "SNCA", "alpha-Syn": "SNCA",
    "A53T alphaS": "SNCA", "alpha-syn ( A53T ) Tg": "SNCA",
    # transcription factors / coactivators / receptors
    "GADD45alpha": "GADD45A", "HNF4gamma": "HNF4G", "HNF4alpha": "HNF4A",
    "ATF6alpha": "ATF6", "RORgammat": "RORC", "mPRalpha": "PAQR7",
    "alpha7nAChR": "CHRNA7",
    # p53 / p63 isoforms
    "Delta133p53alpha": "TP53", "Delta133p53": "TP53", "DeltaNp63": "TP63",
    # antibodies (anti-X) / ligands
    "alpha-CTLA-4": "CTLA4", "alpha-PD-1": "PDCD1", "alpha-PD1": "PDCD1",
    "SDF1alpha": "CXCL12",
    # kinases / subunits
    "JNK2alpha": "MAPK9", "PI3Kalpha": "PIK3CA", "AMPKalpha2": "PRKAA2",
    "SykDeltaMG": "SYK", "Thoc2Delta / Y": "THOC2", "DeltaFVII": "F7",
    # structural / other
    "gamma-taxilin": "TXLNG", "Hsp90beta": "HSP90AB1", "beta-Actin": "ACTB",
    "betaIII-tubulin": "TUBB3", "gammaH2A. X": "H2AX", "AbetaTg": "APP",
    "beta-arrestin 1": "ARRB1", "beta-arrestin 2": "ARRB2", "beta3-AR": "ADRB3",
    # ambiguous: bare IkappaB across the classic IkB family
    "IkappaB": ["NFKBIA", "NFKBIB", "NFKBIE"],
    "beta-arrestin": ["ARRB1", "ARRB2"],
}
for _k, _v in _BATCH.items():
    CURATED_SYNONYMS.setdefault(_k, _v)
del _k, _v

# occ>=5 tier.
_BATCH5 = {
    # p53/p63 isoforms & deletions
    "Delta40p53": "TP53", "p53beta": "TP53", "p53DeltaE5-6": "TP53",
    "mutant Delta133p53alpha R273H": "TP53", "RDeltaA HMGA2": "HMGA2",
    # phosphatases / kinases / their isoforms & constructs
    "PTPmu": "PTPRM", "ChoKalpha": "CHKA", "CK1epsilon": "CSNK1E",
    "PI3Kgamma": "PIK3CG", "p110delta": "PIK3CD", "p110alphanKO": "PIK3CA",
    "p38alphaMCK": "MAPK14", "CaMKIIdelta": "CAMK2D", "PHLPP1alpha": "PHLPP1",
    "Stat3Delta": "STAT3", "beta-STAT3KO": "STAT3", "SykDeltaMG": "SYK",
    "Thoc2Delta": "THOC2",
    # receptors / coactivators
    "IL15Ralpha": "IL15RA", "PDGF-Ralpha": "PDGFRA", "PDGF-Rbeta": "PDGFRB",
    "RIbeta": "PRKAR1B", "LXRalpha": "NR1H3", "sigma-1 receptor": "SIGMAR1",
    "alpha7 nAChR": "CHRNA7", "alpha4-1BB": "TNFRSF9", "mPRalpha": "PAQR7",
    "DeltaIgTrkB": "NTRK2", "DeltaN-Bcl-xL": "BCL2L1",
    # ligands / cytokines / antibodies
    "hIFN-beta": "IFNB1", "IFN-gamma +": "IFNG", "IFNgammahigh": "IFNG",
    "Hrg beta-1": "NRG1", "NRG1beta": "NRG1", "pro-IL-1beta": "IL1B",
    "alpha-PD-L1": "CD274", "alphaPECAM": "PECAM1", "alphaDG": "DAG1",
    # structural / metabolic / other
    "Cx43DeltaCT": "GJA1", "Cln3Deltaex7-8": "CLN3", "Kapbeta2": "TNPO1",
    "NFkappaB1": "NFKB1", "Sen-beta-Gal": "GLB1", "colXIalpha1": "COL11A1",
    "betaIII tubulin": "TUBB3", "beta-III tubulin": "TUBB3",
    "betaIV spectrin": "SPTBN4", "beta-MHC": "MYH7", "beta-cateninWT": "CTNNB1",
    "alphaSMA": "ACTA2", "alphaS": "SNCA", "alphaBC": "CRYAB",
    # amyloid-beta peptide forms -> APP
    "Abetao": "APP", "Abeta 42": "APP", "Abeta11": "APP", "AbetapH": "APP",
    "sAPPalpha": "APP",
    # TGF-beta receptors
    "TGF-betaR1": "TGFBR1", "TGFbetaR2": "TGFBR2", "TGFbeta-1": "TGFB1",
    # ambiguous (coordinate / generic)
    "TGF-betaR": ["TGFBR1", "TGFBR2"], "TNF-alphaR": ["TNFRSF1A", "TNFRSF1B"],
    "alphaCD3": ["CD3D", "CD3E", "CD3G"], "alphaVLA-4": ["ITGA4", "ITGB1"],
    "Galphai1 / 3 DshRNA": ["GNAI1", "GNAI3"],
}
for _k, _v in _BATCH5.items():
    CURATED_SYNONYMS.setdefault(_k, _v)
# alpha-TUB is alpha-tubulin -> same ambiguous set
CURATED_SYNONYMS.setdefault("alpha-TUB", CURATED_SYNONYMS["alpha-tubulin"])
del _k, _v

# Finalized occ>=5 triage: high-frequency values reviewed and intentionally NOT
# mapped because they are not a single HGNC gene. Recorded for auditability and
# surfaced in the report so they read as "reviewed", not "missed".
NON_GENE = {
    "G47Delta": "oncolytic HSV-1 (G47delta)",
    "G47Delta-mIL2": "oncolytic HSV-1 construct",
    "G47Delta-mIL12": "oncolytic HSV-1 construct",
    "Delta-24": "oncolytic adenovirus",
    "Delta-24-RGD": "oncolytic adenovirus",
    "Delta-24-ACT": "oncolytic adenovirus",
    "Delta-24-RGD-H43m": "oncolytic adenovirus",
    "gammaPNA2": "gamma peptide-nucleic-acid reagent",
    "sgammaPNA": "gamma peptide-nucleic-acid reagent",
    "alphaKG": "metabolite (alpha-ketoglutarate)",
    "alpha-KG": "metabolite (alpha-ketoglutarate)",
    "TRIM30alpha": "mouse gene, no human ortholog",
    "alphaLy6G": "mouse Ly6G, no clear human ortholog",
    "beta-elemene": "small-molecule compound",
    "CD28zeta": "CAR signaling construct",
    "13BBzeta": "CAR signaling construct",
    "DeltaB7": "B7 (CD80/CD86) construct, ambiguous",
    "Ad-DeltaB7": "adenoviral B7 construct",
    "Ad-DeltaB7 / 4-1BBL": "adenoviral construct",
    "Ad-DeltaB7 / IL-12 / 4-1BBL": "adenoviral construct",
    "DeltaB7 / 4-1BBL": "construct",
    "Exo-srIkappaB": "IkB super-repressor construct",
    "T-alphaFGL2": "unclear construct",
    "alphaPD": "incomplete anti-PD (PD-1 vs PD-L1 unspecified)",
    "MirgDelta": "unclear",
    "Rho GTPase": "generic family term",
    "beta1": "bare/ambiguous token",
    "beta8-": "bare/ambiguous token",
    "alpha7-": "bare/ambiguous token",
    "beta-": "bare/ambiguous token",
    "-alpha": "bare/ambiguous token",
    "Ealpha": "unclear (MHC E-alpha / mouse)",
    "beta1L": "unclear",
    "InkDelta2 / 3-": "unclear deletion construct",
}
del _sub, _g, _f
CURATED_SYNONYMS.setdefault("IKKalpha-KM", "CHUK")     # kinase-mutant construct
CURATED_SYNONYMS.setdefault("IKKbeta S177A mutant", "IKBKB")   # residue + qualifier
CURATED_SYNONYMS.setdefault("IKKbeta S177 / 181A", "IKBKB")
CURATED_SYNONYMS.setdefault("IKKbeta kinase", "IKBKB")
CURATED_SYNONYMS.setdefault(
    "dominant-negative kinase-dead IKKbeta S177 / 181A mutant", "IKBKB")
# 'integrin alphabeta subunits' names no specific chain -> ambiguous across all
# protein-coding integrin alpha/beta subunit genes.
CURATED_SYNONYMS.setdefault("integrin alphabeta subunits", [
    "ITGA1", "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7",
    "ITGA8", "ITGA9", "ITGA10", "ITGA11", "ITGAD", "ITGAE", "ITGAL", "ITGAM",
    "ITGAV", "ITGAX", "ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6",
    "ITGB7", "ITGB8"])
CURATED_SYNONYMS.setdefault("shPLCbeta1 HA", "PLCB1")   # shRNA + HA-tag construct
CURATED_SYNONYMS.setdefault("opto-PLCbeta", "PLCB1")    # optogenetic PLCb1 construct

_ALPHA_NAMED = {"v": "ITGAV", "iib": "ITGA2B", "l": "ITGAL", "m": "ITGAM",
                "x": "ITGAX", "d": "ITGAD", "e": "ITGAE"}
_FUSED_INTEGRIN = re.compile(r"alpha[a-z0-9]{0,4}beta\s?[0-9]")
_NFKB_SINGLE = {"nfkappab1", "nfkappab2", "nfkappabia", "nfkappabid",
                "nfkappabiz", "nfkappaben", "nfkappabs", "nfkappabib",
                "nfkappabie"}


def _has(s, *subs):
    return any(x in s for x in subs)


def classify_complex(name):
    """Return the complex-family label if `name` (a greek_expanded value)
    denotes a multi-protein complex, else None."""
    s = name.lower()
    if _has(s, "nf-kappab", "nfkappab", "nuclear factor-kappab",
            "nuclear factor kappab", "nf-kappabeta", "nfkappabeta",
            "nf-kappaben"):
        if s in _NFKB_SINGLE:
            return None
        if _has(s, "p65", "p50", "p52", "p105", "p100"):
            return None
        return "NF-kappaB"
    fused = _FUSED_INTEGRIN.search(s.replace(" ", ""))
    if ("integrin" in s and "alpha" in s and "beta" in s) or fused:
        return "integrin heterodimer"
    if _has(s, "gamma-secretase", "gamma secretase"):
        return "gamma-secretase"
    if _has(s, "ck2alpha / beta", "ck2 alpha / beta", "ck2alpha/beta",
            "ck2 alpha/beta"):
        return "CK2 holoenzyme"
    if (_has(s, "ikkalpha / beta", "ikkalpha/beta", "ikkalpha + beta",
             "ikk complex") or ("ikk" in s and "alpha" in s and "beta" in s)):
        return "IKK complex"
    if "ampk" in s and (("alpha" in s and "beta" in s) or "complex" in s):
        return "AMPK complex"
    if _has(s, "ltalpha / beta", "ltalpha/beta", "lt alpha / beta"):
        return "lymphotoxin alpha/beta"
    if "mutsbeta" in s:
        return "MutSbeta"
    if "mutlalpha" in s:
        return "MutLalpha"
    if "destruction complex" in s:
        return "beta-catenin destruction complex"
    if "beta catenin tcf complex" in s:
        return "beta-catenin/TCF complex"
    if (_has(s, "gammadeltatcr", "gammadelta tcr")
            or ("tcr" in s and "alpha" in s and "beta" in s)
            or ("mait" in s and "alphabeta" in s)
            or s == "alphabeta"):
        return "T-cell receptor"
    return None


def integrin_genes(name):
    """Parse the alpha/beta integrin chains named in `name` to HGNC symbols."""
    s = name.lower()
    genes = []
    for m in re.finditer(r"alpha(v|iib|l|m|x|d|e|[0-9]+)", s):
        tok = m.group(1)
        genes.append(_ALPHA_NAMED.get(tok, "ITGA" + tok.upper()))
    for m in re.finditer(r"beta\s?([0-9]+)", s):
        genes.append("ITGB" + m.group(1))
    seen, out = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def tcr_genes(name):
    s = name.lower()
    if "gammadelta" in s or ("gamma" in s and "delta" in s):
        return ["TRG", "TRD"]
    return ["TRA", "TRB"]


def resolve_genes(label, name):
    """Component HGNC symbols for one complex value."""
    if label == "integrin heterodimer":
        return integrin_genes(name)
    if label == "T-cell receptor":
        return tcr_genes(name)
    return list(GENE_SETS.get(label, []))


def load_hgnc_docs(path):
    """Return the list of HGNC gene records."""
    return json.load(open(path, encoding="utf-8"))["response"]["docs"]


# --------------------------------------------------------------------------
# Self-documenting HTML report (generated from the live pipeline results, so it
# never drifts from the JSON outputs).
# --------------------------------------------------------------------------
_HTML_STYLE = """<style>
  body { font: 15px/1.5 -apple-system, Segoe UI, Roboto, sans-serif;
          max-width: 940px; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; }
  h1 { font-size: 1.5rem; }
  h2 { font-size: 1.2rem; margin-top: 2rem; border-bottom: 1px solid #ddd;
        padding-bottom: .3rem; }
  h3 { font-size: 1.02rem; margin-top: 1.5rem; }
  code { background: #f3f3f3; padding: 1px 5px; border-radius: 3px; font-size: .9em; }
  pre { background: #f7f7f7; padding: .7rem 1rem; border-radius: 6px; overflow:auto; }
  table { border-collapse: collapse; margin: .7rem 0; font-size: .9em; }
  th, td { border: 1px solid #ccc; padding: .35rem .7rem; text-align: left; }
  th { background: #f7f7f7; }
  .num { text-align: right; font-variant-numeric: tabular-nums; }
  /* tidy data tables: header band, zebra rows, accented totals */
  table.tidy { border: 1px solid #cfd6dd; }
  table.tidy thead th { background: #eef2f6; border-color: #cfd6dd;
        border-bottom: 2px solid #b9c2cc; }
  table.tidy tbody tr:nth-child(even) td { background: #f6f8fa; }
  table.tidy tbody tr:hover td { background: #eef6ee; }
  table.tidy tr.tot th { background: #eef2f6; border-top: 2px solid #b9c2cc; }
  ol.flow > li { margin: .4rem 0; }
  .headline { background: #eef6ee; border: 1px solid #cde3cd; border-radius: 6px;
              padding: .6rem 1rem; margin: 1rem 0; }
  .headline .big { font-size: 1.6rem; font-weight: 700; }
</style>"""


def _occ_sum(lib):
    return sum(v["occurrences"] for v in lib.values()
               if isinstance(v.get("occurrences"), int))


_MODES = ["key case-insensitive", "curated synonym", "case-sensitive",
          "case-insensitive",
          "hyphen/whitespace", "greek-letter split", "letter-number split",
          "letter-number split (ci)", "prefix-stripped", "IFN expansion",
          "alpha-stripped", "phospho-stripped", "p- stripped",
          "leading-greek strip", "integrin subunit", "separator-deletion",
          "qualifier-stripped", "residue-stripped"]
_MODE_ABBR = {"key case-insensitive": "key-ci", "curated synonym": "curated",
              "case-sensitive": "case-sens.", "case-insensitive": "case-insens.",
              "hyphen/whitespace": "hyphen/ws", "greek-letter split": "gk-split",
              "letter-number split": "ln-split", "letter-number split (ci)": "ln-ci",
              "prefix-stripped": "pfx-strip", "IFN expansion": "ifn-exp",
              "alpha-stripped": "alpha-strip", "phospho-stripped": "phospho-strip",
              "p- stripped": "p-strip", "leading-greek strip": "lead-gk",
              "integrin subunit": "itg-sub",
              "separator-deletion": "sep-del", "qualifier-stripped": "qual-strip",
              "residue-stripped": "res-strip"}

# Full Stage-B cascade, in order, for the strategy section of the report. Each
# entry is (step number, match_mode label, HTML description). Kept in sync with
# match_cascade() and the module docstring's MAPPING STRATEGY section.
STRATEGY_STEPS = [
    ("1", "key case-insensitive",
     "the original <code>full_greek.json</code> key (the pre-expansion surface "
     "form) casefolded is exactly equal to a symbol / alias / prev / name field "
     "&mdash; matched before any glyph expansion or transform "
     "(<code>tau</code>&rarr;MAPT, <code>&beta;-actin</code>&rarr;ACTB). A "
     "deliberate curated override still wins over this for its few keys "
     "(<code>Rho</code>&rarr;RHOA, not the coincidental RHO)"),
    ("2", "curated synonym",
     "hand-curated table for names HGNC carries no field for "
     "(<code>HIF-2alpha</code>&rarr;EPAS1, <code>NF-kappaB p65</code>&rarr;RELA, "
     "and whole families: PKC / PLC / G&alpha; / tubulins / PPAR / &hellip;); a "
     "string value &rarr; single gene, a list &rarr; ambiguous"),
    ("3", "case-sensitive", "exact string equality against the six fields"),
    ("4", "case-insensitive", "casefold equality (on the expanded value)"),
    ("5", "hyphen/whitespace",
     "casefold + <code>-</code>&harr;space (<code>NF-kappaB</code> == "
     "<code>NF kappaB</code>)"),
    ("6", "greek-letter split",
     "split a trailing Greek word, optional hyphen/space "
     "(<code>HIF1alpha</code> / <code>gamma-H2AX</code> &rarr; the bare gene)"),
    ("7", "letter-number split",
     "<code>XXX-&lt;n&gt;&lt;greek&gt;</code> &rarr; "
     "<code>XXX&lt;n&gt; &lt;greek&gt;</code> (<code>HIF-2alpha</code>)"),
    ("8", "letter-number split (ci)",
     "case-insensitive-prefix variant, matched hyphen-folded"),
    ("9", "prefix-stripped",
     "drop a leading <code>si</code>/<code>sh</code>/<code>p</code>/"
     "<code>pro</code> and re-run the cascade on the remainder"),
    ("10", "IFN expansion",
     "<code>IFN-gamma</code> / <code>IFNgamma</code> &rarr; "
     "<code>interferon gamma</code>"),
    ("11", "alpha-stripped",
     "leading <code>alpha</code> = anti-X (<code>alphaPD-1</code>&rarr;PD-1); a "
     "denylist skips intrinsic-alpha forms (<code>alphaSMA</code>&hellip;)"),
    ("12", "phospho-stripped", "drop a leading <code>phospho-</code>, re-cascade"),
    ("13", "p- stripped", "drop a leading <code>p-</code>, re-cascade"),
    ("14", "leading-greek strip",
     "drop a leading Greek word (<code>gammaH2AX</code>&rarr;H2AX); denylist "
     "guards mis-linking remainders"),
    ("15", "integrin subunit",
     "parse a lone integrin chain &rarr; ITGA*/ITGB* "
     "(<code>integrin beta6</code>&rarr;ITGB6)"),
    ("16", "separator-deletion",
     "casefold + drop all hyphens/spaces (<code>HIF1alpha</code> == "
     "<code>HIF-1alpha</code>&rarr;HIF1A)"),
    ("17", "qualifier-stripped",
     "drop leading state words (nuclear/mutant/wild-type) and trailing "
     "molecule/assay words (protein/mRNA/gene/antibody), then re-cascade "
     "(<code>HIF-1alpha protein</code>&rarr;HIF1A)"),
    ("18", "residue-stripped",
     "drop phospho-site / point-mutation notes, then re-cascade "
     "(<code>GSK-3beta S9A</code>, <code>GSK3beta ( Ser9 )</code>&rarr;GSK3B)"),
]


def _mode_split(lib):
    out = {m: 0 for m in _MODES}
    for v in lib.values():
        out[v["match_mode"]] = out.get(v["match_mode"], 0) + 1
    return out


def _by_occ(items):
    return sorted(items, key=lambda kv: -(kv[1]["occurrences"]
                  if isinstance(kv[1]["occurrences"], int) else 0))


def _snippet(lib, mode):
    """A real one-entry JSON snippet for the first (highest-occurrence) entry of
    a given match_mode, for the schema examples."""
    cands = _by_occ([(k, v) for k, v in lib.items()
                     if v.get("match_mode") == mode])
    if not cands:
        return None
    k, v = cands[0]
    return _esc(json.dumps({k: v}, ensure_ascii=False, indent=2))


def _snippet_for(lib, greek_value):
    """One-entry JSON snippet for the entry whose greek_expanded equals
    `greek_value`, for a specific worked example; None if absent."""
    for k, v in lib.items():
        if v.get("greek_expanded") == greek_value:
            return _esc(json.dumps({k: v}, ensure_ascii=False, indent=2))
    return None


def _ex_list(items):
    parts = []
    for k, v in items:
        sym = v["hgnc_symbol"]
        sym = "/".join(sym) if isinstance(sym, list) else sym
        parts.append(f"<code>{_esc(k)}</code>&rarr;{_esc(sym)}")
    return ", ".join(parts)


def render_html(library, n_glyphs, expanded, leftovers,
                complex_lib, single_lib, ambiguous_lib, unmatched_lib,
                cosine_lib=None):
    """Build the full summary HTML from the computed pipeline results."""
    cosine_lib = cosine_lib or {}
    n = len(library)
    P = ['<!doctype html>', '<html lang="en">', '<head>', '<meta charset="utf-8">',
         '<title>greek &mdash; Greek-glyph expansion &amp; HGNC linkage</title>',
         _HTML_STYLE, '</head>', '<body>']

    # --- expansion ---
    P.append("<h1>greek &mdash; Greek-glyph expansion &amp; HGNC linkage</h1>")
    P.append("<p>Generated by <code>greek.py</code>. Expands single-character "
             "Greek glyphs in <code>GENETIC/greek_clean_genetic_ne.tsv</code> "
             "to spelled-out words, then links the entities to HGNC.</p>")
    P.append(f'<div class="headline"><div><span class="big">{n:,}</span> '
             f"entities &middot; Greek glyphs expanded in <strong>{expanded:,}"
             f"</strong> &middot; {n_glyphs} Greek code points mapped.</div></div>")

    P.append("<h2>full_greek.json</h2>")
    P.append("<p>Every entity, keyed by the original <code>clean_genetic_ne</code> "
             "value, with its Greek-expanded form and occurrence count:</p>")
    ex = _by_occ(list(library.items()))[:4]
    P.append("<pre><code>" + _esc(json.dumps(dict(ex), ensure_ascii=False,
             indent=2)) + "</code></pre>")

    if leftovers:
        P.append("<h2>Deliberately left untouched</h2>")
        P.append("<p>Glyphs in the Greek/Coptic blocks that are not "
                 "Greek-alphabet letters (no Greek-word form), reported and "
                 "passed through unchanged:</p>")
        P.append('<table class="tidy"><thead><tr><th>glyph</th><th>code point</th>'
                 '<th class="num">count</th><th>name</th></tr></thead><tbody>')
        for ch, cnt in sorted(leftovers.items(), key=lambda kv: -kv[1]):
            nm = unicodedata.name(ch, "?")
            P.append(f'<tr><td>&#x{ord(ch):X};</td><td>U+{ord(ch):04X}</td>'
                     f'<td class="num">{cnt}</td><td>{_esc(nm.title())}</td></tr>')
        P.append("</tbody></table>")

    # --- full mapping strategy ---
    P.append("<h2>Mapping strategy (greek_expanded &rarr; HGNC)</h2>")
    P.append("<p>Each Greek-expanded value is mapped to HGNC gene symbol(s) from "
             "<code>databases/hgnc_complete_set_2026-05-01.json</code> and routed to one of "
             "four libraries.</p>")
    P.append("<ol class='flow'>")
    P.append("<li><strong>Index HGNC</strong> as <code>value &rarr; {(symbol, "
             "hgnc_id)}</code> over six fields (<code>symbol</code>, "
             "<code>alias_symbol</code>, <code>prev_symbol</code>, <code>name</code>, "
             "<code>alias_name</code>, <code>prev_name</code>; array fields per "
             "element) under four key normalizations: <em>literal</em>, "
             "<em>casefold</em>, <em>hyphen&harr;space fold</em>, and "
             "<em>separator-deletion</em>.</li>")
    P.append("<li><strong>Stage&nbsp;A &mdash; complexes.</strong> A value naming a "
             "multi-protein complex (&ge;2 distinct genes, e.g. NF-kappaB = "
             "NFKB1/RELA) is mapped to its verified subunit genes &rarr; "
             "<code>greek_complex.json</code>. Single subunits and coordinate "
             "paralog mentions are excluded.</li>")
    P.append("<li><strong>Stage&nbsp;B &mdash; single/ambiguous cascade</strong> "
             "(table below): the non-complex remainder is matched by an ordered "
             "cascade whose <strong>first step matches the original "
             "<code>full_greek.json</code> key</strong> (the pre-expansion surface "
             "form) casefolded against the six HGNC fields; later steps key on the "
             "glyph-expanded value. Genes are unioned across the six fields and the "
             "first pass that hits wins, recorded as <code>match_mode</code> "
             "(highest-priority matching field recorded as <code>match_field</code>).</li>")
    P.append("<li><strong>Route</strong> by the winning pass's distinct-gene count: "
             "1 &rarr; <code>greek.json</code> (single); &ge;2 &rarr; "
             "<code>greek_ambiguous.json</code> (list + <code>hgnc_ids</code> + "
             "<code>n_genes</code>); 0 &rarr; <code>unmatched_greek.json</code>. "
             "<strong>Approved-symbol guard:</strong> a &ge;2-gene hit in which "
             "exactly one gene carried the key in its approved <code>symbol</code> "
             "field (the rest only via alias/prev/name) resolves to that gene "
             "&mdash; the string is its official symbol &mdash; rather than going to "
             "ambiguous, and is flagged <code>approved_symbol_guard</code> in "
             "<code>greek.json</code>. Curated and integrin-parse hits bypass it.</li>")
    P.append("</ol>")
    P.append("<h3>Stage B cascade (in order)</h3>")
    P.append('<table class="tidy"><thead><tr><th class="num">#</th>'
             '<th>pass (match_mode)</th><th>what it does</th></tr></thead><tbody>')
    for num, label, desc in STRATEGY_STEPS:
        P.append(f'<tr><td class="num">{num}</td><td><code>{_esc(label)}</code>'
                 f'</td><td>{desc}</td></tr>')
    P.append("</tbody></table>")
    P.append("<p>Pass&nbsp;1 keys on the raw original surface; passes&nbsp;3+ key "
             "on the glyph-expanded value. Passes 9/12/13 re-cascade their stripped "
             "remainder (one level, depth-guarded); passes 17&ndash;18 re-cascade a "
             "cleaned form of the "
             "original value. Every curated and complex-subunit symbol is verified "
             "against the HGNC approved-symbol set at load (a warning lists any that "
             "are not). High-frequency (occ&ge;5) unmatched values that are not a "
             "single gene &mdash; oncolytic viruses, peptide-nucleic-acid reagents, "
             "the metabolite &alpha;-ketoglutarate, mouse-only genes, CAR/viral "
             "constructs &mdash; are recorded in <code>NON_GENE</code> and listed "
             "under &ldquo;Reviewed non-gene&rdquo; below.</p>")

    # --- Stage A: complexes ---
    fam = {}
    for v in complex_lib.values():
        e = fam.setdefault(v["complex"], {"forms": 0, "occ": 0, "genes": set()})
        e["forms"] += 1
        if isinstance(v["occurrences"], int):
            e["occ"] += v["occurrences"]
        e["genes"].update(v["hgnc_symbol"])
    cgenes = sorted({g for v in complex_lib.values() for g in v["hgnc_symbol"]})
    P.append("<h2>Stage A &mdash; multi-protein complexes "
             "(<code>greek_complex.json</code>)</h2>")
    P.append("<p>Values whose expanded form names a multi-protein complex "
             "(an assembly of distinct gene products, e.g. NF-kappaB = "
             "NFKB1/RELA) are mapped to their verified HGNC subunit genes. Single "
             "subunits (IKKbeta=IKBKB) and coordinate paralog mentions "
             "(GSK3alpha/beta) are excluded.</p>")
    P.append(f'<div class="headline"><div><span class="big">{len(complex_lib)}'
             f"</span> values name complexes ({_occ_sum(complex_lib):,} occ), "
             f"{len(fam)} families, {len(cgenes)} distinct HGNC subunit genes."
             "</div></div>")
    P.append('<table class="tidy"><thead><tr><th>complex</th>'
             '<th class="num">forms</th><th class="num">occ.</th>'
             '<th>component genes (HGNC)</th></tr></thead><tbody>')
    for label, e in sorted(fam.items(), key=lambda kv: -kv[1]["occ"]):
        genes = sorted(e["genes"])
        shown = "/".join(genes[:6]) + ("&hellip;" if len(genes) > 6 else "")
        P.append(f'<tr><td>{_esc(label)}</td><td class="num">{e["forms"]}</td>'
                 f'<td class="num">{e["occ"]:,}</td><td>{shown}</td></tr>')
    P.append(f'<tr class="tot"><th>total</th><th class="num">{len(complex_lib)}</th>'
             f'<th class="num">{_occ_sum(complex_lib):,}</th><th></th></tr>'
             '</tbody></table>')

    # --- Stage B: HGNC string match ---
    rest = len(single_lib) + len(ambiguous_lib) + len(unmatched_lib)
    s_mode, a_mode = _mode_split(single_lib), _mode_split(ambiguous_lib)
    P.append("<h2>Stage B &mdash; HGNC string match</h2>")
    P.append("<p>The non-complex remainder is matched by the 18-pass cascade "
             "detailed in <em>Mapping strategy &rarr; Stage B cascade</em> above. "
             "Each pass runs only on what the previous left unmatched; a value's "
             "distinct-gene count routes it to <code>greek.json</code> (single) or "
             "<code>greek_ambiguous.json</code> (multiple). The yields per pass:</p>")
    n_guarded = sum(1 for v in single_lib.values()
                    if v.get("approved_symbol_guard"))
    guard_note = (f" Of the single-gene hits, <strong>{n_guarded}</strong> were "
                  f"resolved by the approved-symbol guard (one approved-"
                  f"<code>symbol</code> match among &ge;2 alias-level candidates)."
                  if n_guarded else "")
    P.append(f'<div class="headline"><div><span class="big">{len(single_lib)}'
             f"</span> values map to a single gene &rarr; <code>greek.json</code>; "
             f"<strong>{len(ambiguous_lib)}</strong> to multiple &rarr; "
             f"<code>greek_ambiguous.json</code>; <strong>{len(unmatched_lib)}"
             f"</strong> remain unmatched (of {rest:,} non-complex values)."
             f"{guard_note}</div></div>")

    P.append("<h3>Output schema (real entries)</h3>")
    snips = [s for s in (_snippet(single_lib, "curated synonym"),
                         _snippet_for(single_lib, "NF-kappaB p65"),
                         _snippet(single_lib, "qualifier-stripped"),
                         _snippet(ambiguous_lib, "case-insensitive")) if s]
    P.append("<pre><code>" + "\n".join(snips) + "</code></pre>")

    # Bucket summary -- a compact 4-column table (no per-mode columns).
    P.append("<h3>Buckets</h3>")
    P.append('<table class="tidy"><thead><tr><th>bucket</th><th>file</th>'
             '<th class="num">values</th><th class="num">occ.</th></tr></thead><tbody>')
    for label, fname, lib in (
            ("single gene", "greek.json", single_lib),
            ("multiple genes", "greek_ambiguous.json", ambiguous_lib),
            ("unmatched", "unmatched_greek.json", unmatched_lib)):
        P.append(f'<tr><td>{label}</td><td><code>{fname}</code></td>'
                 f'<td class="num">{len(lib):,}</td>'
                 f'<td class="num">{_occ_sum(lib):,}</td></tr>')
    P.append('</tbody></table>')

    # Per-pass yield -- one row per matcher, in cascade order (tall, not wide).
    occ_by_mode = {m: 0 for m in _MODES}
    for lib in (single_lib, ambiguous_lib):
        for v in lib.values():
            if isinstance(v.get("occurrences"), int):
                occ_by_mode[v["match_mode"]] += v["occurrences"]
    P.append("<h3>Yield by matching pass</h3>")
    P.append('<table class="tidy"><thead><tr><th class="num">#</th><th>pass</th>'
             '<th class="num">single</th><th class="num">ambig.</th>'
             '<th class="num">values</th><th class="num">occ.</th></tr></thead><tbody>')
    for i, m in enumerate(_MODES, 1):
        sv, av = s_mode[m], a_mode[m]
        P.append(f'<tr><td class="num">{i}</td><td>{_esc(m)}</td>'
                 f'<td class="num">{sv or ""}</td><td class="num">{av or ""}</td>'
                 f'<td class="num">{sv + av}</td>'
                 f'<td class="num">{occ_by_mode[m]:,}</td></tr>')
    tot_s, tot_a = len(single_lib), len(ambiguous_lib)
    P.append(f'<tr class="tot"><td></td><th>total matched</th>'
             f'<th class="num">{tot_s}</th><th class="num">{tot_a}</th>'
             f'<th class="num">{tot_s + tot_a}</th>'
             f'<th class="num">{_occ_sum(single_lib) + _occ_sum(ambiguous_lib):,}</th>'
             '</tr></tbody></table>')

    # NF-kappaB subunit sub-tally: curated subunit mappings (+ their
    # qualifier-stripped variants) that the per-pass table folds into the
    # 'curated synonym' / 'qualifier-stripped' rows.
    _NFKB_GENES = {"RELA", "NFKB1", "NFKB2", "RELB", "REL"}
    nfkb = [(k, v) for k, v in single_lib.items()
            if "kappab" in v["greek_expanded"].lower()
            and v["hgnc_symbol"] in _NFKB_GENES]
    if nfkb:
        nocc = sum(v["occurrences"] for _, v in nfkb
                   if isinstance(v["occurrences"], int))
        P.append("<h3>NF-kappaB subunit sub-tally</h3>")
        P.append("<p>Named NF-kappaB subunits resolved to single genes (a subset "
                 "of the <code>curated synonym</code> and "
                 "<code>qualifier-stripped</code> rows above; the NF-kappaB "
                 "<em>complex</em> itself is handled separately in Stage&nbsp;A):</p>")
        P.append('<table class="tidy"><thead><tr><th>value</th>'
                 '<th>&rarr; gene</th><th>pass</th>'
                 '<th class="num">occ.</th></tr></thead><tbody>')
        for k, v in sorted(nfkb, key=lambda kv: -(kv[1]["occurrences"]
                           if isinstance(kv[1]["occurrences"], int) else 0)):
            P.append(f'<tr><td><code>{_esc(k)}</code></td>'
                     f'<td>{v["hgnc_symbol"]}</td><td>{_esc(v["match_mode"])}</td>'
                     f'<td class="num">{v["occurrences"]}</td></tr>')
        P.append(f'<tr class="tot"><th>{len(nfkb)} values</th><th></th><th></th>'
                 f'<th class="num">{nocc:,}</th></tr></tbody></table>')

    top_single = _ex_list(_by_occ(list(single_lib.items()))[:6])
    if top_single:
        P.append(f"<p><strong>Top single matches:</strong> {top_single}.</p>")
    if ambiguous_lib:
        P.append("<p><strong>Ambiguous (multi-gene) matches:</strong> "
                 + _ex_list(_by_occ(list(ambiguous_lib.items()))) + ".</p>")

    # --- reviewed non-gene (finalized occ>=5 triage) ---
    ng = [(v["greek_expanded"], v["occurrences"], NON_GENE[v["greek_expanded"]])
          for v in unmatched_lib.values() if v["greek_expanded"] in NON_GENE]
    if ng:
        ng.sort(key=lambda r: -r[1] if isinstance(r[1], int) else 0)
        ng_occ = sum(o for _, o, _ in ng if isinstance(o, int))
        P.append("<h2>Reviewed non-gene (occ&ge;5, intentionally unmapped)</h2>")
        P.append("<p>High-frequency unmatched values that were triaged and are "
                 "<strong>not a single HGNC gene</strong> &mdash; recorded so they "
                 "read as reviewed, not missed:</p>")
        P.append('<table class="tidy"><thead><tr><th>value</th>'
                 '<th class="num">occ.</th><th>why not a gene</th></tr>'
                 '</thead><tbody>')
        for ge, o, reason in ng:
            P.append(f'<tr><td><code>{_esc(ge)}</code></td>'
                     f'<td class="num">{o}</td><td>{_esc(reason)}</td></tr>')
        P.append(f'<tr class="tot"><th>{len(ng)} values</th>'
                 f'<th class="num">{ng_occ:,}</th><th></th></tr></tbody></table>')

    # --- Stage C: cosine fallback overlay ---
    if cosine_lib:
        cos_occ = sum(v["occurrences"] for v in cosine_lib.values()
                      if isinstance(v["occurrences"], int))
        P.append("<h2>Stage C &mdash; word-order cosine fallback "
                 "(<code>greek_cosine.json</code>)</h2>")
        P.append("<p>Stage B links by exact equality of a transformed key. This "
                 "fuzzy fallback re-examines the <em>unmatched</em> multi-word "
                 "values and suggests a gene when a value's <strong>word "
                 "multiset</strong> (near-)equals a multi-word HGNC <em>name</em>-"
                 "field surface but the word <strong>order/punctuation</strong> "
                 f"differs (cosine &ge; {COSINE_THRESHOLD}; 1.0 = a pure "
                 "reordering). It is an <strong>advisory overlay</strong>: these "
                 "entities <em>remain</em> in <code>unmatched_greek.json</code> and "
                 "the four-way partition is unchanged &mdash; the overlay only "
                 "records the suggestion and its score for review. Values triaged "
                 "as <code>NON_GENE</code> are skipped.</p>")
        P.append(f'<div class="headline"><div><span class="big">{len(cosine_lib)}'
                 f"</span> unmatched values get a cosine suggestion "
                 f"({cos_occ:,} occ).</div></div>")
        P.append('<table class="tidy"><thead><tr><th>value (greek_expanded)</th>'
                 '<th class="num">occ.</th><th class="num">cosine</th>'
                 '<th>HGNC symbol(s)</th><th>matched name surface</th>'
                 '</tr></thead><tbody>')
        for original, v in _by_occ(list(cosine_lib.items())):
            sym = v["hgnc_symbol"]
            sym = sym if isinstance(sym, str) else "/".join(sym)
            surf = " ; ".join(v["hgnc_value"][:3]) + (
                "&hellip;" if len(v["hgnc_value"]) > 3 else "")
            occ = v["occurrences"]
            occ = f"{occ:,}" if isinstance(occ, int) else _esc(str(occ))
            P.append(f'<tr><td><code>{_esc(v["greek_expanded"])}</code></td>'
                     f'<td class="num">{occ}</td>'
                     f'<td class="num">{v["cosine"]}</td>'
                     f'<td>{_esc(sym)}</td><td>{_esc(surf)}</td></tr>')
        P.append('</tbody></table>')

    # --- partition + reproduce ---
    total = (len(complex_lib) + len(single_lib) + len(ambiguous_lib)
             + len(unmatched_lib))
    P.append("<h2>Partition</h2>")
    P.append("<p>The four HGNC-stage libraries (<code>greek_complex.json</code>, "
             "<code>greek.json</code>, <code>greek_ambiguous.json</code>, "
             "<code>unmatched_greek.json</code>) partition the "
             f"{n:,} entities of <code>full_greek.json</code> exactly: "
             f"{len(complex_lib)} + {len(single_lib)} + {len(ambiguous_lib)} + "
             f"{len(unmatched_lib)} = {total:,}. "
             "(<code>greek_cosine.json</code> is an advisory overlay on the "
             "unmatched set, not a partition member.)</p>")
    P.append("<h2>Reproduce</h2><pre><code>python greek.py</code></pre>")
    P.append("</body></html>")
    return "\n".join(P) + "\n"


def symbol_index(docs):
    """Return {approved_symbol: hgnc_id}."""
    return {d["symbol"]: d.get("hgnc_id") for d in docs if "symbol" in d}


def hsfold(s):
    """Casefold and fold hyphen<->whitespace: a length-preserving 1:1 fold that
    makes 'NF-kappaB' and 'NF kappaB' compare equal (casefold, then '-' -> ' ')."""
    return s.casefold().replace("-", " ")


def delsep(s):
    """Casefold and DELETE separators (hyphens and spaces): a destructive fold so
    'HIF1alpha' and 'HIF-1alpha' compare equal (vs hsfold, which keeps a space)."""
    return s.casefold().replace("-", "").replace(" ", "")


# 24 spelled-out Greek letter words (lambda normalised), longest first.
_GREEK_WORDS = sorted(
    ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
     "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho",
     "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega"],
    key=len, reverse=True)


def split_trailing_greek(value):
    """If `value` ends with a spelled-out Greek word glued to upstream text,
    insert a space before it ('HIF1alpha' -> 'HIF1 alpha'), else return None."""
    low = value.lower()
    for w in _GREEK_WORDS:
        if (low.endswith(w) and len(value) > len(w)
                and value[-len(w) - 1].isalnum()):
            return value[:-len(w)] + " " + value[-len(w):]
    return None


# Upper-case prefix, hyphen, digits, glued Greek word: e.g. HIF-2alpha, IL-1beta.
_LN_GREEK_RE = re.compile(r"^([A-Z]+)-(\d+)(" + "|".join(_GREEK_WORDS) + r")$")
# Case-insensitive prefix variant: e.g. Il-1beta, il-1beta, pGSK-3beta.
_LN_GREEK_CI_RE = re.compile(r"^([a-z]+)-(\d+)(" + "|".join(_GREEK_WORDS) + r")$",
                             re.IGNORECASE)


def regroup_letter_number_greek(value, ci=False):
    """Convert 'XXX-<digits><greek>' to 'XXX<digits> <greek>' -- drop the hyphen
    and space the Greek word ('HIF-2alpha' -> 'HIF2 alpha'), else return None.
    With ci=True the prefix may be any case ('Il-1beta' -> 'Il1 beta')."""
    m = (_LN_GREEK_CI_RE if ci else _LN_GREEK_RE).match(value)
    return f"{m.group(1)}{m.group(2)} {m.group(3)}" if m else None


# Reagent / modification prefixes that sit in front of a gene token, longest
# first: si/sh (knock-down), p (phospho-), pro (pro-form).
_REAGENT_PREFIXES = ["pro", "sh", "si", "p"]


def strip_reagent_prefix(value):
    """Strip a leading reagent/modifier prefix when it precedes an upper-case
    gene-like token ('shHIF-1alpha' -> 'HIF-1alpha', 'pHP1gamma' -> 'HP1gamma');
    return the stripped string, or None if no prefix applies."""
    for p in _REAGENT_PREFIXES:
        if (value.startswith(p) and len(value) > len(p)
                and value[len(p)].isupper()):
            return value[len(p):]
    return None


# 'IFN-<greek>' or 'IFN<greek>' -> the spelled-out interferon name.
_IFN_RE = re.compile(r"^IFN-?(" + "|".join(_GREEK_WORDS) + r")$")


def expand_ifn(value):
    """Expand the IFN abbreviation: 'IFN-gamma'/'IFNgamma' -> 'interferon gamma',
    else return None."""
    m = _IFN_RE.match(value)
    return f"interferon {m.group(1)}" if m else None


# Leading 'alpha' immediately followed by a Latin letter -- the alpha/anti-
# antibody notation, e.g. alphaPD-1, alphaCD40.
_ALPHA_AB_RE = re.compile(r"^alpha([A-Za-z].*)$")
# Post-strip remainders that mis-link (shared by the leading-Greek strips):
# alphaSMA = a-smooth-muscle-actin (ACTA2, not SMN1), alphaGC/Delta... = a
# glycolipid not gene GC, B7 = a CD80/CD86 family label not LRRC23's alias,
# (beta-)arrestin = ARRB1/2 not the visual-arrestin SAG. 'gal'/'ar' are intrinsic-
# Greek abbreviations: beta-gal = beta-galactosidase (GLB1 / lacZ reporter), NOT
# galanin GAL; beta-AR / alpha-AR = adreno(b/a)ceptors (ADRB*/ADRA*), NOT androgen
# receptor AR -- both coincidentally equal an approved symbol once the Greek
# letter is stripped, so they are blocked here.
_STRIP_DENYLIST = {"sma", "gc", "b7", "cat", "syn", "arrestin",
                   "arrestin 1", "arrestin 2", "gal", "ar"}


def strip_leading_alpha(value):
    """Drop a leading 'alpha' that is glued to a Latin letter ('alphaPD-1' ->
    'PD-1'); return None if no such prefix or the remainder is denylisted."""
    m = _ALPHA_AB_RE.match(value)
    if not m or m.group(1).casefold() in _STRIP_DENYLIST:
        return None
    return m.group(1)


# Leading spelled-out Greek word (other than alpha, which has its own pass)
# glued to a Latin letter -- a modifier prefix, e.g. gammaH2AX, DeltaEGFR.
_LEAD_GREEK = [w for w in _GREEK_WORDS if w != "alpha"]
# Greek word, an OPTIONAL hyphen/space, then a Latin token: gammaH2AX, gamma-H2AX.
_LEAD_GREEK_RE = re.compile(r"^(" + "|".join(_LEAD_GREEK) + r")[- ]?([A-Za-z].*)$",
                            re.IGNORECASE)


def strip_leading_greek(value):
    """Drop a leading Greek word (not alpha) and an optional hyphen/space
    ('gammaH2AX'/'gamma-H2AX' -> 'H2AX'); return None if none, the remainder is
    shorter than 2 chars, or it is denylisted."""
    m = _LEAD_GREEK_RE.match(value)
    if not m:
        return None
    rem = m.group(2)
    if len(rem) < 2 or rem.casefold() in _STRIP_DENYLIST:
        return None
    return rem


# Descriptive qualifiers stripped before re-matching the core entity: leading
# state/localization words and trailing molecule-type/assay words. Words that
# occur inside real gene names (factor, transcription, receptor, kinase,
# subunit, complex, family, domain, binding) are deliberately EXCLUDED.
_QUAL_LEAD = {"nuclear", "cytoplasmic", "cytosolic", "membrane", "mutant",
              "wild-type", "wildtype", "wt", "serum", "total", "active",
              "cleaved", "soluble", "recombinant", "endogenous", "secreted",
              "mature"}
_QUAL_TRAIL = {"protein", "proteins", "mrna", "gene", "genes", "antibody",
               "antibodies", "mab", "isoform", "isoforms", "promoter",
               "reporter", "construct", "constructs", "staining", "expression",
               "level", "levels", "mutant", "mutants", "deletion", "fragment",
               "sirna", "shrna", "cdna", "peptide", "target", "targets"}


def strip_qualifiers(value):
    """Remove leading state/localization qualifiers and trailing
    molecule-type/assay qualifiers ('HIF-1alpha protein' -> 'HIF-1alpha',
    'nuclear beta-catenin' -> 'beta-catenin'); return the cleaned string, or
    None if nothing was stripped or nothing remains."""
    toks = value.split()
    changed = False
    while toks and toks[0].lower() in _QUAL_LEAD:
        toks.pop(0)
        changed = True
    while toks and toks[-1].lower() in _QUAL_TRAIL:
        toks.pop()
        changed = True
    return " ".join(toks) if changed and toks else None


def strip_phospho_word(value):
    """Drop a leading 'phospho-' marker, any case ('phospho-tau' -> 'tau'),
    else return None."""
    return value[8:] if value.lower().startswith("phospho-") and len(value) > 8 \
        else None


def strip_p_dash(value):
    """Drop a leading 'p-' phospho marker ('p-tau' -> 'tau'), else return None."""
    return value[2:] if value.startswith("p-") and len(value) > 2 else None


def build_field_index(docs, keyfn=None):
    """Return {field: {key: {(symbol, hgnc_id), ...}}} for the six match fields,
    so a string can be looked up to the gene(s) it occurs in. `keyfn` transforms
    each field value into its lookup key (identity if None; e.g. str.casefold for
    case-insensitive, hsfold for hyphen/whitespace-insensitive lookup)."""
    idx = {f: {} for f in MATCH_FIELDS}
    for d in docs:
        gene = (d.get("symbol"), d.get("hgnc_id"))
        for f in MATCH_FIELDS:
            if f not in d:
                continue
            v = d[f]
            for s in (v if isinstance(v, list) else [v]):
                if isinstance(s, str):
                    idx[f].setdefault(keyfn(s) if keyfn else s,
                                      set()).add(gene)
    return idx


def match_hgnc(value, idx):
    """Exact (case-sensitive) match `value` against every match field.
    Returns (genes, fields): the union of (symbol, hgnc_id) genes it equals and
    the priority-ordered list of fields in which it matched."""
    genes, fields = set(), []
    for f in MATCH_FIELDS:
        hits = idx[f].get(value)
        if hits:
            genes |= hits
            fields.append(f)
    return genes, fields


# Phospho-site / point-mutation annotations to strip off a residue mutant so the
# bare gene remains (e.g. 'GSK-3beta S9A' / 'GSK3beta ( Ser9 )' -> 'GSK3beta').
_AA3 = ("ala|arg|asn|asp|cys|gln|glu|gly|his|ile|leu|lys|met|phe|pro|"
        "ser|thr|trp|tyr|val")
_PMUT_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]?$")
_RESNUM_RE = re.compile(r"(?i)^(" + _AA3 + r")\d+[a-z]?$")
_RESWORD_RE = re.compile(r"(?i)^(serine|threonine|tyrosine)$")
# trailing site suffix: -Ser9 (residue name), -S9A (AA+pos+AA), or -Y216 / -F720
# (AA + >=2 digits). NOT -R1 (one digit, no trailing AA -> usually "receptor 1").
_RES_SUFFIX_RE = re.compile(
    r"(?i)-((" + _AA3 + r")\d+[a-z]?"
    r"|[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]"
    r"|[ACDEFGHIKLMNPQRSTVWY]\d{2,})$")


def strip_residue_annotations(value):
    """Remove phospho-site / point-mutation annotations so the bare gene remains:
    parenthetical site notes ('( Ser9 )'), point mutations ('S9A', 'Y216E'),
    residue names/words ('Ser9', 'serine'), and a trailing residue suffix on a
    token ('GSK3beta-Ser9' -> 'GSK3beta'). Bare numbers are kept (they are often
    isoform numbers, e.g. 'TGFbeta 2'). Returns the cleaned string, or None if
    nothing changed / nothing remains."""
    v = re.sub(r"\([^)]*\)", " ", value)          # drop "( site )" groups
    out = []
    for t in v.split():
        t = _RES_SUFFIX_RE.sub("", t)             # GSK3beta-Ser9 -> GSK3beta
        if not t:
            continue
        if _PMUT_RE.match(t) or _RESNUM_RE.match(t) or _RESWORD_RE.match(t):
            continue
        out.append(t)
    cleaned = " ".join(out).strip()
    return cleaned if cleaned and cleaned != value else None


def _toks(s):
    """Lower-case word tokens for the cosine bag-of-words."""
    return s.casefold().split()


def _cosine(a, b):
    """Cosine of two word-count vectors (Counters); 0 if they share no token."""
    dot = sum(cnt * b.get(t, 0) for t, cnt in a.items())
    if not dot:
        return 0.0
    na = math.sqrt(sum(v * v for v in a.values()))
    nb = math.sqrt(sum(v * v for v in b.values()))
    return dot / (na * nb)


def cosine_fallback(unmatched_lib, docs, threshold=COSINE_THRESHOLD):
    """Fuzzy word-order cosine fallback over the UNMATCHED multi-word values.

    Stage B links by exact equality of a transformed key; this Stage C catches a
    multi-word greek_expanded value whose *word multiset* (near-)equals a
    multi-word HGNC name-field surface but whose word ORDER / punctuation differs
    (cosine >= threshold; 1.0 = a pure reordering, e.g. 'receptor, type II
    activin' vs 'activin receptor type 2'). Bag-of-words is order-blind, so it
    bridges exactly the reshuffles the exact keys miss.

    Returns an ADVISORY overlay library keyed by the original entity. It does NOT
    change the complex/single/ambiguous/unmatched partition: the entries remain in
    unmatched_greek.json, and this overlay merely records the suggested gene(s) and
    the cosine score for review. Values already triaged as NON_GENE are skipped so
    the overlay never contradicts the non-gene record."""
    # multi-word HGNC name surfaces -> word-count vector (the only fields whose
    # surfaces are multi-word phrases worth a bag-of-words comparison).
    candidates = []                                   # (surface, vec, field, sym, hid)
    for d in docs:
        sym = d.get("symbol", "")
        hid = d.get("hgnc_id")
        for field in NAME_FIELDS:
            val = d.get(field)
            if val is None:
                continue
            for s in (val if isinstance(val, list) else [val]):
                if isinstance(s, str) and len(_toks(s)) >= 2:
                    candidates.append((s, Counter(_toks(s)), field, sym, hid))
    # token -> document frequency + postings, so each query only compares against
    # candidates that share one of its two RAREST tokens (keeps it tractable).
    df = Counter()
    postings = defaultdict(list)
    for i, (_, cv, _, _, _) in enumerate(candidates):
        for t in cv:
            df[t] += 1
            postings[t].append(i)

    lib = {}
    for original, entry in unmatched_lib.items():
        value = entry["greek_expanded"]
        if value in NON_GENE:                         # don't contradict the triage
            continue
        itk = _toks(value)
        if len(itk) < 2:
            continue
        iv = Counter(itk)
        cand_ids = set()
        for t in sorted(iv, key=lambda x: df.get(x, 0))[:2]:
            cand_ids.update(postings.get(t, ()))
        best = {}                                     # (sym, hid) -> [cos, fields, surfaces]
        for ci in cand_ids:
            s, cv, field, sym, hid = candidates[ci]
            if s.casefold() == value.casefold():
                continue
            c = _cosine(iv, cv)
            if c >= threshold:
                rec = best.setdefault((sym, hid), [0.0, set(), set()])
                rec[0] = max(rec[0], c)
                rec[1].add(field)
                rec[2].add(s)
        if not best:
            continue
        genes = sorted(best, key=lambda g: g[0])      # by symbol
        syms = [g[0] for g in genes]
        ids = [g[1] for g in genes]
        allfields = set().union(*(best[g][1] for g in genes))
        surfaces = sorted(set().union(*(best[g][2] for g in genes)))
        lib[original] = {
            "greek_expanded": value,
            "occurrences": entry["occurrences"],
            "match_field": next(f for f in NAME_FIELDS if f in allfields),
            "match_mode": "word-order cosine",
            "cosine": round(max(best[g][0] for g in genes), 4),
            "hgnc_symbol": syms[0] if len(syms) == 1 else syms,
            "hgnc_id": ids[0] if len(ids) == 1 else ids,
            "n_genes": len(genes),
            "hgnc_value": surfaces,
        }
    return lib


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("in_path", nargs="?", default=IN_PATH,
                    help=f"input TSV (default: {IN_PATH})")
    ap.add_argument("out_path", nargs="?", default=OUT_PATH,
                    help=f"output JSON (default: {OUT_PATH})")
    ap.add_argument("--hgnc", default=HGNC_PATH,
                    help=f"HGNC complete-set JSON (default: {HGNC_PATH})")
    ap.add_argument("--complex-out", default=COMPLEX_OUT_PATH,
                    help="multi-protein-complex -> HGNC mapping JSON "
                         f"(default: {COMPLEX_OUT_PATH})")
    ap.add_argument("--single-out", default=SINGLE_OUT_PATH,
                    help="single-gene exact-match mapping JSON "
                         f"(default: {SINGLE_OUT_PATH})")
    ap.add_argument("--ambiguous-out", default=AMBIGUOUS_OUT_PATH,
                    help="multi-gene exact-match mapping JSON "
                         f"(default: {AMBIGUOUS_OUT_PATH})")
    ap.add_argument("--unmatched-out", default=UNMATCHED_OUT_PATH,
                    help="entries matched by neither stage "
                         f"(default: {UNMATCHED_OUT_PATH})")
    ap.add_argument("--cosine-out", default=COSINE_OUT_PATH,
                    help="Stage C word-order cosine fallback overlay JSON "
                         f"(default: {COSINE_OUT_PATH})")
    ap.add_argument("--html-out", default=HTML_OUT_PATH,
                    help=f"summary HTML report (default: {HTML_OUT_PATH})")
    args = ap.parse_args()

    table = build_translation_table()

    library = {}
    expanded_count = 0
    leftovers = {}
    dupes = 0
    with open(args.in_path, encoding="utf-8") as fin:
        header = next(fin, "")  # skip the column header
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            if not cols or not cols[0]:
                continue
            original = cols[0]
            occ = cols[1] if len(cols) > 1 else ""
            expanded = original.translate(table)
            if expanded != original:
                expanded_count += 1
            for ch in expanded:
                if ord(ch) not in table and is_greek_block(ch):
                    leftovers[ch] = leftovers.get(ch, 0) + 1
            if original in library:
                dupes += 1
            library[original] = {
                "greek_expanded": expanded,
                "occurrences": int(occ) if occ.isdigit() else occ,
            }

    with open(args.out_path, "w", encoding="utf-8") as fout:
        json.dump(library, fout, ensure_ascii=False, indent=2)
        fout.write("\n")

    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    print(f"Read {len(library)} unique clean_genetic_ne value(s) "
          f"from {args.in_path}")
    print(f"Greek glyphs expanded in {expanded_count} value(s)")
    print(f"Glyphs mapped: {len(table)} Greek code points -> words")
    if dupes:
        print(f"Note: {dupes} duplicate key(s) collapsed")
    print(f"Wrote {args.out_path}")
    if leftovers:
        print(f"\nLeft untouched (not Greek-alphabet letters, "
              f"{len(leftovers)} distinct glyph(s)):")
        for ch in sorted(leftovers, key=lambda c: -leftovers[c]):
            print(f"  {ch}  U+{ord(ch):04X}  x{leftovers[ch]}  "
                  f"{unicodedata.name(ch, '<unnamed>')}")

    # ---- HGNC-dependent stages -------------------------------------------
    try:
        docs = load_hgnc_docs(args.hgnc)
    except FileNotFoundError:
        print(f"\n[skip] HGNC file not found ({args.hgnc}); "
              f"not writing {args.complex_out} / {args.single_out} / "
              f"{args.ambiguous_out} / {args.unmatched_out}")
        return
    sym2id = symbol_index(docs)
    bad_curated = sorted({s for v in CURATED_SYNONYMS.values()
                          for s in ([v] if isinstance(v, str) else v)
                          if s not in sym2id})
    if bad_curated:
        print(f"[warn] curated synonym targets not approved HGNC symbols: "
              f"{', '.join(bad_curated)}")

    # Stage A -- multi-protein complexes (verified vs HGNC) vs the rest.
    complex_lib, rest_lib = {}, {}
    labels, bad, unresolved = {}, {}, []
    for original, entry in library.items():
        label = classify_complex(entry["greek_expanded"])
        if label is None:
            rest_lib[original] = entry
            continue
        verified, ids = [], []
        for g in resolve_genes(label, entry["greek_expanded"]):
            if g in sym2id:
                verified.append(g)
                ids.append(sym2id[g])
            else:
                bad.setdefault(g, set()).add(original)
        if not verified:                      # complex named but no gene resolved
            unresolved.append((original, label))
            rest_lib[original] = entry
            continue
        order = sorted(range(len(verified)), key=lambda i: verified[i])
        verified = [verified[i] for i in order]
        ids = [ids[i] for i in order]
        complex_lib[original] = {
            "greek_expanded": entry["greek_expanded"],
            "occurrences": entry["occurrences"],
            "complex": label,
            "hgnc_symbol": verified,
            "hgnc_ids": ids,
            "n_genes": len(verified),
            "match_mode": "complex→subunits",
        }
        labels[label] = labels.get(label, 0) + 1

    # Stage B -- string-equality match of the remainder against HGNC fields, as
    # a cascade of progressively looser keys: exact case-sensitive, then
    # case-insensitive, then hyphen/whitespace fold. Each match is split
    # single-gene vs multi-gene; match_mode records the pass that hit.
    idx = build_field_index(docs)                       # case-sensitive keys
    idx_ci = build_field_index(docs, keyfn=str.casefold)  # casefolded keys
    idx_hs = build_field_index(docs, keyfn=hsfold)      # casefold + hyphen/space
    idx_del = build_field_index(docs, keyfn=delsep)     # casefold + drop separators
    single_lib, ambiguous_lib, unmatched_lib = {}, {}, {}
    mode_counts = {"key case-insensitive": 0, "curated synonym": 0,
                   "case-sensitive": 0, "case-insensitive": 0,
                   "hyphen/whitespace": 0, "greek-letter split": 0,
                   "letter-number split": 0, "letter-number split (ci)": 0,
                   "prefix-stripped": 0, "IFN expansion": 0,
                   "alpha-stripped": 0, "phospho-stripped": 0,
                   "p- stripped": 0, "leading-greek strip": 0,
                   "integrin subunit": 0,
                   "separator-deletion": 0, "qualifier-stripped": 0,
                   "residue-stripped": 0}

    def match_cascade(value, depth=0, raw_key=None):
        """Run the full Stage-B matcher cascade; return (genes, fields, mode,
        anchor) with genes empty if nothing hit. `anchor` is the subset of `genes`
        that carried the winning key in its approved `symbol` field (vs only an
        alias/prev/name) -- the approved-symbol guard at routing uses it to resolve
        an otherwise-ambiguous tie. The phospho/prefix strips re-run the whole
        cascade on the remainder (so 'p-GSK3beta' -> 'GSK3beta' -> curated GSK3B);
        `depth` guards that one level of stripping (no runaway recursion).

        `raw_key` is the original full_greek.json key (the pre-expansion surface
        form); when given (top level only) it is tried FIRST -- a casefold exact
        match of the raw key against the six HGNC fields -- so a surface that
        already equals an HGNC field links before any glyph expansion / transform."""
        def hit(key, ix, mode):
            """Look `key` up in index `ix`; on a hit return (genes, fields, mode,
            anchor), else None. `anchor` = the genes whose approved `symbol` field
            carried this key (idx[*]['symbol'] is the symbol-field sub-index)."""
            genes, fields = match_hgnc(key, ix)
            if not genes:
                return None
            return genes, fields, mode, ix["symbol"].get(key, set()) & genes

        # STEP 1: the raw key, casefolded, exactly equals an HGNC field (before any
        # expansion / transform). A deliberate CURATED_SYNONYMS override still wins
        # over this for its handful of keys -- e.g. 'Rho' is curated to RHOA (the
        # GTPase), not the coincidental raw-key hit RHO (rhodopsin).
        if depth == 0 and raw_key is not None and value not in CURATED_SYNONYMS:
            r = hit(raw_key.casefold(), idx_ci, "key case-insensitive")
            if r:
                return r
        cur = CURATED_SYNONYMS.get(value)         # hand-curated synonyms win
        if cur:                                   # str -> one gene; list -> ambiguous
            syms = [cur] if isinstance(cur, str) else cur
            genes = {(s, sym2id[s]) for s in syms if s in sym2id}
            if genes:                             # authoritative: no guard (empty anchor)
                return genes, ["curated"], "curated synonym", set()
        r = hit(value, idx, "case-sensitive")
        if r:
            return r
        r = hit(value.casefold(), idx_ci, "case-insensitive")
        if r:
            return r
        r = hit(hsfold(value), idx_hs, "hyphen/whitespace")
        if r:
            return r
        sep = split_trailing_greek(value)         # split trailing Greek letter
        if sep:
            for keyfn, ix in ((None, idx), (str.casefold, idx_ci), (hsfold, idx_hs)):
                r = hit(keyfn(sep) if keyfn else sep, ix, "greek-letter split")
                if r:
                    return r
        conv = regroup_letter_number_greek(value)  # XXX-<n><greek> -> XXX<n> <greek>
        if conv:
            for keyfn, ix in ((None, idx), (str.casefold, idx_ci), (hsfold, idx_hs)):
                r = hit(keyfn(conv) if keyfn else conv, ix, "letter-number split")
                if r:
                    return r
        conv = regroup_letter_number_greek(value, ci=True)   # ci prefix, hyphen-fold
        if conv:
            r = hit(hsfold(conv), idx_hs, "letter-number split (ci)")
            if r:
                return r
        if depth == 0:                            # si/sh/p/pro -> re-cascade remainder
            st = strip_reagent_prefix(value)
            if st:
                g, f, _m, a = match_cascade(st, depth + 1)
                if g:
                    return g, f, "prefix-stripped", a
        conv = expand_ifn(value)                  # IFN<greek> -> interferon <greek>
        if conv:
            r = hit(conv.casefold(), idx_ci, "IFN expansion")
            if r:
                return r
        rem = strip_leading_alpha(value)          # drop leading 'alpha' (anti-X)
        if rem:
            r = hit(rem.casefold(), idx_ci, "alpha-stripped")
            if r:
                return r
        if depth == 0:                            # 'phospho-' -> re-cascade remainder
            rem = strip_phospho_word(value)
            if rem:
                g, f, _m, a = match_cascade(rem, depth + 1)
                if g:
                    return g, f, "phospho-stripped", a
        if depth == 0:                            # 'p-' -> re-cascade remainder
            rem = strip_p_dash(value)
            if rem:
                g, f, _m, a = match_cascade(rem, depth + 1)
                if g:
                    return g, f, "p- stripped", a
        rem = strip_leading_greek(value)          # gammaH2AX -> H2AX
        if rem:
            r = hit(rem.casefold(), idx_ci, "leading-greek strip")
            if r:
                return r
        if "integrin" in value.lower():           # parse integrin chain(s) -> ITGA*/ITGB*
            isyms = [s for s in integrin_genes(value) if s in sym2id]
            if isyms:                             # explicit parse: authoritative, no guard
                return ({(s, sym2id[s]) for s in isyms}, ["integrin"],
                        "integrin subunit", set())
        r = hit(delsep(value), idx_del, "separator-deletion")  # HIF1alpha==HIF-1alpha
        if r:
            return r
        return set(), [], None, set()

    guard_count = 0
    for original, entry in rest_lib.items():
        value = entry["greek_expanded"]
        genes, fields, mode, anchor = match_cascade(value, raw_key=original)
        if not genes:                             # strip descriptive qualifiers, retry
            cleaned = strip_qualifiers(value)
            if cleaned and cleaned != value:
                g, f, _m, a = match_cascade(cleaned)
                if g:
                    genes, fields, mode, anchor = g, f, "qualifier-stripped", a
        if not genes:                             # strip residue/site annotations, retry
            cleaned = strip_residue_annotations(value)
            if cleaned and cleaned != value:
                g, f, _m, a = match_cascade(cleaned)
                if g:
                    genes, fields, mode, anchor = g, f, "residue-stripped", a
        if not genes:
            unmatched_lib[original] = entry
            continue
        # collapse to distinct genes, sorted by symbol
        by_id = {hid: sym for sym, hid in genes}
        # approved-symbol guard: a key that hits >=2 genes is normally held aside as
        # ambiguous, but if exactly ONE of them carried it in its approved `symbol`
        # field (the rest only via alias/prev/name) the string IS that gene's
        # official symbol -- a strong single-gene signal -- so resolve to it.
        guarded = False
        if len(by_id) > 1:
            anchor_ids = {hid for _, hid in anchor}
            if len(anchor_ids) == 1:
                keep = next(iter(anchor_ids))
                genes = {(sym, hid) for sym, hid in genes if hid == keep}
                by_id = {hid: sym for sym, hid in genes}
                guarded = True
                guard_count += 1
        syms = sorted(by_id.values())
        ids = [hid for hid, _ in sorted(by_id.items(), key=lambda kv: kv[1])]
        base = {
            "greek_expanded": value,
            "occurrences": entry["occurrences"],
            "match_field": fields[0],          # highest-priority field hit
            "match_mode": mode,
        }
        if guarded:                            # flag the disambiguation in the JSON
            base["approved_symbol_guard"] = True
        if len(by_id) == 1:
            single_lib[original] = {**base, "hgnc_symbol": syms[0],
                                    "hgnc_id": ids[0]}
        else:
            ambiguous_lib[original] = {**base, "hgnc_symbol": syms,
                                       "hgnc_ids": ids, "n_genes": len(syms)}
        mode_counts[mode] += 1

    # Stage C -- fuzzy word-order cosine fallback over the unmatched leftovers.
    # An advisory OVERLAY: it suggests gene(s) for multi-word unmatched values but
    # leaves the four-way partition untouched (entries stay in unmatched_greek.json).
    cosine_lib = cosine_fallback(unmatched_lib, docs)

    for path, lib in ((args.complex_out, complex_lib),
                      (args.single_out, single_lib),
                      (args.ambiguous_out, ambiguous_lib),
                      (args.unmatched_out, unmatched_lib),
                      (args.cosine_out, cosine_lib)):
        with open(path, "w", encoding="utf-8") as f:
            json.dump(lib, f, ensure_ascii=False, indent=2)
            f.write("\n")

    occ = sum(v["occurrences"] for v in complex_lib.values()
              if isinstance(v["occurrences"], int))
    cgenes = sorted({g for v in complex_lib.values() for g in v["hgnc_symbol"]})
    print(f"\nStage A -- multi-protein complexes: {len(complex_lib)} value(s) "
          f"({occ} occ), {len(cgenes)} subunit genes, {len(labels)} families "
          f"-> {args.complex_out}")
    if unresolved:
        print(f"  complex named but no gene resolved (-> Stage B): "
              f"{', '.join(n for n, _ in unresolved)}")

    print(f"\nStage B -- HGNC string match of the remaining {len(rest_lib)} "
          f"value(s) [key case-insensitive -> curated synonym -> case-sensitive "
          f"-> case-insensitive -> hyphen/whitespace "
          f"-> greek-letter split -> letter-number split -> letter-number "
          f"split (ci) -> prefix-stripped -> IFN expansion -> alpha-stripped "
          f"-> phospho-stripped -> p- stripped -> leading-greek strip "
          f"-> integrin subunit -> separator-deletion; then qualifier-strip / "
          f"residue-strip + re-cascade]:")
    print(f"  single gene  -> {args.single_out}: {len(single_lib)}")
    print(f"  multi  genes -> {args.ambiguous_out}: {len(ambiguous_lib)}")
    print(f"  no match     -> {args.unmatched_out}: {len(unmatched_lib)}")
    print(f"  (of which {guard_count} resolved by the approved-symbol guard "
          f"that would otherwise be ambiguous)")
    print(f"  by mode: key case-insensitive {mode_counts['key case-insensitive']}, "
          f"curated synonym {mode_counts['curated synonym']}, "
          f"case-sensitive {mode_counts['case-sensitive']}, "
          f"case-insensitive {mode_counts['case-insensitive']}, "
          f"hyphen/whitespace {mode_counts['hyphen/whitespace']}, "
          f"greek-letter split {mode_counts['greek-letter split']}, "
          f"letter-number split {mode_counts['letter-number split']}, "
          f"letter-number split (ci) {mode_counts['letter-number split (ci)']}, "
          f"prefix-stripped {mode_counts['prefix-stripped']}, "
          f"IFN expansion {mode_counts['IFN expansion']}, "
          f"alpha-stripped {mode_counts['alpha-stripped']}, "
          f"phospho-stripped {mode_counts['phospho-stripped']}, "
          f"p- stripped {mode_counts['p- stripped']}, "
          f"leading-greek strip {mode_counts['leading-greek strip']}, "
          f"integrin subunit {mode_counts['integrin subunit']}, "
          f"separator-deletion {mode_counts['separator-deletion']}, "
          f"qualifier-stripped {mode_counts['qualifier-stripped']}, "
          f"residue-stripped {mode_counts['residue-stripped']}")
    if bad:
        print(f"  (complex stage) unverified symbols: {', '.join(sorted(bad))}")

    cos_occ = sum(v["occurrences"] for v in cosine_lib.values()
                  if isinstance(v["occurrences"], int))
    print(f"\nStage C -- word-order cosine fallback (advisory overlay, >= "
          f"{COSINE_THRESHOLD}): {len(cosine_lib)} of the {len(unmatched_lib)} "
          f"unmatched value(s) get a suggested gene ({cos_occ} occ) "
          f"-> {args.cosine_out} (partition unchanged)")

    # ---- self-documenting HTML report ------------------------------------
    html = render_html(library, len(table), expanded_count, leftovers,
                       complex_lib, single_lib, ambiguous_lib, unmatched_lib,
                       cosine_lib)
    with open(args.html_out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\nWrote {args.html_out}")


if __name__ == "__main__":
    main()
