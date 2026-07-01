# -*- coding: utf-8 -*-
"""roman.py

End-to-end normalization + HGNC linkage pipeline, unifying the strategies of
four earlier scripts into one:

  STAGE 1  (roman.py)             dash normalization  + Greek/Roman split
  COMPLEX  (priority, runs first) curated complex/heterodimer -> subunit gene(s);
                                  + bare receptor-family name -> family gene set;
                                  + generic NF-kB -> NFKB1 + RELA (p50/p65 dimer),
                                  matching the greek pipeline's NF-kappaB complex
  STAGE 2  (symbols_match.py)     symbol-field linkage (TYPE 1 + TYPE 2)
  STAGE 3  (names_match.py)       descriptive-name-field linkage
  STAGE 3b                        PHOSPHO: strip leading p-/phospho-, re-match 6 fields
  STAGE 3c                        ANTI: strip leading 'anti-' / trailing ' antibody', re-match 6 fields
  STAGE 3c2                       CLEAVED: strip leading 'cleaved-', re-match 6 fields
  STAGE 3d                        PROTEIN: strip ' protein'/'-protein', re-match 6 fields
  STAGE 3e                        TRANSCRIPT: strip ' transcript(s)'/' mRNA(s)', re-match
  STAGE 3f                        VARIANT: strip ' variant(s)'/' mutant(s)'/' mutation(s)', re-match
  STAGE 3g                        MISSENSE: <approved symbol>+[AA]<pos>[AA] -> symbol
  STAGE 3g2                       MISSENSE-SP: <approved symbol> ' ' [AA]<pos>[AA] -> symbol
  STAGE 3h                        INHIBITOR: strip terminal '-i', match symbol fields
  STAGE 3i                        MIRNA: mature miR name -> HGNC precursor gene(s)
  STAGE 3j                        PROMOTER: strip ' promoter' descriptor, match symbol fields
  STAGE 3k                        DEHYPHEN: <UPPER>-<digits> hyphen removed, match symbol fields
  STAGE 3l                        P-PREFIX: strip bare leading 'p', match symbol fields
  STAGE 3m                        WT: strip leading 'WT '/'wild-type ', match symbol fields
  STAGE 3m2                       NUCLEAR: strip leading 'nuclear ', match symbol fields
  STAGE 3m3                       TARGET-GENES: strip ' target genes' -> regulator (tagged separable)
  STAGE 3n                        HISTONE: marks/mutants/variants -> histone gene(s)
  STAGE 3o                        C-ONCOGENE: c-Kit/c-Src/c-Met/c-erbB-N/... spelling -> gene
  STAGE 3p                        GREEK-LETTER: lone Roman letter -> spelled-out Greek word
  STAGE 3q                        DELSEP: all hyphens+spaces removed (PD-L1==PDL1==PD L1), symbol fields
  STAGE 4  (cosine_similarity.py) multi-word word-order cosine match

The unifying idea: STAGES 2 and 3 are the *same* algorithm -- exact equality of
a transformed key, then keep only entities that resolve to a single gene -- run
over different HGNC fields with a looser or tighter key. They share one
`match_pass()`. STAGE 1 produces the input; STAGE 4 is a similarity fallback for
what exact matching leaves behind.

Inputs
------
  sentences/*.json                       BioBERT NER output (GENETIC entities)
  databases/hgnc_complete_set_2026-05-01.json   HGNC complete set (response.docs[])

Outputs (all under GENETIC/)
-------
  clean_genetic_ne.tsv  non_singletons.tsv
  greek_clean_genetic_ne.tsv  roman_clean_genetic_ne.tsv
  roman.json               (all matched entities -> HGNC symbol + match-type)
  roman_ambiguous.json     (all ambiguous entities -> candidate symbols)
  roman_cosine.json   (cosine library; keeps the similarity value)
  unmatched.tsv            (final unmatched values)
  roman.html       (unified strategy + results)

No per-step matched/unmatched TSVs are produced -- the cascade's intermediate
buckets are passed in memory and surface only through roman.json /
roman_ambiguous.json; stale copies from older runs are deleted on start.

NON-GENE triage: prominent values in the final unmatched set that were reviewed
and found NOT to be a single human HGNC gene (reporters/tags, gene-editing & assay
reagents, drugs, engineered variants, circular RNAs, non-protein antigens, a
mouse-only symbol, tokenization fragments, and bare gene-family/pathway/complex
names) are recorded in the NON_GENE dict with a reason and listed in the HTML
report, so prominent leftovers read as deliberately unmapped rather than missed.
It is documentation only and never changes matching; real genes HGNC merely lacks
a surface for (APNG=MPG, Presenilin1=PSEN1, ...) are deliberately excluded as
recall gaps, not non-genes.

Run from anywhere (paths resolve relative to this file)::

    python roman.py
"""

import collections
import glob
import html
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent
SENT_DIR = ROOT / "sentences"
SUMMARY_DIR = ROOT / "GENETIC"
HGNC_JSON = ROOT / "databases" / "hgnc_complete_set_2026-05-01.json"
HTML_OUT = SUMMARY_DIR / "roman.html"

TSV_HEADER = "clean_genetic_ne\toccurrences\tn_source_forms"
LABEL = "GENETIC"


# =======================================================================
# STAGE 1 -- dash normalization + Greek/Roman split   (was roman.py)
# =======================================================================
DASH_VARIANTS = "-‐‑‒–—―−⁃"                 # U+002D 2010 2011 2012 2013 2014 2015 2212 2043
_TO_ASCII = str.maketrans({c: "-" for c in DASH_VARIANTS})
_DESPACE = re.compile(r"\s*-\s*")

GREEK_SYMBOL_RE = re.compile(r"[Ͱ-Ͽἀ-῿]")    # Greek & Coptic + Greek Extended
GREEK_NAMES = [
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
    "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho",
    "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega",
]
GREEK_WORD_RE = re.compile(r"\b(" + "|".join(GREEK_NAMES) + r")\b", re.IGNORECASE)
AMBIGUOUS_SHORT = {"mu", "nu", "xi", "pi", "psi", "eta", "phi", "chi"}


def dash_normalize(text):
    """Unify dash variants to ASCII '-', then strip whitespace around '-'."""
    return _DESPACE.sub("-", (text or "").translate(_TO_ASCII))


def is_greek(name):
    """True if the entity carries a Greek letter as a symbol or a spelled word."""
    if GREEK_SYMBOL_RE.search(name):
        return True
    for m in GREEK_WORD_RE.finditer(name):
        tok = m.group(0)
        if tok.isupper():
            continue
        if tok.lower() in AMBIGUOUS_SHORT and tok != tok.lower():
            continue
        return True
    return False


def stage1():
    """Parse sentences, normalize dashes, split Greek vs Roman. Returns (roman_rows, info)."""
    counts = collections.Counter()
    n_files = 0
    for fpath in sorted(glob.glob(str(SENT_DIR / "*.json"))):
        try:
            rec = json.load(open(fpath, encoding="utf-8"))
        except Exception:
            continue
        n_files += 1
        for sent in rec.get("sentences", []):
            for ent in sent.get("entities", []):
                if ent.get("label") == LABEL:
                    counts[ent.get("text", "")] += 1

    agg = defaultdict(lambda: {"occ": 0, "forms": set()})
    for form, n in counts.items():
        a = agg[dash_normalize(form)]
        a["occ"] += n
        a["forms"].add(form)

    rows = sorted(agg.items(), key=lambda kv: (-kv[1]["occ"], kv[0]))

    def line(key, a):
        return f"{key}\t{a['occ']}\t{len(a['forms'])}"

    def write_rows(path, items):
        body = ("\n".join(line(k, a) for k, a in items) + "\n") if items else ""
        path.write_text(TSV_HEADER + "\n" + body, encoding="utf-8")

    write_rows(SUMMARY_DIR / "clean_genetic_ne.tsv", rows)
    write_rows(SUMMARY_DIR / "non_singletons.tsv",
               [(k, a) for k, a in rows if a["occ"] > 1])

    greek, roman = [], []
    for k, a in rows:
        (greek if is_greek(k) else roman).append(line(k, a))
    (SUMMARY_DIR / "greek_clean_genetic_ne.tsv").write_text(
        "\n".join([TSV_HEADER] + greek) + "\n", encoding="utf-8")
    (SUMMARY_DIR / "roman_clean_genetic_ne.tsv").write_text(
        "\n".join([TSV_HEADER] + roman) + "\n", encoding="utf-8")

    info = {
        "n_files": n_files, "occurrences": sum(counts.values()),
        "unique_before": len(counts), "unique_after": len(agg),
        "greek": len(greek), "roman": len(roman),
    }
    return roman, info


# =======================================================================
# STAGES 2 & 3 -- the shared exact-key, single-gene matcher
#                 (was symbols_match.py + names_match.py)
# =======================================================================
def identity(s):
    return s


def hsfold(s):
    """Case-fold and treat '-' and ' ' as interchangeable.

    Rewriting every '-' to ' ' is a length-preserving 1:1 substitution, so a
    shared key means the strings line up character-for-character, differing only
    at hyphen<->whitespace positions (on top of case). This folds case AND
    hyphen/space surface noise into one key.
    """
    return s.casefold().replace("-", " ")


# a leading phospho marker: 'p-'/'P-' or 'phospho-'/'Phospho-'/'phospho '
_PHOSPHO_RE = re.compile(r"^(?:[Pp]hospho[- ]|[Pp]-)")


def pstrip_hsfold(s):
    """Strip a leading phospho prefix, then hsfold.

    Phospho-protein mentions prepend a phospho marker to a gene name -- 'p-' or
    'phospho-' (e.g. p-AKT, phospho-STAT3) -- which is surface noise for
    identity, so removing it before folding lets them reach AKT1 / STAT3.
    Used for the QUERY only; the HGNC index is built with plain hsfold so that an
    HGNC surface like 'P-protein' is NOT stripped (which would collide a bare
    'protein' query onto it).
    """
    return hsfold(_PHOSPHO_RE.sub("", s, count=1))


# a leading 'anti-' antibody/antagonist marker: 'anti-'/'Anti-'/'anti '
_ANTI_RE = re.compile(r"^[Aa]nti[- ]")
# a trailing ' antibody' reagent descriptor (case-insensitive)
_ANTIBODY_SUF = re.compile(r"(?i) antibody$")


def anti_hsfold(s):
    """Strip a leading 'anti-' marker, and a trailing ' antibody' descriptor too
    -- but the antibody suffix only when the 'anti-' prefix is present, so both
    come off together for the double-decorated 'anti-X antibody' form. Then
    hsfold (anti-PD-1 -> pd 1, anti-EGFR antibody -> egfr).

    'anti-X' / 'anti-X antibody' name a reagent *against* gene X rather than X's
    own product, so this links it to X's gene while the match_mode flags it as
    'anti- stripped' (a target, not the entity itself). Query-only.
    """
    if _ANTI_RE.match(s):
        s = _ANTIBODY_SUF.sub("", _ANTI_RE.sub("", s, count=1))
    return hsfold(s)


# a leading 'cleaved-' / 'Cleaved-' proteolytic-processing marker
_CLEAVED_RE = re.compile(r"^[Cc]leaved-")


def cleaved_strip_hsfold(s):
    """Strip a leading 'cleaved-'/'Cleaved-' marker, then hsfold (cleaved-PARP ->
    PARP, cleaved-caspase 3 -> caspase 3). The cleaved product is still the gene's
    protein, so it folds onto the gene; matched against all six symbol+name
    fields (so 'caspase 3' resolves via the name fields). Query-only.
    """
    return hsfold(_CLEAVED_RE.sub("", s, count=1))


# a ' protein' / '-protein' descriptor suffix (e.g. p53 protein, MGMT-protein)
_PROTEIN_RE = re.compile(r"( protein|-protein)")


def protein_strip_hsfold(s):
    """Strip ' protein'/'-protein' descriptors, then hsfold (e.g. p53 protein -> p53).

    'X protein' denotes X's own gene product, so it folds onto the gene.
    Query-only; the HGNC index keeps its surfaces intact so a genuine name like
    'X binding protein' is not mangled in the index.
    """
    return hsfold(_PROTEIN_RE.sub("", s))


# a ' transcript(s)' / ' mRNA(s)' descriptor (e.g. EGFR mRNA, MYC transcripts,
# EGFR mRNA transcript). The trailing \b keeps 'transcriptase'/'transcriptional'
# intact. Case-sensitive: only the 'mRNA' standard casing occurs in the data.
_TRANSCRIPT_RE = re.compile(r"[ -](?:transcripts?|mRNAs?)\b")


def transcript_strip_hsfold(s):
    """Strip ' transcript(s)' / ' mRNA(s)' descriptors, then hsfold.

    Names the gene's RNA product (e.g. EGFR mRNA -> EGFR), so it folds onto the
    gene. Global removal collapses combos like 'EGFR mRNA transcript' -> EGFR.
    Query-only, so HGNC surfaces are not mangled.
    """
    return hsfold(_TRANSCRIPT_RE.sub("", s))


# a ' variant(s)' / ' mutant(s)' / ' mutation(s)' descriptor on EITHER side
# (suffix 'EGFR mutant', 'TP53 mutation', or prefix 'mutant p53'), with its
# bordering separator. \b keeps 'mutational'/'covariant' intact. Lowercase only.
_VARIANT_RE = re.compile(
    r"[ -](?:variants?|mutants?|mutations?)\b|\b(?:variants?|mutants?|mutations?)[ -]")


def variant_strip_hsfold(s):
    """Strip ' variant(s)'/' mutant(s)'/' mutation(s)' descriptors (either side),
    then hsfold.

    'EGFR mutant', 'mutant p53', 'IDH1-mutant', 'TP53 mutation' name an altered
    form of the gene, so they fold onto it (EGFR, TP53, IDH1). Query-only; HGNC
    index untouched.
    """
    return hsfold(_VARIANT_RE.sub("", s))


# a trailing regulatory-region descriptor; longest alternative first so the
# whole suffix is removed (case-insensitive).
_PROMOTER_RE = re.compile(r"(?i)( gene promoter| promoter region| promoter)$")


def promoter_strip_hsfold(s):
    """Strip a trailing ' promoter'/' promoter region'/' gene promoter', then hsfold.

    'MGMT promoter', 'TERT promoter' name the gene's regulatory region, so they
    fold onto the gene. Query-only; HGNC surfaces are not mangled.
    """
    return hsfold(_PROMOTER_RE.sub("", s, count=1))


# <UPPER>-<digits> where the hyphen is mere punctuation (MMP-9 -> MMP9). Distinct
# from TYPE-2's hyphen->space fold: here the hyphen is *removed* (joined).
_DEHYPHEN_RE = re.compile(r"^[A-Z]+-[0-9]+$")
# drug / compound / cell-line names that coincidentally fit and hit a symbol
DEHYPHEN_EXCLUDE = {
    "YC-1",    # soluble guanylyl cyclase stimulator / HIF-1a inhibitor (drug)
    "ISO-1",   # MIF tautomerase inhibitor (drug)
    "IWR-1",   # Wnt pathway inhibitor (drug)
    "BMS-1",   # BMS- compound code (PD-L1 small molecule)
    "THP-1",   # monocytic leukaemia cell line, not a gene
}


def dehyphen_key(s):
    """For an <UPPER>-<digits> token (not an excluded compound), drop the hyphen
    and casefold (MMP-9 -> 'mmp9') to match the symbol fields; else a sentinel
    that cannot match. Query-only; the index holds casefolded HGNC surfaces, so
    the removal never touches HGNC strings.
    """
    if s not in DEHYPHEN_EXCLUDE and _DEHYPHEN_RE.match(s):
        return s.replace("-", "").casefold()
    return "\x00" + s


# General separator-deletion key. DEHYPHEN above only fires on <UPPER>-<digits>
# and folds the query onto an already-joined HGNC surface; hsfold only swaps
# '-'<->space (so 'PDL1' never reaches 'PD-L1'). delsep drops *all* hyphens AND
# spaces on BOTH sides, so 'PD-L1' / 'PD L1' / 'PDL1' (and 'CTLA-4' / 'CTLA4')
# share one key regardless of which side carries the separator. It runs as a late
# mop-up after the precise stages have had first claim, and is split across fields
# for precision (below).
#
# Precision: a destructive fold makes short keys collide easily. Empirically the
# false positives concentrate in short ALIAS/PREV hits -- coincidental collisions
# with a tiny alias ('Hb'->GSTM1, 'CR3'->CRIPTO3, 'GCB'->NPR2, 'DC-'->DCX) or with
# a non-gene concept/drug ('CAR-T'->CARTPT, 'Nec-1'->PCSK1, 'T-Ag'->LINC01194) --
# while short hits on the approved `symbol` field stay trustworthy (the string IS
# the official symbol modulo separators: 'CD4-'->CD4, 'ATM-'->ATM). So DELSEP runs
# as two passes: the `symbol` field at any length, then the noisier `alias/prev`
# fields with a >=4-char floor. Remaining same-key collisions degrade to
# 'ambiguous' (>=2 genes) rather than a wrong single link.
DELSEP_MIN = 4   # min joined-key length for alias/prev hits (symbol field is exempt)
DELSEP_EXCLUDE = set(DEHYPHEN_EXCLUDE) | {
    "CAR-T", "CAR T",   # CAR T-cell therapy, not the CARTPT neuropeptide ('cart' alias)
    "Nec-1",   # necrostatin-1 (RIPK1 inhibitor drug), not PCSK1's prev 'NEC1'
    "T-Ag", "T Ag",   # SV40 large T antigen, not the LINC01194 'TAG' alias
    "CNI",     # calcineurin inhibitor (drug class), not the NT5C1A 'CNI' alias
    "Cat B", "Cat-B",  # cathepsin B (= CTSB), not TYRP1's catalase-B 'CATB' alias
    "LDL-C",   # LDL cholesterol (a lab measurement), not COG2's prev 'LDLC'
}


def delsep(s):
    """Casefold and DELETE all separators (hyphens and spaces): a destructive fold
    so 'PD-L1' / 'PD L1' / 'PDL1' all share the key 'pdl1'. Used for the HGNC index
    (symmetric with the query keys below)."""
    return s.casefold().replace("-", "").replace(" ", "")


def delsep_key_sym(s):
    """Query key for the symbol-field DELSEP pass: plain delsep, with excluded
    look-alikes sentinel'd. No length floor -- an approved symbol is authoritative
    even when short ('CD4-'->CD4, 'ATM-'->ATM). Query-only; the index uses plain
    delsep, making the fold symmetric for HGNC surfaces."""
    return "\x00" + s if s in DELSEP_EXCLUDE else delsep(s)


def delsep_key_alias(s):
    """Query key for the alias/prev-field DELSEP pass: as delsep_key_sym but also
    sentinels keys shorter than DELSEP_MIN. Short alias/prev collisions are the
    main false-positive source ('Hb'->GSTM1, 'CR3'->CRIPTO3), whereas short
    approved-symbol hits are handled (trustworthily) by the symbol pass."""
    if s in DELSEP_EXCLUDE:
        return "\x00" + s
    k = delsep(s)
    return k if len(k) >= DELSEP_MIN else "\x00" + s


# a bare 'p' prefix (no delimiter) before an upper-case letter: the no-hyphen
# phospho/protein convention (pAKT, pSTAT3, pVHL). Distinct from the PHOSPHO
# stage, which needs a 'p-'/'phospho-' delimiter.
_PPREFIX_RE = re.compile(r"^p[A-Z]")
# constructs / vectors / reporters that coincidentally hit a symbol once de-p'd
PPREFIX_EXCLUDE = {
    "pMIR",    # pMIR-REPORT luciferase vector
    "pCEP4",   # pCEP4 episomal expression vector
    "pORF",    # ORF expression construct
    "pU2",     # U2 snRNA construct (not a phospho-protein)
    "pG8", "pG4", "pE3",   # constructs / non-specific, coincidental hits
}


def pprefix_strip_hsfold(s):
    """Strip a bare leading 'p' before an upper-case letter, then hsfold
    (pAKT -> AKT, pSTAT3 -> STAT3). Captures the no-delimiter phospho/protein
    convention. Excluded constructs/vectors are left intact. Query-only.
    """
    if s not in PPREFIX_EXCLUDE and _PPREFIX_RE.match(s):
        return hsfold(s[1:])
    return hsfold(s)


# a leading wild-type marker: 'WT ' (upper-case) or 'wild-type '/'wild type ' (any case)
_WT_RE = re.compile(r"^(?:WT |[Ww]ild[ -]?type )")


def wt_strip_hsfold(s):
    """Strip a leading wild-type marker ('WT ', 'wild-type ', 'Wild-type '), then
    hsfold (WT IDH1 -> IDH1, wild-type p53 -> p53).

    Denotes the unaltered gene, so it folds onto that gene. Query-only.
    """
    return hsfold(_WT_RE.sub("", s, count=1))


def nuclear_strip_hsfold(s):
    """Strip a leading 'nuclear ' localization prefix, then hsfold
    (nuclear YAP -> YAP). Denotes the gene's product (in the nucleus), so it
    folds onto the gene. Query-only.
    """
    return hsfold(s[8:] if s.startswith("nuclear ") else s)


# 'X target genes' -- the genes REGULATED BY X, not X itself
_TARGETGENES_RE = re.compile(r"(?i) target genes$")


def target_genes_strip_hsfold(s):
    """Strip a trailing ' target genes' suffix, then hsfold (MYC target genes ->
    MYC). NB: 'X target genes' denotes X's downstream targets, not X -- so the
    hit is the REGULATOR, recorded under a distinct match_mode (not a normal
    gene-identity link). Query-only.
    """
    return hsfold(_TARGETGENES_RE.sub("", s))


# spelled-out Greek-letter expansion: a lone Roman letter standing for a Greek
# symbol -> the spelled-out Greek word (e.g. IFN-g -> IFN-gamma, NF-kB ->
# NF-kappaB, TNF-a -> TNF-alpha).
# 'm'/'u' (mu) are deliberately excluded: a lone 'm' is overwhelmingly a 'mouse'
# prefix (mC2, mTOR), which would mis-expand (mC2 -> muC2 -> MUC2).
_GREEK_LET = {"a": "alpha", "b": "beta", "g": "gamma", "y": "gamma", "d": "delta",
              "e": "epsilon", "z": "zeta", "k": "kappa", "l": "lambda",
              "s": "sigma", "t": "tau", "w": "omega"}
# a 'lone' Greek-letter Roman char: not adjacent to lowercase letters (so it sits
# at a word boundary -- after start/hyphen/digit/upper, before end/hyphen/digit/upper)
_GREEK_LET_RE = re.compile(r"(?<![a-z])([abdegklstwyz])(?![a-z])")
# abbreviation expansion so the full HGNC *name* can match (IFN-gamma alone is
# not a surface, but 'interferon gamma' is IFNG). Only before a hyphen/space/
# digit, so it does not fire inside 'ILK', 'ILF3', etc.
_ABBREV = {"IFN": "interferon", "IL": "interleukin",
           "TGF": "transforming growth factor", "TNF": "tumor necrosis factor",
           "IGF": "insulin-like growth factor", "EGF": "epidermal growth factor",
           "PDGF": "platelet-derived growth factor",
           "VEGF": "vascular endothelial growth factor"}
_ABBREV_RE = re.compile(r"(?<![A-Za-z])(IFN|IL|TGF|TNF|IGF|EGF|PDGF|VEGF)(?=[-\s\d])")


def greek_letter_hsfold(s):
    """Expand a lone Roman letter that stands for a Greek symbol to its spelled-
    out Greek word, then hsfold (NF-kB -> NF-kappaB -> NFKB1; TNF-a -> TNF-alpha
    -> TNF; p38a -> p38alpha -> MAPK14). Word-internal letters are untouched, and
    only expansions that actually hit an HGNC surface link (e.g. 'p53' -> 'pi53'
    falls through). Values shorter than 3 chars are skipped. Query-only.
    """
    if len(s) < 3:
        return "\x00" + s
    return hsfold(_GREEK_LET_RE.sub(lambda m: _GREEK_LET[m.group(1)], s))


def greek_abbrev_hsfold(s):
    """As greek_letter_hsfold, but also expands a leading abbreviation to its full
    word (IFN-g -> interferon gamma -> IFNG; IL-6 receptor -> interleukin 6
    receptor -> IL6R). Run only on what the plain Greek-letter pass left over, so
    the abbreviation expansion never overrides a shorter form that already linked
    via an alias (e.g. TNF-a -> 'TNF-alpha' alias, not 'tumor necrosis factor
    alpha' which is no HGNC surface). Query-only.
    """
    if len(s) < 3:
        return "\x00" + s
    t = _GREEK_LET_RE.sub(lambda m: _GREEK_LET[m.group(1)], s)
    t = _ABBREV_RE.sub(lambda m: _ABBREV[m.group(1)], t)
    return hsfold(t)


# classic 'c-' proto-oncogene spellings: HGNC lists some alias forms (KIT 'C-Kit',
# SRC 'c-src', ABL1 'c-ABL'), but hyphen/case variants and descriptor-bearing forms
# miss them ('cKit', 'c-Src kinase', 'c-Abl'). Map the prefix to its symbol; a word
# boundary keeps look-alikes out (e.g. 'ckitCSca-1C' is not c-Kit). c-Met is in the
# map too (the HGF receptor = single gene MET, not an HGNC 'c-Met' surface), plus
# the descriptive 'HGF receptor' / 'MET receptor' forms via a second pattern.
_CONCO_MAP = {"kit": "KIT", "src": "SRC", "abl": "ABL1",
              "fos": "FOS", "jun": "JUN", "myc": "MYC", "ret": "RET",
              "raf": "RAF1", "cbl": "CBL", "mpl": "MPL", "met": "MET",
              "ros": "ROS1", "sis": "PDGFB", "fms": "CSF1R", "yes": "YES1",
              "fgr": "FGR", "fes": "FES", "mos": "MOS", "rel": "REL",
              "mil": "RAF1",   # c-mil = avian c-raf homolog -> RAF1
              "myb": "MYB", "erba": "THRA",   # c-erbA = thyroid receptor (cf. c-erbB = EGFR)
              "crk": "CRK", "maf": "MAF", "ski": "SKI",
              "fps": "FES",    # v-fps / v-fes are the same oncogene -> FES
              "pim": "PIM1", "cot": "MAP3K8",  # c-cot / Tpl2 -> MAP3K8
              "mer": "MERTK", "sea": "MST1R"}  # v-sea -> MST1R (RON)
_CONCO_RE = re.compile(
    r"^c-?(kit|src|abl|fos|jun|myc|ret|raf|cbl|mpl|met|ros|sis|fms|yes|fgr|fes"
    r"|mos|rel|mil|myb|erba|crk|maf|ski|fps|pim|cot|mer|sea)\b", re.I)
# 'HGF receptor' / 'MET receptor' (no leading 'c-') also denote the single gene MET.
_METR_RE = re.compile(r"^(?:hgf|met)[ -]receptors?\b", re.I)
# numbered ets oncogene: c-Ets-1 -> ETS1, c-Ets-2 -> ETS2 (the bare unnumbered
# 'c-Ets' is ambiguous and is handled in the COMPLEX stage via ONCO_AMBIG_MAP).
_CETS_RE = re.compile(r"^c-?ets-?([12])\b", re.I)
# ras family: the isoform prefix sets the gene (c-Ha-ras -> HRAS, c-Ki-ras -> KRAS,
# c-N-ras -> NRAS). Bare 'c-ras' (no isoform) is ambiguous -> ONCO_AMBIG_MAP.
_CRAS_RE = re.compile(r"^c-?(ha|ki|[hkn])-?ras\b", re.I)
_CRAS_GENE = {"ha": "HRAS", "h": "HRAS", "ki": "KRAS", "k": "KRAS", "n": "NRAS"}
# bcl family: numbered c-Bcl-2/3/6 -> BCL2/3/6 (bare 'c-Bcl' is non-standard and
# spans a 20+ gene family, so it is left unmatched).
_CBCL_RE = re.compile(r"^c-?bcl-?([236])\b", re.I)
# erbB family: 'c-erbB-N' is an old synonym for the ErbB receptors, but HGNC lists
# only 'c-ERB-2' (no inner 'b'), so it never matches directly. 1 -> EGFR (ErbB1),
# 2/3/4 -> ERBBn (trailing descriptor 'c-erbB-2 gene/protein' ignored).
_CERBB_RE = re.compile(r"^c-?erbb-?([1-4])\b", re.I)
_CERBB_GENE = {"1": "EGFR", "2": "ERBB2", "3": "ERBB3", "4": "ERBB4"}
# guard: 'c-Jun N-terminal kinase' (JNK) is MAPK8/9/10, NOT the JUN transcription
# factor -- so a c-Jun form mentioning the kinase must NOT map to JUN.
_JNK_RE = re.compile(r"(?i)termin|jnk")


def conco_hsfold(s):
    """Map a leading c-Kit/c-Src/c-Abl/c-Fos/c-Jun/c-Myc/c-Ret/c-Raf/c-Cbl/c-Mpl/
    c-Met spelling -- or 'HGF receptor' / 'MET receptor' -- to its HGNC symbol
    (query-only; no-match sentinel otherwise)."""
    m = _CONCO_RE.match(s)
    if m:
        g = m.group(1).lower()
        if g == "jun" and _JNK_RE.search(s):   # c-Jun N-terminal kinase -> not JUN
            return "\x00" + s
        return hsfold(_CONCO_MAP[g])
    me = _CETS_RE.match(s)
    if me:                                      # c-Ets-1 -> ETS1, c-Ets-2 -> ETS2
        return hsfold("ETS" + me.group(1))
    mr = _CRAS_RE.match(s)
    if mr:                                      # c-Ha-ras -> HRAS, c-Ki-ras -> KRAS, ...
        return hsfold(_CRAS_GENE[mr.group(1).lower()])
    mb = _CBCL_RE.match(s)
    if mb:                                      # c-Bcl-2/3/6 -> BCL2/3/6
        return hsfold("BCL" + mb.group(1))
    mc = _CERBB_RE.match(s)
    if mc:                                      # c-erbB-1 -> EGFR, c-erbB-2/3/4 -> ERBBn
        return hsfold(_CERBB_GENE[mc.group(1)])
    if _METR_RE.match(s):                       # HGF/MET receptor -> MET
        return hsfold("MET")
    return "\x00" + s


# type-1 missense convention: <approved symbol> + [AA]<pos>[AA] (e.g. IDH1R132H).
# Residues restricted to the 20 standard amino-acid single-letter codes.
_AA = "ACDEFGHIKLMNPQRSTVWY"
_MISSENSE_RE = re.compile(rf"^(.+?)[{_AA}]\d+[{_AA}]$")
# curated names that coincidentally fit the pattern but are not gene missense
MISSENSE_EXCLUDE = {
    "KYA1797K",   # KYA-1797K is a Wnt/beta-catenin inhibitor (small molecule)
    "SSTR2A",     # somatostatin receptor 2A gene, not an SST missense
    "PIPK2A",     # PIP4K2A kinase gene, not a PIP missense
}


def missense_strip(s):
    """If s is <prefix> + a valid-AA missense [AA]<digits>[AA], return the
    prefix (to be matched against APPROVED symbols, case-sensitive); else return
    s unchanged. Known small-molecule compounds are left intact so they don't
    link. Query-only; the index is built over approved symbols verbatim.
    """
    if s in MISSENSE_EXCLUDE:
        return s
    m = _MISSENSE_RE.match(s)
    return m.group(1) if m else s


# space-delimited missense: <prefix> ' ' [AA]<pos>[AA]  (e.g. TP53 R248L)
_MISSENSE_SP_RE = re.compile(rf"^(.+) [{_AA}]\d+[{_AA}]$")


def missense_sp_strip(s):
    """Like missense_strip but for the space-delimited form (TP53 R248L -> TP53).

    Same guards: matched against the APPROVED symbol field, valid-AA residues
    only, excluded compounds left intact.
    """
    if s in MISSENSE_EXCLUDE:
        return s
    m = _MISSENSE_SP_RE.match(s)
    return m.group(1) if m else s


def inhibitor_strip(s):
    """Strip a terminal inhibitor 'i' (PARPi -> PARP, BRAFi -> BRAF).

    The '-i' shorthand names a drug/inhibitor AGAINST the gene product (a target,
    like 'anti-'), so it links to the gene but is tagged distinctly. Case-
    sensitive (type-1), query-only; matched against the symbol fields verbatim.
    """
    return s[:-1] if s.endswith("i") else s


def build_index(docs, fields, keyfn):
    """keyfn(name) -> {gene_id: {surface: set(fields)}} + gid->symbol.

    Keeping surface->fields (rather than two flat sets) lets a match be
    explained: which exact HGNC string was hit, in which field.
    """
    index = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    gid_symbol = {}
    for doc in docs:
        gid = doc.get("hgnc_id") or doc.get("symbol")
        gid_symbol[gid] = doc.get("symbol", "")
        for field in fields:
            val = doc.get(field)
            if val is None:
                continue
            for s in set([val] if isinstance(val, str) else val):
                index[keyfn(s)][gid][s].add(field)
    return index, gid_symbol


def _write_block(path, header_line, rows):
    body = ("\n".join(rows) + "\n") if rows else ""
    path.write_text(header_line + "\n" + body, encoding="utf-8")


def match_pass(rows, fields, keyfn, unmatched_tsv, label, mode_override=None,
               index_keyfn=None, approved_guard=False):
    """Bucket rows into matched/ambiguous/unmatched against one field set + key.

    Returns a stats dict that includes the matched/ambiguous JSON *libraries*
    (keyed by the GENETIC named entity) so the caller can consolidate them across
    passes, plus the unmatched rows so it can cascade them onward. No per-step
    matched TSV is written; the unmatched TSV is written only when
    `unmatched_tsv` is given -- intermediate leftovers stay in memory, so only
    the final stage persists its unmatched bucket (to `unmatched.tsv`).

    `mode_override` forces the recorded match_mode (used by the phospho stage,
    where the link is explained by the stripped prefix, not by case/hyphen).
    `index_keyfn` lets the HGNC index use a different (un-stripped) transform
    than the query keyfn -- defaults to keyfn.
    """
    docs = match_pass.docs
    index, gid_symbol = build_index(docs, fields, index_keyfn or keyfn)

    def fields_in_order(fs):
        return "|".join(f for f in fields if f in fs)

    def primary_field(fs):
        return next(f for f in fields if f in fs)

    def classify(entity, sf):
        """(match_field, match_mode) explaining how `entity` hit surfaces `sf`.

        Exact surface present -> case-sensitive; else a casefold-equal surface
        -> case-insensitive; else the only remaining way to share the key is a
        hyphen<->whitespace swap.
        """
        if entity in sf:
            return primary_field(sf[entity]), "case-sensitive"
        ecf = entity.casefold()
        ci = [s for s in sf if s.casefold() == ecf]
        if ci:
            return primary_field(sf[ci[0]]), "case-insensitive"
        return primary_field(sf[next(iter(sf))]), "hyphen vs white space"

    matched, unmatched = [], []
    matched_lib, ambiguous_lib = {}, {}
    field_counts = defaultdict(int)
    mode_counts = defaultdict(int)
    occ = {"total": 0, "matched": 0, "ambiguous": 0, "unmatched": 0}

    for ln in rows:
        name, rest = ln.split("\t", 1)            # rest = "occ\tn_forms"
        n_occ = int(rest.split("\t", 1)[0])
        occ["total"] += n_occ
        genes = index.get(keyfn(name))
        guard_gid = None
        if approved_guard and genes and len(genes) > 1:
            # disambiguate by preferring the gene whose APPROVED symbol (not just
            # an alias) produced the key; only if exactly one such gene exists.
            appr = [g for g in genes
                    if any("symbol" in fs for fs in genes[g].values())]
            if len(appr) == 1:
                guard_gid = appr[0]
        if not genes:
            unmatched.append(ln)
            occ["unmatched"] += n_occ
        elif len(genes) == 1 or guard_gid is not None:
            gid = guard_gid if guard_gid is not None else next(iter(genes))
            sf = genes[gid]                       # {surface: set(fields)}
            f = fields_in_order(set().union(*sf.values()))
            surf = "|".join(sorted(sf))
            matched.append(f"{name}\t{rest}\t{f}\t{gid}"
                           f"\t{gid_symbol.get(gid, '')}\t{surf}")
            field_counts[f] += 1
            occ["matched"] += n_occ
            mfield, mmode = classify(name, sf)
            if mode_override:
                mmode = mode_override
            mode_counts[mmode] += 1
            matched_lib[name] = {                 # entity -> symbol + annotations
                "hgnc_symbol": gid_symbol.get(gid, ""),
                "occurrences": n_occ,
                "match_field": mfield,
                "match_mode": mmode,
            }
        else:
            occ["ambiguous"] += n_occ
            gids = sorted(genes)
            combined = defaultdict(set)            # surface -> fields, across genes
            for g in genes:
                for s, fs in genes[g].items():
                    combined[s] |= fs
            mfield, mmode = classify(name, combined)
            if mode_override:
                mmode = mode_override
            ambiguous_lib[name] = {                # entity -> candidate symbols + annotations
                "hgnc_symbol": [gid_symbol.get(g, "") for g in gids],
                "occurrences": n_occ,
                "match_field": mfield,
                "match_mode": mmode,
                "n_genes": len(gids),
                "hgnc_ids": gids,
            }

    if unmatched_tsv is not None:                 # only the final stage persists
        _write_block(unmatched_tsv, TSV_HEADER, unmatched)

    return {
        "label": label, "fields": fields,
        "n_total": len(rows), "n_matched": len(matched),
        "n_ambiguous": len(ambiguous_lib), "n_unmatched": len(unmatched),
        "by_field": dict(field_counts), "by_mode": dict(mode_counts), "occ": occ,
        "sample": matched[:10], "unmatched_rows": unmatched,
        "matched_lib": matched_lib, "ambiguous_lib": ambiguous_lib,
    }


# =======================================================================
# STAGE 4 -- multi-word word-order cosine        (was cosine_similarity.py)
# =======================================================================
COSINE_THRESHOLD = 0.99
NAME_FIELDS = ["name", "alias_name", "prev_name"]


def _toks(s):
    return s.casefold().split()


def _cosine(a, b):
    dot = sum(a[t] * b.get(t, 0) for t in a)
    if not dot:
        return 0.0
    na = math.sqrt(sum(v * v for v in a.values()))
    nb = math.sqrt(sum(v * v for v in b.values()))
    return dot / (na * nb)


def stage4_cosine(rows, docs):
    """Multi-word value vs multi-word HGNC name value, word-cosine >= threshold."""
    candidates = []                                   # (surface, vec, field, gid, sym)
    for d in docs:
        gid = d.get("hgnc_id") or d.get("symbol")
        sym = d.get("symbol", "")
        for field in NAME_FIELDS:
            val = d.get(field)
            if val is None:
                continue
            for s in ([val] if isinstance(val, str) else val):
                tk = _toks(s)
                if len(tk) >= 2:
                    candidates.append((s, Counter(tk), field, gid, sym))
    df = Counter()
    postings = defaultdict(list)
    for i, (_, cv, _, _, _) in enumerate(candidates):
        for t in cv:
            df[t] += 1
            postings[t].append(i)

    out = []
    n_inputs_mw = 0
    for ln in rows:
        name, rest = ln.split("\t", 1)
        occ = int(rest.split("\t", 1)[0])
        itk = _toks(name)
        if len(itk) < 2:
            continue
        n_inputs_mw += 1
        iv = Counter(itk)
        cand_ids = set()
        for t in sorted(iv, key=lambda x: df.get(x, 0))[:2]:
            cand_ids.update(postings.get(t, ()))
        per_surface = {}
        for ci in cand_ids:
            s, cv, field, gid, sym = candidates[ci]
            if s.casefold() == name.casefold():
                continue
            c = _cosine(iv, cv)
            if c >= COSINE_THRESHOLD:
                rec = per_surface.setdefault(s, [c, set(), set()])
                rec[0] = max(rec[0], c)
                rec[1].add(field)
                rec[2].add((gid, sym))
        for s, (c, fields, idsym) in per_surface.items():
            fstr = "|".join(f for f in NAME_FIELDS if f in fields)
            ids = ";".join(sorted({g for g, _ in idsym}))
            syms = "|".join(sorted({y for _, y in idsym}))
            out.append({"input": name, "hgnc": s, "cos": c, "field": fstr,
                        "ids": ids, "syms": syms, "occ": occ})
    out.sort(key=lambda r: (-r["occ"], r["input"].casefold(), r["hgnc"].casefold()))

    # Aggregate pairs per GENETIC entity into a JSON library, same shape as the
    # exact-step libraries plus the retained cosine value. A handful of entities
    # match several HGNC surfaces of the *same* gene; those collapse to one entry
    # (max cosine, joined surfaces, highest-priority field).
    agg = {}
    for r in out:
        a = agg.setdefault(r["input"], {"syms": set(), "fields": set(),
                                        "surfaces": set(), "cos": 0.0,
                                        "occ": r["occ"]})
        a["syms"].update(r["syms"].split("|"))
        a["fields"].update(r["field"].split("|"))
        a["surfaces"].add(r["hgnc"])
        a["cos"] = max(a["cos"], r["cos"])

    lib = {}
    for e, a in agg.items():
        primary = next((f for f in NAME_FIELDS if f in a["fields"]),
                       sorted(a["fields"])[0])
        lib[e] = {
            "hgnc_symbol": "|".join(sorted(a["syms"])),
            "occurrences": a["occ"],
            "match_field": primary,
            "match_mode": "word-order cosine",
            "cosine": round(a["cos"], 4),
            "hgnc_value": "|".join(sorted(a["surfaces"])),
        }

    (SUMMARY_DIR / "roman_cosine.json").write_text(
        json.dumps(lib, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"n_candidates": len(candidates), "n_inputs_mw": n_inputs_mw,
            "rows": out, "n_entities": len(lib)}


# =======================================================================
# MIRNA stage -- map mature microRNA names to HGNC precursor gene(s)
# =======================================================================
# HGNC tracks miRNA *precursor* genes (MIR21, MIR124-1/2/3); the corpus uses
# *mature* names (miR-21, miR-124, miR-34c-5p). We canonicalise both to a mature
# key and map; a mature miR from several genomic loci is genuinely 1-to-many
# (ambiguous). Bespoke because the generic case/hyphen key cannot bridge
# 'miR-21' <-> 'hsa-mir-21' / 'MIR21'.
MIR_FIELDS = ["symbol", "alias_symbol", "prev_symbol"]
_MIR_KEYRE = re.compile(r"^((?:mir|let)-\d+[a-z]*)(?:-\d+)?$")
# compound / family / reagent forms that are not a single mappable miRNA
_MIR_EXCLUDE = re.compile(
    r"[/~&]| and |cluster|family|combo|sponge|\binh\b|-/-|^mir-?nc$|mirfp", re.I)


def _mir_mature(s):
    """Canonical mature key ('mir-<id>' / 'let-7<id>') for a single miRNA, else None.

    Handles both the corpus form (miR-34c-5p, miRNA-34a, hsa-miR-21) and the HGNC
    forms (MIR21, hsa-mir-124-1) -> strips species prefix, mature arm (-5p/-3p),
    and the precursor copy number, lower-cases, and normalises 'mir-'/'let-7'.
    """
    if _MIR_EXCLUDE.search(s):
        return None
    t = re.sub(r"^hsa[-_]", "", s.strip().lower())
    t = t.replace("microrna", "mir").replace("mirna", "mir")
    t = re.sub(r"[-_ ]*(?:3p|5p)\b", "", t)
    t = t.rstrip("*").strip()
    t = re.sub(r"^mir\s*-?\s*", "mir-", t)
    t = re.sub(r"^let\s*-?\s*7", "let-7", t)
    m = _MIR_KEYRE.match(t)
    return m.group(1) if m else None


def mirna_pass(rows, docs, unmatched_tsv, label):
    """Map mature miR names to HGNC precursor gene(s); 1 locus -> matched,
    >=2 loci -> ambiguous. Returns a stats dict shaped like match_pass."""
    key2gf = defaultdict(lambda: defaultdict(set))      # mature key -> {gid: fields}
    gid_symbol = {}
    for d in docs:
        sym = d.get("symbol", "")
        if not sym.startswith("MIR"):
            continue
        gid = d.get("hgnc_id") or sym
        gid_symbol[gid] = sym
        for field in MIR_FIELDS:
            v = d.get(field)
            for s in ([v] if isinstance(v, str) else (v or [])):
                k = _mir_mature(s)
                if k:
                    key2gf[k][gid].add(field)

    def primary(fs):
        return next(f for f in MIR_FIELDS if f in fs)

    matched_lib, ambiguous_lib, unmatched = {}, {}, []
    field_counts = defaultdict(int)
    occ = {"total": 0, "matched": 0, "ambiguous": 0, "unmatched": 0}
    sample = []
    for ln in rows:
        name, rest = ln.split("\t", 1)
        n_occ = int(rest.split("\t", 1)[0])
        occ["total"] += n_occ
        k = _mir_mature(name)
        genes = key2gf.get(k) if k else None
        if not genes:
            unmatched.append(ln)
            occ["unmatched"] += n_occ
        elif len(genes) == 1:
            gid = next(iter(genes))
            f = primary(genes[gid])
            matched_lib[name] = {"hgnc_symbol": gid_symbol[gid], "occurrences": n_occ,
                                 "match_field": f, "match_mode": "miR mature→gene"}
            field_counts[f] += 1
            occ["matched"] += n_occ
            if len(sample) < 10:
                sample.append(f"{name}\t{rest}\t{f}\t{gid}\t{gid_symbol[gid]}")
        else:
            gids = sorted(genes)
            allf = primary({fld for fs in genes.values() for fld in fs})
            ambiguous_lib[name] = {"hgnc_symbol": [gid_symbol[g] for g in gids],
                                   "occurrences": n_occ, "match_field": allf,
                                   "match_mode": "miR mature→gene",
                                   "n_genes": len(gids), "hgnc_ids": gids}
            occ["ambiguous"] += n_occ

    if unmatched_tsv is not None:
        _write_block(unmatched_tsv, TSV_HEADER, unmatched)
    return {"label": label, "fields": MIR_FIELDS,
            "n_total": len(rows), "n_matched": len(matched_lib),
            "n_ambiguous": len(ambiguous_lib), "n_unmatched": len(unmatched),
            "by_field": dict(field_counts),
            "by_mode": {"miR mature→gene": len(matched_lib)},
            "occ": occ, "sample": sample, "unmatched_rows": unmatched,
            "matched_lib": matched_lib, "ambiguous_lib": ambiguous_lib}


# =======================================================================
# HISTONE stage -- map histone marks / mutants / variants to HGNC histone genes
# =======================================================================
# HGNC tracks histone *genes* (clusters): the H3 protein is encoded by ~23 genes,
# H4 by 16, etc. A mark like H3K27me3 or a mutant like H3.3K27M is on the histone
# PROTEIN, so it maps to the histone variant's gene(s) -- a single gene for the
# specific variants (H2A.X -> H2AX), a small set for H3.3 (H3-3A/B), or the whole
# family for a generic mark (H3K27me3 -> all H3 genes, ambiguous).
_HIST_ENZ = re.compile(
    r"deacetylas|demethylas|methyltransferas|acetyltransferas|acetyl transferas|"
    r"methyl transferas|monomethyltransferas|acetylase|methylase|kinase|inhibitor", re.I)
_HIST_PREFIX = re.compile(
    r"^(?:phospho|cit|acetyl(?:ated)?|anti[- ]?phospho|anti|g|wt|oncohistones?|core)[- ]+",
    re.I)
# non-histone look-alikes that slip through (e.g. histamine H3 receptor + 5-HT)
_HIST_EXCLUDE = {"H3-5HT"}


def _histone_variant(s):
    """Return the histone variant key ('H3', 'H3.3', 'H2A.X', ...) for a mark /
    mutant / variant string, else None. Cell-line look-alikes (H358, H460, H19)
    are rejected: a histone token must be followed by a residue letter, a
    separator, a known variant, or end -- not a bare multi-digit number.
    """
    if s in _HIST_EXCLUDE or _HIST_ENZ.search(s):
        return None
    t = _HIST_PREFIX.sub("", s.lower().strip())
    m = re.match(r"^histones?[ -]*", t)      # only treat a bare number as Hn when
    if m:                                    # it actually followed 'histone'
        t = t[m.end():]
        if t[:1].isdigit():
            t = "h" + t                      # 'histone 3' -> '3' -> 'h3'
    if re.match(r"(?:macroh2a|mh2a)", t): return "macroH2A"
    if re.match(r"h2a[.\- ]*z", t): return "H2A.Z"
    if re.match(r"h2a[.\- ]*x", t): return "H2A.X"
    if re.match(r"(?:cenp[.\- ]*a|cenh3)", t): return "CENP-A"
    if re.match(r"h3f3", t): return "H3.3"   # H3F3A/B are the H3.3 genes
    if re.match(r"h3[.\- ]*3", t): return "H3.3"
    if re.match(r"h3[.\- ]*4", t): return "H3.4"
    if re.match(r"h3[.\- ]*5", t): return "H3.5"
    if re.match(r"h3(?:[.\- ]|[krsgtqa]|f3|$)", t): return "H3"
    if re.match(r"h4(?:[.\- ]|[krs]|$)", t): return "H4"
    if re.match(r"h2b", t): return "H2B"
    if re.match(r"h2a", t): return "H2A"
    if re.match(r"h1(?:[.\- ]|[krs]\d|f|$)", t): return "H1"
    return None


def histone_pass(rows, docs, unmatched_tsv, label):
    """Map histone marks/mutants/variants to HGNC histone gene(s); a variant with
    one gene -> matched, a family (or H3.3/H2A.Z) -> ambiguous candidate set."""
    grp = defaultdict(set)
    sym2id = {}
    for d in docs:
        sym = d.get("symbol")
        sym2id[sym] = d.get("hgnc_id") or sym
        for g in (d.get("gene_group") or []):
            grp[g].add(sym)
    # generic H3 marks sit on canonical H3 only: the clustered H3C* genes
    # (H3.1/H3.2) + H3.3 (H3-3A/B); not the centromeric/testis variants.
    canonical_h3 = grp["H3 histones"] - {"CENPA", "H3-4", "H3-5", "H3-7",
                                         "H3Y1", "H3Y2"}
    VAR = {
        "macroH2A": {"MACROH2A1", "MACROH2A2"}, "H2A.Z": {"H2AZ1", "H2AZ2"},
        "H2A.X": {"H2AX"}, "CENP-A": {"CENPA"},
        "H3.3": {"H3-3A", "H3-3B"}, "H3.4": {"H3-4"}, "H3.5": {"H3-5"},
        "H3": canonical_h3, "H4": grp["H4 histones"],
        "H2B": grp["H2B histones"], "H2A": grp["H2A histones"],
        "H1": grp["H1 histones"],
    }
    matched_lib, ambiguous_lib, unmatched = {}, {}, []
    field_counts = defaultdict(int)
    occ = {"total": 0, "matched": 0, "ambiguous": 0, "unmatched": 0}
    sample = []
    for ln in rows:
        name, rest = ln.split("\t", 1)
        n_occ = int(rest.split("\t", 1)[0])
        occ["total"] += n_occ
        v = _histone_variant(name)
        genes = VAR.get(v) if v else None
        if not genes:
            unmatched.append(ln)
            occ["unmatched"] += n_occ
        elif len(genes) == 1:
            sym = next(iter(genes))
            gid = sym2id.get(sym, sym)
            matched_lib[name] = {"hgnc_symbol": sym, "occurrences": n_occ,
                                 "match_field": v, "match_mode": "histone"}
            field_counts[v] += 1
            occ["matched"] += n_occ
            if len(sample) < 10:
                sample.append(f"{name}\t{rest}\t{v}\t{gid}\t{sym}")
        else:
            syms = sorted(genes)
            ambiguous_lib[name] = {"hgnc_symbol": syms, "occurrences": n_occ,
                                   "match_field": v, "match_mode": "histone",
                                   "n_genes": len(syms),
                                   "hgnc_ids": [sym2id.get(s, s) for s in syms]}
            occ["ambiguous"] += n_occ
    if unmatched_tsv is not None:
        _write_block(unmatched_tsv, TSV_HEADER, unmatched)
    return {"label": label, "fields": ["gene_group"],
            "n_total": len(rows), "n_matched": len(matched_lib),
            "n_ambiguous": len(ambiguous_lib), "n_unmatched": len(unmatched),
            "by_field": dict(field_counts),
            "by_mode": {"histone": len(matched_lib)},
            "occ": occ, "sample": sample, "unmatched_rows": unmatched,
            "matched_lib": matched_lib, "ambiguous_lib": ambiguous_lib}


# =======================================================================
# COMPLEX stage -- map protein-complex / heterodimer names to subunit gene(s)
# =======================================================================
# Some entities name an assembled complex that HGNC catalogs only as its subunit
# genes -- e.g. IL-12 (the IL-12p70 heterodimer) = IL12A (p35) + IL12B (p40);
# there is no single 'IL12' gene. Such names map to the subunit set (ambiguous),
# while a subunit-specific form (IL-12 p40) resolves to one gene.
COMPLEX_MAP = {
    "il12":    ("IL12A", "IL12B"),   # IL-12 (p70) heterodimer = p35 + p40
    "il12p70": ("IL12A", "IL12B"),
    "il12p35": ("IL12A",),           # p35 subunit
    "il12p40": ("IL12B",),           # p40 subunit
    "il23":    ("IL12B", "IL23A"),   # IL-23 = p19 (IL23A) + p40 (IL12B)
    "il23p19": ("IL23A",),
    "il23p40": ("IL12B",),
    "il27":    ("EBI3", "IL27"),     # IL-27 = p28 (IL27) + EBI3
    "il27p28": ("IL27",),
    "ap1":     ("FOS", "FOSB", "JUN", "JUND"),  # AP-1 = Fos/Jun family dimer
}
# Paralog receptor families: a bare family name (no specific member) is genuinely
# ambiguous -- it is ONE-OF the family's genes, not an assembled complex, so it
# carries a distinct tag from the heterodimers above. e.g. 'VEGF receptor' /
# 'VEGFRs' = VEGFR-1/2/3 = FLT1 / KDR / FLT4. A specific member (VEGFR2, KDR)
# still links to its single gene upstream and never reaches this map.
_ERBB = ("EGFR", "ERBB2", "ERBB3", "ERBB4")     # EGFR/ERBB (HER) family members
_TRK = ("NTRK1", "NTRK2", "NTRK3")              # Trk (TrkA/B/C) high-affinity receptors
_NTR = ("NGFR", "NTRK1", "NTRK2", "NTRK3")      # neurotrophin receptors = Trk + p75 (NGFR)
FAMILY_MAP = {
    "vegfreceptor":  ("FLT1", "KDR", "FLT4"),   # VEGF receptor  (VEGFR-1/2/3)
    "vegfreceptors": ("FLT1", "KDR", "FLT4"),   # VEGF receptors
    "vegfrs":        ("FLT1", "KDR", "FLT4"),   # VEGFRs
    "pdgfreceptor":  ("PDGFRA", "PDGFRB"),      # PDGF receptor  (PDGFR-alpha/beta)
    "pdgfreceptors": ("PDGFRA", "PDGFRB"),      # PDGF receptors
    "pdgfrs":        ("PDGFRA", "PDGFRB"),      # PDGFRs (bare PDGFR -> PDGFRB upstream)
    "fgfreceptor":   ("FGFR1", "FGFR2", "FGFR3", "FGFR4"),  # FGF receptor (FGFR1-4)
    "fgfreceptors":  ("FGFR1", "FGFR2", "FGFR3", "FGFR4"),  # FGF receptors
    "fgfr":          ("FGFR1", "FGFR2", "FGFR3", "FGFR4"),  # bare FGFR (unmatched upstream)
    "fgfrs":         ("FGFR1", "FGFR2", "FGFR3", "FGFR4"),  # FGFRs
    # EGFR/ERBB (HER) family = EGFR (ErbB1) + ERBB2/3/4. Only the family-level
    # surfaces route here; a specific member (EGFR, HER2/ERBB2, c-erbB-2) and the
    # HGNC-sanctioned singulars (ErbB/ERBB -> EGFR via prev_symbol) link upstream.
    "egfrfamily":         _ERBB, "egfrfamilymembers": _ERBB,  # EGFR family (members)
    "erbbreceptor":       _ERBB, "erbbreceptors":     _ERBB,  # (er)bB receptor(s)
    "erbbproteins":       _ERBB, "erbboncoproteins":  _ERBB,  # ERBB proteins/oncoproteins
    "erbbfamilyreceptors": _ERBB, "erbbs":            _ERBB,  # erbB family receptors / ErbBs
    "her":                _ERBB, "herfamily":         _ERBB,  # HER (HER family)
    "herreceptors":       _ERBB,                              # HER receptors
    # NGF/TRK neurotrophin-receptor family. The Trk (high-affinity) receptors are
    # NTRK1/2/3 (TrkA/B/C); specific members (TrkA->NTRK1 ...) link upstream. Bare
    # 'TRK'/'Trk' otherwise mis-links to [TPM3, NTRK1] via a fusion alias, so it is
    # overridden here. Generic 'neurotrophin receptor' also covers p75 (NGFR).
    "trk":           _TRK, "trks":          _TRK, "trkr":   _TRK,  # TRK / Trk / TRKs / TrkR
    "trkreceptor":   _TRK, "trkreceptors":  _TRK,                  # Trk receptor(s)
    "trktyrosinekinase": _TRK,                                     # Trk tyrosine kinase
    "ntrk":          _TRK, "ntrks":         _TRK, "ntrk13": _TRK,  # NTRK / NTRKs / NTRK1-3
    "ntrkreceptor":  _TRK, "ntrkreceptors": _TRK,                  # NTRK receptor(s)
    "ntrkfusion":    _TRK, "ntrkfusions":   _TRK,                  # NTRK fusion(s)
    "ntrkfusionproteins": _TRK,                                    # NTRK fusion proteins
    "neurotrophinreceptor":  _NTR, "neurotrophinreceptors": _NTR,  # neurotrophin receptor(s)
}
# proto-oncogene aliases that denote >1 cellular gene, so a bare (unnumbered) form
# is ambiguous. e.g. the 'ets' oncogene has two cellular homologs -- bare 'c-Ets'
# is ETS1 or ETS2 (a numbered 'c-Ets-1'/'c-Ets-2' resolves to one gene in the
# C-ONCOGENE stage). Tagged distinctly from the receptor families above.
ONCO_AMBIG_MAP = {
    "cets": ("ETS1", "ETS2"),            # c-Ets (unnumbered) -> ETS1 / ETS2
    "cras": ("HRAS", "KRAS", "NRAS"),    # c-ras (no isoform) -> HRAS / KRAS / NRAS
}
_COMPLEX_DESC = re.compile(r"(?i) (?:mrna|gene|protein|cytokine)$")


def _complex_key(s):
    """Normalise a complex name: drop a trailing form descriptor, casefold, and
    remove hyphens/spaces ('IL-12' / 'IL-12 mRNA' / 'IL12' -> 'il12')."""
    return re.sub(r"[- ]", "", _COMPLEX_DESC.sub("", s.casefold()))


# NF-kB transcription factor: like the greek pipeline's NF-kappaB complex, the
# generic NF-kB name denotes the canonical p50/p65 dimer NFKB1 + RELA -- not a
# single gene -- so it maps to that subunit set (ambiguous). In the Roman set the
# 'k' stands for kappa, so both the 'k' spelling (NF-kB, NFkB) and the spelled
# 'kappa' form (NF-kappaB) are recognised. Specific genes are NOT the complex and
# must keep linking to themselves: NFKB1 / NFKB2 (the p50/p52 precursor genes) and
# the IkB-family inhibitors NFKBIA/IB/IE/IZ/ID, plus the named-subunit forms
# (NF-kB p65 = RELA, p50 = NFKB1, ...) -- all excluded below, mirroring greek's
# `_NFKB_SINGLE` + p65/p50 guards.
NFKB_GENES = ("NFKB1", "RELA")             # canonical p50/p65 dimer
_NFKB_STEM = re.compile(r"nfk(?:appab|b)")  # the NF-kB stem in a separator-free key
_NFKB_SUBUNIT = re.compile(r"p(?:65|50|52|105|100)")


def nfkb_complex_genes(name):
    """If `name` denotes the generic NF-kB complex, return (NFKB1, RELA); else None.

    Works on the separator-free `_complex_key` form so 'NF-kB', 'NFkB', 'NF-KB',
    'NF-kappaB', 'NF-kB protein/gene', 'p-NF-kB', 'NFkB complex' and 'NF-kB target
    genes' all qualify, while a specific subunit gene (the stem followed by a digit
    -> NFKB1/2, or by 'i' -> the NFKBI* inhibitors) or a named subunit (p65/p50/...)
    is rejected so it still links to its own gene."""
    k = _complex_key(name)
    m = _NFKB_STEM.search(k)
    if not m:
        return None
    after = k[m.end():]
    if after[:1].isdigit():        # nfkb1, nfkb2 (and shnfkb2, nfkb12 ...): specific genes
        return None
    if after[:1] == "i":           # nfkbia / nfkbib / nfkbie / nfkbiz / nfkbid: IkB inhibitors
        return None
    if _NFKB_SUBUNIT.search(k):     # NF-kB p65 / p50 / ...: a named single subunit
        return None
    return NFKB_GENES


def complex_pass(rows, docs, unmatched_tsv, label):
    """Map curated complex/heterodimer names to HGNC subunit gene(s); one subunit
    -> matched, the full complex -> ambiguous subunit set."""
    sym2id = {d["symbol"]: (d.get("hgnc_id") or d["symbol"])
              for d in docs if d.get("symbol")}
    matched_lib, ambiguous_lib, unmatched = {}, {}, []
    field_counts = defaultdict(int)
    occ = {"total": 0, "matched": 0, "ambiguous": 0, "unmatched": 0}
    sample = []
    for ln in rows:
        name, rest = ln.split("\t", 1)
        n_occ = int(rest.split("\t", 1)[0])
        occ["total"] += n_occ
        key = _complex_key(name)
        genes = COMPLEX_MAP.get(key)
        field, mode = "complex", "complex→subunits"
        if genes is None and key in FAMILY_MAP:
            genes, field, mode = FAMILY_MAP[key], "family", "receptor family (ambiguous)"
        elif genes is None and key in ONCO_AMBIG_MAP:
            genes, field, mode = ONCO_AMBIG_MAP[key], "oncogene", "c-oncogene→genes (ambiguous)"
        elif genes is None:
            nfkb = nfkb_complex_genes(name)         # generic NF-kB -> NFKB1 + RELA dimer
            if nfkb:
                genes, field, mode = nfkb, "complex", "complex→subunits"
        if not genes:
            unmatched.append(ln)
            occ["unmatched"] += n_occ
        elif len(genes) == 1:
            sym = genes[0]
            gid = sym2id.get(sym, sym)
            matched_lib[name] = {"hgnc_symbol": sym, "occurrences": n_occ,
                                 "match_field": field, "match_mode": mode}
            field_counts[field] += 1
            occ["matched"] += n_occ
            if len(sample) < 10:
                sample.append(f"{name}\t{rest}\t{field}\t{gid}\t{sym}")
        else:
            syms = sorted(genes)
            ambiguous_lib[name] = {"hgnc_symbol": syms, "occurrences": n_occ,
                                   "match_field": field,
                                   "match_mode": mode,
                                   "n_genes": len(syms),
                                   "hgnc_ids": [sym2id.get(s, s) for s in syms]}
            occ["ambiguous"] += n_occ
    if unmatched_tsv is not None:
        _write_block(unmatched_tsv, TSV_HEADER, unmatched)
    return {"label": label, "fields": ["complex"],
            "n_total": len(rows), "n_matched": len(matched_lib),
            "n_ambiguous": len(ambiguous_lib), "n_unmatched": len(unmatched),
            "by_field": dict(field_counts),
            "by_mode": {"complex→subunits": len(matched_lib)},
            "occ": occ, "sample": sample, "unmatched_rows": unmatched,
            "matched_lib": matched_lib, "ambiguous_lib": ambiguous_lib}


# =======================================================================
# manual correction map -- last-resort override for entities the cascade links
# wrongly (typically an alias collision the approved-symbol guard can't catch,
# because only one -- wrong -- gene matches). entity -> approved symbol.
# =======================================================================
CORRECTIONS = {
    "pYAP": "YAP1",              # phospho-YAP; auto-links to YY1AP1 via the alias 'YAP'
    "nuclear YAP": "YAP1",       # nuclear YAP; same 'YAP' alias collision -> YY1AP1
    "YAP target genes": "YAP1",  # YAP target-set regulator; same 'YAP' collision
}


# =======================================================================
# NON-GENE triage -- high-frequency final-unmatched values that were reviewed
# and found NOT to be a single human HGNC gene, each with the reason. This is a
# DOCUMENTATION record only (it does not change matching): the report lists these
# so prominent leftovers read as deliberately unmapped, not missed. Values that
# ARE a real gene HGNC merely lacks a surface for (APNG=MPG, IF1=ATP5IF1,
# Presenilin1=PSEN1, LC3B=MAP1LC3B, actin->ACTB) are intentionally EXCLUDED -- they
# are recall gaps, not non-genes. Every entry below is verified absent from all
# HGNC symbol/alias/prev/name fields and present in the final unmatched set.
# =======================================================================
NON_GENE = {
    # -- reporter / tag / foreign proteins (not human genes) --
    "GFP": "green fluorescent protein (Aequorea reporter), not a human gene",
    "EGFP": "enhanced GFP reporter",
    "eGFP": "enhanced GFP reporter",
    "YFP": "yellow fluorescent protein reporter",
    "luciferase": "firefly/Renilla luciferase reporter enzyme, not a human gene",
    "luciferase reporter": "luciferase reporter construct",
    "FLAG": "FLAG epitope tag, not a gene",
    "Cas9": "Streptococcus pyogenes CRISPR nuclease, not a human gene",
    # -- gene-editing / knock-down / assay reagents & techniques --
    "sgRNA": "single-guide RNA reagent (CRISPR), not a gene",
    "sgRNAs": "single-guide RNA reagents (CRISPR), not a gene",
    "shRNA": "short-hairpin RNA reagent, not a gene",
    "sh": "short-hairpin (shRNA) prefix fragment, not a gene",
    "shCtrl": "non-targeting shRNA control reagent, not a gene",
    "CRISPR": "gene-editing technique, not a gene",
    "TUNEL": "TdT dUTP nick-end labeling apoptosis assay, not a gene",
    "CCK8": "Cell Counting Kit-8 viability assay reagent, not the CCK gene",
    "CCK-8": "Cell Counting Kit-8 viability assay reagent, not the CCK gene",
    # -- drugs / small molecules --
    "TMZ": "temozolomide (alkylating chemotherapeutic), not a gene",
    "RSL3": "ferroptosis inducer (GPX4 inhibitor), a small molecule",
    "erastin": "ferroptosis inducer (system xc- inhibitor), a small molecule",
    "JQ1": "BET-bromodomain inhibitor, a small molecule",
    "Torin1": "ATP-competitive mTOR inhibitor, a small molecule",
    "MEKi": "MEK inhibitor (drug class), not a gene",
    # -- engineered / variant constructs --
    "EGFRvIII": "EGFR variant III in-frame deletion mutant, not a distinct HGNC gene",
    # -- circular RNAs (named after a host gene; not HGNC entries) --
    "circKIF4A": "circular RNA (host gene KIF4A), not an HGNC gene symbol",
    "circRNF10": "circular RNA (host gene RNF10), not an HGNC gene symbol",
    "circRPPH1": "circular RNA (host gene RPPH1), not an HGNC gene symbol",
    # -- biomaterials / non-protein antigens --
    "ctDNA": "circulating tumor DNA (a biomaterial), not a gene",
    "GD2": "disialoganglioside GD2 (a glycolipid antigen), not a gene",
    # -- mouse-only symbol --
    "Trp53": "mouse Tp53 symbol (human ortholog is TP53)",
    # -- tokenization fragments / bare tokens --
    "-": "bare hyphen, a tokenization artifact",
    "AS1": "orphan '-AS1' antisense-RNA suffix fragment, not a gene on its own",
    # -- gene families / pathways / complexes (not a single gene) --
    "MAPK": "MAP-kinase family/pathway name (MAPK1/3/...), not a single gene",
    "MEK": "MEK family (MAP2K1/MAP2K2), not a single gene",
    "ERK1 / 2": "coordinate ERK1/2 mention (MAPK3 + MAPK1), not a single gene",
    "Wnt": "WNT ligand family/pathway, not a single gene",
    "WNT": "WNT ligand family/pathway, not a single gene",
    "Notch": "NOTCH receptor family (NOTCH1-4), not a single gene",
    "Hippo": "Hippo signaling pathway, not a single gene",
    "Hh": "Hedgehog pathway (SHH/IHH/DHH), not a single gene",
    "JAK": "Janus-kinase family (JAK1/2/3, TYK2), not a single gene",
    "RAS": "RAS family (HRAS/KRAS/NRAS), not a single gene",
    "Ras": "RAS family (HRAS/KRAS/NRAS), not a single gene",
    "IDH": "isocitrate-dehydrogenase family (IDH1/2/3), not a single gene",
    "PDGF": "PDGF ligand family (PDGFA-D), not a single gene",
    "BMP": "bone-morphogenetic-protein ligand family, not a single gene",
    "MMP": "matrix-metalloproteinase family, not a single gene",
    "caspase": "caspase family (CASP1-14), not a single gene",
    "caspase 3 / 7": "coordinate caspase-3/7 mention (CASP3 + CASP7), not a single gene",
    "histone": "histone family (many genes), not a single gene",
    "collagen": "collagen family (COL* genes), not a single gene",
    "integrin": "integrin family (a heterodimer of ITGA*/ITGB* chains), not a single gene",
    "interferon": "interferon family (IFNA*/IFNB1/IFNG/...), not a single gene",
    "PRC2": "Polycomb Repressive Complex 2 (EZH2/EED/SUZ12/...), a complex",
    "mTORC1": "mTOR complex 1 (MTOR/RPTOR/...), a complex",
    "mTORC2": "mTOR complex 2 (MTOR/RICTOR/...), a complex",
    "CK2": "casein kinase 2 holoenzyme (CSNK2A1/A2/B), a complex",
    "PKA": "protein kinase A holoenzyme (PRKACA/B/G + regulatory), a complex",
    "HSP70": "HSP70 chaperone family (HSPA1A/...), not a single gene",
    "SOD": "superoxide-dismutase family (SOD1/2/3), not a single gene",
    "LDH": "lactate-dehydrogenase family (LDHA/B/C) / serum marker, not a single gene",
}


def apply_corrections(roman_lib, roman_amb, unmatched_rows, approved, name_occ):
    """Force CORRECTIONS entities onto their specified gene, overriding whatever
    the cascade did. The original match_mode is preserved (so a semantic tag like
    'regulatory target-set (not the gene)' survives the correction); the override
    is flagged via match_field='manual'. Returns the (possibly reduced) unmatched
    rows."""
    drop = set()
    for entity, sym in CORRECTIONS.items():
        if sym not in approved:           # guard against typos in the map
            continue
        if entity in roman_lib:
            occ, mode = roman_lib[entity]["occurrences"], roman_lib[entity]["match_mode"]
        elif entity in roman_amb:
            occ, mode = roman_amb[entity]["occurrences"], roman_amb[entity]["match_mode"]
        else:
            occ, mode = name_occ.get(entity), "manual correction"
        if occ is None:                   # entity not in this corpus -> skip
            continue
        roman_amb.pop(entity, None)
        roman_lib[entity] = {"hgnc_symbol": sym, "occurrences": occ,
                             "match_field": "manual", "match_mode": mode}
        drop.add(entity)
    if not drop:
        return unmatched_rows
    return [r for r in unmatched_rows if r.split("\t", 1)[0] not in drop]


# =======================================================================
# driver
# =======================================================================
def main():
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
    docs = json.loads(HGNC_JSON.read_text(encoding="utf-8"))["response"]["docs"]
    match_pass.docs = docs

    roman_rows, s1 = stage1()

    # Intermediate stages keep their matched/unmatched buckets in memory (the
    # leftovers cascade onward); only the final stage persists unmatched.tsv. The
    # matched data lives entirely in the consolidated roman.json below. Any stale
    # per-step TSVs from older runs are deleted so they cannot mislead.
    for stale in ("types_1_matched.tsv", "types_2_matched.tsv", "names_matched.tsv",
                  "unmatched_type_1.tsv", "unmatched_type_2.tsv"):
        (SUMMARY_DIR / stale).unlink(missing_ok=True)

    SYM = ["symbol", "alias_symbol", "prev_symbol"]
    ALL_FIELDS = SYM + NAME_FIELDS
    # COMPLEX runs FIRST so the curated complex/heterodimer map is authoritative
    # (e.g. IL-23 -> IL12B/IL23A, not the partial alias match IL23A; IL-27 -> the
    # right subunits, not stale aliases). Its leftovers feed the normal cascade.
    cx = complex_pass(roman_rows, docs, None,
                      "COMPLEX - heterodimer/complex -> subunit gene(s) [priority]")
    t1 = match_pass(cx["unmatched_rows"], SYM, identity, None,
                    "TYPE 1 - symbol, case-sensitive")
    t2 = match_pass(t1["unmatched_rows"], SYM, hsfold, None,
                    "TYPE 2 - symbol, case- & hyphen/space-insensitive")
    nm = match_pass(t2["unmatched_rows"], NAME_FIELDS, hsfold, None,
                    "NAME - descriptive name, case- & hyphen/space-insensitive")
    # PHOSPHO: strip a leading phospho prefix ('p-' or 'phospho-') and re-match
    # the leftovers against all six symbol+name fields.
    ph = match_pass(nm["unmatched_rows"], ALL_FIELDS, pstrip_hsfold, None,
                    "PHOSPHO - leading p-/phospho- stripped, symbol+name fields",
                    mode_override="phospho stripped", index_keyfn=hsfold)
    # ANTI: strip a leading 'anti-' prefix and re-match. 'anti-X' is the antibody
    # AGAINST gene X (a target, not the entity itself) -- tagged 'anti- stripped'.
    # Final exact stage.
    an = match_pass(ph["unmatched_rows"], ALL_FIELDS, anti_hsfold, None,
                    "ANTI - leading 'anti-' / trailing ' antibody' stripped, symbol+name fields",
                    mode_override="anti- stripped", index_keyfn=hsfold)
    # CLEAVED: strip a leading 'cleaved-' proteolytic marker; the cleaved product
    # is still the gene's protein, so re-match against all six fields.
    cl = match_pass(an["unmatched_rows"], ALL_FIELDS, cleaved_strip_hsfold, None,
                    "CLEAVED - leading 'cleaved-' stripped, symbol+name fields",
                    mode_override="cleaved- stripped", index_keyfn=hsfold)
    # PROTEIN: strip a ' protein'/'-protein' descriptor and re-match. These name
    # the gene's own product, so they fold into normal gene-identity matches.
    pr = match_pass(cl["unmatched_rows"], ALL_FIELDS, protein_strip_hsfold, None,
                    "PROTEIN - ' protein'/'-protein' stripped, symbol+name fields",
                    mode_override="protein stripped", index_keyfn=hsfold)
    # TRANSCRIPT: strip ' transcript(s)' / ' mRNA(s)' descriptors and re-match.
    # Names the gene's RNA product -> folds onto the gene.
    tr = match_pass(pr["unmatched_rows"], ALL_FIELDS, transcript_strip_hsfold, None,
                    "TRANSCRIPT - ' transcript(s)'/' mRNA(s)' stripped, symbol+name fields",
                    mode_override="transcript/mRNA stripped", index_keyfn=hsfold)
    # VARIANT: strip ' variant(s)'/' mutant(s)' descriptors (either side) and
    # re-match -- an altered form of the gene folds onto it.
    vr = match_pass(tr["unmatched_rows"], ALL_FIELDS, variant_strip_hsfold, None,
                    "VARIANT - ' variant(s)'/' mutant(s)'/' mutation(s)' stripped, symbol+name fields",
                    mode_override="variant/mutant stripped", index_keyfn=hsfold)
    # MISSENSE: <approved symbol> + a valid-AA missense [AA]<pos>[AA] folds onto
    # the gene (IDH1R132H -> IDH1). Approved 'symbol' field only, case-sensitive
    # (type-1); small-molecule compounds excluded. Final exact stage.
    ms = match_pass(vr["unmatched_rows"], ["symbol"], missense_strip, None,
                    "MISSENSE - <approved symbol>+[AA]<pos>[AA] stripped, symbol field",
                    mode_override="missense stripped", index_keyfn=identity)
    # MISSENSE-SP: the space-delimited form '<approved symbol> [AA]<pos>[AA]'
    # (TP53 R248L -> TP53). Same approved-symbol + valid-AA guards as MISSENSE.
    ms2 = match_pass(ms["unmatched_rows"], ["symbol"], missense_sp_strip, None,
                     "MISSENSE-SP - <approved symbol> ' ' [AA]<pos>[AA] stripped, symbol field",
                     mode_override="missense stripped", index_keyfn=identity)
    # INHIBITOR: strip a terminal '-i' inhibitor marker (PARPi -> PARP) and match
    # case-sensitively against the symbol fields. The '-i' is a target/inhibitor
    # (like 'anti-'), tagged distinctly. Final exact stage.
    ih = match_pass(ms2["unmatched_rows"], SYM, inhibitor_strip, None,
                    "INHIBITOR - terminal '-i' stripped, symbol fields",
                    mode_override="inhibitor-i stripped", index_keyfn=identity)
    # MIRNA: map mature microRNA names to HGNC precursor gene(s); multi-locus
    # mature miRs go to the ambiguous library.
    mr = mirna_pass(ih["unmatched_rows"], docs, None,
                    "MIRNA - mature miR name -> precursor gene(s)")
    # PROMOTER: strip a trailing ' promoter'/' promoter region'/' gene promoter'
    # regulatory descriptor and re-match the symbol fields. Final stage.
    pm = match_pass(mr["unmatched_rows"], SYM, promoter_strip_hsfold, None,
                    "PROMOTER - ' promoter' descriptor stripped, symbol fields",
                    mode_override="promoter stripped", index_keyfn=hsfold)
    # DEHYPHEN: <UPPER>-<digits> tokens whose hyphen is mere punctuation -- remove
    # it and match the symbol fields (MMP-9 -> MMP9). Drugs/cell-lines excluded.
    # Final stage.
    dh = match_pass(pm["unmatched_rows"], SYM, dehyphen_key, None,
                    "DEHYPHEN - <UPPER>-<digits> hyphen removed, symbol fields",
                    mode_override="hyphen-removed", index_keyfn=str.casefold,
                    approved_guard=True)
    # P-PREFIX: strip a bare leading 'p' before an upper-case letter (no-delimiter
    # phospho/protein form, pAKT -> AKT) and match the symbol fields. Constructs/
    # vectors excluded; approved-symbol guard for disambiguation. Final stage.
    pp = match_pass(dh["unmatched_rows"], SYM, pprefix_strip_hsfold, None,
                    "P-PREFIX - bare leading 'p' stripped, symbol fields",
                    mode_override="p-prefix stripped", index_keyfn=hsfold,
                    approved_guard=True)
    # WT: strip a leading 'WT ' (wild-type) prefix and match the symbol fields
    # (WT IDH1 -> IDH1). Final stage.
    wt = match_pass(pp["unmatched_rows"], SYM, wt_strip_hsfold, None,
                    "WT - leading 'WT ' stripped, symbol fields",
                    mode_override="WT stripped", index_keyfn=hsfold)
    # NUCLEAR: strip a leading 'nuclear ' localization prefix and match the
    # symbol fields (nuclear YAP -> YAP).
    nu = match_pass(wt["unmatched_rows"], SYM, nuclear_strip_hsfold, None,
                    "NUCLEAR - leading 'nuclear ' stripped, symbol fields",
                    mode_override="nuclear stripped", index_keyfn=hsfold)
    # TARGET-GENES: 'X target genes' -> the regulator X, but tagged distinctly --
    # it denotes X's downstream targets, NOT X, so it is kept separable.
    tg = match_pass(nu["unmatched_rows"], SYM, target_genes_strip_hsfold, None,
                    "TARGET-GENES - ' target genes' stripped -> regulator, symbol fields",
                    mode_override="regulatory target-set (not the gene)", index_keyfn=hsfold)
    # HISTONE: map histone marks/mutants/variants to HGNC histone gene(s); most
    # are families -> ambiguous, a few variants -> a single gene.
    hs = histone_pass(tg["unmatched_rows"], docs, None,
                      "HISTONE - marks/mutants/variants -> histone gene(s)")
    # GREEK-LETTER: expand a lone Roman letter standing for a Greek symbol to its
    # spelled-out word and re-match all six fields (NF-kB -> NF-kappaB -> NFKB1).
    # c-X proto-oncogene spellings (c-Kit/c-Src/c-Abl/.../c-Met/c-erbB-N) and the
    # HGF/MET receptor -> their HGNC gene; numbered/bare ambiguous forms handled too.
    co = match_pass(hs["unmatched_rows"], ["symbol"], conco_hsfold, None,
                    "C-ONCOGENE - c-Kit/c-Src/c-Met/c-erbB-N/... spelling -> gene (symbol)",
                    mode_override="c-oncogene→gene alias", index_keyfn=hsfold,
                    approved_guard=True)
    gx = match_pass(co["unmatched_rows"], ALL_FIELDS, greek_letter_hsfold, None,
                    "GREEK-LETTER - lone Roman letter -> spelled-out Greek word, symbol+name fields",
                    mode_override="greek-letter spelled out", index_keyfn=hsfold,
                    approved_guard=True)
    # second pass on the leftovers: also expand a leading abbreviation to its full
    # word so the HGNC *name* can match (IFN-g -> interferon gamma -> IFNG). Runs
    # after the plain pass so it never clobbers a shorter alias link (TNF-a).
    gx2 = match_pass(gx["unmatched_rows"], ALL_FIELDS, greek_abbrev_hsfold, None,
                     "GREEK-LETTER (abbrev) - leading abbreviation -> full word + Greek letter, symbol+name fields",
                     mode_override="greek-letter spelled out", index_keyfn=hsfold,
                     approved_guard=True)
    # DELSEP: a final separator-deletion mop-up. Drop ALL hyphens and spaces on
    # both sides so 'PDL1' <-> 'PD-L1' <-> 'PD L1' unify (the no-separator surfaces
    # hsfold and the narrow DEHYPHEN both miss). Split for precision: the approved
    # `symbol` field at any length (authoritative), then the noisier alias/prev
    # fields with a >=4-char floor. Runs last so every precise stage gets first claim.
    ds = match_pass(gx2["unmatched_rows"], ["symbol"], delsep_key_sym, None,
                    "DELSEP - all separators removed, symbol field",
                    mode_override="separator-removed", index_keyfn=delsep)
    ds2 = match_pass(ds["unmatched_rows"], ["alias_symbol", "prev_symbol"],
                     delsep_key_alias, None,
                     "DELSEP (alias) - all separators removed, alias/prev fields (>=4-char key)",
                     mode_override="separator-removed", index_keyfn=delsep)

    # consolidate the matched / ambiguous JSON libraries across all passes
    # (COMPLEX first, then the cascade). Entities are unique across passes, so the
    # merges are collision-free; match_field/match_mode record how each linked.
    roman_lib, roman_amb = {}, {}
    for st in (cx, t1, t2, nm, ph, an, cl, pr, tr, vr, ms, ms2, ih, mr, pm, dh, pp, wt, nu, tg, hs, co, gx, gx2, ds, ds2):
        roman_lib.update(st["matched_lib"])
        roman_amb.update(st["ambiguous_lib"])

    # manual correction map -- override mis-links, then persist the final leftovers
    approved_symbols = {d.get("symbol") for d in docs}
    name_occ = {r.split("\t", 1)[0]: int(r.split("\t")[1]) for r in roman_rows}
    final_unmatched = apply_corrections(roman_lib, roman_amb,
                                        ds2["unmatched_rows"], approved_symbols, name_occ)
    _write_block(SUMMARY_DIR / "unmatched.tsv", TSV_HEADER, final_unmatched)
    (SUMMARY_DIR / "roman.json").write_text(
        json.dumps(roman_lib, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (SUMMARY_DIR / "roman_ambiguous.json").write_text(
        json.dumps(roman_amb, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    cos = stage4_cosine(final_unmatched, docs)

    # NON-GENE triage: of the documented non-gene values, list those that really
    # are in the final unmatched set (verified, not contradicting any link), with
    # occurrences, sorted by frequency -- a reviewed-leftovers record for the report.
    final_unmatched_vals = {r.split("\t", 1)[0] for r in final_unmatched}
    nongene_rows = sorted(
        ((v, name_occ.get(v, 0), reason) for v, reason in NON_GENE.items()
         if v in final_unmatched_vals),
        key=lambda r: -r[1])

    roman_occ = sum(int(r.split("\t")[1]) for r in roman_rows)
    passes = [cx, t1, t2, nm, ph, an, cl, pr, tr, vr, ms, ms2, ih, mr, pm, dh, pp, wt, nu, tg, hs, co, gx, gx2, ds, ds2]
    write_html(s1, passes, cos, roman_occ, nongene_rows)
    report(s1, passes, cos, roman_occ, nongene_rows)
    print(f"consolidated : roman.json ({len(roman_lib):,} matched), "
          f"roman_ambiguous.json ({len(roman_amb):,} ambiguous); "
          f"final unmatched -> unmatched.tsv")


def report(s1, passes, cos, roman_occ, nongene_rows=()):
    print(f"STAGE 1  files={s1['n_files']:,}  occ={s1['occurrences']:,}  "
          f"unique {s1['unique_before']:,}->{s1['unique_after']:,}  "
          f"greek={s1['greek']:,} roman={s1['roman']:,}")
    linked = 0
    for st in passes:
        o = st["occ"]
        linked += o["matched"] + o["ambiguous"]
        print(f"{st['label']}")
        print(f"   in={st['n_total']:,}  matched={st['n_matched']:,}  "
              f"ambiguous={st['n_ambiguous']:,}  unmatched={st['n_unmatched']:,}")
    print(f"STAGE 4  cosine pairs={len(cos['rows']):,} "
          f"(over {cos['n_inputs_mw']:,} multi-word leftovers)")
    print(f"cumulative exact-linked occ (matched+ambiguous): {linked:,} "
          f"({100*linked/roman_occ:.1f}% of {roman_occ:,} Roman occ)")
    if nongene_rows:
        ng_occ = sum(o for _, o, _ in nongene_rows)
        print(f"NON-GENE triage : {len(nongene_rows)} reviewed non-gene leftover(s) "
              f"documented ({ng_occ:,} occ) -- e.g. "
              f"{', '.join(v for v, _, _ in nongene_rows[:5])}")
    print(f"written : {HTML_OUT}")


def write_html(s1, passes, cos, roman_occ, nongene_rows=()):
    def esc(s):
        return html.escape(s)

    linked = sum(st["occ"]["matched"] + st["occ"]["ambiguous"] for st in passes)

    nongene_section = ""
    if nongene_rows:
        ng_occ = sum(o for _, o, _ in nongene_rows)
        ng_body = "".join(
            f"<tr><td><code>{esc(v)}</code></td><td class='num'>{o:,}</td>"
            f"<td>{esc(reason)}</td></tr>" for v, o, reason in nongene_rows)
        nongene_section = f"""
<h2>Reviewed non-gene leftovers</h2>
<p>Prominent values in the final <code>unmatched.tsv</code> that were reviewed and
found <strong>not to be a single human HGNC gene</strong> &mdash; reporter / tag
proteins (GFP, luciferase, Cas9), gene-editing &amp; assay reagents (sgRNA, shRNA,
TUNEL), drugs (TMZ, JQ1, erastin), engineered variants (EGFRvIII), circular RNAs,
non-protein antigens (GD2), a mouse-only symbol (Trp53), tokenization fragments,
and bare gene-family / pathway / complex names (MAPK, Wnt, RAS, mTORC1, &hellip;).
Documented here so these leftovers read as <em>deliberately</em> unmapped, not
missed. Real genes that HGNC merely lacks a surface for (e.g. APNG=MPG,
Presenilin1=PSEN1) are deliberately excluded &mdash; they are recall gaps, not
non-genes.</p>
<table>
  <tr><th>value</th><th class="num">occ</th><th>why not a single gene</th></tr>
  {ng_body}
  <tr><td><strong>{len(nongene_rows)} values</strong></td>
      <td class="num"><strong>{ng_occ:,}</strong></td><td></td></tr>
</table>"""

    pass_sections = ""
    for st in passes:
        o = st["occ"]
        by_field = " &nbsp; ".join(
            f"<code>{esc(f)}</code>&nbsp;{st['by_field'].get(f, 0):,}"
            for f in st["fields"])
        by_mode = " &nbsp; ".join(
            f"<code>{esc(m)}</code>&nbsp;{c:,}"
            for m, c in sorted(st["by_mode"].items(), key=lambda kv: -kv[1])) or "&mdash;"
        sample = "".join(
            "<tr>" + "".join(f"<td>{esc(c)}</td>" for c in r.split("\t")[:6])
            + "</tr>" for r in st["sample"])
        pass_sections += f"""
<h3>{esc(st['label'])}</h3>
<table>
  <tr><th>Bucket</th><th class="num">Rows</th><th class="num">Occurrences</th></tr>
  <tr><td><strong>Matched (1 gene)</strong></td>
      <td class="num"><strong>{st['n_matched']:,}</strong></td>
      <td class="num">{o['matched']:,}</td></tr>
  <tr><td>Ambiguous (&ge;2 genes)</td><td class="num">{st['n_ambiguous']:,}</td>
      <td class="num">{o['ambiguous']:,}</td></tr>
  <tr><td>Unmatched</td><td class="num">{st['n_unmatched']:,}</td>
      <td class="num">{o['unmatched']:,}</td></tr>
  <tr><td>Pass input</td><td class="num">{st['n_total']:,}</td>
      <td class="num">{o['total']:,}</td></tr>
</table>
<p>Matched by field: {by_field}<br>Matched by mode: {by_mode}<br>
Contributes {st['n_matched']:,} matched / {st['n_ambiguous']:,} ambiguous
entries to the consolidated <code>roman.json</code> /
<code>roman_ambiguous.json</code>.</p>
<details><summary>sample matches</summary><table>
  <tr><th>value</th><th>occ</th><th>n_forms</th><th>field</th><th>hgnc_id</th><th>symbol</th></tr>
  {sample}</table></details>
"""

    cos_rows = "".join(
        f"<tr><td>{esc(r['input'])}</td><td>{esc(r['hgnc'])}</td>"
        f"<td class='num'>{r['cos']:.4f}</td><td><code>{esc(r['field'])}</code></td>"
        f"<td>{esc(r['ids'])}</td><td>{esc(r['syms'])}</td>"
        f"<td class='num'>{r['occ']:,}</td></tr>" for r in cos["rows"])

    doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>roman &mdash; unified normalization &amp; HGNC linkage</title>
<style>
  body {{ font: 15px/1.5 -apple-system, Segoe UI, Roboto, sans-serif;
          max-width: 940px; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; }}
  h1 {{ font-size: 1.5rem; }}
  h2 {{ font-size: 1.2rem; margin-top: 2rem; border-bottom: 1px solid #ddd;
        padding-bottom: .3rem; }}
  h3 {{ font-size: 1.02rem; margin-top: 1.5rem; }}
  code {{ background: #f3f3f3; padding: 1px 5px; border-radius: 3px; font-size: .9em; }}
  table {{ border-collapse: collapse; margin: .7rem 0; font-size: .9em; }}
  th, td {{ border: 1px solid #ccc; padding: .35rem .7rem; text-align: left; }}
  th {{ background: #f7f7f7; }}
  .num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  details {{ margin: .5rem 0; }} summary {{ cursor: pointer; color: #357; }}
  ol.flow > li {{ margin: .4rem 0; }}
</style>
</head>
<body>
<h1>roman &mdash; unified normalization &amp; HGNC linkage</h1>
<p>One pipeline unifying the strategies of <code>roman.py</code>,
<code>symbols_match.py</code>, <code>names_match.py</code> and
<code>cosine_similarity.py</code>: from raw BioBERT GENETIC spans to
gene-level HGNC links.</p>

<h2>Unified strategy</h2>
<p>The pipeline is a <strong>cascade of progressively looser matchers</strong>.
Every stage hands its <em>unmatched</em> leftovers to the next, so each entity is
linked by the strictest rule that succeeds, and effort concentrates on what
remains. Throughout, an entity is only &ldquo;linked&rdquo; when it resolves to a
<strong>single</strong> gene; strings that map to two or more genes are held
aside as <em>ambiguous</em> rather than guessed.</p>
<ol class="flow">
  <li><strong>STAGE 1 &mdash; surface normalization &amp; routing.</strong> Parse
      every <code>sentences/*.json</code>, collect GENETIC spans, and unify dash
      noise: all dash/hyphen variants &rarr; ASCII <code>-</code>, whitespace
      around <code>-</code> removed. Identical normalized forms are aggregated
      (with occurrence counts) into <code>clean_genetic_ne.tsv</code>. Entities
      are then routed: those carrying a <strong>Greek letter</strong> (a Unicode
      symbol, or a spelled-out letter name as a whole word &mdash; with all-caps
      and short ambiguous tokens like <code>Xi/Pi</code> excluded) go to the
      Greek file; the rest, the <strong>Roman</strong> set, feed the linker.</li>
  <li><strong>Exact single-gene linkage (a four-pass cascade).</strong> One
      shared routine matches a value to HGNC by exact equality of a
      <em>transformed key</em>, indexing the relevant fields and bucketing by
      how many distinct genes the key hits (1 = matched, &ge;2 = ambiguous,
      0 = unmatched). It is applied as a sequence of passes with progressively
      looser keys / different fields (plus one bespoke microRNA mapping):
    <ul>
      <li><strong>TYPE&nbsp;1</strong> &mdash; <em>symbol</em> fields
          (<code>symbol</code>, <code>alias_symbol</code>,
          <code>prev_symbol</code>), <strong>case-sensitive</strong> identity
          key. The strict benchmark.</li>
      <li><strong>TYPE&nbsp;2</strong> &mdash; same symbol fields, key folds
          <strong>case</strong> and <strong>hyphen&harr;whitespace</strong>
          (<code>casefold</code> then rewrite <code>-</code>&rarr;space; a
          length-preserving 1:1 fold). Recovers casing and hyphen/space
          variants (<code>Akt</code>&rarr;AKT, <code>PD L1</code>&rarr;PD-L1).</li>
      <li><strong>NAME</strong> &mdash; the <em>descriptive name</em> fields
          (<code>name</code>, <code>alias_name</code>, <code>prev_name</code>)
          with the same case + hyphen/space key. Catches entities written out as
          a full gene name (<code>caspase-3</code>&rarr;caspase&nbsp;3).</li>
      <li><strong>PHOSPHO</strong> &mdash; strips a leading phospho prefix
          (<code>p-</code> or <code>phospho-</code>) from the remaining
          leftovers, then re-matches the stripped form against <em>all six</em>
          symbol+name fields with the case/hyphen-folded key. Recovers
          phospho-protein mentions (<code>p-AKT</code>&rarr;AKT1,
          <code>phospho-STAT3</code>&rarr;STAT3); recorded with
          <code>match_mode</code>&nbsp;=&nbsp;<code>phospho stripped</code>.</li>
      <li><strong>ANTI</strong> &mdash; strips a leading <code>anti-</code>
          prefix, and a trailing <code>&nbsp;antibody</code> too when
          <code>anti-</code> is present (both come off together for
          <code>anti-X antibody</code>), then re-matches against all six fields.
          Recovers antibody / antagonist mentions (<code>anti-PD-1</code>&rarr;PDCD1,
          <code>anti-EGFR antibody</code>&rarr;EGFR); recorded with
          <code>match_mode</code>&nbsp;=&nbsp;<code>anti- stripped</code>.
          Note <code>anti-X</code> denotes a reagent <em>against</em> gene X
          &mdash; a target, not X's own product &mdash; so this mode is kept
          distinct.</li>
      <li><strong>CLEAVED</strong> &mdash; strips a leading <code>cleaved-</code>
          proteolytic-processing marker and re-matches against all six fields;
          the cleaved product is still the gene's protein, so it folds onto the
          gene (<code>cleaved-PARP</code>&rarr;PARP, <code>cleaved-caspase 3</code>
          &rarr;CASP3 via the name fields); recorded with
          <code>match_mode</code>&nbsp;=&nbsp;<code>cleaved- stripped</code>.</li>
      <li><strong>PROTEIN</strong> &mdash; strips a <code>&nbsp;protein</code> /
          <code>-protein</code> descriptor and re-matches against all six fields.
          Recovers &ldquo;<em>X protein</em>&rdquo; mentions
          (<code>p53 protein</code>&rarr;TP53, <code>MGMT protein</code>&rarr;MGMT);
          recorded with <code>match_mode</code>&nbsp;=&nbsp;<code>protein stripped</code>.
          These name the gene's own product, so they fold in with the
          gene-identity matches.</li>
      <li><strong>TRANSCRIPT</strong> &mdash; strips RNA-product descriptors
          <code>&nbsp;transcript(s)</code> and <code>&nbsp;mRNA(s)</code>
          (collapsing combos like <code>EGFR mRNA transcript</code>) and
          re-matches against all six fields (<code>EGFR mRNA</code>&rarr;EGFR,
          <code>MYC transcripts</code>&rarr;MYC); recorded with
          <code>match_mode</code>&nbsp;=&nbsp;<code>transcript/mRNA stripped</code>.
          A trailing word boundary keeps <code>transcriptase</code> /
          <code>transcriptional</code> intact.</li>
      <li><strong>VARIANT</strong> &mdash; strips <code>&nbsp;variant(s)</code> /
          <code>&nbsp;mutant(s)</code> / <code>&nbsp;mutation(s)</code> descriptors
          on <em>either side</em> (suffix <code>EGFR mutant</code>,
          <code>TP53 mutation</code>, or prefix <code>mutant p53</code>) and
          re-matches against all six fields (&rarr;EGFR, &rarr;TP53); recorded
          with <code>match_mode</code>&nbsp;=&nbsp;<code>variant/mutant stripped</code>.
          An altered form of a gene folds onto that gene.</li>
      <li><strong>MISSENSE</strong> &mdash; recognises the type-1 missense
          convention <em>&lt;approved symbol&gt; + [AA]&lt;pos&gt;[AA]</em>
          (<code>IDH1R132H</code>&rarr;IDH1, <code>BRAFV600E</code>&rarr;BRAF) and
          links it to that gene. Requires an <strong>approved</strong>
          <code>symbol</code> (case-sensitive, not an alias) and <strong>valid
          amino-acid</strong> residues; recorded as
          <code>missense stripped</code>. Curated small-molecule compounds
          (e.g. <code>KYA1797K</code>) are excluded. The space-delimited form
          (<code>TP53 R248L</code>&rarr;TP53) is handled by a sibling pass with
          the same guards and the same <code>missense stripped</code> tag.</li>
      <li><strong>INHIBITOR</strong> &mdash; strips a terminal inhibitor
          <code>i</code> (<code>PARPi</code>&rarr;PARP, <code>BRAFi</code>&rarr;BRAF)
          and matches case-sensitively against the symbol fields; recorded as
          <code>inhibitor-i stripped</code>. Like <code>anti-</code>, the
          &ldquo;-i&rdquo; denotes a drug <em>against</em> the gene product &mdash;
          a target, not the entity itself &mdash; so the mode is kept distinct.</li>
      <li><strong>MIRNA</strong> (bespoke) &mdash; maps mature microRNA names to
          HGNC <em>precursor</em> genes. HGNC tracks precursors
          (<code>MIR21</code>, <code>MIR124-1/2/3</code>) while the corpus uses
          mature names (<code>miR-21</code>, <code>miR-34c-5p</code>); both are
          canonicalised to a mature key (strip <code>hsa-</code>, unify
          <code>miRNA/microRNA&rarr;miR</code>, strip arm <code>-5p/-3p</code> and
          the precursor copy number) and matched against the
          <code>hsa-mir-&hellip;</code> aliases / <code>MIR</code> symbols. One
          locus &rarr; matched (<code>miR-21</code>&rarr;MIR21); several loci
          &rarr; <em>ambiguous</em> candidate set
          (<code>miR-124</code>&rarr;MIR124-1/2/3). Compound / family / reagent
          forms (<code>miR-29a/b/c</code>, clusters, <code>miRFP&hellip;</code>)
          are excluded. Tagged <code>miR mature&rarr;gene</code>.</li>
      <li><strong>PROMOTER</strong> &mdash; strips a trailing
          <code>&nbsp;promoter</code> / <code>&nbsp;promoter region</code> /
          <code>&nbsp;gene promoter</code> regulatory descriptor
          (case-insensitive) and re-matches against the symbol fields with the
          case/hyphen-folded key (<code>MGMT promoter</code>&rarr;MGMT,
          <code>TERT promoter</code>&rarr;TERT, <code>c-jun promoter</code>&rarr;JUN);
          recorded with <code>match_mode</code>&nbsp;=&nbsp;<code>promoter stripped</code>.
          Names the gene's regulatory region, so it folds in with the
          gene-identity matches.</li>
      <li><strong>DEHYPHEN</strong> &mdash; for <code>&lt;UPPER&gt;-&lt;digits&gt;</code>
          tokens (<code>^[A-Z]+-[0-9]+$</code>) the hyphen is mere punctuation, so
          it is <em>removed</em> (not folded to a space) and the joined form
          matched against the symbol fields (<code>MMP-9</code>&rarr;MMP9,
          <code>ICAM-1</code>&rarr;ICAM1, <code>IGF-1</code>&rarr;IGF1); recorded
          with <code>match_mode</code>&nbsp;=&nbsp;<code>hyphen-removed</code>.
          A curated exclusion list drops drug / cell-line look-alikes
          (<code>YC-1</code>, <code>ISO-1</code>, <code>IWR-1</code>,
          <code>THP-1</code>, &hellip;). An <em>approved-symbol guard</em>
          resolves an otherwise-ambiguous hit to the gene whose approved symbol
          (not merely an alias) matched, when exactly one such gene exists
          (<code>FGF-2</code>&rarr;FGF2, <code>XBP-1</code>&rarr;XBP1).</li>
      <li><strong>P-PREFIX</strong> &mdash; strips a bare leading <code>p</code>
          before an upper-case letter (the no-delimiter phospho/protein form
          <code>pAKT</code>&rarr;AKT, <code>pSTAT3</code>&rarr;STAT3,
          <code>pVHL</code>&rarr;VHL) and matches the symbol fields with the
          folded key; recorded as <code>p-prefix stripped</code>. Distinct from
          PHOSPHO, which needs a <code>p-</code>/<code>phospho-</code> delimiter.
          A curated exclusion list drops constructs / vectors
          (<code>pMIR</code>, <code>pCEP4</code>, <code>pORF</code>, &hellip;) and
          the approved-symbol guard disambiguates (<code>pMET</code>&rarr;MET,
          <code>pTF</code>&rarr;TF).</li>
      <li><strong>WT</strong> &mdash; strips a leading wild-type marker
          (<code>WT&nbsp;</code> or <code>wild-type&nbsp;</code>) and matches the
          symbol fields with the folded key (<code>WT IDH1</code>&rarr;IDH1,
          <code>wild-type p53</code>&rarr;TP53). Denotes
          the unaltered gene, so it folds in with the gene-identity matches;
          recorded with <code>match_mode</code>&nbsp;=&nbsp;<code>WT stripped</code>.</li>
      <li><strong>NUCLEAR</strong> &mdash; strips a leading <code>nuclear&nbsp;</code>
          localization prefix and matches the symbol fields with the folded key
          (<code>nuclear YAP</code>&rarr;YAP, <code>nuclear p53</code>&rarr;TP53).
          Denotes the gene's product, so it folds in with the gene-identity
          matches; recorded with
          <code>match_mode</code>&nbsp;=&nbsp;<code>nuclear stripped</code>.</li>
      <li><strong>TARGET-GENES</strong> &mdash; strips a trailing
          <code>&nbsp;target genes</code> suffix and matches the regulator
          against the symbol fields (<code>MYC target genes</code>&rarr;MYC).
          <em>Caveat:</em> <code>X target genes</code> denotes the genes
          <em>regulated by</em> X, not X itself &mdash; so this is the
          <strong>regulator</strong>, recorded under a distinct
          <code>match_mode</code>&nbsp;=&nbsp;<code>regulatory target-set (not the
          gene)</code> to keep it separable from true gene-identity matches.</li>
      <li><strong>HISTONE</strong> (bespoke) &mdash; maps histone marks, mutants
          and variants to HGNC histone <em>genes</em>. HGNC encodes each histone
          protein in many genes (H3&nbsp;in 23, H4&nbsp;in 16&hellip;), so a mark
          like <code>H3K27me3</code> or mutant <code>H3.3K27M</code> &mdash; which
          is on the <em>protein</em> &mdash; maps to the variant's gene(s):
          variant-specific ones resolve to a single gene
          (<code>H2A.X</code>&rarr;H2AX, <code>CENP-A</code>&rarr;CENPA), H3.3
          forms to <code>H3-3A/B</code>, and a generic family mark to the
          <em>canonical</em> set (<code>H3K27me3</code>&rarr; the 17 canonical H3
          genes &mdash; <code>H3C*</code> + <code>H3-3A/B</code>, excluding the
          centromeric/testis/primate variants CENPA/H3-4/H3-5/H3-7/H3Y1/H3Y2
          &mdash; <em>ambiguous</em>).
          Cell-line look-alikes (<code>H358</code>, <code>H460</code>,
          <code>H19</code>) and histone-modifying enzymes are rejected. Tagged
          <code>histone</code>.</li>
      <li><strong>C-ONCOGENE</strong> &mdash; classic <code>c-</code> proto-oncogene
          spellings whose hyphen/case variants and descriptor-bearing forms miss
          HGNC's own alias entries: <code>c-Kit</code>&nbsp;/&nbsp;<code>cKit</code>&rarr;KIT,
          <code>c-Src kinase</code>&rarr;SRC, <code>c-Abl</code>&rarr;ABL1,
          <code>cFos</code>&rarr;FOS, <code>cJun</code>&rarr;JUN,
          <code>cMyc</code>&rarr;MYC, <code>Cret</code>&rarr;RET,
          <code>c-Raf</code>&rarr;RAF1, <code>c-Cbl</code>&rarr;CBL,
          <code>c-Mpl</code>&rarr;MPL, <code>c-Met</code>&rarr;MET,
          <code>cROS</code>&rarr;ROS1, <code>c-Fgr</code>&rarr;FGR,
          <code>c-Fes</code>&rarr;FES, <code>c-Yes</code>&rarr;YES1,
          <code>c-Sis</code>&rarr;PDGFB, <code>c-Fms</code>&rarr;CSF1R,
          <code>c-Mos</code>&rarr;MOS, <code>c-Rel</code>&rarr;REL,
          <code>c-Mil</code>&rarr;RAF1, <code>c-Myb</code>&rarr;MYB,
          <code>c-ErbA</code>&rarr;THRA (cf. c-erbB=EGFR), <code>c-Crk</code>&rarr;CRK,
          <code>c-Maf</code>&rarr;MAF, <code>c-Ski</code>&rarr;SKI,
          <code>c-Fps</code>&rarr;FES, <code>c-Pim</code>&rarr;PIM1,
          <code>c-Cot</code>&rarr;MAP3K8, <code>c-Mer</code>&rarr;MERTK,
          <code>c-Sea</code>&rarr;MST1R. The HGF receptor
          is the single gene MET, so the descriptive <code>HGF receptor</code>&nbsp;/
          <code>MET receptor</code> forms also map to <code>MET</code> (trailing
          descriptors such as <code>c-Met protein/kinase/promoter</code> ignored).
          A word boundary excludes look-alikes (<code>ckitCSca-1C</code> is not
          c-Kit, <code>c-Metis</code> is not c-Met, <code>c-MycUHRF1</code> is a
          fusion), and the map is a curated allowlist &mdash; the corpus is full of
          unrelated <code>c-</code>/<code>C</code> tokens (cyclin, CAR, circRNA,
          CREB, CRISPR, &hellip;) that must <em>not</em> map to a gene. A guard also
          keeps <code>c-Jun N-terminal kinase</code> (JNK = MAPK8/9/10) from
          mis-mapping to the <code>JUN</code> transcription factor.
          Tagged <code>c-oncogene&rarr;gene alias</code>. Numbered oncogene families
          resolve per isoform &mdash; <code>c-erbB-1/2/3/4</code>&rarr;EGFR/ERBB2/3/4
          (HGNC lists only <code>c-ERB-2</code>, no inner &lsquo;b&rsquo;),
          <code>c-Ets-1/2</code>&rarr;ETS1/2,
          <code>c-Ha/Ki/N-ras</code>&rarr;HRAS/KRAS/NRAS,
          <code>c-Bcl-2/3/6</code>&rarr;BCL2/3/6 &mdash; while a bare
          <code>c-Ets</code> or <code>c-ras</code> is ambiguous, routed by the
          COMPLEX stage to its gene set ([ETS1, ETS2] / [HRAS, KRAS, NRAS]) and
          tagged <code>c-oncogene&rarr;genes (ambiguous)</code>.</li>
      <li><strong>GREEK-LETTER</strong> &mdash; expands a lone Roman letter that
          stands for a Greek symbol to its <em>spelled-out Greek word</em>
          (<code>a</code>&rarr;alpha, <code>b</code>&rarr;beta,
          <code>g</code>&rarr;gamma, <code>k</code>&rarr;kappa&hellip;) <em>and</em>
          a leading abbreviation to its full word (<code>IFN</code>&rarr;interferon,
          <code>IL</code>&rarr;interleukin,
          <code>TGF</code>&rarr;transforming growth factor,
          <code>TNF</code>&rarr;tumor necrosis factor,
          <code>IGF</code>&rarr;insulin-like growth factor,
          <code>EGF</code>&rarr;epidermal growth factor,
          <code>PDGF</code>&rarr;platelet-derived growth factor,
          <code>VEGF</code>&rarr;vascular endothelial growth factor), then re-matches against all six
          fields (<code>NF-kB</code>&rarr;NF-kappaB&rarr;NFKB1,
          <code>TNF-a</code>&rarr;TNF, <code>b-catenin</code>&rarr;CTNNB1,
          <code>IFN-g</code>&rarr;interferon-gamma&rarr;IFNG,
          <code>IL-6 receptor</code>&rarr;IL6R via the <em>name</em> field). A
          &ldquo;lone&rdquo; letter is one not adjacent to lowercase letters, so
          word-internal letters are untouched; only expansions that actually hit
          an HGNC surface link (so <code>p53</code>&rarr;<code>pi53</code> just
          falls through), <code>m</code>/<code>u</code> (mu) is excluded to avoid
          the &lsquo;mouse&rsquo; prefix (<code>mC2</code>), and the abbreviation
          only expands before a hyphen/space/digit (so it leaves <code>ILK</code>
          alone). Tagged <code>greek-letter spelled out</code>.</li>
      <li><strong>DELSEP</strong> &mdash; a final separator-deletion mop-up that
          removes <em>all</em> hyphens and spaces from both the query and the HGNC
          surfaces (a destructive fold, vs <strong>TYPE&nbsp;2</strong>'s
          <code>-</code>&harr;space swap and <strong>DEHYPHEN</strong>'s narrow
          <code>^[A-Z]+-[0-9]+$</code> case), so <code>PD-L1</code>,
          <code>PD L1</code> and <code>PDL1</code> share one key
          (<code>pdl1</code>&rarr;CD274) regardless of which side carries the
          separator &mdash; the no-separator surfaces the looser folds miss.
          Matched against the symbol fields and tagged
          <code>separator-removed</code>. It runs last so every precise stage gets
          first claim; the drug / cell-line exclusion list is carried over from
          DEHYPHEN, the approved-symbol guard resolves the alias-only ties this
          looser key can create, and any remaining key collision degrades to
          <em>ambiguous</em> (&ge;2 genes) rather than a wrong single link.</li>
    </ul></li>
  <li><strong>COMPLEX (bespoke, runs first with priority).</strong> Maps curated
      protein-complex / heterodimer names to their HGNC <em>subunit</em> genes.
      Some entities name an assembled complex that HGNC catalogs only as subunit
      genes &mdash; e.g. <code>IL-12</code> (the IL-12p70 heterodimer) has no
      single gene, only <code>IL12A</code>&nbsp;(p35) + <code>IL12B</code>&nbsp;(p40).
      It runs before the cascade so the curated map is authoritative (otherwise
      <code>IL-23</code> would alias-match <code>IL23A</code> alone, and
      <code>IL-27</code> a stale junk set). The full complex maps to the subunit
      set (<em>ambiguous</em>: <code>IL-12</code>&rarr;[IL12A, IL12B],
      <code>IL-23</code>&rarr;[IL12B, IL23A], <code>IL-27</code>&rarr;[EBI3, IL27],
      <code>AP-1</code>&rarr;[FOS, FOSB, JUN, JUND]); a subunit-specific form
      resolves to one gene (<code>IL-12 p40</code>&rarr;IL12B). Tagged
      <code>complex&rarr;subunits</code>. The same stage also routes a bare
      <em>receptor-family</em> name &mdash; which is <em>one-of</em> a paralog set,
      not an assembled complex &mdash; to the family's genes with a distinct
      <code>receptor family (ambiguous)</code> tag:
      <code>VEGF receptor</code>&nbsp;/&nbsp;<code>VEGFRs</code>&rarr;[FLT1, KDR, FLT4]
      (VEGFR-1/2/3),
      <code>PDGF receptor</code>&nbsp;/&nbsp;<code>PDGFRs</code>&rarr;[PDGFRA, PDGFRB],
      <code>FGFR</code>&nbsp;/&nbsp;<code>FGFRs</code>&rarr;[FGFR1&ndash;4],
      <code>EGFR family</code>&nbsp;/&nbsp;<code>ERBB receptors</code>&nbsp;/&nbsp;<code>HER</code>&rarr;[EGFR, ERBB2, ERBB3, ERBB4],
      <code>TRK</code>&nbsp;/&nbsp;<code>NTRK</code>&nbsp;/&nbsp;<code>Trk receptors</code>&rarr;[NTRK1, NTRK2, NTRK3]
      (with <code>neurotrophin receptor</code>&rarr;[NGFR, NTRK1, NTRK2, NTRK3], adding p75).
      A specific member (<code>VEGFR2</code>&rarr;KDR, <code>PDGFR&beta;</code>&rarr;PDGFRB,
      <code>HER2</code>&rarr;ERBB2, <code>TrkA</code>&rarr;NTRK1) and HGNC-sanctioned
      singulars (<code>ErbB</code>&rarr;EGFR via prev_symbol) still link to one gene
      upstream; bare <code>TRK</code> is overridden here from its [TPM3, NTRK1]
      fusion-alias artifact to the true Trk family.
      The same stage also recognises the <strong>NF-kB transcription factor</strong>
      &mdash; mirroring the greek pipeline's NF-kappaB complex &mdash; mapping the
      generic name to its canonical p50/p65 dimer
      <code>NF-kB</code>&nbsp;/&nbsp;<code>NFkB</code>&nbsp;/&nbsp;<code>NF-kappaB</code>&rarr;[NFKB1, RELA]
      (the Roman <code>k</code> stands for kappa, so both spellings qualify; the
      <code>&nbsp;protein</code>/<code>&nbsp;gene</code> descriptor, a leading
      <code>p-</code>, and trailing <code>complex</code>/<code>target genes</code>
      forms come along too). Specific genes are <em>not</em> the complex and keep
      linking to themselves &mdash; the precursor genes <code>NFKB1</code>/<code>NFKB2</code>,
      the I&kappa;B inhibitors <code>NFKBIA</code>/<code>NFKBIB</code>/<code>NFKBIE</code>/<code>NFKBIZ</code>,
      and named subunits (<code>NF-kB p65</code>) &mdash; matching greek's
      <code>_NFKB_SINGLE</code> + p65/p50 guards. Tagged <code>complex&rarr;subunits</code>.</li>
  <li><strong>Word-order cosine fallback.</strong> For the
      multi-word values still unmatched, compare against multi-word HGNC name
      values by <strong>word-level term-frequency cosine</strong>
      (case-insensitive) and keep pairs &ge;{COSINE_THRESHOLD}. At these short
      lengths that threshold is reachable only at 1.0, so it surfaces pure
      <em>word-order variants</em> (<code>aurora kinase A</code> vs
      <code>Aurora A kinase</code>) that exact keys cannot align.</li>
</ol>

<p><strong>Matched output: one consolidated JSON library.</strong> All three
exact steps write into a single dictionary, <code>roman.json</code> &mdash;
<em>keys</em> are the GENETIC named entities, <em>values</em> are the linked HGNC
approved symbols, each carrying annotations: the entity's
<code>occurrences</code>, and the <strong>match-type</strong> of the value
&mdash; which HGNC field linked it
(<code>symbol</code>/<code>alias_symbol</code>/<code>prev_symbol</code> for the
symbol steps, <code>name</code>/<code>alias_name</code>/<code>prev_name</code>
for the name/phospho steps) and the mode that explains the link
(<code>case-sensitive</code>, <code>case-insensitive</code>,
<code>hyphen vs white space</code>, <code>phospho stripped</code>,
<code>anti- stripped</code>, <code>cleaved- stripped</code>, <code>protein stripped</code>,
<code>transcript/mRNA stripped</code>, <code>variant/mutant stripped</code>,
<code>missense stripped</code>, <code>inhibitor-i stripped</code>,
<code>miR mature&rarr;gene</code>, <code>promoter stripped</code>,
<code>hyphen-removed</code>, <code>p-prefix stripped</code>,
<code>WT stripped</code>, <code>nuclear stripped</code>,
<code>regulatory target-set (not the gene)</code>, <code>histone</code>,
<code>complex&rarr;subunits</code>, <code>receptor family (ambiguous)</code>,
<code>c-oncogene&rarr;gene alias</code>,
<code>c-oncogene&rarr;genes (ambiguous)</code>,
<code>greek-letter spelled out</code>, or
<code>manual correction</code>).
Because each entity is linked by the
strictest step that succeeds and leftovers flow forward, the keys are unique
across steps, so the merge is collision-free and the field/mode annotations tell
you exactly how (and effectively at which step) each entity matched. Example:</p>
<pre><code>// roman.json
{{ "EGFR":      {{ "hgnc_symbol": "EGFR",  "occurrences": 2328,
                 "match_field": "symbol", "match_mode": "case-sensitive" }},
  "caspase-3": {{ "hgnc_symbol": "CASP3", "occurrences": 662,
                 "match_field": "name",   "match_mode": "hyphen vs white space" }} }}</code></pre>
<p>For the exact steps the mode is decided per entity: an exact surface hit is
<code>case-sensitive</code>; otherwise a casefold-equal surface is
<code>case-insensitive</code>; otherwise the only way the key can coincide is a
hyphen&harr;whitespace swap.</p>
<p><strong>Ambiguous output: one consolidated JSON library.</strong> The
ambiguous buckets &mdash; entities whose key hits <em>&ge;2 different genes</em>
&mdash; are likewise merged into <code>roman_ambiguous.json</code>, in the same
shape, except the value's <code>hgnc_symbol</code> is the <em>list</em> of
candidate symbols (with <code>n_genes</code> and the parallel
<code>hgnc_ids</code>); <code>match_field</code>/<code>match_mode</code> are
computed over the surfaces of all candidate genes. Example:</p>
<pre><code>// roman_ambiguous.json
{{ "ERK": {{ "hgnc_symbol": ["EPHB2", "MAPK1"], "occurrences": 714,
           "match_field": "alias_symbol", "match_mode": "case-sensitive",
           "n_genes": 2, "hgnc_ids": ["HGNC:3393", "HGNC:6871"] }} }}</code></pre>
<p>The word-order cosine fallback keeps its own library
<code>roman_cosine.json</code> (its values additionally carry the retained
<code>cosine</code> similarity and the matched <code>hgnc_value</code>). Every
value the whole pipeline still cannot link is written to
<code>unmatched.tsv</code>.</p>
<p><strong>No per-step TSVs.</strong> The cascade's intermediate buckets &mdash;
the per-step matched tables (<code>types_1_matched.tsv</code>,
<code>types_2_matched.tsv</code>, <code>names_matched.tsv</code>) and the
intermediate leftovers (<code>unmatched_type_1.tsv</code>,
<code>unmatched_type_2.tsv</code>) &mdash; are no longer written; they are held
in memory and surface only through the consolidated <code>roman.json</code> /
<code>roman_ambiguous.json</code>. Only <code>unmatched.tsv</code> (the final
leftover) persists as a TSV. Stale copies of those per-step files from earlier
runs are deleted when the pipeline starts.</p>
<p><strong>Manual correction map.</strong> A small curated map
(entity&nbsp;&rarr;&nbsp;approved symbol) is applied last to override the few
links the cascade gets wrong &mdash; typically a single-gene <em>alias
collision</em> the approved-symbol guard cannot catch (only the wrong gene
matches). For example <code>pYAP</code> auto-links to YY1AP1 via its alias
<code>YAP</code>, but is corrected to <strong>YAP1</strong>. A correction changes
only the gene: the entity's original <code>match_mode</code> is <em>preserved</em>
(so a semantic tag like <code>regulatory target-set (not the gene)</code> on
<code>YAP target genes</code> survives), and the override is flagged via
<code>match_field</code>&nbsp;=&nbsp;<code>manual</code>; a previously-unmatched
corrected entity gets <code>match_mode</code>&nbsp;=&nbsp;<code>manual correction</code>.</p>

<h2>Results</h2>
<p>STAGE 1: read <strong>{s1['n_files']:,}</strong> files,
<strong>{s1['occurrences']:,}</strong> GENETIC occurrences,
{s1['unique_before']:,}&rarr;{s1['unique_after']:,} unique after dash
normalization; split <strong>{s1['greek']:,}</strong> Greek /
<strong>{s1['roman']:,}</strong> Roman entities.</p>
{pass_sections}

<h3>Word-order cosine fallback (&ge; {COSINE_THRESHOLD})</h3>
<p><strong>{len(cos['rows'])}</strong> pair(s) over {cos['n_inputs_mw']:,}
multi-word leftovers, written as the JSON library
<code>roman_cosine.json</code> ({cos['n_entities']} entit(y/ies), each with
its retained <code>cosine</code> value):</p>
<table>
  <tr><th>input_value</th><th>hgnc_value</th><th class="num">cosine</th>
      <th>field</th><th>hgnc_id</th><th>symbol</th><th class="num">occ</th></tr>
  {cos_rows or '<tr><td colspan="7"><em>none</em></td></tr>'}
</table>

<h2>Cumulative coverage</h2>
<p>Exact stages (TYPE&nbsp;1 + TYPE&nbsp;2 + NAME) link
<strong>{linked:,}</strong> of the <strong>{roman_occ:,}</strong> Roman GENETIC
occurrences (<strong>{100*linked/roman_occ:.1f}%</strong>), counting
matched&nbsp;+&nbsp;ambiguous; STAGE&nbsp;4 adds word-order variants on top of
the remainder.</p>
{nongene_section}

</body>
</html>
"""
    HTML_OUT.write_text(doc, encoding="utf-8")


if __name__ == "__main__":
    main()
