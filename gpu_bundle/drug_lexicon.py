#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""drug_lexicon.py -- recognise drug names the BioBERT chemical model misses.

WHY THIS EXISTS
---------------
The chemical NER behaves like a BC5CDR-trained model: excellent on small molecules, poor on
biologics. Measured on the lung_adeno corpus, twelve named monoclonal antibodies appear 680
times in the candidate sentences and only 181 of those mentions (27%) were tagged at all --
and 18 of them came back as GENETIC, because antibodies are routinely written by their target.

Two cheap recognisers close most of that gap, and neither needs a model retrained:

  INN STEMS      The WHO's International Nonproprietary Name scheme puts the substance class in
                 the suffix, so "-mab" (and the 2021 replacements "-tug", "-bart", "-mig", plus
                 "-cept" for fusion proteins) identifies an antibody therapeutic with almost no
                 false positives. This is what catches agents no dictionary carries yet:
                 inetetamab and tinurilimab appear 62 times in this corpus.

  NCI THESAURUS  37,619 drug concepts and 140,783 names, including every antibody above. It is
                 the best oncology drug vocabulary that can be downloaded as a flat file, and it
                 carries a stable code per concept, so a match can be normalised later instead
                 of being thrown away for lacking a ChEBI id.

The two are complementary: the dictionary knows brand names, codes and spellings the regex
cannot guess ("MK-3475", "Keytruda"), while the regex covers what the dictionary has not caught
up with. Matches carry their provenance so either can be audited or switched off.

Build the dictionary once (needs network), then it is a local file::

    python drug_lexicon.py --build                     # -> databases/ncit_drugs.json
    python drug_lexicon.py --test "treated with pembrolizumab plus carboplatin"
"""

import argparse
import io
import json
import os
import re
import time
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LEXICON = ROOT / "databases" / "ncit_drugs.json"
NCIT_URL = os.environ.get(
    "NCIT_FLAT_URL", "https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Thesaurus.FLAT.zip")

# NCIt semantic types that denote a substance one could administer.
#
# "Amino Acid, Peptide, or Protein" is deliberately NOT here, though it is where the antibodies
# also sit: NCIt types every gene product with it, so including it made the lexicon tag KRAS,
# EGFR, TP53 and eleven other symbols as drugs -- the gene/drug confusion this file exists to
# reduce, reintroduced from the other side. Measured over the whole thesaurus, dropping it costs
# nothing and fixes everything: 15/15 test antibodies still resolve (they are Pharmacologic
# Substances and Immunologic Factors too) and all 14 gene-symbol false positives disappear.
# "Immunologic Factor" is excluded for the same reason: it is where NCIt keeps the cytokines,
# chemokines and receptors, so including it dragged in IL-6, CTLA4, CXCL8, CCL20, PSCA and ITGAL
# -- genes the gene model had already tagged correctly. Measured: dropping it takes gene-ish
# names from 11 to 2 and still resolves 15/15 test antibodies, which carry
# "Pharmacologic Substance" as well.
DRUG_SEMANTIC_TYPES = {
    "Pharmacologic Substance", "Antibiotic", "Organic Chemical",
}

# INN stems for antibodies and fusion proteins. Small-molecule stems (-tinib, -parib, ...) are
# deliberately absent: the model already finds those well -- osimertinib and gefitinib are among
# the corpus' best-covered drugs -- so adding them would buy noise, not recall.
INN_STEMS = ("mab", "tug", "bart", "mig", "cept")
INN_RE = re.compile(r"\b[a-z][a-z-]{4,}(?:" + "|".join(INN_STEMS) + r")\b", re.I)

# Names too generic to tag on sight. NCIt lists them as substances and they are, but in prose
# they are far more often ordinary words -- "lead to", "in the light of", "the gold standard".
STOP_NAMES = {
    "lead", "gold", "iron", "water", "oxygen", "air", "salt", "sugar", "alcohol", "acid",
    "protein", "gas", "oil", "tin", "carbon", "light", "cream", "factor", "agent", "drug",
    "control", "vehicle", "saline", "serum", "sodium", "calcium", "potassium", "silver",
    "copper", "zinc", "nitrogen", "hydrogen", "helium", "ice", "food", "milk", "honey",
    "tobacco", "coal", "dust", "smoke", "ash", "clay", "resin", "wax", "starch", "fiber",
    "gel", "foam", "powder", "tablet", "capsule", "solution", "suspension", "ointment",
    # class nouns: true substances, useless as nodes, and they swamp the real drug names
    "antibody", "antibodies", "antigen", "antigens", "cytokine", "cytokines", "chemokine",
    "chemokines", "complement", "conjugate", "lipid", "lipids", "peptide", "hormone", "enzyme",
    "medicine", "vaccine", "placebo", "chemotherapy", "il-6", "cxcl10",
    # ordinary words NCIt also lists as (or as synonyms of) substances
    "same", "impact", "stability", "chip", "shape", "spleen", "adrenal", "target", "bridge",
    "focus", "matrix", "vision", "balance", "motion", "reason", "spirit", "novel", "select",
    # surfaced by auditing what the lexicon adds to this corpus: assay names, anatomy and
    # collective nouns that are not agents anyone administered
    "thca", "rtca", "stomach", "regulatory t cell", "drug treatment", "medications",
    "alkaloid", "corticosteroid", "chemotherapies", "immunotherapy",
}
MIN_NAME_LEN = 4          # shorter surfaces are ambiguous far more often than they are useful
MAX_NGRAM = 5             # longest dictionary name, in tokens, worth scanning for

_TOKEN_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9\-']*")


def build(dest=LEXICON, url=NCIT_URL, verbose=True):
    """Download the NCIt flat file and distil it to {name: code} for drug concepts."""
    req = urllib.request.Request(url, headers={"User-Agent": "relationship-graphs/1.0"})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=300) as r:
        blob = r.read()
    rows = zipfile.ZipFile(io.BytesIO(blob)).read("Thesaurus.txt").decode("utf-8", "replace")
    names, labels, concepts = {}, {}, 0
    for line in rows.splitlines():
        p = line.split("\t")
        if len(p) < 8:
            continue
        if not ({s.strip() for s in p[7].split("|")} & DRUG_SEMANTIC_TYPES):
            continue
        concepts += 1
        code = p[0].strip()
        # NCIt lists the preferred term first; keep it so a match can be reported
        # with a human-readable label, the way a ChEBI match is
        pref = " ".join(p[3].split("|")[0].split())
        if pref:
            labels.setdefault(code, pref)
        for syn in p[3].split("|"):
            n = " ".join(syn.split()).lower()
            if len(n) < MIN_NAME_LEN or n in STOP_NAMES:
                continue
            if not any(c.isalpha() for c in n):
                continue
            if len(n.split()) > MAX_NGRAM:
                continue
            # NCIt forms with two or more capitals are acronyms (SAMe, ChIP, TXA2). Matching
            # those case-insensitively fires on ordinary words, so remember the exact form and
            # demand it; single-capital names (Bevacizumab) stay case-insensitive.
            form = " ".join(syn.split())
            cased = form if sum(1 for c in form if c.isupper()) >= 2 else None
            names.setdefault(n, [code, cased] if cased else code)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps({"source": url, "built": time.strftime("%Y-%m-%d"),
                                "concepts": concepts, "names": names, "labels": labels}),
                    encoding="utf-8")
    if verbose:
        print(f"  {concepts:,} drug concepts -> {len(names):,} names "
              f"({dest.stat().st_size/1e6:.1f} MB) in {time.time()-t0:.0f}s -> {dest}")
    return dest


class DrugLexicon:
    """Dictionary + INN-stem matcher over raw sentence text.

    Longest match wins and matches never overlap, so "carboplatin" inside "carboplatin-paclitaxel"
    yields the longer name when NCIt knows it and the two components separately when it does not.
    """

    def __init__(self, path=LEXICON, use_inn=True):
        self.names, self.labels, self.use_inn, self.source = {}, {}, use_inn, None
        p = Path(path)
        if p.exists():
            try:
                d = json.loads(p.read_text(encoding="utf-8"))
                self.names = d.get("names") or {}
                self.labels = d.get("labels") or {}
                self.source = d.get("built")
            except Exception:
                self.names, self.labels = {}, {}

    @property
    def ready(self):
        return bool(self.names) or self.use_inn

    def lookup(self, surface):
        """(code, label) for an exact drug name, honouring the acronym casing rule."""
        v = self.names.get(" ".join((surface or "").split()).lower())
        if not v:
            return None
        code, cased = (v, None) if isinstance(v, str) else (v[0], v[1])
        if cased and " ".join((surface or "").split()) != cased:
            return None
        return code, self.labels.get(code, cased or surface)

    def find(self, text):
        """[(start, end, surface, code|None, origin)] -- non-overlapping, longest first."""
        if not text:
            return []
        toks = [(m.start(), m.end(), m.group(0)) for m in _TOKEN_RE.finditer(text)]
        hits, taken = [], []

        def free(a, b):
            return all(b <= s or a >= e for s, e, *_ in taken)

        if self.names:
            for i in range(len(toks)):
                for n in range(min(MAX_NGRAM, len(toks) - i), 0, -1):
                    a, b = toks[i][0], toks[i + n - 1][1]
                    if not free(a, b):
                        continue
                    phrase = " ".join(t[2] for t in toks[i:i + n]).lower()
                    v = self.names.get(phrase)
                    if not v:
                        continue
                    code, cased = (v, None) if isinstance(v, str) else (v[0], v[1])
                    if cased and text[a:b] != cased:
                        continue                       # acronym: only the exact casing counts
                    taken.append((a, b))
                    hits.append((a, b, text[a:b], code, "ncit"))
                    break
        if self.use_inn:
            for m in INN_RE.finditer(text):
                a, b = m.start(), m.end()
                if free(a, b):
                    taken.append((a, b))
                    v = self.names.get(m.group(0).lower())
                    hits.append((a, b, text[a:b], v if isinstance(v, str) else (v[0] if v else None), "inn"))
        hits.sort(key=lambda h: h[0])
        return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", action="store_true", help="download NCIt and write the lexicon")
    ap.add_argument("--dest", default=str(LEXICON))
    ap.add_argument("--test", metavar="SENTENCE", help="show what the lexicon finds in a sentence")
    args = ap.parse_args()
    if args.build:
        build(Path(args.dest))
    if args.test:
        lex = DrugLexicon(args.dest)
        print(f"lexicon: {len(lex.names):,} names (built {lex.source})")
        for a, b, s, code, origin in lex.find(args.test):
            print(f"  [{a:4}:{b:<4}] {s:28} {origin:5} {code or '-'}")
    if not args.build and not args.test:
        ap.print_help()


if __name__ == "__main__":
    main()
