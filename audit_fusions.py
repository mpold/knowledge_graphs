#!/usr/bin/env python3
"""audit_fusions.py -- regression test for the ChimerDB fusion attribute on gene nodes.

`fus` / `fpart` / `fn` / `fseq` decide how the graph reads a gene with no drug pocket: a fusion
junction exists in no normal cell, so it is a selectivity handle even where inhibition has
nothing to bind. That claim is only as good as the attribute behind it, and the attribute rests
on a database that will be re-released. This is what catches a silent regression -- the fusion
counterpart to audit_drug_targets.py, kept separate because the checks and the sources are
disjoint.

WHAT IT CHECKS (1-5 fail the run, 6-7 are drift reports)

  1. staging          ChimerKB4.tsv is present and carries the columns the loader reads. The
                      file is OPTIONAL, so an absent one SKIPS rather than fails -- but a
                      present file with renamed columns must fail loudly, because the loader
                      would silently return {} and every node would go quietly unflagged.
  2. anchors          Genes whose side is not a judgement call. TMPRSS2 must be 5' and only 5':
                      it donates the promoter and keeps none of its protein, which is the whole
                      reason TMPRSS2-ERG is treated as an ERG lesion. NRG1 and CD74 must be 3'
                      and 5'. If a release flips a side, every inference built on it inverts.
  3. negative         PRDM14, USP39, SOX2, MCL1 must NOT be flagged. They are the genes the
                      write-up calls pocketless AND junctionless; if ChimerKB starts listing
                      them, the argument in candidate_targets.md section 4a changes and should
                      be re-read rather than silently drifting.
  4. seq never flags  The flagged set must come from ChimerKB or from the small curated
                      FUSION_SUPPLEMENT, never from ChimerSeq, which is algorithmic -- TP53 in
                      65 rows, CTNNB1 in 50, USP39 in 40 --
                      and no recurrence threshold rescues it: >=3 samples still admits USP39
                      while losing PAX3, FOXO1, YAP1 and MYB.
  5. payload          Every gene node carries the four fields; every non-gene node carries them
                      empty. A missing field is how a template edit breaks the HTML silently.
  6. curation drift   Genes candidate_targets.md section 6a lists that ChimerKB does not, and
                      vice versa. Neither source is a superset of the other, so this is drift to
                      read, not an error to fix.
  7. promoter swaps   Flagged genes whose partners include an immunoglobulin or TCR locus. These
                      are promoter substitutions, not chimeric proteins -- MYC and BCL6 arrive
                      this way -- and the distinction is invisible in the flag itself.

Run::  python audit_fusions.py [--data-root kaggle_working] [--graph PATH] [--quiet]
Exit code 0 when checks 1-5 pass (or staging is absent), 1 otherwise.
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

import high_confidence_g as HC

# side assignments that are structural, not statistical: if one of these flips, something is wrong
ANCHORS = {
    "TMPRSS2": "5p",     # donates the promoter, keeps none of its protein
    "CD74": "5p",        # CD74-NRG1, CD74-ROS1
    "NRG1": "3p",        # the ligand domain the fusion is named for
    "EWSR1": "both",     # 33 partners, both orientations
    "ERG": "both",
    "PAX3": "both",
}
# pocketless AND junctionless -- the pairing candidate_targets.md 4a rests on
NEGATIVE = ("PRDM14", "USP39", "SOX2", "MCL1")
# section 6a, for the drift report
CURATED_6A = ("EWSR1 FLI1 ERG TMPRSS2 PML RUNX1 PAX3 FOXO1 FOXO6 MEIS1 HOXA9 HMGA2 MYB RELA "
              "PPARG NRG1 CD74 YAP1 DDIT3 POU5F1 TP63 CREBBP EP300 MAP3K8 SMAD3 CRLF2 PVT1 "
              "NBEA CD44 FOS RARB RARG CTNNB1").split()
# Immunoglobulin and TCR loci, matched as SYMBOLS rather than by prefix. A prefix test is what
# this project keeps auditing out of other people's data: "TRA" also opens TRAF1, TRAK1 and
# TRAPPC10, none of which is a receptor locus, and the first draft of this check duly reported
# ALK as a promoter swap on the strength of TRAF1.
IG_TR = re.compile(r"^(IG[HKL]|TR[ABDG])([VDJC])?[\d@-]*$")
NEEDED_COLUMNS = ("H_gene", "T_gene", "Fusion_pair")


def gene_nodes(path):
    s = Path(path).read_text(encoding="utf-8")
    i = s.index("const DATA=")
    j = s.index("\n", i)
    payload = json.loads(s[i + len("const DATA="):j].rstrip().rstrip(";"))
    return payload["nodes"]


def kb_columns(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return next(csv.reader(fh, delimiter="\t"), [])


def newest_graph(root):
    """The dated copy stage 3 writes; the name carries the run date, so take the newest."""
    hits = sorted(Path(root).glob("*_M.html"), key=lambda p: p.stat().st_mtime, reverse=True)
    return hits[0] if hits else None


def audit(data_root, graph, quiet=False):
    HC.set_data_root(data_root)
    fails = 0

    # -- 1. staging -------------------------------------------------------------------------
    if not Path(HC.CHIMER_KB).exists():
        print(f"[1] ChimerKB4.tsv absent under {Path(data_root)/'databases'} -- the attribute is "
              f"optional, so nothing to audit.\n    Stage it (see PLACE_DATABASES_HERE.md) and "
              f"convert with chimerdb_to_tsv.py.")
        return 0
    cols = kb_columns(HC.CHIMER_KB)
    missing = [c for c in NEEDED_COLUMNS if c not in cols]
    print(f"[1] staging ............................ {len(cols)} columns"
          + (f"  MISSING {missing}" if missing else ""))
    fails += bool(missing)
    if missing:                       # nothing below can mean anything
        print("\nFAIL: the loader would return {} and every node would go quietly unflagged")
        return 1

    fusion = HC._fusion_partners()    # the real loader, not a reimplementation
    print(f"    {len(fusion):,} genes carry a fusion record")

    # -- 2. anchors -------------------------------------------------------------------------
    bad = {g: (fusion.get(g, ("absent",))[0], want) for g, want in ANCHORS.items()
           if fusion.get(g, ("absent",))[0] != want}
    print(f"[2] anchor sides ....................... {len(ANCHORS) - len(bad)}/{len(ANCHORS)} as expected"
          + ("" if not bad else "  " + ", ".join(f"{g}: got {got}, want {want}"
                                                 for g, (got, want) in bad.items())))
    fails += bool(bad)

    # -- 3. negative controls ---------------------------------------------------------------
    hit = [g for g in NEGATIVE if g in fusion]
    print(f"[3] negative controls unflagged ........ {len(NEGATIVE) - len(hit)}/{len(NEGATIVE)}"
          + (f"  NOW FLAGGED: {hit}" if hit else ""))
    fails += bool(hit)

    # -- 4. ChimerSeq must not set flags ----------------------------------------------------
    kb_genes = set()
    with open(HC.CHIMER_KB, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for g in ((r.get("H_gene") or "").strip(), (r.get("T_gene") or "").strip()):
                if g:
                    kb_genes.add(g)
    # the curated supplement is allowed to add genes ChimerKB lacks; ChimerSeq is not
    supp_genes = {g for pair in HC.FUSION_SUPPLEMENT for g in pair[:2]}
    extra = sorted(set(fusion) - kb_genes - supp_genes)
    print(f"[4] flags come from ChimerKB or the curated supplement "
          f"({len(supp_genes)} genes) ... {'yes' if not extra else 'NO: ' + str(extra[:8])}")
    fails += bool(extra)

    # -- 5. payload integrity ---------------------------------------------------------------
    if graph:
        nodes = gene_nodes(graph)
        genes = [n for n in nodes if n.get("kind") == "gene"]
        lack = [n["id"] for n in genes if any(k not in n for k in ("fus", "fpart", "fn", "fseq"))]
        dirty = [n["id"] for n in nodes if n.get("kind") != "gene" and n.get("fus")]
        flagged = [n for n in genes if n.get("fus")]
        print(f"[5] payload ............................ {len(genes)} gene nodes, "
              f"{len(flagged)} flagged"
              + (f"  MISSING FIELDS: {lack[:5]}" if lack else "")
              + (f"  NON-GENE FLAGGED: {dirty[:5]}" if dirty else ""))
        fails += bool(lack or dirty)
    else:
        print("[5] payload ............................ no graph HTML found, skipped")
        genes, flagged = [], []

    # -- 6. drift against the hand curation --------------------------------------------------
    absent = [g for g in CURATED_6A if g not in fusion]
    print()
    print(f"[6] curated in section 6a, carried by neither ChimerKB nor the supplement: "
          f"{len(absent)}")
    if absent:
        print("      " + ", ".join(absent))
    if genes and not quiet:
        extra_flagged = sorted(n["id"] for n in flagged if n["id"] not in CURATED_6A)
        print(f"    flagged in the graph but not in section 6a: {len(extra_flagged)}")
        print("      " + ", ".join(extra_flagged[:24]) + (" ..." if len(extra_flagged) > 24 else ""))

    # -- 7. promoter swaps hiding inside the flag ---------------------------------------------
    swaps = sorted(g for g, v in fusion.items()
                   if any(IG_TR.match(p) for p in v[1]))
    print(f"\n[7] flagged through an Ig/TCR partner (promoter substitution, not a fusion "
          f"protein): {len(swaps)}")
    if swaps and not quiet:
        for g in swaps[:12]:
            ig = [p for p in fusion[g][1] if IG_TR.match(p)]
            print(f"      {g:9} {fusion[g][0]:5} {len(fusion[g][1]):>3} partners incl. {', '.join(ig[:3])}")

    print(f"\n{'FAIL' if fails else 'PASS'}: checks 1-5 {'found problems' if fails else 'are clean'}")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=str(HC.DATA_ROOT),
                    help="pipeline output tree holding databases/ (default: kaggle_working)")
    ap.add_argument("--graph", help="graph HTML to check the payload of "
                                    "(default: the newest *_M.html beside this script)")
    ap.add_argument("--quiet", action="store_true", help="counts only, no per-gene listings")
    args = ap.parse_args()
    graph = args.graph or newest_graph(Path(__file__).resolve().parent)
    return audit(args.data_root, graph, args.quiet)


if __name__ == "__main__":
    sys.exit(main())
