#!/usr/bin/env python3
"""audit_drug_targets.py -- check the gene <- chemical cross-link for FALSE TARGETS.

The graph colours a gene node when a corpus chemical targets it (green = DGIdb approved
anti-neoplastic, amber = approved, pink = neither). That colour is only as good as the
cross-link behind it, and an audit of the reference run found three ways a gene got a
colour it had not earned. Each has a fix in the pipeline; this script is the regression
test for all three, so a future DGIdb or ChEBI release cannot quietly reopen one.

WHAT IT CHECKS (1-3 fail the run, 4-5 are review output)

  1. non-HGNC gene keys
     target_pharm.py used to keep an unresolvable DGIdb gene claim verbatim
     (`map_gene(gene) or {gene}`), writing entities like "NULL", "F7R" and "GAPDHL17"
     into the cross-link as if they were genes. Expect 0.

  2. identifier-namespace artifacts
     GuideToPharmacology carries its target as an "NCBIGENE:<id>" claim, but the id is
     GtoPdb's OWN target id wearing an NCBI label. DGIdb resolves most of them through
     the GtoPdb relation and lands on the right gene -- imatinib's claim is 1923 while
     ABL1 is entrez 25, tazemetostat's is 2654 while EZH2 is 2146. Where it read the id
     AS an NCBI gene id, the gene it produced is whatever holds that number: A1BG
     "targeted" by 119 drugs, and four separate ALK/FLT3 inhibitors all landing on HBEGF,
     whose entrez id is exactly the claim, 1839. The equality IS the signature, and
     load_dgidb drops those rows; this re-derives the signature independently and checks
     none survived. Expect 0.

  3. non-drug chemicals
     Metabolites and NER accidents that acquired a target: ADP reached BIRC2 through a
     consecutive run of GtoPdb ids, and L-cysteine reached HAL because the one-letter NER
     surface "C" normalizes onto it. Expect 0.

  4. surface-only matches (informational)
     A chemical whose corpus SURFACE, not its own ChEBI label, matched the drug name --
     "dasatinib monohydrate" through "dasatinib", "PD-0332991" through palbociclib. This
     is the mechanism working as intended; it is printed so a new pattern gets noticed.

  5. targets with no binding-sourced row (informational)
     DGIdb also ingests clinical-context sources (MyCancerGenome, TALC,
     ClearityFoundation, CancerCommons) whose claim is "this drug is used where this gene
     is altered", recorded with interaction_type "inhibitor" so the type gate cannot see
     them. Those that a drug's pharmacology CONTRADICTS are dropped by name in
     high_confidence_g.FALSE_TARGETS. What is left is listed here, and the ones already
     reviewed carry their verdict in REVIEWED below -- so a re-run after a database
     update shows only what is NEW to judge.

The audit deliberately re-derives its checks from the databases rather than trusting the
pipeline's own output: it imports FALSE_TARGETS and _drug_targets from high_confidence_g
(so it audits the cross-link the graph actually draws, blocklist applied) but re-reads
HGNC and the DGIdb TSV itself.

Run::  python audit_drug_targets.py [--data-root kaggle_working] [--quiet]
Exit code 0 when checks 1-3 are clean, 1 otherwise.
"""
import argparse
import collections
import csv
import json
import re
import sys
from pathlib import Path

import high_confidence_g as HC

# Chemicals that are metabolites, cofactors or single-letter NER accidents. Kept here rather
# than imported from target_pharm.py: an audit that reuses the blocklist it is auditing can
# only ever agree with itself. Matched against the ChEBI label.
NON_DRUG = re.compile(
    r"^(ADP|AMP|ATP|GTP|NAD\+?|acetyl-CoA|L-cysteine|cysteine|L-serine|serine|threonine|"
    r"tyrosine|glycine|cytosine|thymine|uracil|adenine|guanine|glucose|lactate|pyruvate|"
    r"glutathione|cholesterol|water|dimethyl sulfoxide|ethanol|"
    r"(iron|copper|zinc|silicon|platinum|calcium|magnesium|potassium|sodium|carbon|"
    r"nitrogen|oxygen|selenium|plutonium)( atom| dication| cation)?|"
    r"(charm|down|up|strange|top|bottom) quark)$", re.I)

# DGIdb sources that assert BINDING (an assay or a curated pharmacology record). Everything
# else in the file -- MyCancerGenome, MyCancerGenomeClinicalTrial, TALC,
# ClearityFoundationClinicalTrial, CancerCommons -- asserts clinical context instead.
PHARM_SOURCES = {"ChEMBL", "GuideToPharmacology", "DTC", "TTD", "DrugBank"}

# Check 5 findings already reviewed by hand, with the verdict. A pair here is NOT re-flagged;
# anything else the check turns up is new and needs judging. Keyed (gene, chemical label
# casefolded). A pair whose verdict is that the claim is WRONG belongs in
# high_confidence_g.FALSE_TARGETS instead -- this table is for the ones that stay.
REVIEWED = {
    ("JAK2", "ruxolitinib"): "true -- ruxolitinib is a JAK1/2 inhibitor, only MCG carries it",
    ("JAK1", "ruxolitinib"): "true -- the other half of the JAK1/2 pair",
    ("TYMS", "5-fluorouracil"): "true -- FdUMP inhibits thymidylate synthase",
    ("NFKB1", "bortezomib"): "pathway-level -- the proteasome is the target; NF-kB is downstream",
    ("PSMD5", "bortezomib"): "complex-level -- 19S assembly chaperone, not the catalytic site",
    ("PSMD5", "carfilzomib"): "complex-level -- as above",
    ("PSMD9", "bortezomib"): "complex-level -- 19S assembly chaperone",
    ("PSMD9", "carfilzomib"): "complex-level -- as above",
    ("PSME3", "bortezomib"): "complex-level -- PA28gamma activator, not the catalytic beta subunit",
    ("PSME3", "carfilzomib"): "complex-level -- as above",
    ("FGFR4", "bgj-398"): "true -- infigratinib is pan-FGFR, weaker on FGFR4",
    ("FGFR4", "dovitinib"): "true -- pan-FGFR/VEGFR",
    ("BIRC2", "lcl161"): "true -- LCL161 is a SMAC mimetic; cIAP1 is the primary target",
    ("XIAP", "lcl161"): "true -- same agent, XIAP is the secondary target",
    ("HSPB1", "apatorsen"): "true by mechanism -- antisense to HSPB1; knockdown, not binding",
    ("BIRC5", "sepantronium bromide"): "mechanism -- YM155 suppresses survivin transcription",
    ("NOTCH1", "ro4929097"): "pathway-level -- a gamma-secretase inhibitor; NOTCH1 is substrate",
    ("RICTOR", "dactolisib"): "complex-level -- RICTOR is an mTORC2 subunit, not the binding site",
}


def norm_drug(s):
    """DGIdb drug name -> match key. Mirrors target_pharm._norm_drug."""
    s = (s or "").strip().lower()
    if not s or s == "null" or "+" in s:
        return ""
    s = re.sub(r"\s+\d+[a-z]+$", "", s)
    s = re.sub(r"-[a-z]{4}$", "", s)
    return s.strip()


def load_hgnc(path):
    """(valid symbols, symbol -> entrez id)."""
    docs = json.loads(Path(path).read_text(encoding="utf-8"))["response"]["docs"]
    syms = {d["symbol"] for d in docs if d.get("symbol")}
    entrez = {d["symbol"]: str(d["entrez_id"]) for d in docs
              if d.get("symbol") and d.get("entrez_id")}
    return syms, entrez


def load_dgidb_rows(path, entrez):
    """(drug key, gene) -> rows, plus the gene -> drug keys whose only support is a
    namespace artifact (a GtoPdb NCBIGENE claim equal to the resolved gene's entrez id)."""
    rows = collections.defaultdict(list)
    artifact = collections.defaultdict(set)
    if not Path(path).exists():
        return rows, artifact
    with open(path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            itype = (row.get("interaction_type") or "").strip()
            if not itype or itype.upper() == "NULL":      # the gate load_dgidb applies
                continue
            gene = (row.get("gene_name") or "").strip()
            keys = {norm_drug(row.get("drug_name")), norm_drug(row.get("drug_claim_name"))} - {""}
            for k in keys:
                rows[(k, gene)].append(row)
            if row.get("interaction_source_db_name") == "GuideToPharmacology":
                m = re.match(r"NCBIGENE:(\d+)$", (row.get("gene_claim_name") or "").upper())
                if m and entrez.get(gene) == m.group(1):
                    artifact[gene] |= keys
    return rows, artifact


def corpus_surfaces(data_root):
    """chebi_id -> every corpus surface that normalized onto it."""
    out = collections.defaultdict(set)
    for name in ("chemical.json", "chemical_ambiguous.json"):
        p = Path(data_root) / "CHEMICAL" / name
        if not p.exists():
            continue
        for surface, e in json.loads(p.read_text(encoding="utf-8")).items():
            cid = e["chebi_id"]
            for c in ([cid] if isinstance(cid, str) else cid):
                out[c].add(surface)
    return out


def audit(data_root, quiet=False):
    HC.set_data_root(data_root)
    root = Path(data_root)
    c2t = json.loads(HC.TARGET_FILE.read_text(encoding="utf-8"))
    syms, entrez = load_hgnc(root / "databases" / "hgnc_complete_set_2026-05-01.json")
    rows, artifact = load_dgidb_rows(root / "databases" / "interactions.tsv", entrez)
    surf = corpus_surfaces(root)
    # the cross-link as the GRAPH sees it: in-corpus genes, FALSE_TARGETS already applied
    tgt, chems_by_gene, tsrc, tcat = HC._drug_targets()

    fails = 0
    print(f"cross-link: {len(c2t):,} gene keys, {sum(1 for v in c2t.values() if v.get('in_corpus_GENETIC')):,} "
          f"in corpus; {len(tgt):,} genes carry a target after FALSE_TARGETS")

    # -- 1. gene keys HGNC does not know -------------------------------------------------
    bad = sorted(g for g in c2t if g not in syms)
    print(f"\n[1] non-HGNC gene keys ................. {len(bad)}"
          + (f"  {bad}" if bad else ""))
    fails += bool(bad)

    # -- 2/3/4/5 walk the mappings the graph would draw -----------------------------------
    leaks, nondrug, surf_only, clin_only = [], [], [], []
    for gene, chems in sorted(chems_by_gene.items()):
        entry = c2t.get(gene) or {}
        by_label = {(c.get("chebi_label") or ""): c for c in entry.get("chemicals") or []}
        no_binding = []
        for label in chems:
            c = by_label.get(label)
            if c is None:
                continue
            keys = {norm_drug(label)} | {norm_drug(x) for x in
                                         surf.get(c.get("chebi_id"), set())
                                         | set(c.get("chemical_surfaces") or [])}
            keys -= {""}
            if keys & artifact.get(gene, set()):
                leaks.append((gene, label))
            if NON_DRUG.match(label.strip()):
                nondrug.append((gene, label))
            hit_by = [k for k in keys if rows.get((k, gene))]
            if hit_by and norm_drug(label) not in hit_by:
                surf_only.append((gene, label, sorted(hit_by)))
            srcs = {(r.get("interaction_source_db_name") or "?")
                    for k in keys for r in rows.get((k, gene), [])}
            if srcs and not (srcs & PHARM_SOURCES):
                no_binding.append(label)
        if no_binding and len(no_binding) == len(chems):
            clin_only.append((gene, no_binding))

    print(f"[2] namespace-artifact mappings ........ {len(leaks)}"
          + (f"  {leaks}" if leaks else ""))
    fails += bool(leaks)
    print(f"[3] non-drug chemicals as targets ...... {len(nondrug)}"
          + (f"  {nondrug}" if nondrug else ""))
    fails += bool(nondrug)

    print(f"\n[4] matched through a corpus surface, not the chemical's own label: {len(surf_only)}"
          "  (expected: salt/hydrate and code names)")
    if not quiet:
        for gene, label, keys in surf_only:
            print(f"      {gene:9} <- {label[:46]:46} via {keys}")

    new = [(g, cs) for g, cs in clin_only
           if any((g, c.casefold()) not in REVIEWED for c in cs)]
    print(f"\n[5] targets with no binding-sourced DGIdb row: {len(clin_only)}"
          f"  ({len(clin_only) - len(new)} already reviewed, {len(new)} new)")
    for gene, cs in clin_only:
        for c in cs:
            verdict = REVIEWED.get((gene, c.casefold()), "** NEW -- review this **")
            print(f"      {gene:9} <- {c[:34]:34} {verdict}")

    print(f"\n{'FAIL' if fails else 'PASS'}: checks 1-3 {'found problems' if fails else 'are clean'}")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=str(HC.DATA_ROOT),
                    help="pipeline output tree holding CHEMICAL/ and databases/ "
                         "(default: kaggle_working next to this script)")
    ap.add_argument("--quiet", action="store_true",
                    help="skip the per-mapping listing under check 4")
    args = ap.parse_args()
    return audit(args.data_root, args.quiet)


if __name__ == "__main__":
    sys.exit(main())
