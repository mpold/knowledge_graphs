#!/usr/bin/env python3
"""pubmed_fusions.py -- ask PubMed whether a gene forms FUSIONS or REARRANGEMENTS.

The graph says which genes a tumour depends on and which carry drugs; it says nothing about
whether a gene is a fusion partner, which changes what "undruggable" means. A fusion junction
exists in no normal cell, so it is a selectivity handle even for a protein with no pocket --
which is why HMGA2, MEIS1, HOXA9 and PAX3 read differently from PRDM14 despite sitting in the
same "no chemistry" list. This script is the reproducible version of that lookup.

TWO QUERIES

  tier 1 (default)  Does this gene form fusions at all?

      ("GENE"[tiab])
      AND ("Oncogene Proteins, Fusion"[MeSH] OR "Gene Fusion"[MeSH]
           OR "Translocation, Genetic"[MeSH] OR "Gene Rearrangement"[MeSH]
           OR fusion[tiab] OR translocation[tiab] OR rearrangement[tiab] OR chimeric[tiab])
      NOT ("Spinal Fusion"[MeSH] OR "Cell Fusion"[MeSH] OR "Membrane Fusion"[MeSH])

  tier 2 (--pair)   Is THIS pair a described fusion?

      "G1-G2"[tiab] OR "G2-G1"[tiab]
      OR ("G1"[tiab] AND "G2"[tiab] AND (fusion[tiab] OR translocation[tiab] OR chimeric[tiab]))

WHY IT IS SHAPED THIS WAY

  [tiab], not [All Fields]   otherwise affiliations, grant text and keyword soup match.
  MeSH OR'd with free text   MeSH is precise but lags 6-18 months, so free text carries the
                             recent literature and MeSH carries the well-indexed older work.
  the NOT block              "fusion" in PubMed is dominated by SPINAL fusion, CELL fusion,
                             MEMBRANE fusion and image fusion. Without it any gene expressed in
                             bone or studied in virology looks like a fusion partner.
  the third tier-2 clause    most case reports name both partners in prose and never hyphenate.

READ THE SHARE, NOT THE COUNT

  A raw count tracks how well studied a gene is, not whether it fuses -- TP53 scores 1766 with no
  fusions at all, above HMGA2's 317, which has several. So each gene is also queried bare and the
  fusion count reported as a FRACTION of its whole literature. Measured 2026-08-20:

      EWSR1  1733/2089  83.0%     the fusion gene par excellence
      PAX3    626/2226  28.1%     PAX3::FOXO1, ::NCOA1/2, ::MAML3
      HOXA9   301/1139  26.4%     NUP98::HOXA9
      MEIS1   218/ 890  24.5%     MEIS1::NCOA2
      HMGA2   317/1987  16.0%     ::LPP, ::NFIB, ::RAD51B, ::WIF1
      TP53   1766/31043  5.7%     no fusions -- this is the background level
      CTNNB1  242/4535   5.3%     no fusions
      PRDM16   47/1022   4.6%     HAS fusions, and still reads as background (see below)
      PRDM14    9/ 205   4.4%     no fusions

  Above ~15% means look; at background, look anyway if you have a reason. PRDM16 is the warning:
  RPN1::, ETV6::, CBFA2T3:: and NUP98::PRDM16 are all real, but a thousand papers on its role in
  brown-fat thermogenesis dilute the share to noise. A gene famous for something ELSE will hide
  its fusions here, which is why this is a screen and --pair is the verdict.

FOUR TRAPS, ALL MET WHILE BUILDING THIS

  * `::` breaks the request. A long OR string containing "%3A%3A" returns HTTP 500 from
    E-utilities, and PubMed does not index the HGVS-style notation anyway. Hyphens only.
  * Hyphens are tokenized: "PAX3-FOXO1" is searched as the phrase `pax3 foxo1`, so word ORDER
    is what matters, not punctuation. Both orders are cheap and authors are inconsistent.
  * Legacy aliases are mandatory or the older, often foundational literature is missed:
    PAX3-FKHR (=FOXO1), C11orf95-RELA (=ZFTA-RELA), EVI1 (=MECOM), MTGR1 (=CBFA2T2).
    Pass them with --alias; a few well-known ones are built in.
  * Residual noise survives: "phage FUSION PROTEIN-targeted liposomes" matched PRDM14 and no
    MeSH exclusion removes that sense cleanly. --strict drops the free-text arm for near-perfect
    precision at a real cost in recall.

A COUNT IS PUBLICATIONS, NOT CASES. PubMed tells you a fusion exists and is studied, never how
often it occurs -- for frequency use Mitelman (`mitelman-db.prod.MolClinGene` in BigQuery, which
is case-level) or COSMIC.

Run::
    python pubmed_fusions.py HMGA2 MEIS1 HOXA9 PAX3
    python pubmed_fusions.py --pair PAX3 FOXO1 --alias FKHR
    python pubmed_fusions.py --from-graph oncogene_addiction_2026_08_19_M.html --blue --max 40
    python pubmed_fusions.py PRDM14 --titles --strict

NCBI allows 3 requests/second without a key and 10 with one; set NCBI_API_KEY (or --api-key)
and the pacing follows suit.
"""
import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

FUSION_MESH = ('"Oncogene Proteins, Fusion"[MeSH] OR "Gene Fusion"[MeSH] '
               'OR "Translocation, Genetic"[MeSH] OR "Gene Rearrangement"[MeSH]')
FUSION_TIAB = ('fusion[tiab] OR translocation[tiab] OR rearrangement[tiab] OR chimeric[tiab]')
# the senses of "fusion" that have nothing to do with genes
NOT_BLOCK = '"Spinal Fusion"[MeSH] OR "Cell Fusion"[MeSH] OR "Membrane Fusion"[MeSH]'

# legacy names that still carry a chunk of each fusion's literature
BUILTIN_ALIASES = {
    "FOXO1": ["FKHR"],
    "RELA": ["C11orf95"],
    "MECOM": ["EVI1"],
    "CBFA2T2": ["MTGR1"],
    "RUNX1T1": ["ETO"],
    "KMT2A": ["MLL"],
    "NCOA2": ["TIF2"],
}


def aliases_for(gene, extra):
    return [gene] + BUILTIN_ALIASES.get(gene.upper(), []) + list(extra or [])


def tier1_term(gene, extra=None, strict=False):
    names = " OR ".join(f'"{n}"[tiab]' for n in aliases_for(gene, extra))
    concept = FUSION_MESH if strict else f"{FUSION_MESH} OR {FUSION_TIAB}"
    return f"({names}) AND ({concept}) NOT ({NOT_BLOCK})"


def baseline_term(gene, extra=None):
    """Everything PubMed holds on the gene -- the denominator for the fusion count."""
    return "(%s)" % _tiab(aliases_for(gene, extra))


def _tiab(names):
    return " OR ".join('"%s"[tiab]' % n for n in names)


def tier2_term(g1, g2, extra=None):
    a1, a2 = aliases_for(g1, extra), aliases_for(g2, extra)
    pairs = ['"%s-%s"[tiab]' % (x, y) for x in a1 for y in a2]
    pairs += ['"%s-%s"[tiab]' % (y, x) for x in a1 for y in a2]
    # the prose clause: both symbols plus a fusion word, for reports that never hyphenate
    prose = ("((%s) AND (%s) AND (fusion[tiab] OR translocation[tiab] OR chimeric[tiab]))"
             % (_tiab(a1), _tiab(a2)))
    return " OR ".join(pairs) + " OR " + prose


def eutils(endpoint, params, api_key=None, pause=0.34):
    if api_key:
        params = dict(params, api_key=api_key)
    url = f"{EUTILS}/{endpoint}.fcgi?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read().decode("utf-8"))
    except Exception as e:                       # a 500 here is usually an over-long term
        print(f"  ! {endpoint} failed ({type(e).__name__}: {e})", file=sys.stderr)
        data = {}
    time.sleep(pause)
    return data


def search(term, retmax, api_key, pause):
    d = eutils("esearch", {"db": "pubmed", "retmode": "json", "retmax": retmax,
                           "term": term}, api_key, pause)
    res = d.get("esearchresult", {})
    return int(res.get("count", 0)), res.get("idlist", [])


def titles(pmids, api_key, pause):
    if not pmids:
        return []
    d = eutils("esummary", {"db": "pubmed", "retmode": "json",
                            "id": ",".join(pmids)}, api_key, pause)
    r = d.get("result", {})
    return [(p, (r.get(p) or {}).get("title", ""),
             (r.get(p) or {}).get("pubdate", "")[:4]) for p in pmids]


def blue_genes(path):
    """Gene nodes the graph draws as undrugged (target == 0), read straight from the payload."""
    s = open(path, encoding="utf-8").read()
    i = s.index("const DATA=")
    j = s.index("\n", i)
    payload = json.loads(s[i + len("const DATA="):j].rstrip().rstrip(";"))
    return sorted(n["id"] for n in payload["nodes"]
                  if n.get("kind") == "gene" and not n.get("target"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genes", nargs="*", help="gene symbols to screen (tier 1)")
    ap.add_argument("--pair", nargs=2, metavar=("GENE1", "GENE2"),
                    help="ask whether this specific pair is a described fusion (tier 2)")
    ap.add_argument("--from-graph", metavar="HTML",
                    help="read gene symbols out of a graph payload instead of the command line")
    ap.add_argument("--blue", action="store_true",
                    help="with --from-graph: only the undrugged (blue) gene nodes")
    ap.add_argument("--max", type=int, default=25,
                    help="cap on how many genes to screen from --from-graph (default 25); "
                         "the whole blue set is 352 queries, which is ~2 minutes unkeyed")
    ap.add_argument("--alias", action="append", default=[],
                    help="extra legacy name to OR in (repeatable), e.g. --alias FKHR")
    ap.add_argument("--strict", action="store_true",
                    help="MeSH arm only: near-perfect precision, real loss of recall")
    ap.add_argument("--titles", action="store_true", help="print the top PMIDs with titles")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of a table")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                    help="NCBI API key (or set NCBI_API_KEY); raises the rate limit to 10/s")
    args = ap.parse_args()
    pause = 0.11 if args.api_key else 0.34       # NCBI: 10/s with a key, 3/s without

    if args.pair:
        g1, g2 = args.pair
        term = tier2_term(g1, g2, args.alias)
        n, ids = search(term, 10 if args.titles else 0, args.api_key, pause)
        print(f"{g1}::{g2}  {n} publications")
        for p, t, y in titles(ids, args.api_key, pause) if args.titles else []:
            print(f"    {p}  {y}  {t[:96]}")
        return 0

    genes = list(args.genes)
    if args.from_graph:
        g = blue_genes(args.from_graph)
        if not args.blue:
            print("(--from-graph currently reads the blue set; pass --blue to make that explicit)",
                  file=sys.stderr)
        if len(g) > args.max:
            print(f"note: {len(g)} genes available, screening the first {args.max} "
                  f"(raise with --max)", file=sys.stderr)
        genes += g[:args.max]
    if not genes:
        ap.error("give gene symbols, --pair, or --from-graph")

    out = []
    if not args.json:
        print(f"{'gene':12}{'fusion':>7}{'total':>8}{'share':>8}   (share = fusion-query hits"
              f" as a fraction of all papers on the gene)")
    for gene in genes:
        n, ids = search(tier1_term(gene, args.alias, args.strict),
                        5 if args.titles else 0, args.api_key, pause)
        tot, _ = search(baseline_term(gene, args.alias), 0, args.api_key, pause)
        share = (n / tot) if tot else 0.0
        rec = {"gene": gene, "fusion": n, "total": tot, "share": round(share, 4)}
        if args.titles:
            rec["top"] = [{"pmid": p, "year": y, "title": t}
                          for p, t, y in titles(ids, args.api_key, pause)]
        out.append(rec)
        if not args.json:
            print(f"{gene:12}{n:>7}{tot:>8}{share*100:>7.1f}%")
            for e in rec.get("top", []):
                print(f"    {e['pmid']}  {e['year']}  {e['title'][:96]}")
    if args.json:
        print(json.dumps(out, indent=2))
    if not args.json and len(out) > 1:
        print()
        print("Read the SHARE, not the count: a raw count tracks how well studied a gene is, not whether it fuses.")
        print("TP53 scores 1766 fusion-query hits with no fusions at all, above HMGA2's 317, which has several.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
