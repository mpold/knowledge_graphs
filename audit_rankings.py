#!/usr/bin/env python3
"""audit_rankings.py -- regression test for rankings.py and the page it writes.

A ranking is the most dangerous artefact this project produces. A wrong edge is one wrong edge; a
wrong ranking is a reading order, and whoever reads it will start at the top. Four of the errors in
candidate_targets.md appendix A were ranking errors -- PRDM14 first on one publication, USP39 and
MASTL nominated as tractable while pan-essential -- and each was fluent enough to survive several
re-readings. This is what catches the next one.

It is the third audit in the repo and deliberately disjoint from the other two: audit_drug_targets.py
checks what colours a node, audit_fusions.py checks the fusion attribute, and this checks what the
score does with them.

WHAT IT CHECKS (1-6 fail the run, 7-9 are drift reports)

  1. staging          The cache rankings.py collects into is present and shaped as the loader
                      expects. Absent cache SKIPS rather than fails -- it is re-fetchable -- but a
                      present cache with missing keys must fail loudly, because every check below
                      would otherwise pass vacuously on empty data.
  2. anchors          Verdicts that are not judgement calls. TP53 must read `lesion` on its mutation
                      rate; HSF1 must read `block lesion`, because its 8q amplification moves with
                      PTP4A3 and calling it a lesion is the error the arm-level test was written to
                      stop; FASN must NOT read `lesion`, or the page contradicts sections 7-9 of
                      candidate_targets.md; RELA must not read `ligand-bound`, MDM2 must; USP39 and
                      MASTL must be unranked as pan-essential.
  3. arithmetic       Every row's score is recomputed from the terms its own tooltip states. This
                      is the check that matters: the tooltip is what a reader trusts, and a score
                      that no longer matches its explanation is worse than a wrong score, because
                      it looks reasoned.
  4. invariants       Ranks are a dense permutation of 1..N in both columns, ordered by descending
                      score; only blue genes carry one; pan-essential genes carry none; no score is
                      negative. And the trap that produced them: any row whose lesion term is zero
                      must have sub-threshold evidence in EVERY component -- an earlier scorer read
                      the fallback prose "max amp 3%, mut 1%, del 0%" and awarded the amplification
                      bonus to genes with no lesion at all, which put H2AX second.
  5. prose            The numbers written into the page's narrative, checked against the cache they
                      came from. Hard-coded prose is how this project's documents go stale: the page
                      claimed 352 blue genes when there were 351, 49 of 60 block calls when the arm
                      test said 35, and CTNNB1's "eight ligand entries" survived the chain-length
                      filter that had reduced them to two.
  6. payload          The rendered HTML: one row per gene, ten cells each, ranks unique, no
                      unsubstituted %(name)s left by a template edit.
  7. distributions    Docking grades, Q1 and Q3 tallies. Drift to read, not an error.
  8. both columns     The top of the candidate ranking and of the neglect ranking, with the docking
                      grade that moved each one.
  9. suppressors      Loss-of-function-shaped genes inside the top 20. The score cannot tell an
                      oncogene from a tumour suppressor, and APC and STK11 rank in the twenties on
                      alterations whose therapy would be restoration, not inhibition.

Run::  python audit_rankings.py [--data-root kaggle_working] [--cache DIR] [--page PATH] [--quiet]
Exit code 0 when checks 1-6 pass (or the cache is absent), 1 otherwise.
"""
import argparse
import collections
import json
import re
import sys
from pathlib import Path

import rankings as R

# --- 2. anchors ------------------------------------------------------------------------------
# (gene, field, predicate, why it is not a judgement call)
ANCHORS = [
    ("TP53", "v1", lambda v: v == "lesion",
     "53% mutated in colorectal, 40x the cohort median"),
    ("KRAS", "v1", lambda v: v == "lesion", "37% mutated in colorectal"),
    ("HSF1", "v1", lambda v: v == "block lesion",
     "8q24.3 amplification shares 95% of its samples with PTP4A3"),
    ("FASN", "v1", lambda v: v != "lesion",
     "6.1% mutated against a 1.3% cohort median is passenger burden, not a lesion"),
    ("H2AX", "v1", lambda v: v == "none found", "no amplification, no outlier mutation, no fusion"),
    ("MDM2", "vd", lambda v: v == "ligand-bound", "80 entries hold a drug-like ligand"),
    ("RELA", "vd", lambda v: v != "ligand-bound",
     "its 22 ligand entries are 14-3-3 sigma holding a RelA phosphopeptide"),
    ("MYB", "vd", lambda v: v == "model only", "no experimental structure of any kind"),
    ("USP39", "rank", lambda v: v is None, "pan-essential: appendix A2"),
    ("MASTL", "rank", lambda v: v is None, "pan-essential in 97% of lines: appendix A2"),
]

# --- 5. prose claims -------------------------------------------------------------------------
# (label, regex over the page's text with one capturing group, value it must equal)
def claims(d, rows):
    arm, pdb, les = d['arm'], d['pdb'], d['lesions']
    blue = [r for r in rows if r['colour'] == 'blue']

    def ov(g, coh):
        return round(arm.get(g, {}).get(coh, {}).get('overlap', 0) * 100)

    def mut(g, coh='colorectal'):
        # the sentence quotes one cohort, not each gene's best -- TP53 is 53% in colorectal and
        # 58% in lung adeno, and reading the wrong one turns a correct page into a failing check
        return round(les.get(g, {}).get(coh, {}).get('mut', 0.0))

    return [
        ("gene nodes", r"All (\d+) gene nodes", len(rows)),
        ("blue genes", r"orders the (\d+) undrugged", len(blue)),
        ("disqualified", r"(\d+) blue genes are disqualified", sum(1 for r in blue if not r['rank'])),
        ("ranked", r"leaving (\d+) ranked", sum(1 for r in blue if r['rank'])),
        ("block calls", r"(\d+) of \d+ amplification calls are block",
         sum(1 for r in rows if r['ev']['block'])),
        ("amplification calls", r"\d+ of (\d+) amplification calls are block",
         sum(1 for r in rows if r['ev']['amp'] >= R.AMP_CALL)),
        ("Q2 assessed", r"(\d+) of \d+ genes.{0,40}carry a verdict",
         sum(1 for r in rows if r['v2'] == 'assessed')),
        ("HSF1/PTP4A3 overlap", r"HSF1 \(8q24\.3\) shares (\d+)% of its amplified", ov('HSF1', 'breast')),
        ("PVT1/MYC overlap", r"PVT1 shares (\d+)% with MYC", ov('PVT1', 'breast')),
        ("CDK12/ERBB2 overlap", r"CDK12 (\d+)% with ERBB2", ov('CDK12', 'breast')),
        ("RELA peptide entries", r"RELA scored (\d+) ligand-bound entries", pdb['RELA']['pep']),
        ("RELA after the filter", r"takes RELA to (\d+)", pdb['RELA']['lig']),
        ("MDM2 ligand entries", r"leaves MDM2 \((\d+)\)", pdb['MDM2']['lig']),
        ("MCL1 ligand entries", r"MCL1 \((\d+)\)", pdb['MCL1']['lig']),
        ("BCL6 ligand entries", r"BCL6 \((\d+)\) untouched", pdb['BCL6']['lig']),
        ("FASN colorectal mutation", r"FASN reads ([\d.]+)% in colorectal",
         les['FASN']['colorectal']['mut']),
        ("NBEA mutation", r"NBEA ([\d.]+)%", les['NBEA']['colorectal']['mut']),
        ("colorectal median", r"cohort median of ([\d.]+)%", round(d['medmut']['colorectal'], 1)),
        ("TP53 mutation", r"TP53 \((\d+)%\)", mut('TP53')),
        ("APC mutation", r"APC \((\d+)%\)", mut('APC')),
        ("KRAS mutation", r"KRAS \((\d+)%\)", mut('KRAS')),
        ("DepMap cell lines", r"across ([\d,]+) cell lines", f"{d['n_lines']:,}"),
    ]


def page_text(path):
    """Tags out, whitespace collapsed, entities that carry numbers restored."""
    s = Path(path).read_text(encoding="utf-8")
    s = s.split('<table id="t">')[0]                 # narrative only, never the rows
    s = re.sub(r"<[^>]+>", " ", s)
    s = s.replace("&mdash;", "-").replace("&ndash;", "-").replace("&ge;", ">=")
    return re.sub(r"\s+", " ", s)


def tooltip_math(why):
    """Recompute a score from the terms its own tooltip states: (a +1; b +0.6; ...) x 0.8 for ...

    Returns None when the tooltip is not an arithmetic one (pan-essential, or a drugged gene)."""
    m = re.match(r"\((.*)\) x ([\d.]+) for ", why or "")
    if not m:
        return None
    terms = re.findall(r"([+-][\d.]+)(?:;|$)", m.group(1) + ";")
    return round(max(0.0, sum(float(t) for t in terms)) * float(m.group(2)), 2)


def audit(data_root, cache, page, quiet=False):
    fails = 0

    # -- 1. staging ---------------------------------------------------------------------------
    need = {"gene_entrez": (), "lesions": ("breast",), "armlevel": (), "pdb": ("acc", "pdb", "lig")}
    missing = [n for n in need if not (Path(cache) / f"{n}.json").exists()]
    if len(missing) == len(need):
        print(f"[1] no cache under {cache} -- it is re-fetchable, so nothing to audit.\n"
              f"    Build it with: python rankings.py")
        return 0
    if missing:
        print(f"[1] staging ............................ MISSING {missing}")
        print("\nFAIL: a partial cache would make every check below pass on empty data")
        return 1
    d = R.load_data(data_root=data_root, cache=cache, offline=True, quiet=True)
    shape = []
    for g, v in list(d['pdb'].items())[:50]:
        if any(k not in v for k in need['pdb']):
            shape.append(g)
    print(f"[1] staging ............................ {len(d['lesions'])} genes with lesion data, "
          f"{len(d['arm'])} amplification tests, {len(d['pdb'])} PDB records"
          + (f"  MALFORMED: {shape[:5]}" if shape else ""))
    fails += bool(shape)

    rows = R.build_rows(d)
    by = {r['g']: r for r in rows}

    # -- 2. anchors ---------------------------------------------------------------------------
    bad = [(g, f, by.get(g, {}).get(f), why) for g, f, ok, why in ANCHORS
           if g in by and not ok(by[g].get(f))]
    absent = [g for g, *_ in ANCHORS if g not in by]
    print(f"[2] anchor verdicts .................... {len(ANCHORS) - len(bad) - len(absent)}"
          f"/{len(ANCHORS)} as expected"
          + (f"  ABSENT FROM GRAPH: {absent}" if absent else ""))
    for g, f, got, why in bad:
        print(f"      {g} {f} = {got!r} -- {why}")
    fails += bool(bad)

    # -- 3. the score matches its own explanation ---------------------------------------------
    off = []
    for r in rows:
        for sk, wk in (('score', 'why'), ('nscore', 'nwhy')):
            want = tooltip_math(r[wk])
            if want is not None and r[sk] is not None and abs(want - r[sk]) > 0.011:
                off.append((r['g'], sk, r[sk], want))
    print(f"[3] score = the sum its tooltip states . {2 * len(rows) - len(off)}/{2 * len(rows)}"
          + ("" if not off else "  MISMATCH: " + ", ".join(
              f"{g}.{k} says {w} computes {s}" for g, k, s, w in off[:4])))
    fails += bool(off)

    # -- 4. invariants ------------------------------------------------------------------------
    problems = []
    blue = [r for r in rows if r['colour'] == 'blue']
    for key, rk in (('score', 'rank'), ('nscore', 'nrank')):
        have = sorted(r[rk] for r in rows if r[rk])
        if have != list(range(1, len(have) + 1)):
            problems.append(f"{rk} is not a dense 1..N permutation")
        order = sorted([r for r in rows if r[rk]], key=lambda r: r[rk])
        if any(a[key] < b[key] for a, b in zip(order, order[1:])):
            problems.append(f"{rk} is not ordered by descending {key}")
        if any(r['colour'] != 'blue' for r in order):
            problems.append(f"a non-blue gene carries a {rk}")
        if any((r[key] or 0) < 0 for r in rows):
            problems.append(f"{key} went negative")
    pan = [r['g'] for r in blue if r['dep'] == 'common' and r['rank']]
    if pan:
        problems.append(f"pan-essential genes carry a rank: {pan[:5]}")
    # the fallback-prose trap: no lesion term means no lesion evidence, in every component
    leak = [r['g'] for r in blue if r['score'] is not None
            and 'no lesion found' in (r['why'] or '')
            and (r['ev']['amp'] >= R.AMP_CALL or r['ev']['mut'] or r['ev']['del'] >= R.DEL_CALL
                 or r['ev']['fn'])]
    if leak:
        problems.append(f"scored 'no lesion found' but has lesion evidence: {leak[:5]}")
    silent = [r['g'] for r in blue if r['score'] is not None
              and 'no lesion found' not in (r['why'] or '')
              and not (r['ev']['amp'] >= R.AMP_CALL or r['ev']['mut']
                       or r['ev']['del'] >= R.DEL_CALL or r['ev']['fn'])]
    if silent:
        problems.append(f"scored a lesion term with no lesion evidence: {silent[:5]}")
    print(f"[4] ranking invariants ................. {'all hold' if not problems else 'BROKEN'}")
    for p in problems:
        print(f"      {p}")
    fails += bool(problems)

    # -- 5. the prose against the data --------------------------------------------------------
    if page and Path(page).exists():
        text = page_text(page)
        wrong, unfound = [], []
        for label, rx, want in claims(d, rows):
            m = re.search(rx, text)
            if not m:
                unfound.append(label)
            elif m.group(1).rstrip('0').rstrip('.') != str(want).rstrip('0').rstrip('.'):
                wrong.append((label, m.group(1), want))
        n = len(claims(d, rows))
        print(f"[5] prose claims vs the cache .......... {n - len(wrong) - len(unfound)}/{n} verified"
              + (f"  UNREADABLE: {unfound}" if unfound else ""))
        for label, got, want in wrong:
            print(f"      {label}: page says {got}, data says {want}")
        fails += bool(wrong or unfound)
    else:
        print("[5] prose claims ....................... no page found, skipped")
        text = None

    # -- 6. payload -------------------------------------------------------------------------
    if page and Path(page).exists():
        s = Path(page).read_text(encoding="utf-8")
        tr = re.findall(r"<tr data-g=.*?</tr>", s, re.S)
        cells = collections.Counter(len(re.findall(r"<td", t)) for t in tr)
        heads = len(re.findall(r"<th[ >]", s.split("<tbody>")[0]))
        left = re.findall(r"%\((\w+)\)s", s)
        ranks = [int(m) for m in re.findall(r'class="num rk" data-sort="(\d+)"', s) if m != '99999']
        bad6 = []
        if len(tr) != len(rows):
            bad6.append(f"{len(tr)} rows rendered for {len(rows)} genes")
        if list(cells) != [heads]:
            bad6.append(f"cells per row {dict(cells)} against {heads} headers")
        if left:
            bad6.append(f"unsubstituted placeholders: {sorted(set(left))}")
        if len(set(ranks)) != len(ranks):
            bad6.append("duplicate ranks in the rendered table")
        print(f"[6] payload ............................ {len(tr)} rows x {heads} columns, "
              f"{len(ranks)} ranked" + ("" if not bad6 else "  " + "; ".join(bad6)))
        fails += bool(bad6)
    else:
        print("[6] payload ............................ no page found, skipped")

    # -- 7-9. drift ---------------------------------------------------------------------------
    print()
    dk = collections.Counter(r['vd'] for r in rows)
    print("[7] docking grades ..... " + ", ".join(f"{k} {v}" for k, v in dk.most_common()))
    print("    Q1 ................. " + ", ".join(
        f"{k} {v}" for k, v in collections.Counter(r['v1'] for r in rows).most_common()))
    print("    Q3 ................. " + ", ".join(
        f"{k} {v}" for k, v in collections.Counter(r['v3'] for r in rows).most_common()))

    if not quiet:
        top = sorted([r for r in rows if r['rank']], key=lambda r: r['rank'])[:10]
        neg = sorted([r for r in rows if r['nrank']], key=lambda r: r['nrank'])[:10]
        print("\n[8] best candidates      " + ", ".join(f"{r['g']} {r['score']:g}" for r in top))
        print("    most neglected       " + ", ".join(f"{r['g']} {r['nscore']:g}" for r in neg))
        print("    in both top tens:    " + ", ".join(sorted(
            {r['g'] for r in top} & {r['g'] for r in neg}) or ["none"]))

        lof = [r for r in rows if r.get('lof') and r['rank'] and r['rank'] <= 20]
        print(f"\n[9] loss-of-function shape inside the top 20: {len(lof)}")
        for r in sorted(lof, key=lambda r: r['rank']):
            print(f"      #{r['rank']:<3} {r['g']:9} {r['d1'][:64]}")

    print(f"\n{'FAIL' if fails else 'PASS'}: checks 1-6 "
          f"{'found problems' if fails else 'are clean'}")
    return 1 if fails else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=str(R.DEFAULT_ROOT),
                    help="pipeline output tree holding databases/ (default: kaggle_working)")
    ap.add_argument("--cache", default=str(R.DEFAULT_CACHE),
                    help="the cache rankings.py collects into (default: rankings_cache)")
    ap.add_argument("--page", default=str(R.DEFAULT_OUT),
                    help="the rendered HTML to check (default: ../Q1_3.html)")
    ap.add_argument("--quiet", action="store_true", help="counts only, no listings")
    args = ap.parse_args()
    return audit(args.data_root, args.cache, args.page, args.quiet)


if __name__ == "__main__":
    sys.exit(main())
