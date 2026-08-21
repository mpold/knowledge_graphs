#!/usr/bin/env python3
"""rankings.py -- put every gene in the graph through the three-question test and rank the
undrugged ones. Writes Q1_3.html.

THE THREE QUESTIONS. A dependency reported in a paper is not yet a target. Three things have to
hold, and they fail independently:

    Q1  Is there a LESION that predicts the dependency?   -> can the trial be designed?
    Q2  Can the cell obtain the product another way?      -> does inhibition achieve anything?
    Q3  Does normal tissue need the same enzyme?          -> is there a therapeutic window?

The point of the page is that these are not equally answerable from staged data, and pretending
otherwise is how a target list fills with fluent nonsense. Q1 comes from cBioPortal plus the graph's
own fusion attribute, Q3 from DepMap, and **Q2 from nothing** -- bypass needs isotope tracing or
rescue experiments, so 415 of 424 genes read "not assessed" rather than carrying an invented answer.
A fourth column, docking readiness, is not one of the questions but decides whether any answer can
be acted on with a small molecule.

WHAT IT PRODUCES

    Q1_3.html   one row per gene node: the four verdicts with their evidence, a candidate rank and
                a neglect rank over the blue (undrugged) genes, sortable, filterable.

TWO RANKINGS, ONE FORMULA

    score = (strongest lesion + publications + addiction class + docking) x dependency

Dependency MULTIPLIES because it is a gate: without something to inhibit, a lesion and a literature
are worth nothing -- TP53 has 34 publications and the clearest lesion in cancer, and 0.4% of lines
depend on it. Docking ADDS because a pocket is not a gate: an apo structure is a fine docking input
and an unsolved protein can still be drugged once somebody tries.

The two columns differ only in the docking table. `candidate` rewards an occupied pocket
(ligand-bound +1.5 ... model only 0) and answers "what is most actionable". `neglect` reverses it
and centres on zero (model only +1.5 ... ligand-bound -1.5) and answers "where has the chemistry
never been tried". Centring matters: with a merely reordered table MDM2 topped both columns, since
its lesion and its nine publications swamp everything else -- and MDM2 has nutlins in trials, so
calling it the field's most overlooked target is absurd. A filled pocket has to cost, not merely
fail to earn. The publication term is deliberately NOT inverted in the neglect column: publications
measure attention too, but flipping them floats the single-paper curiosities to the top, which is
the PRDM14 error this project already made once.

FIVE TRAPS THIS SCRIPT EXISTS TO AVOID. Each cost a wrong answer before it was found:

  1. Flat mutation cuts call every LARGE gene mutated, because passenger burden scales with coding
     length. FASN reads 6.1% in colorectal and NBEA 10.6% against a cohort median of 1.3%. A
     mutation counts only as a cohort OUTLIER -- >=3% and >=10x that cohort's median across the
     genes in the graph. TP53 (53%), APC (65%) and KRAS (37%) clear it by an order of magnitude.
  2. Amplifications move in BLOCKS. If >=80% of a gene's amplified samples also amplify a larger
     same-arm neighbour, the segment is what is selected and which gene inside it drives the
     selection is not resolvable from copy number. HSF1 shares 95% of its amplified breast samples
     with PTP4A3; before this test it was called a lesion on 10.6% amplification.
  3. PDB ligand counts are inflated by PEPTIDES. RELA scored 22 ligand-bound entries and every one
     was a 14-3-3 sigma structure holding a RelA phosphopeptide -- the pocket was someone else's.
     The entry must hold the protein as a chain of >=50 residues; that takes RELA to 0 and leaves
     MDM2, MCL1 and BCL6 untouched.
  4. Scoring by matching strings in the evidence prose. An earlier version read the fallback line
     "max amp 3%, mut 1%, del 0%" and awarded lesion-free genes the amplification bonus, which put
     H2AX second. q1() returns structured evidence and the scorer reads numbers, never text.
  5. One publication is not corroboration. It earns nothing in either column.

STAGING. Everything except the two web APIs is read from <data-root>/databases:

    depmap_dependency.tsv                 built by depmap_to_tsv.py
    CRISPRInferredCommonEssentials.csv    DepMap, staged as-is
    hgnc_complete_set_*.json              symbol -> Entrez id and UniProt accession

and the graph HTML itself supplies the gene list, the addiction class, the fusion attribute and the
per-node publication counts.

NETWORK. Two APIs, both anonymous and both slow: cBioPortal (three TCGA PanCancer cohorts, copy
number and mutations) and RCSB PDB (two searches per gene). A full collection is 15-25 minutes, so
every stage caches its JSON under --cache and is skipped when the cache is present. Use --refresh to
force one back, or --offline to fail loudly instead of reaching for the network.

Run::
    python rankings.py                                  # cached where possible
    python rankings.py --refresh pdb                    # re-collect one stage
    python rankings.py --refresh all --out Q1_3.html
    python rankings.py --offline --graph oncogene_addiction_2026_08_20_M.html
"""
import argparse
import collections
import csv
import html
import json
import re
import statistics
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_ROOT = ROOT / "kaggle_working"
DEFAULT_CACHE = ROOT / "rankings_cache"
DEFAULT_OUT = ROOT.parent / "Q1_3.html"

CBIO = "https://www.cbioportal.org/api"
RCSB = "https://search.rcsb.org/rcsbsearch/v2/query?json="
# Three PanCancer Atlas cohorts, not pan-cancer. A gene lesioned only in leukaemia or sarcoma reads
# "none found" here -- NPM1, RUNX1 and CRLF2 are exactly that case, and the page says so.
STUDIES = [("brca_tcga_pan_can_atlas_2018", "breast"),
           ("coadread_tcga_pan_can_atlas_2018", "colorectal"),
           ("luad_tcga_pan_can_atlas_2018", "lung adeno")]
COHORT_LABEL = {"breast": "breast", "colorectal": "colorectal", "lung adeno": "lung adenocarcinoma"}

# Thresholds. Kept together because every one of them was moved at least once, and moving one
# silently is how the page starts disagreeing with candidate_targets.md.
AMP_CALL = 5.0        # % of a cohort with a GISTIC amplification before it counts at all
DEL_CALL = 5.0
MUT_FLOOR = 3.0       # a mutation must clear this AND 10x the cohort median
MUT_MULT = 10.0
WEAK = 2.0            # below the call thresholds but above this reads "weak"
BLOCK_OVERLAP = 0.80  # share of amplified samples shared with a bigger same-arm neighbour
MIN_CHAIN = 50        # residues: below this the "protein" in a PDB entry is a peptide of it
LIGAND_DA = 300       # a ligand this heavy is drug-like enough to prove a pocket was occupied

DK_ACT = {'ligand-bound': 1.5, 'few ligands': 0.75, 'apo only': 0.25, 'model only': 0.0,
          'unknown': 0.0}
DK_NEG = {'model only': 1.5, 'apo only': 0.5, 'few ligands': -0.5, 'ligand-bound': -1.5,
          'unknown': 0.0}

# The nine genes examined by hand in candidate_targets.md sections 7-9. Q2 has no database, so this
# is the whole of what can honestly be said about bypass; every other gene reads "not assessed".
HAND_Q2 = {
    'FASN': 'yes - lipid uptake (LPL co-silencing in this corpus)',
    'GLS': 'yes - glutamine scavenging, macropinocytosis',
    'SLC2A1': 'yes - other transporters, contextual',
    'PKM': 'yes - PKM1 isoform',
    'HK1': 'yes - HK2 covers it',
    'ACACA': 'yes - lipid uptake, as FASN',
    'HSPD1': 'NO - no second mitochondrial chaperonin',
    'HSF1': 'partial - HSF2 does not compensate',
    'XBP1': 'yes - ATF6 and PERK arms',
}


# ---------------------------------------------------------------- staging ---------------------
def newest_graph(root):
    """The dated copy stage 3 writes; the name carries the run date, so take the newest."""
    hits = sorted(Path(root).glob("*_M.html"), key=lambda p: p.stat().st_mtime, reverse=True)
    return hits[0] if hits else None


def read_payload(path):
    """The graph embeds its nodes and edges as one JSON line: `const DATA={...};`"""
    s = Path(path).read_text(encoding="utf-8")
    i = s.index("const DATA=")
    j = s.index("\n", i)
    return json.loads(s[i + len("const DATA="):j].rstrip().rstrip(";"))


def hgnc_docs(db):
    hits = sorted(Path(db).glob("hgnc_complete_set_*.json"))
    if not hits:
        raise SystemExit(f"no hgnc_complete_set_*.json in {db} -- see PLACE_DATABASES_HERE.md")
    return json.loads(hits[-1].read_text(encoding="utf-8"))["response"]["docs"]


def cached(cache, name, refresh, offline, build, note=""):
    """Stage results live as JSON under --cache. A stage runs only when its file is missing or
    named in --refresh, because a full collection is 15-25 minutes of somebody else's API."""
    p = Path(cache) / f"{name}.json"
    if p.exists() and not refresh:
        return json.loads(p.read_text(encoding="utf-8"))
    if offline:
        raise SystemExit(f"--offline and {p} is absent; run without --offline to collect it")
    print(f"[collect] {name}{note}", file=sys.stderr)
    data = build()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data), encoding="utf-8")
    return data


# ---------------------------------------------------------------- collectors ------------------
def post(url, body, tries=3, timeout=300):
    req = urllib.request.Request(url, data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"}, method="POST")
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except Exception as e:                                   # noqa: BLE001 - report and go on
            if attempt == tries - 1:
                print(f"  ! {type(e).__name__} {e}", file=sys.stderr)
                return []
            time.sleep(5)


def get_json(url, tries=3, timeout=120):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return {} if r.status == 204 else json.loads(r.read().decode())
        except Exception:                                        # noqa: BLE001
            if attempt == tries - 1:
                return None
            time.sleep(4)


def collect_lesions(entrez, quiet=False):
    """Per gene and cohort: % of samples amplified, deeply deleted, and mutated.

    Mutations are counted per SAMPLE, not per row -- a gene with three mutations in one tumour is
    one mutated sample. The cohort size is taken as the largest per-gene CNA row count rather than
    from a sample-list call, which is one fewer request per batch and agrees with the study page."""
    ids = sorted(set(entrez.values()))
    amp, dele, mut, nsam = collections.Counter(), collections.Counter(), collections.Counter(), {}
    for sid, lab in STUDIES:
        n = None
        for k in range(0, len(ids), 100):
            batch = ids[k:k + 100]
            cna = post(f"{CBIO}/molecular-profiles/{sid}_gistic/molecular-data/fetch"
                       f"?projection=SUMMARY",
                       {"entrezGeneIds": batch, "sampleListId": f"{sid}_all"})
            per = collections.Counter()
            for r in cna:
                per[r['entrezGeneId']] += 1
                if r.get('value') == 2:
                    amp[(sid, r['entrezGeneId'])] += 1
                if r.get('value') == -2:
                    dele[(sid, r['entrezGeneId'])] += 1
            if per:
                n = max(per.values())
            mu = post(f"{CBIO}/molecular-profiles/{sid}_mutations/mutations/fetch"
                      f"?projection=SUMMARY",
                      {"entrezGeneIds": batch, "sampleListId": f"{sid}_all"})
            seen = collections.defaultdict(set)
            for r in mu:
                seen[r['entrezGeneId']].add(r['sampleId'])
            for g, ss in seen.items():
                mut[(sid, g)] += len(ss)
            if not quiet:
                print(f"    {lab}: {min(k + 100, len(ids))}/{len(ids)} genes", file=sys.stderr)
        nsam[sid] = n
    out = {}
    for g, e in entrez.items():
        rec = {}
        for sid, lab in STUDIES:
            n = nsam.get(sid) or 0
            if not n:
                continue
            rec[lab] = {"n": n,
                        "amp": round(100 * amp[(sid, e)] / n, 1),
                        "del": round(100 * dele[(sid, e)] / n, 1),
                        "mut": round(100 * mut[(sid, e)] / n, 1)}
        out[g] = rec
    return out


def collect_armlevel(entrez, lesions, arm_of, quiet=False):
    """Is an amplification the gene's own, or the arm's?

    For every gene called amplified, compare its amplified SAMPLE SET with that of each larger
    amplified gene on the same arm. High overlap means the segment is what is selected. Only a
    BIGGER neighbour can explain a gene's amplification -- otherwise MYC would be explained by every
    passenger it drags along, and the test would report the driver as the passenger."""
    cand = {g for g, r in lesions.items() if any(v['amp'] >= AMP_CALL for v in r.values())}
    ids = sorted({entrez[g] for g in cand if g in entrez})
    rev = {v: k for k, v in entrez.items()}
    out = {}
    for sid, lab in STUDIES:
        ampsam = collections.defaultdict(set)
        for k in range(0, len(ids), 100):
            for r in post(f"{CBIO}/molecular-profiles/{sid}_gistic/molecular-data/fetch"
                          f"?projection=SUMMARY",
                          {"entrezGeneIds": ids[k:k + 100], "sampleListId": f"{sid}_all"}):
                if r.get('value') == 2:
                    ampsam[rev.get(r['entrezGeneId'])].add(r['sampleId'])
        if not quiet:
            print(f"    {lab}: {len(ampsam)} amplified genes", file=sys.stderr)
        for g, S in ampsam.items():
            if not g or g not in lesions or lesions[g].get(lab, {}).get('amp', 0) < AMP_CALL:
                continue
            a = arm_of.get(g)
            best = (0.0, None)
            for h, T in ampsam.items():
                if h == g or not h or arm_of.get(h) != a or len(T) <= len(S):
                    continue
                ov = len(S & T) / len(S)
                if ov > best[0]:
                    best = (ov, h)
            out.setdefault(g, {})[lab] = {"amp": lesions[g][lab]['amp'], "arm": a,
                                          "overlap": round(best[0], 2), "partner": best[1],
                                          "n_amp": len(S)}
    return out


def _acc_node(acc):
    return {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers."
                     "reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": acc}}


def _rcsb(nodes, return_type):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
         "return_type": return_type,
         "request_options": {"paginate": {"start": 0, "rows": 3000}}}
    return get_json(RCSB + urllib.parse.quote(json.dumps(q)))


def collect_pdb(genes, acc_of, quiet=False):
    """Per gene: entries holding it as a real chain, and how many of those hold a drug-like ligand.

    Two searches, deliberately. The chain search is restricted to polymer entities of at least
    MIN_CHAIN residues; the ligand search is not restricted at all, and the two are INTERSECTED. The
    difference between them -- entries with a ligand but no real chain -- is kept as `pep`, because
    it is the interesting failure: those are structures of somebody else's protein holding a peptide
    of this one, and reporting the count in the row is what makes the correction visible."""
    out = {}
    for i, g in enumerate(genes, 1):
        a = acc_of.get(g)
        if not a:
            out[g] = {"acc": None, "pdb": None, "lig": None, "pep": None}
            continue
        ch = _rcsb([_acc_node(a),
                    {"type": "terminal", "service": "text", "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "greater_or_equal", "value": MIN_CHAIN}}], "polymer_entity")
        time.sleep(0.15)
        lg = _rcsb([_acc_node(a),
                    {"type": "terminal", "service": "text", "parameters": {
                        "attribute": "chem_comp.formula_weight",
                        "operator": "greater", "value": LIGAND_DA}}], "entry") if ch is not None \
            else None
        time.sleep(0.15)
        if ch is None or lg is None:
            out[g] = {"acc": a, "pdb": None, "lig": None, "pep": None}
        else:
            chains = {x['identifier'].split('_')[0] for x in ch.get('result_set', [])}
            ligs = {x['identifier'] for x in lg.get('result_set', [])}
            out[g] = {"acc": a, "pdb": len(chains), "lig": len(chains & ligs),
                      "pep": len(ligs - chains)}
        if not quiet and i % 40 == 0:
            print(f"    {i}/{len(genes)} genes", file=sys.stderr)
    return out


# ---------------------------------------------------------------- the questions ---------------
def median_mutation(lesions):
    """One median per cohort, over the genes in this graph -- the baseline the outlier rule uses."""
    out = {}
    for _, lab in STUDIES:
        vals = [r[lab]['mut'] for r in lesions.values() if lab in r]
        out[lab] = statistics.median(vals) if vals else 0.0
    return out


def q1(g, d):
    """Lesion. Returns (verdict, evidence prose, structured evidence).

    The structured third value is what the scorer reads. An earlier version scored by matching
    strings in the prose and duly awarded the amplification bonus to the fallback line
    "max amp 3%, mut 1%, del 0%", which put H2AX second in the ranking."""
    r = d['lesions'].get(g, {})
    amp = max((v['amp'] for v in r.values()), default=0.0)
    dl = max((v['del'] for v in r.values()), default=0.0)
    mut, mutcoh = 0.0, None
    for coh, v in r.items():
        bar = max(MUT_FLOOR, MUT_MULT * d['medmut'].get(coh, 0.0))
        if v['mut'] >= bar and v['mut'] > mut:
            mut, mutcoh = v['mut'], coh
    mutmax = max((v['mut'] for v in r.values()), default=0.0)
    node = d['node'].get(g, {})
    fus, fn = node.get('fus') or '', node.get('fn') or 0

    bits, block = [], None
    if amp >= AMP_CALL:
        a = d['arm'].get(g, {})
        coh = max(a, key=lambda c: a[c]['amp'], default=None) if a else None
        if coh and a[coh]['overlap'] >= BLOCK_OVERLAP:
            block = (a[coh]['partner'], int(a[coh]['overlap'] * 100), a[coh]['arm'])
            bits.append(f"amp {amp}% <i>in a {block[2]} block with {block[0]} "
                        f"({block[1]}% overlap)</i>")
        else:
            bits.append(f"amp {amp}%")
    if mut:
        bits.append(f"mut {mut}% ({mutcoh}, outlier)")
    if dl >= DEL_CALL:
        bits.append(f"del {dl}%")
    if fus:
        side = {'5p': '5&prime;', '3p': '3&prime;', 'both': '5&prime;+3&prime;'}.get(fus, fus)
        bits.append(f"fusion {side}, {fn} partner" + ("s" if fn != 1 else ""))

    solid = [b for b in bits if 'block with' not in b]
    if solid:
        verdict = 'lesion'
    elif bits:
        verdict = 'block lesion'        # the segment is amplified; the gene is not resolved
    else:
        verdict = 'weak' if (amp >= WEAK or mutmax >= WEAK or dl >= WEAK) else 'none found'
    detail = ", ".join(bits) if bits else f"max amp {amp}%, mut {mutmax}%, del {dl}%"
    return verdict, detail, {'amp': amp, 'block': bool(block), 'mut': mut, 'del': dl, 'fn': fn}


def q2(g):
    """Bypass. There is no database for this, and inventing 415 answers is the failure mode the
    whole page is arranged against."""
    if g in HAND_Q2:
        return 'assessed', HAND_Q2[g]
    return 'not assessed', 'no bypass evidence derivable from the staged data'


def q3(g, d):
    """Normal-tissue need, proxied by pan-essentiality. A proxy, not a measurement: cell lines have
    no microenvironment, so for CTLA4, CD274, TIGIT, CD40 and lineage markers a null means nothing."""
    r = d['dep'].get(g)
    if not r:
        return 'not screened', 'gene absent from the DepMap matrix'
    c, f = r['class'], float(r['frac_dep']) * 100
    if c == 'common' or g in d['essential']:
        return 'yes', (f"pan-essential ({f:.1f}% of lines"
                       + ("; on DepMap's common-essential list)" if g in d['essential'] else ")"))
    if c == 'selective':
        return 'not established', f"selective ({f:.1f}%) - normal-tissue need not testable here"
    if c == 'rare':
        return 'not established', f"rare ({f:.1f}%)"
    return 'no dependency', f"no dependent line ({f:.1f}%)"


def dock(g, d):
    """Is there something to dock into? Entries holding the gene as a real chain, and how many of
    those hold a ligand over LIGAND_DA.

    A proxy in both directions: the ligand filter admits detergents, lipids and nucleotides
    (CTNNB1's are ADP, Mg and Na from co-crystallised partners, not binders in its groove), and an
    apo-only gene may never have been screened rather than lacking a pocket."""
    r = d['pdb'].get(g) or {}
    n, l, pep = r.get('pdb'), r.get('lig'), r.get('pep') or 0
    tail = (f"; {pep} further ligand entr{'y' if pep == 1 else 'ies'} hold only a peptide of it"
            if pep else "")
    if n is None:
        return 'unknown', 'no UniProt accession mapped'
    if n == 0:
        return 'model only', (f"no entry contains it as a chain{tail}; "
                              f"AlphaFold model AF-{r.get('acc')}-F1")
    if (l or 0) >= 10:
        return 'ligand-bound', f"{n} PDB entries, {l} with a ligand >{LIGAND_DA} Da{tail}"
    if (l or 0) >= 1:
        return 'few ligands', f"{n} PDB entries, only {l} with a ligand >{LIGAND_DA} Da{tail}"
    return 'apo only', f"{n} PDB entries, none with a drug-like ligand{tail}"


# ---------------------------------------------------------------- the score -------------------
def score(r, dk=DK_ACT):
    """(strongest lesion + publications + class + docking) x dependency.

    Lesion components are NOT summed -- several kinds of lesion do not make a gene several times
    more selectable, and summing hid a 34-partner fusion gene and a one-record curiosity behind the
    same word."""
    dep, f = r['dep'], r['depfrac']
    if dep == 'common':
        return None, 'pan-essential: no therapeutic window'
    if dep == 'selective' and 10 <= f <= 60:
        m, mw = 1.0, f'dependency {f:.0f}% (selective, in band)'
    elif dep == 'selective':
        m, mw = 0.8, f'dependency {f:.0f}% (selective, edge of band)'
    elif dep == 'rare':
        m, mw = 0.45, f'dependency {f:.0f}% (rare)'
    elif dep == '':
        m, mw = 0.35, 'not screened in DepMap'
    else:
        m, mw = 0.10, f'no dependent line ({f:.1f}%)'

    e, comp = r['ev'], []
    if e['mut']:
        comp.append((f"mutation outlier {e['mut']}%", 2.0))
    if e['del'] >= DEL_CALL:
        comp.append((f"deep deletion {e['del']}%", 2.0))
    if e['amp'] >= AMP_CALL:
        comp.append((f"amplification {e['amp']}%" + (' (arm-level block)' if e['block'] else ''),
                     1.0 if e['block'] else 2.0))
    fn = e['fn']
    if fn >= 5:
        comp.append((f'fusion, {fn} partners', 2.0))
    elif fn >= 2:
        comp.append((f'fusion, {fn} partners', 1.25))
    elif fn == 1:
        # CTNNB1::PLAG1 is one record in a salivary adenoma and had been scoring exactly as much
        # as EWSR1's 34 partners.
        comp.append(('fusion, a single partner', 0.5))
    lab, L = max(comp, key=lambda c: c[1]) if comp else ('no lesion found', 0.0)

    pb = r['pubs']                                  # one paper is not corroboration
    p_ = 1.5 if pb >= 5 else (1.0 if pb >= 3 else (0.6 if pb == 2 else 0.0))
    c = 0.5 if r['acat'] == 'driver' else (0.25 if r['acat'] == 'noa' else 0.0)
    k = dk.get(r['vd'], 0.0)
    base = max(0.0, L + p_ + c + k)   # DK_NEG can push a well-worked gene negative, and a negative
                                      # score would invert the multiplier, so the floor is zero
    parts = [f"{lab} +{L:g}", f"{pb} publication{'s' if pb != 1 else ''} +{p_:g}"]
    if c:
        parts.append(f"{r['acat']} +{c:g}")
    parts.append(f"{r['vd']} {k:+g}")               # signed: the neglect table goes negative
    return round(base * m, 2), f"({'; '.join(parts)}) x {m:g} for {mw}"


def lof_shaped(r, d):
    """Loss-of-function shape: a mutation outlier or a deep deletion, no amplification, no fusion.
    Tumour suppressors look like this, and for them the therapy is RESTORATION, so a high rank means
    the wrong thing. A caution, not a verdict: PIK3CA is mutated without amplification and is an
    oncogene."""
    les = d['lesions'].get(r['g'], {})
    amp = max((v['amp'] for v in les.values()), default=0.0)
    dl = max((v['del'] for v in les.values()), default=0.0)
    return (('outlier' in (r['d1'] or '')) or dl >= WEAK) and amp < AMP_CALL and not r['fus']


def build_rows(d):
    rows = []
    for g in sorted(d['node']):
        n = d['node'][g]
        v1, d1, ev = q1(g, d)
        v2, d2 = q2(g)
        v3, d3 = q3(g, d)
        vd, dd = dock(g, d)
        colour = 'green' if (n['target'] and n['tcat'] == 'green') else \
                 ('pink' if n['target'] else 'blue')
        rows.append(dict(g=g, colour=colour, acat=n.get('acat') or '', fus=n.get('fus') or '',
                         fn=n.get('fn') or 0, ev=ev, sents=n['sent95'],
                         pubs=len(d['pubs'].get(g, ())),
                         dep=(d['dep'].get(g, {}) or {}).get('class', ''),
                         depfrac=float((d['dep'].get(g, {}) or {}).get('frac_dep', 0)) * 100,
                         v1=v1, d1=d1, v2=v2, d2=d2, v3=v3, d3=d3, vd=vd, dd=dd))
    for r in rows:
        if r['colour'] != 'blue':
            r['rank'] = r['nrank'] = r['score'] = r['nscore'] = None
            r['why'] = r['nwhy'] = 'has a drug: not an undrugged candidate'
            r['lof'] = False
        else:
            r['score'], r['why'] = score(r, DK_ACT)
            r['nscore'], r['nwhy'] = score(r, DK_NEG)
            r['rank'] = r['nrank'] = None
            r['lof'] = lof_shaped(r, d)
    for key, rk in (('score', 'rank'), ('nscore', 'nrank')):
        ordered = sorted([r for r in rows if r['colour'] == 'blue' and r[key] is not None],
                         key=lambda r: (-r[key], -r['pubs'], -r['depfrac'], r['g']))
        for i, r in enumerate(ordered, 1):
            r[rk] = i
    return rows


# ---------------------------------------------------------------- render ----------------------
BADGE = {'ligand-bound': 'ok', 'few ligands': 'mid', 'apo only': 'mid', 'model only': 'no',
         'unknown': 'na',
         'lesion': 'ok', 'block lesion': 'mid', 'weak': 'mid', 'none found': 'no', 'yes': 'no',
         'not established': 'mid', 'no dependency': 'mid', 'not screened': 'na',
         'assessed': 'ok', 'not assessed': 'na'}

STYLE = """
 body{font:14px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:#1c2330;margin:0;background:#fff}
 .wrap{max-width:1180px;margin:0 auto;padding:24px 18px 60px}
 h1{font-size:22px;margin:0 0 4px} h2{font-size:17px;margin:26px 0 8px}
 .lede{color:#5b6677;margin:0 0 18px}
 .nrk{color:#7b8798;font-size:12px}
 .q{border:1px solid #cdd5e0;border-radius:8px;padding:12px 14px;margin:10px 0;background:#fbfcfe}
 .q h3{margin:0 0 6px;font-size:15px}
 .src{font-size:12px;color:#5b6677;border-left:3px solid #cdd5e0;padding-left:9px;margin-top:8px}
 .warn{border-left-color:#d59a2e;background:#fdf7e8;padding:8px 10px;border-radius:0 6px 6px 0}
 table{border-collapse:collapse;width:100%;margin-top:10px}
 th,td{border-bottom:1px solid #e7ebf1;padding:5px 7px;vertical-align:top;text-align:left}
 th{position:sticky;top:0;background:#f4f7fa;font-size:12px;cursor:pointer;white-space:nowrap}
 td.g{font-weight:600;white-space:nowrap} td.num{text-align:right;white-space:nowrap}
 .sub{font-size:11px;color:#5b6677}
 .b{display:inline-block;font-size:11px;padding:1px 6px;border-radius:10px;border:1px solid}
 .b.ok{background:#e8f4ec;border-color:#1b7837;color:#145a28}
 .b.mid{background:#fdf3e3;border-color:#d59a2e;color:#7a5a10}
 .b.no{background:#fdeaee;border-color:#b3243b;color:#8a1b2e}
 .b.na{background:#eef2f7;border-color:#cdd5e0;color:#5b6677}
 .dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:5px;vertical-align:-1px;border:1px solid #999}
 .dot.green{background:#1b7837;border-color:#145a28} .dot.pink{background:#c2185b;border-color:#7a0f3a}
 .dot.blue{background:#cfe3ff;border-color:#2b6cb0}
 #f{padding:6px 9px;border:1px solid #cdd5e0;border-radius:6px;width:230px;font-size:13px}
 .counts{font-size:12px;color:#5b6677;margin-left:10px}
 td.rk{font-weight:600;color:#1c2330;min-width:118px}
 td.rk .sub{font-weight:400}
 .lof{font-size:11px;color:#8a1b2e;margin-top:2px}
"""

SCRIPT = """
const tb=document.querySelector('#t tbody'), f=document.getElementById('f'), cnt=document.getElementById('cnt');
function upd(){const q=f.value.trim().toLowerCase();let n=0;
 [].forEach.call(tb.rows,function(r){const s=!q||r.dataset.g.indexOf(q)>=0;r.style.display=s?'':'none';if(s)n++;});
 cnt.textContent=n+' of '+tb.rows.length+' genes';}
f.addEventListener('input',upd);upd();
[].forEach.call(document.querySelectorAll('#t th'),function(th,k){
 th.addEventListener('click',function(){
  const rows=[].slice.call(tb.rows), asc=th.dataset.asc!=='1';
  rows.sort(function(a,b){
    const ca=a.cells[k], cb=b.cells[k];
    // an explicit data-sort wins: the rank cell reads "#12", which parseFloat cannot order
    if(ca.dataset.sort!==undefined&&cb.dataset.sort!==undefined){
      const na=+ca.dataset.sort, nb=+cb.dataset.sort; return asc?na-nb:nb-na;}
    const x=ca.innerText.trim(), y=cb.innerText.trim();
    const nx=parseFloat(x), ny=parseFloat(y);
    if(!isNaN(nx)&&!isNaN(ny)) return asc?nx-ny:ny-nx;
    return asc?x.localeCompare(y):y.localeCompare(x);});
  th.dataset.asc=asc?'1':'0';
  rows.forEach(function(r){tb.appendChild(r);});});});
"""

DOC = """<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Q1-Q3 for every gene in the graph</title>
<style>%(style)s</style></head><body><div class="wrap">
<h1>Q1&ndash;Q3 for every gene in the graph</h1>
<p class="lede">All %(n_genes)d gene nodes of <code>%(graph)s</code> put through
the three-question test. <b>The three questions are not equally answerable, and this page does not
pretend otherwise</b> &mdash; each column states what its answer was derived from, and Q2 is
<i>unanswered for %(q2_no)d of the %(n_genes)d genes</i>, because nothing in the staged data can
answer it.</p>

<h2>The three questions</h2>

<div class="q"><h3>Q1 &middot; Is there a lesion that predicts the dependency?</h3>
A companion-diagnostic question, not a biological one: <i>if the dependency is real, what do you
assay in a biopsy to find the patients it applies to?</i> A <b>lesion</b> is a discrete, testable
alteration <b>present in the tumour and absent from normal tissue</b> &mdash; a point mutation, an
amplification, a deep deletion, a fusion. That asymmetry is what makes it a selection tool rather
than a description. A dependency without a lesion is not false, it is <b>unprospectable</b>: the
trial cannot be designed. EGFR and HSF1 are each required by about 21%% of cell lines, but a
mutation names which 21%% for EGFR and nothing names it for HSF1.
<div class="src"><b>Answered here from:</b> cBioPortal, TCGA PanCancer Atlas, three cohorts &mdash;
%(cohorts)s &mdash; taking each gene's <i>best</i> cohort, plus ChimerDB fusion status carried on
the graph node.
<b>lesion</b> = amplification &ge;%(amp_call)g%%, deep deletion &ge;%(del_call)g%%, a curated
fusion, or a mutation rate that is a <i>cohort outlier</i> &mdash; at least %(mut_floor)g%% and at
least %(mut_mult)g&times; that cohort's median across these %(n_genes)d genes. The outlier rule
matters: a flat cut calls every <i>large</i> gene mutated, because passenger burden scales with
coding length. FASN reads 6.1%% in colorectal and NBEA 10.6%% against a cohort median of 1.3%%;
TP53 (53%%), APC (65%%) and KRAS (37%%) clear the bar by an order of magnitude. Under a flat 3%% cut
FASN was called a lesion, which contradicts the hand analysis in
<code>candidate_targets.md</code> &mdash; the rule was changed, not the conclusion.
<b>weak</b> = above roughly half those thresholds. <b>none found</b> = below them.
<b>Arm-level amplification is detected and marked.</b> For every gene called amplified, the
sample-level calls are compared with those of larger-amplified genes on the same chromosome arm. If
&ge;%(block_pct)d%% of a gene's amplified samples also amplify that neighbour, the <i>segment</i> is
what is being selected and the verdict becomes <b>block lesion</b>: a real lesion in the region, but
which gene inside it the selection is for cannot be resolved from copy number. <b>%(n_block)d of
%(n_amp)d amplification calls are block calls</b> &mdash; HSF1 (8q24.3) shares 95%% of its amplified
breast samples with PTP4A3, PVT1 shares 98%% with MYC, CDK12 98%% with ERBB2. Note what the test
does <i>not</i> claim: MYC and ERBB2 are the largest amplifications on their arms so they keep a
plain <b>lesion</b> call, while FGFR1 reads <b>block lesion</b> even though its 8p11 amplicon is a
recognised breast-cancer event &mdash; the method finds blocks, not drivers. And a curated fusion
with <i>one</i> partner is much weaker evidence than EWSR1's 34: TP53 has two rare ChimerKB records
and is not a fusion-driven cancer, so the partner count is printed for every gene rather than hidden
behind the word "fusion".
<b>Limits:</b> three cohorts, not pan-cancer, so a gene lesioned only in leukaemia or sarcoma reads
"none found" here &mdash; NPM1, RUNX1 and CRLF2 are exactly that case. And a recurrent lesion is
<i>necessary but not sufficient</i>: it shows a marker exists, not that the marker predicts the
dependency. That last step needs a lesion-versus-dependency cross-tabulation, which is not run
here.</div></div>

<div class="q"><h3>Q2 &middot; Can the cell obtain the product another way?</h3>
Bypass. Blocking synthesis achieves nothing if the cell can scavenge the product, use a paralogue,
or switch to a parallel branch. This is what kills most metabolic targets: FASN dependence is
rescued by lipid uptake, glutaminase by scavenged glutamine, one UPR arm by the other two.
<div class="src warn"><b>Not derivable from the staged data.</b> Answering it needs isotope tracing,
rescue experiments or media-composition screens &mdash; none of them a database lookup.
<b>%(q2_yes)d of %(n_genes)d genes</b> carry a verdict, the ones examined by hand in sections
7&ndash;9 of <code>candidate_targets.md</code>; the other %(q2_no)d read <i>not assessed</i>.
Filling them from memory would have produced %(q2_no)d fluent, unverifiable answers, which is the
exact failure this project logged nine times.</div></div>

<div class="q"><h3>Q3 &middot; Does normal tissue need the same enzyme?</h3>
The therapeutic window. If healthy cells need the gene as much as the tumour does, inhibiting it
buys toxicity rather than selectivity.
<div class="src"><b>Answered here from:</b> DepMap CRISPR dependency across %(n_lines)s cell lines,
plus DepMap's own inferred common-essential list. <b>yes</b> = pan-essential, i.e. nearly every line
dies without it, the strongest available proxy for "normal cells need it too".
<b>not established</b> = a selective or rare dependency; cell lines cannot show normal-tissue need
either way. <b>no dependency</b> = no line depends on it.
<b>Limits:</b> a proxy, not a measurement of normal tissue. Cell lines have no microenvironment, no
immune system and no differentiation state, so for immune and stromal genes (CTLA4, CD274, TIGIT,
CD40) and lineage markers (AFP, UPK2, RAG2) a null here means nothing at all. The direct sources
would be gnomAD constraint, GTEx expression breadth and germline-disease databases; none is
staged.</div></div>

<div class="q"><h3>Docking readiness &mdash; is there something to dock into?</h3>
Not one of the three questions, but the one that decides whether any of them can be acted on with a
small molecule. A dependency with a lesion and a therapeutic window is still nothing to a chemist
without a pocket.
<div class="src"><b>Answered here from:</b> the RCSB PDB, searched by each gene's UniProt accession
&mdash; how many entries exist, and how many hold a <b>ligand over %(ligand_da)d Da</b>, a proxy for
a pocket somebody has already occupied with a drug-like molecule. <b>ligand-bound</b> = 10 or more
such entries, <b>few ligands</b> = 1&ndash;9, <b>apo only</b> = structures but no drug-like ligand,
<b>model only</b> = no experimental structure at all.
<b>The entry must hold the protein as a real chain</b> (&ge;%(min_chain)d residues), not as a docked
peptide. Without that condition the count answers a different question &mdash; <i>does this
accession appear near a ligand</i> &mdash; and RELA scored 22 ligand-bound entries of which every
one was a <b>14-3-3 sigma</b> structure holding a RelA phosphopeptide, with the compounds covalently
glued to 14-3-3. The pocket was someone else's. Requiring a real chain takes RELA to <b>0</b> and
leaves MDM2 (80), MCL1 (95) and BCL6 (139) untouched. Where a gene still has entries that hold only
a peptide of it, the count is reported in its row.
<b>Limits, in both directions.</b> The %(ligand_da)d Da filter admits detergents, lipids and
nucleotides: CTNNB1's two ligand entries are ADP and Mg from co-crystallised partners, not
binders in its groove, so a low count deserves a look at what the ligands actually are. And an
<i>apo only</i> gene may simply never have been screened rather than lacking a pocket.
<b>model only</b> means no experimental coordinates &mdash; an AlphaFold model exists for
essentially every human protein, so it is a statement about evidence, not about whether the protein
folds.</div></div>

<h2>The ranking</h2>
<div class="q"><h3>Candidate rank &mdash; blue nodes only</h3>
The <b>rank</b> column orders the %(n_blue)d undrugged (blue) genes as inhibition candidates. Green
and pink genes already have a drug and read &ndash;.
<div class="src"><b>Score = (Q1 + evidence + class + docking) &times; dependency.</b> Dependency
<i>multiplies</i> rather than adds, because without something to inhibit a lesion and a literature
are worth nothing: TP53 has 34 publications and the clearest lesion in cancer, and 0.4%% of lines
depend on it. Multipliers: selective and in a 10&ndash;60%% band &times;1.0, selective outside it
&times;0.8, rare &times;0.45, not screened &times;0.35, no dependent line &times;0.10.
Added terms: the <b>strongest single piece of lesion evidence</b>, not the label &mdash; mutation
outlier, deep deletion or a plain amplification +2, an arm-level block amplification +1, and a
fusion weighted by how often it has been seen: <b>&ge;5 partners +2</b>, 2&ndash;4 +1.25, <b>a
single partner +0.5</b>. That last weight matters: CTNNB1::PLAG1 is one record in a salivary
adenoma and had been scoring exactly as much as EWSR1's 34 partners. Components are not summed
&mdash; several kinds of lesion do not make a gene several times more selectable. Corpus support is
counted in <i>publications</i>, steeply: &ge;5 +1.5, 3&ndash;4 +1.0, two +0.6, and <b>a single paper
earns nothing</b>, because one paper is not corroboration &mdash; the PRDM14 lesson;
curated class <b>driver</b> +0.5, <b>noa</b> +0.25; and <b>docking readiness</b> &mdash;
<b>ligand-bound +1.5</b>, few ligands +0.75, apo only +0.25, model only 0. Docking <i>adds</i>
rather than multiplying, and the asymmetry is deliberate: dependency is a gate, because without
something to inhibit there is nothing to do, but a pocket is not &mdash; an apo structure is a
perfectly good docking input and a protein with no structure at all can still be drugged once
somebody tries. Evidence of an occupied pocket is a bonus; its absence is not a veto.
<b>%(n_disq)d blue genes are disqualified outright</b> as pan-essential: no therapeutic window at
any rank, leaving %(n_ranked)d ranked.
<b>The neglect column is the same score with the docking term reversed and centred on zero</b>
&mdash; model only <b>+1.5</b>, apo only +0.5, few ligands &minus;0.5, <b>ligand-bound
&minus;1.5</b> (floored at zero) &mdash; and it answers a different question: not <i>what is most
actionable</i> but <i>where has the chemistry never been tried</i>.
Sort by it to find genes carrying a lesion and a selective dependency that no one has yet put a
ligand into. Everything else in the formula is unchanged, and the publication term deliberately is
not inverted: publications measure attention too, but flipping them would rank the single-paper
curiosities first, which is the PRDM14 error this project already made once. Corroboration is
required either way; only the <i>chemistry</i> attention is reversed. A gene high on both columns
is well supported and under-served; high on neglect but low on rank usually means the biology is
thin as well as the chemistry. The penalty rather than a mere absence of bonus is what makes the
column work: a protein whose pocket has already been filled is not neglected, and without the
minus sign <b>MDM2 topped both columns</b> &mdash; it has nutlins in trials and is blue here only
because this corpus carries no chemical edge for it.
<b>What it cannot do.</b> Q2 is absent from the formula &mdash; bypass is unassessed for %(q2_no)d
genes &mdash; so a well-ranked gene can still be rescued by a route nobody has tested; that is how
FASN would have scored well before its uptake bypass was known. And <b>the score cannot tell an
oncogene from a tumour suppressor</b>: APC and STK11 rank in the twenties on real lesions and real
dependencies, but their alterations are loss-of-function and the therapy would be restoration, not
inhibition. Genes whose lesion has that shape &mdash; mutation outlier or deep deletion, no
amplification, no fusion &mdash; carry a red caution under the rank. The caution has false
positives (PIK3CA is mutated without amplification and is an oncogene) and false negatives.</div>
</div>

<h2>How the answers came out</h2>
<p class="lede">Q1 &mdash; %(c1)s. &nbsp;&nbsp; Q3 &mdash; %(c3)s. &nbsp;&nbsp; Q2 &mdash;
%(q2_yes)d assessed, %(q2_no)d not.
<br>A gene answering <b>lesion / no / no</b> is what a target looks like. Nothing in this table
answers all three, because Q2 is unanswered almost everywhere.</p>

<p><input id="f" placeholder="filter by gene symbol&hellip;" autocomplete="off">
<span class="counts" id="cnt"></span></p>
<table id="t"><thead><tr>
<th>rank</th><th title="same score with the docking term reversed">neglect</th><th>gene</th>
<th>drug status</th><th>addiction class</th><th>sents/pubs</th>
<th>Q1 &middot; lesion</th><th>Q2 &middot; bypass</th><th>Q3 &middot; normal need</th>
<th>docking readiness</th>
</tr></thead><tbody>
%(body)s
</tbody></table>
<script>%(script)s</script>
</div></body></html>"""


def badge(v):
    return '<span class="b %s">%s</span>' % (BADGE.get(v, 'na'), html.escape(v))


def render(rows, d, out):
    body = []
    for r in rows:
        body.append(
            '<tr data-g="%s"><td class="num rk" data-sort="%d">%s<div class="sub">%s</div></td>'
            '<td class="num nrk" data-sort="%d" title="%s">%s</td>'
            '<td class="g">%s</td>'
            '<td><span class="dot %s"></span>%s</td><td>%s</td>'
            '<td class="num">%d<span class="sub">/%dp</span></td>'
            '<td>%s<div class="sub">%s</div></td>'
            '<td>%s<div class="sub">%s</div></td>'
            '<td>%s<div class="sub">%s</div></td>'
            '<td>%s<div class="sub">%s</div></td></tr>' % (
                html.escape(r['g'].lower()),
                r['rank'] or 99999,
                ('#%d' % r['rank']) if r['rank'] else '&ndash;',
                html.escape((('score %.2f - ' % r['score']) if r['score'] is not None else '')
                            + (r['why'] or ''))
                + ('<div class="lof">loss-of-function shape &mdash; may be a suppressor</div>'
                   if r['lof'] else ''),
                r['nrank'] or 99999,
                html.escape((('neglect score %.2f - ' % r['nscore'])
                             if r['nscore'] is not None else '') + (r['nwhy'] or '')),
                ('#%d' % r['nrank']) if r['nrank'] else '&ndash;',
                html.escape(r['g']),
                r['colour'], r['colour'], html.escape(r['acat'] or '-'),
                r['sents'], r['pubs'],
                badge(r['v1']), r['d1'],
                badge(r['v2']), html.escape(r['d2']),
                badge(r['v3']), html.escape(r['d3']),
                badge(r['vd']), html.escape(r['dd'])))

    c1 = collections.Counter(r['v1'] for r in rows)
    c3 = collections.Counter(r['v3'] for r in rows)
    n = {lab: v[lab]['n'] for v in d['lesions'].values() for _, lab in STUDIES if lab in v}
    blue = [r for r in rows if r['colour'] == 'blue']
    # Counts the prose used to hard-code. Every one of them has already been wrong once: the page
    # claimed 352 blue genes when there were 351, and 49 of 60 block calls when the arm test had
    # moved on to 35.
    fields = dict(
        style=STYLE, script=SCRIPT, body="\n".join(body),
        graph=Path(d['graph']).name,
        n_genes=len(rows), n_blue=len(blue),
        n_ranked=sum(1 for r in blue if r['rank']),
        n_disq=sum(1 for r in blue if r['score'] is None),
        q2_yes=sum(1 for r in rows if r['v2'] == 'assessed'),
        q2_no=sum(1 for r in rows if r['v2'] != 'assessed'),
        n_amp=sum(1 for r in rows if r['ev']['amp'] >= AMP_CALL),
        n_block=sum(1 for r in rows if r['ev']['block']),
        n_lines=f"{d['n_lines']:,}",
        cohorts=", ".join(f"{COHORT_LABEL[lab]} (n={n[lab]})" for _, lab in STUDIES if lab in n),
        amp_call=AMP_CALL, del_call=DEL_CALL, mut_floor=MUT_FLOOR, mut_mult=MUT_MULT,
        block_pct=int(BLOCK_OVERLAP * 100), min_chain=MIN_CHAIN, ligand_da=LIGAND_DA,
        c1=", ".join("%s %d" % kv for kv in c1.most_common()),
        c3=", ".join("%s %d" % kv for kv in c3.most_common()))
    doc = DOC % fields
    Path(out).write_text(doc, encoding="utf-8")
    return doc, fields


# ---------------------------------------------------------------- main ------------------------
def load_data(graph=None, data_root=DEFAULT_ROOT, cache=DEFAULT_CACHE, refresh=(),
              offline=False, quiet=False):
    """Everything the questions are answered from, in one dict. Split out of main() so that
    audit_rankings.py re-runs the real loader instead of a copy of it that can drift."""
    graph = graph or newest_graph(ROOT)
    if not graph:
        raise SystemExit(f"no *_M.html found in {ROOT}; pass --graph")
    db = Path(data_root) / "databases"
    payload = read_payload(graph)
    node = {n['id']: n for n in payload['nodes'] if n.get('kind') == 'gene'}
    pubs = collections.defaultdict(set)
    for e in payload['edges']:
        for x in e['sents']:
            pubs[e['from']].add(x['pmid'])
            pubs[e['to']].add(x['pmid'])
    print(f"{Path(graph).name}: {len(node)} gene nodes")

    docs = hgnc_docs(db)
    entrez_all = {d['symbol']: int(d['entrez_id']) for d in docs
                  if d.get('symbol') and d.get('entrez_id')}
    acc_all = {d['symbol']: (d.get('uniprot_ids') or [None])[0] for d in docs if d.get('symbol')}
    arm_all = {}
    for d0 in docs:
        s, loc = d0.get('symbol'), d0.get('location') or ''
        mt = re.match(r'^(\d+|X|Y)([pq])', loc)
        if s and mt:
            arm_all[s] = mt.group(1) + mt.group(2)

    dep = {r['gene']: r for r in csv.DictReader(
        open(db / "depmap_dependency.tsv", encoding="utf-8"), delimiter="\t")}
    with open(db / "CRISPRInferredCommonEssentials.csv", encoding="utf-8") as fh:
        essential = {r[0].split(" (")[0] for r in csv.reader(fh) if r and r[0] != "Essentials"}
    n_lines = max((int(r['n_lines']) for r in dep.values()), default=0)

    ref = set(refresh or ())
    rf = lambda k: ("all" in ref) or (k in ref)                              # noqa: E731
    cache = Path(cache)
    entrez = cached(cache, "gene_entrez", rf("entrez"), False,
                    lambda: {g: entrez_all[g] for g in node if g in entrez_all})
    lesions = cached(cache, "lesions", rf("lesions"), offline,
                     lambda: collect_lesions(entrez, quiet),
                     f" -- cBioPortal, {len(STUDIES)} cohorts x {len(entrez)} genes")
    arm = cached(cache, "armlevel", rf("armlevel"), offline,
                 lambda: collect_armlevel(entrez, lesions, arm_all, quiet),
                 " -- sample-level co-amplification")
    pdb = cached(cache, "pdb", rf("pdb"), offline,
                 lambda: collect_pdb(sorted(node), acc_all, quiet),
                 f" -- RCSB, 2 searches x {len(node)} genes")

    return dict(node=node, pubs=pubs, dep=dep, essential=essential, lesions=lesions, arm=arm,
                pdb=pdb, medmut=median_mutation(lesions), graph=str(graph), n_lines=n_lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph", help="graph HTML to read (default: the newest *_M.html beside this "
                                    "script)")
    ap.add_argument("--data-root", default=str(DEFAULT_ROOT),
                    help="pipeline output tree holding databases/ (default: kaggle_working)")
    ap.add_argument("--cache", default=str(DEFAULT_CACHE),
                    help="where collected API results are kept (default: rankings_cache)")
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="HTML to write (default: ../Q1_3.html)")
    ap.add_argument("--refresh", action="append", default=[],
                    choices=["entrez", "lesions", "armlevel", "pdb", "all"],
                    help="re-collect a stage even if it is cached; repeatable")
    ap.add_argument("--offline", action="store_true",
                    help="never touch the network; fail if a stage is not cached")
    ap.add_argument("--quiet", action="store_true", help="no per-batch progress on stderr")
    args = ap.parse_args()

    d = load_data(args.graph, args.data_root, args.cache, args.refresh, args.offline, args.quiet)
    rows = build_rows(d)
    doc, f = render(rows, d, args.out)

    ranked = sorted([r for r in rows if r['rank']], key=lambda r: r['rank'])
    negl = sorted([r for r in rows if r['nrank']], key=lambda r: r['nrank'])
    print(f"{f['n_genes']} rows | {f['n_blue']} blue, {f['n_ranked']} ranked, "
          f"{f['n_disq']} disqualified as pan-essential")
    print(f"  Q1 {f['c1']}")
    print(f"  Q3 {f['c3']}")
    print(f"  {f['n_block']} of {f['n_amp']} amplification calls are arm-level blocks")
    print("  best candidates: " + ", ".join(f"{r['g']} ({r['score']:g})" for r in ranked[:10]))
    print("  most neglected:  " + ", ".join(f"{r['g']} ({r['nscore']:g})" for r in negl[:10]))
    print(f"wrote {args.out} ({len(doc):,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
