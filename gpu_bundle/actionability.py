#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""actionability.py -- AMP/ASCO/CAP actionability of the corpus variants, from CIViC.

Reads the gene -> variants table of step 4b (VARIANT/variants.json: tmVar 3 variants from NCBI
PubTator3) and annotates every variant CIViC knows with its clinical evidence and an AMP/ASCO/CAP
tier (Li et al., J Mol Diagn 2017;19:4-23). CIViC (https://civicdb.org) is CC0: its nightly TSV
exports are downloaded into databases/civic/ and reused for CACHE_DAYS days (--refresh refetches).

HOW A VARIANT IS MATCHED. By NCBI Gene id (CIViC's entrez_id, mapped to the HGNC symbol the graph
uses) and HGVS: CIViC's variant name read as a protein change ("T790M" -> p.T790M), its aliases,
and the p. / c. parts of its hgvs_descriptions with the transcript dropped ("NP_005219.2:
p.Thr790Met" -> p.T790M, "NM_005228.4:c.2369C>T" -> c.2369C>T). CIViC files its evidence under
molecular profiles; a profile of one variant applies to it directly, and a profile combining
several ("BCR::ABL1 Fusion AND ABL1 T315I", "EML4::ALK Fusion AND ALK G1202R") applies to each of its
variants IN THAT CONTEXT -- that is where most TKI-resistance evidence sits -- and every item and
tier keeps the profile's name so the context is never lost. Category variants ("Exon 19 Deletion",
"Amplification") match nothing.

HOW THE TIER IS SET, per variant, best first:
  1. a CIViC ASSERTION with an AMP category: that category, as curated ("Tier I - Level A");
  2. otherwise the best accepted EVIDENCE ITEM that supports a therapeutic (Predictive), diagnostic
     or prognostic significance, of somatic or unstated origin -- AMP tiers classify somatic
     variants by exactly these three kinds of significance. CIViC's evidence levels are read as
     AMP's: A (validated, e.g. FDA / guidelines) -> Tier I, Level A; B (clinical evidence) -> Tier
     I, Level B; C (case studies, other tumour types) -> Tier II, Level C; D (preclinical) and E
     (inferential) -> Tier II, Level D.
  3. otherwise POPULATION FREQUENCY decides between the two tiers AMP gives variants without
     actionable evidence: Tier IV (benign or likely benign) when the variant's own allele reaches
     AF_BENIGN (1%) in gnomAD, overall or in any population (AMP: "observed at significant allele
     frequency in the general or specific subpopulation databases"); Tier III (unknown clinical
     significance) otherwise. Frequencies come from gnomAD via MyVariant.info, looked up by the
     rsID tmVar gave the variant; an rsID can name several alleles, so the one whose protein or
     cDNA change is the variant's is used. Tier III records why: its gnomAD frequency, "not in
     gnomAD", or "no rsID" (frequency not checked -- about half the variants).
An AMP tier is a statement about a variant IN A TUMOUR TYPE: Level A/B evidence in another tumour
type is Level C. A tier here is therefore the best one CIViC supports in the disease of that
evidence, and every tier is listed with its disease and therapies. AMP defines evidence levels
(A-D) for Tiers I and II only; Tiers III and IV have none, so none is invented here.

Inputs : VARIANT/variants.json                 step 4b
         databases/hgnc_complete_set_*.json     NCBI Gene id -> HGNC symbol
         databases/civic/nightly-*.tsv           downloaded here (internet), CC0
         myvariant.info                          gnomAD allele frequencies by rsID (internet), cached
                                                 in VARIANT/gnomad_af.json
Outputs: VARIANT/actionability.json             gene -> hgvs -> {tier, level, source, items, ...}
                                                 for EVERY corpus variant

Run::  python actionability.py [--root DIR] [--refresh]
"""

import argparse
import collections
import csv
import datetime
import glob
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CIVIC_URL = "https://civicdb.org/downloads/nightly/nightly-{}.tsv"
CIVIC_FILES = ("VariantSummaries", "ClinicalEvidenceSummaries", "AssertionSummaries",
               "MolecularProfileSummaries")
CIVIC_DIR = "databases/civic"
CACHE_DAYS = 7
UA = {"User-Agent": "relationship-graphs/1.0"}
OUT = "VARIANT/actionability.json"
MAX_ITEMS = 25          # evidence items kept per variant (best level first); counts stay complete

csv.field_size_limit(10 ** 8)
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*"}
LEVEL_TIER = {"A": ("I", "A"), "B": ("I", "B"), "C": ("II", "C"), "D": ("II", "D"), "E": ("II", "D")}
TIER_TYPES = {"Predictive", "Diagnostic", "Prognostic"}
GERMLINE = {"Rare Germline", "Common Germline"}
AMP_CAT = re.compile(r"Tier\s+(I{1,3}|IV)\s*-\s*Level\s+([A-D])", re.I)
MYVARIANT = "https://myvariant.info/v1/query"
AF_FILE = "VARIANT/gnomad_af.json"
AF_BENIGN = 0.01        # Tier IV at or above this allele frequency (overall or any population)
AF_BATCH = 500          # rsIDs per MyVariant.info request (it takes up to 1000)


def _anns(hit):
    a = (hit.get("snpeff") or {}).get("ann") or []
    return a if isinstance(a, list) else [a]


# gnomAD groups left out of the population maximum, as gnomAD's own grpmax leaves them out: small
# or bottlenecked (Amish, Ashkenazi Jewish, Finnish, Middle Eastern) or mixed ("other"/"remaining"),
# where one founder variant can reach a frequency it has nowhere else
POP_SKIP = ("ami", "asj", "fin", "mid", "oth", "remaining")


def _afs(hit):
    """{source:key -> AF} of one allele: gnomAD exome and genome, overall and per population
    (the sex-split af_*_female / af_*_male fields are left out)."""
    out = {}
    for src in ("gnomad_exome", "gnomad_genome"):
        af = (hit.get(src) or {}).get("af") or {}
        for k, v in af.items():
            if k.startswith("af") and not k.endswith(("_female", "_male")) and isinstance(v, (int, float)):
                out[f"{src[7:]}:{k}"] = v
    return out


def popmax(afs):
    """The highest AF of an allele over overall and the grpmax populations; None if not in gnomAD."""
    vals = [v for k, v in afs.items() if k.split(":", 1)[1] == "af" or k.split(":", 1)[1][3:] not in POP_SKIP]
    return max(vals) if vals else None


def fetch_af(rsids, cache_path):
    """rsID -> [allele, ...] with each allele's protein / cDNA changes and gnomAD popmax; cached, so
    a re-run asks only for new rsIDs. An rsID MyVariant.info does not know is cached as []."""
    cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}
    todo = sorted(r for r in rsids if r not in cache)
    if todo:
        print(f"  gnomAD frequencies: {len(todo):,} rsID(s) to look up, {len(cache):,} cached", flush=True)
    for i in range(0, len(todo), AF_BATCH):
        batch = todo[i:i + AF_BATCH]
        data = urllib.parse.urlencode({
            "q": ",".join(batch), "scopes": "dbsnp.rsid",          # up to 10 alleles per rsID (default)
            "fields": "gnomad_exome.af,gnomad_genome.af,snpeff.ann.hgvs_p,snpeff.ann.hgvs_c,snpeff.ann.genename"}).encode()
        req = urllib.request.Request(MYVARIANT, data=data, headers=UA)
        for k in range(4):
            try:
                with urllib.request.urlopen(req, timeout=300) as r:
                    hits = json.load(r)
                break
            except (urllib.error.URLError, TimeoutError, ValueError):
                if k == 3:
                    raise
                time.sleep(2 ** (k + 1))
        got = collections.defaultdict(list)
        for h in hits:
            if h.get("notfound"):
                continue
            anns = _anns(h)
            got[h["query"]].append({
                "id": h.get("_id", ""), "afs": _afs(h),
                "p": sorted({p1(a["hgvs_p"]) for a in anns if a.get("hgvs_p")}),
                "c": sorted({a["hgvs_c"] for a in anns if a.get("hgvs_c")}),
                "genes": sorted({a["genename"] for a in anns if a.get("genename")})})
        for r in batch:
            cache[r] = got.get(r, [])
        cache_path.write_text(json.dumps(cache, ensure_ascii=False) + "\n", encoding="utf-8")
        time.sleep(0.5)
    return cache


def allele_for(hgvs, gene, alleles):
    """The allele of an rsID that is this variant: its protein or cDNA change matches. An rsID
    with a single allele is taken as it is; several alleles and no match -> None."""
    hit = [a for a in alleles if hgvs in a["p"] or hgvs in a["c"]]
    if hit:
        return max(hit, key=lambda a: (popmax(a["afs"]) is not None, popmax(a["afs"]) or 0))
    return alleles[0] if len(alleles) == 1 else None


def population_tier(x, gene, af_cache):
    """(tier, record) for a variant with no Tier I/II evidence: IV if common, else III."""
    rs = x.get("rsid") or ""
    if not rs:
        return "III", {"af": None, "af_basis": "no rsID: frequency not checked"}
    a = allele_for(x["hgvs"], gene, af_cache.get(rs, []))
    if a is None:
        why = "rsID unknown to MyVariant.info" if not af_cache.get(rs) else "no allele of the rsID matches"
        return "III", {"af": None, "af_basis": f"{rs}: {why}"}
    af = popmax(a["afs"])
    if af is None:
        return "III", {"af": None, "af_basis": f"{rs} ({a['id']}): not in gnomAD"}
    basis = f"gnomAD grpmax {af:.4g} ({rs}, {a['id']})"
    return ("IV" if af >= AF_BENIGN else "III"), {"af": af, "af_basis": basis}


def download(civic_dir, refresh):
    """The CIViC nightly TSVs, fetched unless a copy younger than CACHE_DAYS is on disk."""
    civic_dir.mkdir(parents=True, exist_ok=True)
    for name in CIVIC_FILES:
        f = civic_dir / f"nightly-{name}.tsv"
        fresh = f.exists() and f.stat().st_size > 0 and (time.time() - f.stat().st_mtime) < CACHE_DAYS * 86400
        if fresh and not refresh:
            continue
        req = urllib.request.Request(CIVIC_URL.format(name), headers=UA)
        with urllib.request.urlopen(req, timeout=300) as r:
            data = r.read()
        if b"\t" not in data[:2000]:
            raise SystemExit(f"CIViC {name}: unexpected answer (not a TSV)")
        tmp = f.with_suffix(".part")
        tmp.write_bytes(data)
        tmp.replace(f)
        print(f"  CIViC {name}: {len(data):,} bytes -> {f}")
    return {n: list(csv.DictReader(open(civic_dir / f"nightly-{n}.tsv", encoding="utf-8", newline=""),
                                   delimiter="\t")) for n in CIVIC_FILES}


def p1(three):
    """'p.Thr790Met' -> 'p.T790M' (simple substitutions; anything else unchanged)."""
    m = re.fullmatch(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=)\)?", three)
    if m and m.group(1) in AA3 and (m.group(3) == "=" or m.group(3) in AA3):
        alt = "=" if m.group(3) == "=" else AA3[m.group(3)]
        return f"p.{AA3[m.group(1)]}{m.group(2)}{alt}"
    return three


def civic_keys(v):
    """HGVS strings a CIViC variant answers to, in the forms tmVar writes."""
    keys = set()
    name = (v.get("variant") or "").strip()
    if re.fullmatch(r"[A-Z*]\d+[A-Z*]", name):            # T790M, R248*
        keys.add("p." + name)
    for a in (v.get("variant_aliases") or "").split(","):
        a = a.strip()
        if re.fullmatch(r"[A-Z]\d+[A-Z*]", a):
            keys.add("p." + a)
    for h in (v.get("hgvs_descriptions") or "").split(","):
        h = h.strip().split(":", 1)[-1]
        if h.startswith("p."):
            keys.add(h)
            keys.add(p1(h))
        elif h.startswith("c."):
            keys.add(h)
    return keys


def hgnc_by_entrez(root):
    hits = sorted(glob.glob(str(root / "databases" / "hgnc_complete_set_*.json")))
    if not hits:
        raise SystemExit(f"no databases/hgnc_complete_set_*.json under {root}")
    docs = json.load(open(hits[-1], encoding="utf-8"))["response"]["docs"]
    return {str(d["entrez_id"]): d["symbol"] for d in docs if d.get("entrez_id")}


def tier_rank(t, lv):
    """Sort key: Tier I-A best."""
    return ({"I": 0, "II": 1, "III": 2, "IV": 3}.get(t, 9), lv or "Z")


def annotate(tables, genes, symbol_of):
    """gene -> hgvs -> annotation, for every corpus variant CIViC knows."""
    V, E, A, M = (tables[n] for n in CIVIC_FILES)
    # (gene symbol, hgvs) -> the molecular profiles naming that variant, alone or with others
    var_by_id = {}
    for v in V:
        sym = symbol_of.get((v.get("entrez_id") or "").strip())
        if sym:
            var_by_id[v["variant_id"].strip()] = (sym, v)
    mp_name, index = {}, collections.defaultdict(set)
    first_var = {}                                   # profile -> its variant row (for the CIViC link)
    for m in M:
        mp = m["molecular_profile_id"].strip()
        mp_name[mp] = m["name"]
        for vid in (m.get("variant_ids") or "").split(","):
            hit = var_by_id.get(vid.strip())
            if not hit:
                continue
            sym, v = hit
            for k in civic_keys(v):
                index[(sym, k)].add(mp)
                first_var.setdefault((sym, k), v)
    ev_by_mp, as_by_mp = collections.defaultdict(list), collections.defaultdict(list)
    for e in E:
        if e.get("evidence_status", "accepted") == "accepted":
            ev_by_mp[e["molecular_profile_id"]].append(e)
    for a in A:
        as_by_mp[a["molecular_profile_id"]].append(a)

    out = collections.defaultdict(dict)
    for g, vs in genes.items():
        for x in vs:
            mps = sorted(index.get((g, x["hgvs"]), set()))
            if not mps:
                continue
            items, asserts, counts = [], [], collections.Counter()
            for mp in mps:
                for e in ev_by_mp.get(mp, []):
                    counts[e["evidence_level"]] += 1
                    items.append({"profile": mp_name.get(mp, e["molecular_profile"]),
                                  "level": e["evidence_level"], "type": e["evidence_type"],
                                  "direction": e["evidence_direction"], "significance": e["significance"],
                                  "disease": e["disease"], "therapies": e["therapies"],
                                  "origin": e["variant_origin"], "rating": e.get("rating") or "",
                                  "id": e["evidence_id"], "url": e["evidence_civic_url"]})
                for a in as_by_mp.get(mp, []):
                    asserts.append({"profile": mp_name.get(mp, a["molecular_profile"]),
                                    "amp": a["amp_category"], "type": a["assertion_type"],
                                    "direction": a["assertion_direction"], "significance": a["significance"],
                                    "disease": a["disease"], "therapies": a["therapies"],
                                    "fda": a.get("regulatory_approval") or "", "nccn": a.get("nccn_guideline") or "",
                                    "id": a["assertion_id"], "url": a["assertion_civic_url"]})
            # 1. curated assertion tier
            best = None
            for a in asserts:
                m = AMP_CAT.search(a["amp"] or "")
                if m and a["direction"] == "Supports":
                    t = (m.group(1).upper(), m.group(2).upper(), "assertion", a)
                    if best is None or tier_rank(t[0], t[1]) < tier_rank(best[0], best[1]):
                        best = t
            # 2. else the best qualifying evidence item
            if best is None:
                for it in items:
                    if (it["type"] in TIER_TYPES and it["direction"] == "Supports"
                            and it["origin"] not in GERMLINE and it["level"] in LEVEL_TIER):
                        t, lv = LEVEL_TIER[it["level"]]
                        if best is None or tier_rank(t, lv) < tier_rank(best[0], best[1]):
                            best = (t, lv, "evidence", it)
            items.sort(key=lambda it: (it["level"], -int(it["rating"] or 0), it["type"], it["disease"]))
            v0 = first_var[(g, x["hgvs"])]
            rec = {"civic_variant": v0["variant"], "civic_url": v0["variant_civic_url"],
                   "evidence_levels": dict(sorted(counts.items())), "n_evidence": sum(counts.values()),
                   "assertions": asserts, "evidence": items[:MAX_ITEMS]}
            if best:
                rec.update({"tier": best[0], "level": best[1], "source": best[2],
                            "basis": {k: best[3].get(k, "") for k in
                                      ("profile", "disease", "therapies", "significance", "type", "id", "url")}})
            out[g][x["hgvs"]] = rec
    return dict(sorted(out.items()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT),
                    help="run tree holding VARIANT/ and databases/ (default: next to this script)")
    ap.add_argument("--refresh", action="store_true", help=f"download CIViC even if a copy younger "
                                                         f"than {CACHE_DAYS} days is on disk")
    args = ap.parse_args()
    root = Path(args.root).resolve()
    vfile = root / "VARIANT" / "variants.json"
    if not vfile.exists():
        raise SystemExit(f"no {vfile} -- run variants.py (step 4b) first")
    genes = json.loads(vfile.read_text(encoding="utf-8"))["genes"]
    tables = download(root / CIVIC_DIR, args.refresh)
    ann = annotate(tables, genes, hgnc_by_entrez(root))
    n_civic = sum(len(v) for v in ann.values())

    # 3. every variant without Tier I/II evidence: Tier IV if common in gnomAD, else Tier III
    af_cache = fetch_af({x["rsid"] for vs in genes.values() for x in vs if x.get("rsid")}, root / AF_FILE)
    for g, vs in genes.items():
        for x in vs:
            rec = ann.setdefault(g, {}).setdefault(x["hgvs"], {})
            if rec.get("tier") in ("I", "II"):
                continue
            t, pop = population_tier(x, g, af_cache)
            rec.update({"tier": t, "level": None, "source": "population", **pop})
    ann = dict(sorted(ann.items()))

    n_var = sum(len(v) for v in genes.values())
    label = lambda r: f"{r['tier']}-{r['level']}" if r.get("level") else r["tier"]
    order = ("I-A", "I-B", "II-C", "II-D", "III", "IV")
    tiers = collections.Counter(label(r) for v in ann.values() for r in v.values() if r.get("tier"))
    src = collections.Counter(r["source"] for v in ann.values() for r in v.values() if r.get("tier"))
    why3 = collections.Counter(("no rsID" if "no rsID" in r["af_basis"] else "not in gnomAD" if "not in gnomAD" in r["af_basis"]
                                else "rare in gnomAD" if r.get("af") is not None else "rsID not resolved")
                               for v in ann.values() for r in v.values() if r.get("tier") == "III")
    out = {"source": "CIViC nightly (CC0), https://civicdb.org; gnomAD via https://myvariant.info",
           "built": datetime.date.today().isoformat(),
           "civic_files": {n: len(tables[n]) for n in CIVIC_FILES},
           "tier_rule": "CIViC assertion AMP category if any, else best accepted supporting "
                        "Predictive/Diagnostic/Prognostic evidence of non-germline origin: "
                        "A->I-A, B->I-B, C->II-C, D/E->II-D; tier holds in the evidence's disease. "
                        f"Without such evidence: Tier IV if the variant's allele has a gnomAD popmax "
                        f">= {AF_BENIGN}, else Tier III. Tiers III and IV have no levels.",
           "genes": ann}
    (root / OUT).write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"{n_var:,} corpus variants, {n_civic:,} in CIViC; AMP tiers: "
          + ", ".join(f"Tier {k} {tiers[k]:,}" for k in order if tiers[k])
          + f" (Tier I/II from an assertion {src['assertion']:,}, from evidence {src['evidence']:,})")
    print("  Tier III because: " + ", ".join(f"{k} {n:,}" for k, n in why3.most_common()))
    common = sorted(((r["af"], g, h) for g, v in ann.items() for h, r in v.items() if r.get("tier") == "IV"), reverse=True)
    print("  Tier IV, most common: " + ", ".join(f"{g} {h} ({af:.2g})" for af, g, h in common[:10]))
    top = sorted(((tier_rank(r["tier"], r["level"]), -r["n_evidence"], g, h) for g, v in ann.items()
                  for h, r in v.items() if r.get("tier") in ("I", "II")))[:12]
    print("best supported: " + ", ".join(f"{g} {h} ({ann[g][h]['tier']}-{ann[g][h]['level']}, "
                                         f"{ann[g][h]['n_evidence']} items)" for _, _, g, h in top))
    print(f"-> {root / OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
