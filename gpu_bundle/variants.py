#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""variants.py -- sequence variants per gene, from NCBI PubTator3 (tmVar 3).

The NER of sentences.py tags genes, diseases and chemicals but not sequence variants, and
most variants sit in passages it never selects (methods, tables, figure legends). PubTator3
has already run tmVar 3 over PubMed and the PMC open-access subset: every variant it finds
carries its HGVS form (it normalizes "T790M" and "Thr790Met" to p.T790M), the NCBI Gene id of
the gene it belongs to, and often an rsID. This step reads those annotations for every corpus
document and writes a gene -> variants table that high_confidence_g.py attaches to gene nodes.

  1. PMC full text: /publications/pmc_export/biocjson?pmcids=... in batches of BATCH.
  2. Not in PubTator3's PMC set (no open-access full text, typically the GROBID/PDF files):
     the PubMed abstract instead, /publications/export/biocjson?pmids=... The PMID comes from
     the JATS <article-id pub-id-type="pmid">, else from NCBI's PMC ID converter.
  3. Kept: annotations of type Variant whose HGVS is p., c., g. or r. and whose gene maps to an
     approved HGNC symbol (via HGNC entrez_id). Sites without a change ("Y356") and rsID-only
     mentions are left out, and so are non-human genes, as their NCBI Gene ids are not in HGNC.
     tmVar's own species tag is NOT used: it follows the paper's model organism, so EGFR T790M
     in a mouse xenograft comes tagged 10090 while its gene id is still human EGFR (1956).

Every fetched document is appended to VARIANT/pubtator3_variants.jsonl as it arrives, so an
interrupted run resumes where it stopped; --refresh refetches everything. Needs internet;
NCBI asks for at most 3 requests a second, and this stays under that.

Inputs : experimental_ner/*.xml        document ids (PMC accessions) and JATS PMIDs
         databases/hgnc_complete_set_*.json   NCBI Gene id -> HGNC symbol
Outputs: VARIANT/pubtator3_variants.jsonl     one line per document (the resumable cache)
         VARIANT/variants.json                gene symbol -> variants, read by high_confidence_g.py

Run::  python variants.py [--root DIR] [--xml-dir DIR] [--refresh]
"""

import argparse
import collections
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
RAW_DIR = "experimental_ner"
OUT_DIR = "VARIANT"
CACHE = "pubtator3_variants.jsonl"
TABLE = "variants.json"
API = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications"
IDCONV = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
UA = {"User-Agent": "relationship-graphs/1.0"}
BATCH = 20          # ids per request: larger batches came back with documents missing
PAUSE = 0.4         # seconds between requests (NCBI: at most 3 a second)
LEVELS = ("p", "c", "g", "r")
# PubTator3 gene ids that are wrong for this corpus, by NCBI Gene id -> the right one. SLTM's
# alias is "Met", so the MET kinase-domain mutations (D1228N, Y1230C, H1094Y ...) of MET-amplified
# NSCLC papers land on SLTM.
GENE_FIX = {"79811": "4233"}      # SLTM -> MET
# a zero-padded position ("p.C0157S", "p.P0013C") is a sample or patient code, not a variant
PADDED = re.compile(r"^[pcgr]\.\(?[A-Za-z*]*0\d")
PMID_RE = re.compile(rb'<article-id pub-id-type="pmid">\s*(\d+)\s*</article-id>')


def get_json(url, tries=4):
    """GET url as JSON; None on HTTP 400/404 (PubTator3's answer for ids it does not hold)."""
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            if k == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if k == tries - 1:
                raise
        time.sleep(2 ** (k + 1))


def docs_of(payload):
    """The BioC documents of a PubTator3 biocjson answer."""
    if not payload:
        return []
    if isinstance(payload, dict):
        return payload.get("PubTator3") or payload.get("documents") or []
    return payload


def pmc_of(doc):
    v = str(doc.get("pmcid") or "")
    return ("PMC" + v.replace("PMC", "")) if v else None


def variants_in(doc):
    """[(text, infons), ...] of every Variant annotation in one BioC document."""
    out = []
    for ps in doc.get("passages", []):
        for a in ps.get("annotations", []):
            inf = a.get("infons") or {}
            if inf.get("type") == "Variant":
                out.append((a.get("text") or "", inf))
    return out


def slim(text, inf):
    """The fields of one annotation this step keeps."""
    gids = inf.get("gene_ids") or ([inf["gene_id"]] if inf.get("gene_id") is not None else [])
    return {"text": text, "hgvs": inf.get("hgvs") or "", "subtype": inf.get("subtype") or "",
            "gene_ids": [str(g) for g in gids], "rsid": inf.get("rsid") or "",
            "species": str(inf.get("species_id") or "")}


def jats_pmids(xml_dir):
    """PMC accession -> PMID from the JATS headers (GROBID TEI files carry none)."""
    out = {}
    for f in xml_dir.glob("*.xml"):
        with f.open("rb") as fh:
            m = PMID_RE.search(fh.read(65536))
        if m:
            out[f.name.split(".")[0]] = m.group(1).decode()
    return out


def idconv_pmids(pmcs):
    """PMC accession -> PMID from NCBI's PMC ID converter, for documents without one."""
    out = {}
    pmcs = sorted(pmcs)
    for i in range(0, len(pmcs), 100):
        q = urllib.parse.urlencode({"ids": ",".join(pmcs[i:i + 100]), "format": "json",
                                    "tool": "relationship-graphs"})
        d = get_json(f"{IDCONV}?{q}") or {}
        for r in d.get("records", []):
            if r.get("pmcid") and r.get("pmid"):
                out[r["pmcid"]] = str(r["pmid"])
        time.sleep(PAUSE)
    return out


def hgnc_by_entrez(root):
    """NCBI Gene id -> approved HGNC symbol."""
    hits = sorted(glob.glob(str(root / "databases" / "hgnc_complete_set_*.json")))
    if not hits:
        raise SystemExit(f"no databases/hgnc_complete_set_*.json under {root}")
    docs = json.load(open(hits[-1], encoding="utf-8"))["response"]["docs"]
    return {str(d["entrez_id"]): d["symbol"] for d in docs
            if d.get("entrez_id") and d.get("status", "Approved") == "Approved"}


def load_cache(path):
    done = {}
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                done[r["pmc"]] = r
    return done


def fetch(ids, cache_path, done, pmid_of):
    """Fetch every id not yet in `done`, appending each result to the cache."""
    todo = [p for p in ids if p not in done]
    if not todo:
        return
    print(f"  PubTator3 full text: {len(todo):,} document(s) to fetch, {len(done):,} cached", flush=True)
    with cache_path.open("a", encoding="utf-8") as out:
        def keep(rec):
            done[rec["pmc"]] = rec
            out.write(json.dumps(rec, ensure_ascii=False) + "\n")

        missing = []
        for i in range(0, len(todo), BATCH):
            batch = todo[i:i + BATCH]
            got = {}
            for doc in docs_of(get_json(f"{API}/pmc_export/biocjson?pmcids={','.join(batch)}")):
                p = pmc_of(doc)
                if p in batch:
                    got[p] = [slim(t, inf) for t, inf in variants_in(doc)]
            for p in batch:
                if p in got:
                    keep({"pmc": p, "src": "pmc", "variants": got[p]})
                else:
                    missing.append(p)
            time.sleep(PAUSE)
            if (i // BATCH) % 25 == 24:
                print(f"    {i + len(batch):,}/{len(todo):,} asked, {len(missing):,} not in the PMC set",
                      flush=True)
        if not missing:
            return
        need = [p for p in missing if p not in pmid_of]
        if need:
            pmid_of.update(idconv_pmids(need))
        with_pmid = [p for p in missing if p in pmid_of]
        print(f"  PubMed abstracts for the {len(missing):,} not in the PMC set "
              f"({len(with_pmid):,} have a PMID)", flush=True)
        for p in missing:
            if p not in pmid_of:
                keep({"pmc": p, "src": "none", "variants": []})
        by_pmid = {pmid_of[p]: p for p in with_pmid}
        pmids = sorted(by_pmid)
        for i in range(0, len(pmids), BATCH):
            batch = pmids[i:i + BATCH]
            got = {}
            for doc in docs_of(get_json(f"{API}/export/biocjson?pmids={','.join(batch)}")):
                q = str(doc.get("pmid") or doc.get("id") or "")
                if q in by_pmid:
                    got[q] = [slim(t, inf) for t, inf in variants_in(doc)]
            for q in batch:
                keep({"pmc": by_pmid[q], "src": "abstract" if q in got else "none",
                      "pmid": q, "variants": got.get(q, [])})
            time.sleep(PAUSE)


def level_of(hgvs):
    """'p' | 'c' | 'g' | 'r' for an HGVS string, else None."""
    m = re.match(r"([pcgr])\.", hgvs or "")
    return m.group(1) if m else None


# HGVS variant types, in the order summaries list them. Protein changes split the substitution by
# its consequence; DNA and RNA keep HGVS's own classes.
P_TYPES = ("missense", "nonsense", "synonymous", "unspecified (X)", "deletion", "insertion",
           "duplication", "delins", "frameshift", "extension", "other")
NA_TYPES = ("substitution", "deletion", "insertion", "duplication", "delins", "inversion", "other")
_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
        "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
        "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
        "Xaa": "X", "Sec": "U", "Pyl": "O"}
# one residue replaced -- by one residue, or by several alternatives ("p.E255K,V": K or V)
_P_SUB = re.compile(r"^p\.\(?([A-Z][a-z]{2}|[A-Z*])?(\d+)((?:[A-Z][a-z]{2}|[A-Z*=?])(?:,(?:[A-Z][a-z]{2}|[A-Z*]))*)\)?$")


def hgvs_type(hgvs):
    """The HGVS variant type of one p./c./g./r. string, tolerant of tmVar's shorthand
    ("p.E19,746del", "c.503,504insY,AY", "p.V14fsX", "p.G719X")."""
    h = hgvs or ""
    lv = level_of(h)
    if not lv:
        return "other"
    body = h[2:]
    low = body.lower()
    if "fs" in low and lv == "p":
        return "frameshift"
    if "del" in low and "ins" in low:
        return "delins"
    if "ext" in low and lv == "p":
        return "extension"
    for word, name in (("dup", "duplication"), ("ins", "insertion"), ("del", "deletion"),
                       ("inv", "inversion")):
        if word in low and (lv != "p" or word != "inv"):
            return name
    if lv != "p":
        return "substitution" if ">" in body else "other"
    m = _P_SUB.match(h)
    if not m:
        return "other"
    ref = _AA3.get(m.group(1), m.group(1)) if m.group(1) else None
    alts = [_AA3.get(x, x) for x in m.group(3).split(",")]
    if all(a == "=" or (ref and a == ref) for a in alts):
        return "synonymous"
    if all(a == "*" for a in alts):
        return "nonsense"
    # X is "any residue" in current HGVS but a stop in older papers ("p.R216X"): left open
    if any(a in ("X", "?") for a in alts):
        return "unspecified (X)"
    return "missense"


def build_table(done, ids, symbol_of):
    """gene symbol -> [variant, ...], most-reported first; plus the run's counts."""
    agg = collections.defaultdict(lambda: {"papers": set(), "mentions": 0, "surfaces": collections.Counter(),
                                           "rsid": "", "subtype": ""})
    drop = collections.Counter()
    for p in ids:
        rec = done.get(p)
        if not rec:
            continue
        for v in rec["variants"]:
            lv = level_of(v["hgvs"])
            if not lv:
                drop["no p./c./g./r. HGVS"] += 1
                continue
            if PADDED.match(v["hgvs"]):
                drop["zero-padded position (a sample code)"] += 1
                continue
            syms = sorted({symbol_of[g] for g in (GENE_FIX.get(x, x) for x in v["gene_ids"])
                           if g in symbol_of})
            if not syms:
                drop["gene not an HGNC symbol"] += 1
                continue
            for s in syms:
                a = agg[(s, v["hgvs"])]
                a["papers"].add(p)
                a["mentions"] += 1
                a["surfaces"][v["text"]] += 1
                a["rsid"] = a["rsid"] or v["rsid"]
                a["subtype"] = a["subtype"] or v["subtype"]
    genes = collections.defaultdict(list)
    for (s, h), a in agg.items():
        genes[s].append({"hgvs": h, "level": level_of(h), "type": hgvs_type(h), "papers": sorted(a["papers"]),
                         "mentions": a["mentions"], "rsid": a["rsid"], "subtype": a["subtype"],
                         "surfaces": [w for w, _ in a["surfaces"].most_common(5)]})
    for v in genes.values():
        v.sort(key=lambda x: (-len(x["papers"]), -x["mentions"], x["hgvs"]))
    return dict(sorted(genes.items())), drop


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT),
                    help="run tree holding databases/ and receiving VARIANT/ (default: next to this script)")
    ap.add_argument("--xml-dir", default=None, help=f"corpus (default: <root>/{RAW_DIR})")
    ap.add_argument("--refresh", action="store_true", help="ignore the cache and refetch every document")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    xml_dir = Path(args.xml_dir).resolve() if args.xml_dir else root / RAW_DIR
    ids = sorted({f.name.split(".")[0] for f in xml_dir.glob("*.xml") if f.name.startswith("PMC")})
    if not ids:
        raise SystemExit(f"no PMC*.xml in {xml_dir}")
    out_dir = root / OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_path = out_dir / CACHE
    if args.refresh and cache_path.exists():
        cache_path.unlink()
    symbol_of = hgnc_by_entrez(root)
    done = load_cache(cache_path)
    fetch(ids, cache_path, done, jats_pmids(xml_dir))

    genes, drop = build_table(done, ids, symbol_of)
    src = collections.Counter(done[p]["src"] for p in ids if p in done)
    n_var = sum(len(v) for v in genes.values())
    lv = collections.Counter(x["level"] for v in genes.values() for x in v)
    table = {"source": "NCBI PubTator3 (tmVar 3)", "built": datetime.date.today().isoformat(),
             "documents": {"corpus": len(ids), "full_text": src["pmc"], "abstract_only": src["abstract"],
                           "not_in_pubtator3": src["none"]},
             "genes": genes}
    (out_dir / TABLE).write_text(json.dumps(table, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")

    print(f"{len(ids):,} documents: full text {src['pmc']:,}, abstract only {src['abstract']:,}, "
          f"not in PubTator3 {src['none']:,}")
    print(f"{n_var:,} distinct variant(s) on {len(genes):,} gene(s): "
          + ", ".join(f"{k}. {lv[k]:,}" for k in LEVELS))
    for k in LEVELS:
        ty = collections.Counter(x["type"] for v in genes.values() for x in v if x["level"] == k)
        if ty:
            print(f"  {k}. by type: " + ", ".join(f"{t} {ty[t]:,}" for t in (P_TYPES if k == "p" else NA_TYPES)
                                               if ty[t]))
    if drop:
        print("left out (mentions): " + ", ".join(f"{k} {n:,}" for k, n in drop.most_common()))
    top = sorted(((len(x["papers"]), g, x["hgvs"]) for g, v in genes.items() for x in v), reverse=True)[:10]
    print("most reported: " + ", ".join(f"{g} {h} ({n} papers)" for n, g, h in top))
    print(f"-> {out_dir / TABLE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
