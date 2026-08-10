#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pub_years.py -- resolve the publication year of every article the pipeline cites.

The product is databases/pmc_years.json, a PMC accession -> year map that relationships.py
and high_confidence_g.py both read offline. It is built from two sources, in this order:

  1. pmids/pmid_pmc_ids.tsv -- the stage-1 PubMed table, which already carries a `year`
     column beside each `pmc_id`. Local, free, and on the reference corpus it covers
     100% of the articles the graph cites.
  2. NCBI E-utilities (esummary, db=pmc) -- a network top-up for whatever the table is
     missing. Skip it with --no-fetch to build the cache entirely offline.

An earlier version of this file opened by asserting that "the local corpus stores only the
PMC accession; it has no publication date", and went to NCBI for every article on that
basis. That was wrong twice over: the stage-1 table above has the years, and 14,825 of the
16,778 corpus documents are PMC JATS XML carrying <article-meta><pub-date> directly. (The
other 1,953 are GROBID TEI converted from PDFs, which genuinely have neither -- that is
presumably where the claim came from.) The cost of believing it was that the cache only
ever held whatever the SOURCES list happened to reach, and on this corpus that left 3,974
of 7,764 graph publications undated.

It also writes a "year" field (int, or null when no source has a date) onto every triple
of the files given to --annotate.

NB the two sources can legitimately disagree by a year, on 1.6% of the reference corpus:
a paper e-published in December 2023 inside the 2024 print collection has two publication
years, and PubMed and esummary do not always pick the same one. The table wins by default
because it is the one the rest of stage 1 already reported against.

SOURCES -- why this reads more than genetic_genetic.json
--------------------------------------------------------
It used to harvest accessions from TRIPLES/genetic_genetic.json alone, which is the
GENE-GENE split. Once the graph moved to --nodes all it drew on the normalized multi-type
triples instead, a much wider set: a paper contributing only a gene-disease or a chemical
relation never appears in the gene-gene file, so it was never looked up and its sentences
came out undated. In one lung corpus that was 6,540 of 10,662 articles (61%) -- and since
the year slider drops undated evidence as soon as the range narrows, a quarter of the
graph's edges silently vanished on the first nudge. So the default is now every triples
file the pipeline produces; missing ones are skipped, and the cache means re-running only
fetches what is new.

Run::  python pub_years.py                       # all triples files under <root>/TRIPLES
       python pub_years.py --triples a.json ...  # explicit sources
       python pub_years.py --no-annotate         # only refresh the cache
       python pub_years.py --no-fetch --no-annotate   # build the cache offline, table only
"""

import argparse
import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# every triples file the downstream steps read, widest first; absent ones are skipped
# relation_extraction.py writes the single "triples_re_" normalized file; the other two
# normalized splits come from triples.py and carry no "_re_" (they are pre-scoring). Naming
# all three "triples_re_" meant two never matched a real file, were dropped by the exists()
# filter below, and their accessions went unfetched.
SOURCES = ["TRIPLES/triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json",
           "TRIPLES/triples_GENETIC_DISEASE_normalized.json",
           "TRIPLES/triples_GENETIC_normalized.json",
           "TRIPLES/triples_re.json",
           "TRIPLES/genetic_genetic.json"]
# annotated in place (the others are inputs to later steps that stamp their own year field)
ANNOTATE = ["TRIPLES/genetic_genetic.json"]
GG = ROOT / "TRIPLES" / "genetic_genetic.json"
CACHE = ROOT / "databases" / "pmc_years.json"
# stage-1 PubMed table (pmid, pmc_id, source_publication, issn, journal_impact_factor, year).
# It sits at the PROJECT root, not under the data root -- which is exactly why stage 3, whose
# every path resolves under --data-root, could not see it and went to the network instead.
PMID_TSV = "pmids/pmid_pmc_ids.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
BATCH = 200           # accessions per esummary request
PAUSE = 0.34          # seconds between requests (NCBI: <= 3/sec without an API key)
TOOL, EMAIL = "normalization", "your-email@example.com"

# NCBI E-utilities intermittently returns transient 5xx/429 errors under load; retry those
# (and network blips) with exponential backoff rather than aborting the whole run.
RETRIES = 5
RETRY_STATUS = {429, 500, 502, 503, 504}


def accession(pmid):
    """Bare PMC accession, dropping any '.grobid.tei' / '.xml' suffix."""
    return pmid.split(".")[0]


def find_tsv(explicit, root):
    """Locate the stage-1 PubMed table, or None.

    It lives beside the pipeline rather than inside the data root, so try the obvious
    places: the data root itself, its parent (the project dir), and the same two relative
    to this script -- which covers both the gpu_bundle/ and kaggle_working/ layouts."""
    if explicit:
        p = Path(explicit)
        return p if p.exists() else None
    for base in (root, root.parent, ROOT, ROOT.parent):
        p = base / PMID_TSV
        if p.exists():
            return p
    return None


def tsv_years(path):
    """PMC accession -> year from the stage-1 table; rows without both are skipped."""
    out = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pmc = (row.get("pmc_id") or "").strip()
            yr = (row.get("year") or "").strip()
            if not pmc or not yr.isdigit():
                continue
            out[pmc if pmc.upper().startswith("PMC") else "PMC" + pmc] = int(yr)
    return out


def year_from(rec):
    for key in ("pubdate", "epubdate", "printpubdate"):
        m = re.search(r"\b(\d{4})\b", rec.get(key, "") or "")
        if m:
            return int(m.group(1))
    return None


def esummary(ids):
    """POST one esummary request, retrying transient NCBI failures with backoff."""
    data = urllib.parse.urlencode(
        {"db": "pmc", "id": ids, "retmode": "json",
         "tool": TOOL, "email": EMAIL}).encode()
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(EUTILS, data=data, timeout=60) as r:
                return json.load(r).get("result", {})
        except urllib.error.HTTPError as e:
            if e.code not in RETRY_STATUS or attempt == RETRIES:
                raise
            reason = f"HTTP {e.code}"
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == RETRIES:
                raise
            reason = str(getattr(e, "reason", e))
        backoff = PAUSE * 2 ** attempt        # 0.68, 1.36, 2.72, 5.44 s ...
        print(f"  {reason}; retry {attempt}/{RETRIES - 1} in {backoff:.1f}s")
        time.sleep(backoff)


def fetch_years(accessions, cache=None):
    """Map bare PMC accession -> year (int or None) via NCBI esummary, in batches.

    When `cache` is given it is updated and flushed to disk after each batch so a later
    failure never discards the accessions already fetched this run."""
    out = cache if cache is not None else {}
    for i in range(0, len(accessions), BATCH):
        chunk = accessions[i:i + BATCH]
        ids = ",".join(a[3:] for a in chunk)        # strip "PMC" -> numeric uid
        res = esummary(ids)
        for uid in res.get("uids", []):
            out["PMC" + uid] = year_from(res[uid])
        print(f"  fetched {min(i + BATCH, len(accessions)):,}/{len(accessions):,}")
        if cache is not None:
            CACHE.parent.mkdir(parents=True, exist_ok=True)
            CACHE.write_text(json.dumps(cache, indent=0) + "\n", encoding="utf-8")
        time.sleep(PAUSE)
    return out


def load_triples(path, chunk=1 << 20):
    """Yield each triple from a pipeline JSON (a bare list, or {'triples': [...]}).

    A generator, not a list: SOURCES resolves to files of 300 MB and more, and parsing one
    of those whole costs several GB of resident dicts -- enough to thrash a 16 GB box, which
    is the whole reason the harvest below only ever needed the pmids. A bare top-level list
    (what every pipeline file actually is) is walked one element at a time with raw_decode(),
    so peak memory is one triple plus the read buffer. This is safe only because the elements
    are objects: a half-read object can never parse as a complete shorter one, so a failed
    decode always means "read more", never "silently truncated". The documented
    {'triples': [...]} wrapper cannot be walked this way and falls back to a whole-file parse.
    """
    with path.open(encoding="utf-8") as fh:
        buf = fh.read(chunk)
        if buf.lstrip()[:1] != "[":                      # wrapped form -- parse it whole
            data = json.loads(buf + fh.read())
            yield from (data if isinstance(data, list) else data.get("triples", []))
            return
        buf = buf.lstrip()[1:]
        dec, eof = json.JSONDecoder(), False
        while True:
            buf = buf.lstrip()
            if not buf and not eof:
                more = fh.read(chunk)
                eof, buf = not more, buf + more
                continue
            if not buf or buf[0] == "]":                 # end of the array (or of the file)
                return
            if buf[0] == ",":
                buf = buf[1:]
                continue
            while True:
                try:
                    obj, end = dec.raw_decode(buf)
                    break
                except ValueError:                       # partial element -- refill and retry
                    if eof:
                        raise
                    more = fh.read(chunk)
                    eof, buf = not more, buf + more
            yield obj
            buf = buf[end:]


def main():
    global CACHE
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT),
                    help="pipeline tree holding TRIPLES/ and databases/ (default: next to this script)")
    ap.add_argument("--triples", nargs="*", default=None,
                    help=f"files to harvest accessions from (default: {' '.join(SOURCES)})")
    ap.add_argument("--annotate", nargs="*", default=None,
                    help=f"files to write the per-triple year field onto (default: {' '.join(ANNOTATE)})")
    ap.add_argument("--no-annotate", action="store_true", help="only refresh the cache")
    ap.add_argument("--pmid-tsv", default=None,
                    help=f"stage-1 PubMed table supplying the years (default: {PMID_TSV} "
                         f"found near the data root or this script)")
    ap.add_argument("--no-tsv", action="store_true",
                    help="ignore the stage-1 table and resolve everything through NCBI")
    ap.add_argument("--no-fetch", action="store_true",
                    help="never call NCBI; build the cache from the table alone (offline)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    CACHE = root / "databases" / "pmc_years.json"
    srcs = ([Path(p) for p in args.triples] if args.triples is not None
            else [root / s for s in SOURCES])
    present = [p for p in srcs if p.exists()]
    for p in srcs:
        if not p.exists():
            print(f"  (skipping absent {p.name})")
    if not present:
        raise SystemExit(f"no triples files found under {root/'TRIPLES'}")

    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    needed, ntrip = set(), 0
    for p in present:
        acc, nrows = set(), 0                  # counted while streaming; the file is never held
        for t in load_triples(p):
            nrows += 1
            if t.get("pmid"):
                acc.add(accession(t["pmid"]))
        ntrip += nrows
        print(f"  {p.name}: {nrows:,} triples, {len(acc):,} articles "
              f"({len(acc - set(cache)):,} not yet cached)")
        needed |= acc

    print(f"{ntrip:,} triples; {len(needed):,} unique articles")

    # 1) the stage-1 table, free and offline. It fills only accessions the cache has no year
    # for -- an absent key, or a null left by an earlier run where NCBI had no date. Existing
    # years are never overwritten, so a rebuild cannot silently move dates under the graph.
    tsv_path = None if args.no_tsv else find_tsv(args.pmid_tsv, root)
    if tsv_path:
        table = tsv_years(tsv_path)
        filled = [a for a in needed if cache.get(a) is None and a in table]
        for a in filled:
            cache[a] = table[a]
        print(f"  {tsv_path.name}: {len(table):,} accessions with a year; "
              f"filled {len(filled):,} the cache was missing")
        if filled:
            CACHE.parent.mkdir(parents=True, exist_ok=True)
            CACHE.write_text(json.dumps(cache, indent=0) + "\n", encoding="utf-8")
    elif not args.no_tsv:
        print(f"  (no {PMID_TSV} found -- falling back to NCBI for everything)")

    # 2) NCBI, for whatever the table could not supply
    missing = sorted(a for a in needed if a not in cache)
    if missing and args.no_fetch:
        print(f"  --no-fetch: leaving {len(missing):,} accessions unresolved")
    elif missing:
        print(f"  {len(missing):,} still unresolved -> NCBI")
        fetch_years(missing, cache)        # updates + flushes `cache` to disk per batch
        print(f"cache -> {CACHE} ({len(cache):,} accessions)")

    covered = sum(1 for a in needed if cache.get(a))
    print(f"years known for {covered:,}/{len(needed):,} articles "
          f"({len(needed) - covered:,} with no date in either source)")

    targets = ([] if args.no_annotate else
               [Path(p) for p in args.annotate] if args.annotate is not None
               else [root / a for a in ANNOTATE])
    for path in targets:
        if not path.exists():
            print(f"  (skipping absent {path.name})")
            continue
        data = list(load_triples(path))    # annotate rewrites the file, so this one is held
        resolved = 0
        for t in data:
            t["year"] = cache.get(accession(t["pmid"])) if t.get("pmid") else None
            resolved += t["year"] is not None
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"annotated {resolved:,}/{len(data):,} triples in {path.name} "
              f"({len(data) - resolved:,} unresolved)")

    years = [y for a in needed if (y := cache.get(a))]
    if years:
        from collections import Counter
        print(f"year range: {min(years)}-{max(years)}")
        print("by year:", dict(sorted(Counter(years).items())))


if __name__ == "__main__":
    main()
