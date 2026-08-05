#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pub_years.py -- resolve the publication year of every article the pipeline cites.

The local corpus stores only the PMC accession ("pmid", e.g. "PMC10006201"); it has no
publication date. This script resolves each unique accession to a year through NCBI's
public E-utilities (esummary, db=pmc) and caches the lookups in databases/pmc_years.json,
which is the real product: relationships.py and high_confidence_g.py both read the cache
offline. It also writes a "year" field (int, or null when NCBI returns no date) onto every
triple of the files given to --annotate.

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
"""

import argparse
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# every triples file the downstream steps read, widest first; absent ones are skipped
SOURCES = ["TRIPLES/triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json",
           "TRIPLES/triples_re_GENETIC_DISEASE_normalized.json",
           "TRIPLES/triples_re_GENETIC_normalized.json",
           "TRIPLES/triples_re.json",
           "TRIPLES/genetic_genetic.json"]
# annotated in place (the others are inputs to later steps that stamp their own year field)
ANNOTATE = ["TRIPLES/genetic_genetic.json"]
GG = ROOT / "TRIPLES" / "genetic_genetic.json"
CACHE = ROOT / "databases" / "pmc_years.json"

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


def load_triples(path):
    """Triple list from a pipeline JSON (a bare list, or {'triples': [...]})."""
    data = json.loads(path.read_text(encoding="utf-8"))
    return data if isinstance(data, list) else data.get("triples", [])


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
        rows = load_triples(p)
        acc = {accession(t["pmid"]) for t in rows if t.get("pmid")}
        ntrip += len(rows)
        print(f"  {p.name}: {len(rows):,} triples, {len(acc):,} articles "
              f"({len(acc - set(cache)):,} not yet cached)")
        needed |= acc

    missing = sorted(needed - set(cache))
    print(f"{ntrip:,} triples; {len(needed):,} unique articles; "
          f"{len(missing):,} to fetch ({len(needed) - len(missing):,} cached)")

    if missing:
        fetch_years(missing, cache)        # updates + flushes `cache` to disk per batch
        print(f"cache -> {CACHE} ({len(cache):,} accessions)")

    covered = sum(1 for a in needed if cache.get(a))
    print(f"years known for {covered:,}/{len(needed):,} articles "
          f"({len(needed) - covered:,} without a date at NCBI)")

    targets = ([] if args.no_annotate else
               [Path(p) for p in args.annotate] if args.annotate is not None
               else [root / a for a in ANNOTATE])
    for path in targets:
        if not path.exists():
            print(f"  (skipping absent {path.name})")
            continue
        data = load_triples(path)
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
