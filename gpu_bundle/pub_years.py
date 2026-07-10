#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pub_years.py -- annotate each triple in TRIPLES/genetic_genetic.json with the
publication year of its source article.

The local corpus stores only the PMC accession ("pmid", e.g. "PMC10006201"); it has no
publication date. This script resolves each unique accession to a year through NCBI's
public E-utilities (esummary, db=pmc), caches the lookups in databases/pmc_years.json
(so it is only fetched once and relationships.py can reuse it offline), and writes a
"year" field (int, or null when NCBI returns no date) onto every triple.

Run::  python pub_years.py
"""

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
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


def main():
    data = json.loads(GG.read_text(encoding="utf-8"))
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}

    needed = sorted({accession(t["pmid"]) for t in data})
    missing = [a for a in needed if a not in cache]
    print(f"{len(data):,} triples; {len(needed):,} unique articles; "
          f"{len(missing):,} to fetch ({len(needed) - len(missing):,} cached)")

    if missing:
        fetch_years(missing, cache)        # updates + flushes `cache` to disk per batch
        print(f"cache -> {CACHE} ({len(cache):,} accessions)")

    resolved = 0
    for t in data:
        t["year"] = cache.get(accession(t["pmid"]))
        resolved += t["year"] is not None
    GG.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    years = [t["year"] for t in data if t["year"]]
    print(f"annotated {resolved:,}/{len(data):,} triples with a year "
          f"({len(data) - resolved:,} unresolved)")
    if years:
        print(f"year range: {min(years)}-{max(years)}")
        from collections import Counter
        print("by year:", dict(sorted(Counter(years).items())))


if __name__ == "__main__":
    main()
