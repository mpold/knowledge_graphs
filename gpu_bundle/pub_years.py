#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""pub_years.py -- collect the publication year of every corpus document.

The product is databases/pmc_years.json, a PMC accession -> year map that relationships.py,
triples_per_year.py and high_confidence_g.py read. It is read straight off the corpus: stage 1
(pub_year_xml.py, step 6c) stamps every experimental_ner/ XML with a processing instruction

    <?pub-year 2023?>

in its prolog, so this step only harvests those stamps -- offline, no PubMed table, no NCBI.
It runs before sentences.py, so the cache exists before any step that reads it.

A corpus file without a stamp was not produced by the current stage 1 (or 6c was skipped):
the step writes the years it did find, lists the unstamped files and exits non-zero, so the
orchestrator reports the gap. Fix it at the source with ``python pub_year_xml.py``.

When experimental_ner/ is not present (a run reusing a prebuilt sentences/), an uploaded
databases/pmc_years.json is kept as-is.

Run::  python pub_years.py [--root DIR] [--xml-dir DIR]
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW_DIR = "experimental_ner"
CACHE = "databases/pmc_years.json"
STAMP = re.compile(rb"<\?pub-year\s+(\d{4})\s*\?>")   # written by pub_year_xml.py (stage 1, step 6c)
HEAD = 4096   # the stamp sits right after the XML declaration


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(ROOT),
                    help="run tree holding experimental_ner/ and databases/ (default: next to this script)")
    ap.add_argument("--xml-dir", default=None, help=f"stamped corpus (default: <root>/{RAW_DIR})")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    xml_dir = Path(args.xml_dir) if args.xml_dir else root / RAW_DIR
    cache = root / CACHE
    if not xml_dir.is_dir():
        if cache.exists():
            print(f"  no {xml_dir} -- keeping the uploaded {CACHE}")
            return 0
        raise SystemExit(f"no corpus at {xml_dir} and no {CACHE} to fall back on")

    years, unstamped = {}, []
    for path in sorted(xml_dir.glob("*.xml")):
        with path.open("rb") as fh:
            m = STAMP.search(fh.read(HEAD))
        if m:
            years[path.name.split(".")[0]] = int(m.group(1))   # "PMC123.grobid.tei.xml" -> "PMC123"
        else:
            unstamped.append(path.name)

    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(years, indent=0, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{len(years):,} dated documents -> {cache}")
    if years:
        print(f"year range: {min(years.values())}-{max(years.values())}")
        print("by year:", dict(sorted(Counter(years.values()).items())))
    if unstamped:
        print(f"ERROR: {len(unstamped):,} XML carry no <?pub-year?> stamp (run pub_year_xml.py "
              f"in stage 1): {', '.join(unstamped[:10])}{' ...' if len(unstamped) > 10 else ''}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
