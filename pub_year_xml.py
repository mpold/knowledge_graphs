#!/usr/bin/env python3
"""pub_year_xml.py -- stamp every corpus XML with its publication year.

Step 6c of the Step 1 publications pipeline: runs after the corpus in
``gpu_bundle/experimental_ner/`` is final (step 6 rebuilds it, 6b re-seeds the
archive hits) and before anything reads it. Each file gets one processing
instruction in its prolog::

    <?xml version="1.0" ?><?pub-year 2023?><!DOCTYPE ...

so the year travels WITH the document. Stage 2 (gpu_bundle/pub_years.py) reads it
back offline; nothing downstream needs the PubMed table or the network to date a
paper. A processing instruction rather than an element because it is valid in both
formats the corpus holds (PMC JATS and GROBID TEI), leaves the document tree
untouched, and every XML parser ignores it.

Year sources, first hit wins:

  1. pmids/pmid_pmc_ids.tsv -- the stage-1 PubMed table (``year`` beside ``pmc_id``).
     Every corpus document came from this table, so it normally covers all of them.
  2. the XML itself: JATS ``<article-meta><pub-date>...<year>``, or the TEI
     ``<publicationStmt><date when="YYYY...">`` GROBID sometimes extracts.

The table wins because it is the year the rest of stage 1 reports against; JATS and
PubMed disagree by a year on ~10% of papers (e-pub in December, print issue in
January). Re-running is safe: an existing stamp is replaced, and a file whose stamp
is already right is not rewritten.

Exits non-zero, listing the files, if any document is left without a year -- the
point of this step is that stage 2 receives a FULLY dated corpus.

Usage
-----
    python pub_year_xml.py              # stamp gpu_bundle/experimental_ner/*.xml
    python pub_year_xml.py --check      # report only, write nothing
    python pub_year_xml.py --xml-dir D --pmid-tsv F
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter

BASE    = os.path.dirname(os.path.abspath(__file__))
XML_DIR = os.path.join(BASE, "gpu_bundle", "experimental_ner")
IN_TSV  = os.path.join(BASE, "pmids", "pmid_pmc_ids.tsv")

STAMP      = re.compile(rb"<\?pub-year\s+(\d{4})\s*\?>")      # also read by gpu_bundle/pub_years.py
DECL       = re.compile(rb"^(\xef\xbb\xbf)?\s*<\?xml[^>]*\?>")  # optional BOM + XML declaration
JATS_YEAR  = re.compile(rb"<pub-date\b[^>]*>.*?<year>\s*(\d{4})\s*</year>", re.S)
TEI_YEAR   = re.compile(rb"<publicationStmt>.*?<date\b[^>]*\bwhen=\"(\d{4})", re.S)
TEI_SUFFIX = ".grobid.tei.xml"


def table_years(path):
    """PMC accession -> year from the stage-1 table; rows without both are skipped."""
    out = {}
    with open(path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pmc = (row.get("pmc_id") or "").strip()
            yr = (row.get("year") or "").strip()
            if pmc and yr.isdigit():
                out[pmc if pmc.upper().startswith("PMC") else "PMC" + pmc] = int(yr)
    return out


def own_year(data, name):
    """The year the document states itself, or None."""
    m = (TEI_YEAR if name.endswith(TEI_SUFFIX) else JATS_YEAR).search(data)
    return int(m.group(1)) if m else None


def stamped(data, year):
    """`data` with exactly one <?pub-year?> stamp, placed right after the declaration."""
    body = STAMP.sub(b"", data, count=1)
    pi = b"<?pub-year %d?>" % year
    m = DECL.match(body)
    return body[:m.end()] + pi + body[m.end():] if m else pi + body


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xml-dir", default=XML_DIR, help="corpus to stamp (default: %(default)s)")
    ap.add_argument("--pmid-tsv", default=IN_TSV, help="stage-1 PubMed table (default: %(default)s)")
    ap.add_argument("--check", action="store_true", help="report coverage only; write nothing")
    args = ap.parse_args()

    if not os.path.isdir(args.xml_dir):
        sys.exit("error: no corpus at %s" % args.xml_dir)
    table = table_years(args.pmid_tsv) if os.path.isfile(args.pmid_tsv) else {}
    print("[pub_year] %s: %d accessions with a year" % (os.path.basename(args.pmid_tsv), len(table))
          if table else "[pub_year] no %s -- falling back to the dates inside the XML" % args.pmid_tsv)

    names = sorted(n for n in os.listdir(args.xml_dir) if n.endswith(".xml"))
    source, written, disagree, undated, years = Counter(), 0, 0, [], Counter()
    for name in names:
        path = os.path.join(args.xml_dir, name)
        with open(path, "rb") as fh:
            data = fh.read()
        own = own_year(data, name)
        year = table.get(name.split(".")[0])
        if year is not None:
            source["table"] += 1
            disagree += own is not None and own != year
        elif own is not None:
            year = own
            source["xml"] += 1
        else:
            undated.append(name)
            continue
        years[year] += 1
        if args.check:
            continue
        new = stamped(data, year)
        if new != data:
            with open(path, "wb") as fh:
                fh.write(new)
            written += 1

    print("[pub_year] %d XML: %d dated from the table, %d from the XML itself, %d undated"
          % (len(names), source["table"], source["xml"], len(undated)))
    if disagree:
        print("[pub_year] %d documents state a different year than PubMed (table kept)" % disagree)
    if years:
        print("[pub_year] year range %d-%d" % (min(years), max(years)))
    print("[pub_year] %s" % ("--check: nothing written" if args.check
                             else "%d files stamped, %d already current" % (written, len(names) - len(undated) - written)))
    if undated:
        print("[pub_year] ERROR: no publication year for %d file(s):" % len(undated), file=sys.stderr)
        for n in undated[:20]:
            print("    " + n, file=sys.stderr)
        if len(undated) > 20:
            print("    ... and %d more" % (len(undated) - 20), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
