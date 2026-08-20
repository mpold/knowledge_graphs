#!/usr/bin/env python3
"""chimerdb_to_tsv.py -- convert the staged ChimerDB .xlsx releases to TSV.

ChimerDB ships .xlsx; every other database this project stages is a flat file the pipeline can
read with the standard library (hgnc_complete_set_*.json, chebi.json, interactions.tsv). Adding
openpyxl as a runtime dependency for one file would be the wrong trade, so the workbooks are
converted ONCE, here, and the pipeline only ever sees TSV -- diffable, greppable, and readable
with `csv` like the DGIdb interactions file it will sit beside.

An .xlsx is a zip of XML, so no third-party reader is needed:
  xl/sharedStrings.xml   every distinct string in the workbook, referenced by index
  xl/worksheets/*.xml    rows of <c> cells; t="s" means the value is a shared-string index,
                         t="inlineStr" means the text is inline, no t at all means a number
Cells that are empty are OMITTED from the XML rather than written blank, so the column of each
cell is recovered from its A1-style reference -- skipping that step silently shifts every row
with a gap in it, which in ChimerSeq is most of them.

The sheet is streamed with iterparse: ChimerSeq4.xlsx is 15 MB compressed and expands to
several hundred MB of XML, which a DOM parse will not survive.

Run::
    python chimerdb_to_tsv.py                     # convert every ChimerDB xlsx under databases/
    python chimerdb_to_tsv.py --data-root kaggle_working --force
    python chimerdb_to_tsv.py databases/ChimerKB4.xlsx --head 3
"""
import argparse
import csv
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parent
DEFAULT_ROOT = ROOT / "kaggle_working"
CHIMER_FILES = ("ChimerKB4.xlsx", "ChimerSeq4.xlsx", "ChimerPub4.xlsx", "Recurrent_table.xlsx")


def _tag(el):
    """Element tag without its namespace."""
    return el.tag.rsplit("}", 1)[-1]


def shared_strings(zf):
    """The workbook's string table. <si> entries may be split across <r> runs, and the runs
    have to be concatenated or the cell reads as only its first formatting fragment."""
    try:
        raw = zf.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    out = []
    for _, el in ET.iterparse(_bytes_io(raw), events=("end",)):
        if _tag(el) == "si":
            out.append("".join(t.text or "" for t in el.iter() if _tag(t) == "t"))
            el.clear()
    return out


def _bytes_io(b):
    import io
    return io.BytesIO(b)


def col_index(ref):
    """'BC12' -> 54. The letters are base-26 with no zero, hence the -1 per digit."""
    n = 0
    for ch in ref:
        if ch.isdigit():
            break
        n = n * 26 + (ord(ch.upper()) - 64)
    return n - 1


def first_sheet(zf):
    names = [n for n in zf.namelist() if n.startswith("xl/worksheets/sheet")]
    if not names:
        raise SystemExit("no worksheet found in the workbook")
    return sorted(names)[0]


def rows(path):
    """Yield each row of the first worksheet as a list of strings."""
    with zipfile.ZipFile(path) as zf:
        strings = shared_strings(zf)
        with zf.open(first_sheet(zf)) as fh:
            row, width = [], 0
            for event, el in ET.iterparse(fh, events=("end",)):
                tag = _tag(el)
                if tag == "c":
                    ref = el.get("r") or ""
                    idx = col_index(ref) if ref else len(row)
                    ctype = el.get("t")
                    val = ""
                    if ctype == "s":                       # shared-string index
                        v = el.find("./{*}v")
                        if v is not None and v.text and v.text.isdigit():
                            i = int(v.text)
                            val = strings[i] if i < len(strings) else ""
                    elif ctype == "inlineStr":
                        val = "".join(t.text or "" for t in el.iter() if _tag(t) == "t")
                    else:                                   # number, or a formula's cached value
                        v = el.find("./{*}v")
                        val = (v.text or "") if v is not None else ""
                    while len(row) <= idx:                  # empty cells are omitted, not blank
                        row.append("")
                    row[idx] = val.replace("\t", " ").replace("\n", " ").strip()
                    el.clear()
                elif tag == "row":
                    width = max(width, len(row))
                    while len(row) < width:                 # ragged tail -> pad to the header
                        row.append("")
                    yield row
                    row = []
                    el.clear()


def convert(src, dst, head=0):
    n = 0
    with open(dst, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        for r in rows(src):
            w.writerow(r)
            n += 1
            if head and n <= head:
                print("    " + " | ".join(c[:22] for c in r[:8]))
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", help="specific .xlsx files (default: the ChimerDB set "
                                            "found under <data-root>/databases)")
    ap.add_argument("--data-root", default=str(DEFAULT_ROOT),
                    help="pipeline output tree holding databases/ (default: kaggle_working)")
    ap.add_argument("--force", action="store_true", help="re-convert even if the TSV is newer")
    ap.add_argument("--head", type=int, default=0, help="echo the first N rows of each file")
    args = ap.parse_args()

    db = Path(args.data_root) / "databases"
    targets = [Path(f) for f in args.files] or [db / n for n in CHIMER_FILES if (db / n).exists()]
    if not targets:
        raise SystemExit(f"no ChimerDB .xlsx found in {db} -- see databases/PLACE_DATABASES_HERE.md")
    for src in targets:
        dst = src.with_suffix(".tsv")
        if dst.exists() and not args.force and dst.stat().st_mtime >= src.stat().st_mtime:
            print(f"{src.name}: up to date -> {dst.name}")
            continue
        n = convert(src, dst, args.head)
        print(f"{src.name}: {n:,} rows -> {dst.name} ({dst.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
