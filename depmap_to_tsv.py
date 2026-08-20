#!/usr/bin/env python3
"""depmap_to_tsv.py -- reduce the DepMap CRISPR matrix to one row per gene.

The graph can say a tumour DEPENDS on a gene, on the strength of what a paper wrote. DepMap can
say how OFTEN, across ~1,100 cell lines, and -- the part that matters more -- whether the
dependency is SELECTIVE or whether every cell has it. Those are different therapeutic
propositions and the write-up currently conflates them: PRC1, EXOSC10, MARS1 and the proteasome
subunits sit in its "non-oncogene addiction" list beside genuinely selective dependencies, with
nothing to separate a therapeutic window from its absence.

CRISPRGeneDependency.csv is a 421 MB matrix -- rows are cell lines, columns are ~18,000 genes,
values are the PROBABILITY that the line depends on the gene. Nothing downstream should parse
that on every run, so it is reduced once, here, to a per-gene summary the pipeline reads like
any other flat file:

    gene  n_lines  n_dep  frac_dep  mean_prob  class

`class` is the split the analysis actually needs:

    common     >=90% of lines depend on it. No therapeutic window -- these are the ribosome,
               the proteasome, the spliceosome. A "dependency" here is not a target.
    selective  5-90%. The interesting band: some tumours need it and others do not, which is
               what "addiction" is supposed to mean.
    rare       0.5-5%. Real but narrow; a handful of lines.
    none       <0.5%. Either no dependency, or -- the trap -- no line in the panel represents
               the context where it would show. PRDM14's biology is germ-cell and ESC-like, and
               that is barely in DepMap, so a null here is weak evidence of absence.

The matrix is streamed row by row: 421 MB will not fit comfortably in memory alongside a
per-gene tally, and pandas is not a dependency of this project.

CAVEAT WORTH CARRYING FORWARD. Cell lines are not tumours -- no microenvironment, no immune
system, no differentiation state. For the immune and stromal genes in this corpus (CTLA4, CD274,
TIGIT, VSIR, CD40) and for lineage markers (AFP, UPK2, RAG2) a DepMap null means nothing at all.

Run::
    python depmap_to_tsv.py                       # <data-root>/databases/CRISPRGeneDependency.csv
    python depmap_to_tsv.py --data-root kaggle_working --force
    python depmap_to_tsv.py --threshold 0.5
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_ROOT = ROOT / "kaggle_working"
SRC_NAME = "CRISPRGeneDependency.csv"
OUT_NAME = "depmap_dependency.tsv"
# DepMap's own convention: a probability above this is called a dependency
DEP_THRESHOLD = 0.5
COMMON, SELECTIVE, RARE = 0.90, 0.05, 0.005


def gene_symbol(col):
    """'PRDM14 (63978)' -> 'PRDM14'. The header carries the Entrez id in parentheses."""
    return col.split(" (", 1)[0].strip()


def reduce_matrix(src, threshold=DEP_THRESHOLD, progress=True):
    """Stream the matrix, tallying per gene: lines scored, lines dependent, probability sum."""
    with open(src, encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        genes = [gene_symbol(c) for c in header[1:]]      # column 0 is the model/cell-line id
        n = len(genes)
        scored = [0] * n
        dep = [0] * n
        total = [0.0] * n
        lines = 0
        for row in rd:
            lines += 1
            if progress and lines % 200 == 0:
                print(f"    {lines} cell lines ...", file=sys.stderr)
            for i, v in enumerate(row[1:]):
                if not v:                                  # not screened in this line
                    continue
                try:
                    p = float(v)
                except ValueError:
                    continue
                scored[i] += 1
                total[i] += p
                if p >= threshold:
                    dep[i] += 1
    return genes, scored, dep, total, lines


def classify(frac):
    if frac >= COMMON:
        return "common"
    if frac >= SELECTIVE:
        return "selective"
    if frac >= RARE:
        return "rare"
    return "none"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=str(DEFAULT_ROOT),
                    help="pipeline output tree holding databases/ (default: kaggle_working)")
    ap.add_argument("--threshold", type=float, default=DEP_THRESHOLD,
                    help=f"probability at or above which a line counts as dependent "
                         f"(default {DEP_THRESHOLD}, DepMap's own convention)")
    ap.add_argument("--force", action="store_true", help="rebuild even if the TSV is newer")
    ap.add_argument("--quiet", action="store_true", help="no progress on stderr")
    args = ap.parse_args()

    db = Path(args.data_root) / "databases"
    src, dst = db / SRC_NAME, db / OUT_NAME
    if not src.exists():
        raise SystemExit(f"{src} not found. It is the DepMap CRISPR dependency matrix "
                         f"(~421 MB); see databases/PLACE_DATABASES_HERE.md")
    if dst.exists() and not args.force and dst.stat().st_mtime >= src.stat().st_mtime:
        print(f"{dst.name} is up to date")
        return 0

    genes, scored, dep, total, lines = reduce_matrix(src, args.threshold, not args.quiet)
    with open(dst, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "n_lines", "n_dep", "frac_dep", "mean_prob", "class"])
        counts = {}
        for i, g in enumerate(genes):
            if not scored[i]:
                continue
            frac = dep[i] / scored[i]
            cls = classify(frac)
            counts[cls] = counts.get(cls, 0) + 1
            w.writerow([g, scored[i], dep[i], f"{frac:.4f}",
                        f"{total[i]/scored[i]:.4f}", cls])
    print(f"{lines:,} cell lines x {len(genes):,} genes -> {dst.name}")
    print("  " + ", ".join(f"{k} {v:,}" for k, v in
                           sorted(counts.items(), key=lambda kv: -kv[1])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
