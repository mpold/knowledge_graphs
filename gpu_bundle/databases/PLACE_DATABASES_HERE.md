# Ontology databases (place here)

These files are **required by step 2** (`gpu.py`) but are **not shipped in the repo** —
they are 30 MB – 500 MB each (ChEBI alone exceeds GitHub's 100 MB/file limit) and are
`.gitignore`d. Put them in this directory before running step 2:

| File | Ontology | Used by |
|------|----------|---------|
| `hgnc_complete_set_2026-05-01.json` | HGNC gene symbols | `roman.py`, `greek.py`, `controls.py`, `nonchemical.py`, `target_pharm.py` |
| `mondo-clingen.json` | MONDO disease ontology | `disease.py` |
| `chebi.json` | ChEBI chemical ontology | `chemical.py`, `target_pharm.py` |

Optional:

- `pmc_years.json` — produced by step 15 (`pub_years.py`); drop it in only to run the
  publication-year lookup with the internet off.

- `interactions.tsv` — DGIdb interactions TSV (open drug–gene targets; typed interactions),
  read by `chemical.py` to recover drugs the NER misses **and** ChEBI does not carry
  (e.g. bevacizumab). If it is absent the drug set is simply empty — no error.

- `ChimerKB4.xlsx`, `ChimerSeq4.xlsx`, `Recurrent_table.xlsx` — ChimerDB 4.0 fusion catalogue,
  from https://www.kobic.re.kr/chimerdb/download. Stage-3 only, so these belong under the DATA
  ROOT's `databases/` (e.g. `kaggle_working/databases/`), not here. Convert them once with
  `python chimerdb_to_tsv.py`, which writes `.tsv` beside each workbook using the standard
  library only — nothing downstream reads .xlsx and nothing needs openpyxl.

  ChimerKB (3,138 rows) is manually curated; ChimerSeq (132,979 rows over 122 cancer types) is
  called from TCGA RNA-seq and is algorithmic, so it carries false positives: it lists TP53,
  CTNNB1 and USP39 as fusion partners, none of which fuse in any meaningful sense. Gate on
  ChimerKB, and treat ChimerSeq as supporting evidence with a frame and read-count threshold.
  Only the 32,172 `In-Frame` ChimerSeq rows can produce a chimeric PROTEIN at all; the other
  100,807 are out-of-frame or UTR-CDS.

- `CRISPRGeneDependency.csv` — DepMap CRISPR dependency matrix (~421 MB, cell lines x ~18,000
  genes, values are the probability that the line depends on the gene). Stage-3 only, so it
  belongs under the DATA ROOT's `databases/`. Reduce it once with `python depmap_to_tsv.py`,
  which writes the per-gene `depmap_dependency.tsv` the pipeline actually reads; nothing
  downstream parses the matrix. Available from DepMap's downloads page, or from the release
  article on figshare (24Q4 Public is article 27993248, file 51064631).

On Kaggle, upload these alongside the scripts as part of the dataset. See the repo README
and `step_2_triples.html` for the full step-2 setup.


