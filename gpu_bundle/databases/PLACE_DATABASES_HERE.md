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

- `pmc_years.json` — produced by step 14 (`pub_years.py`); drop it in only to run the
  publication-year lookup with the internet off.

- 'interactions.tsv' - DGIdb interactions TSV (open drug-gene targets; typed interactions)

On Kaggle, upload these alongside the scripts as part of the dataset. See the repo README
and `step_2_triples.html` for the full step-2 setup.
