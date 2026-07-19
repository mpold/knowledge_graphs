# cancer_knowledge_graph

A biomedical relation-extraction pipeline: from a single **PubMed query** to an interactive
**gene–gene relationship graph** for a disease/chemical context, in three stages.

| Stage | What | Where it runs | Entry point |
|------:|------|---------------|-------------|
| **1** | Publications → full-text NER corpus | local (network + Docker/GROBID) | 7 root scripts, orchestrated by `run_pipeline.py` |
| **2** | NER corpus → normalized, model-scored relation **triples** | GPU (Kaggle or local) | `gpu_bundle/gpu.py` (16-step chain) |
| **3** | Triples → high-confidence gene–gene **graph** | local | `high_confidence.py` |

Each stage hands off to the next **by files**. Rendered walk-throughs of every stage ship with
this bundle: [`step_1_publications.html`](step_1_publications.html),
[`step_2_triples.html`](step_2_triples.html), [`step_3_graph.html`](step_3_graph.html), and a
consolidated [`requirements.html`](requirements.html).

---

## Quick start (all three stages, locally)

```bash
# 1. install the light local deps (stage 1 + stage 3)
python -m pip install -r requirements.txt

# 2. (stage 2 only) GPU deps — normally you run stage 2 on Kaggle instead; see below
python -m pip install -r gpu_bundle/requirements.txt

# 3. run the whole pipeline
python run_pipeline.py --query "pancreatic cancer"
```

`run_pipeline.py` runs the stages in order and aborts on the first hard failure. Run a
subset with `--steps` (e.g. `--steps 1`, `--steps 2,3`).

> **Realistically, stage 2 runs on Kaggle**, not your laptop — it needs a CUDA GPU and the
> large ontology databases. The typical flow is: **stage 1 locally → stage 2 on Kaggle →
> download `kaggle_working.zip`, unzip it here → stage 3 locally.** See below.

---

## Requirements at a glance

| | Stage 1 | Stage 2 | Stage 3 |
|--|--|--|--|
| Python | 3.9+ | 3.9+ | 3.9+ |
| Packages | `requests` | `torch`, `transformers`, `datasets<4`, `numpy`, `lxml` | *stdlib only* |
| Hardware | any | **CUDA GPU** (CPU = very slow) | any |
| Network | NCBI / OpenAlex / CrossRef / PMC | HuggingFace (models + BigBIO), NCBI | optional (graph CDN) |
| Extra | **Docker + GROBID** (for `grobid_xml.py`) | ontology DB files (below) | — |

Full details per script are in the `step_*.html` docs and `requirements.html`.

---

## The three stages

### Stage 1 — publications (local)
Seven scripts in the bundle root, run in order by `run_pipeline.py`:
`pubmed_query.py` → `high_impact_xml.py` → `xml_structure.py` → `ncbi_pdf.py` →
`grobid_xml.py` → `named_entity_xml.py` → `pre_ner_xml_structure.py`.

- Input: a PubMed query (read from **STDIN** by `pubmed_query.py`; `run_pipeline.py --query`
  pipes it in). The query used for this project is:
  ```
  "non-small cell lung cancer"[Title/Abstract] NOT "small cell lung carcinoma"[Title/Abstract]
  ```
  > **Watch the hyphen.** Do *not* write the exclusion as `NOT "small cell lung cancer"`:
  > PubMed splits `non-small` into `non` + `small`, so the phrase `"small cell lung cancer"`
  > is a token-substring of every `"non-small cell lung cancer"` record and the `NOT`
  > excludes all of them → **0 hits**. Excluding `"small cell lung carcinoma"` (carcinoma,
  > not cancer) avoids the trap and returns the intended set. A query that matches 0 records
  > now aborts step 1 fast with an explanation instead of hanging.
- **Impact percentile prompt:** when `step_1_orchestrator.py` runs step 2 it prompts
  `Publication impact percentile (decimal between 0 and 1):` on its own line right after
  the query, and passes the entered value to `high_impact_xml.py` via the `PERCENTILE`
  env var — the only channel that script reads it from. It selects articles whose journal
  impact factor is at or above that percentile (e.g. `0.90` → top 10%). A blank line falls
  back to any inherited `PERCENTILE` env var, or the built-in `0.90` default.
- Reaches NCBI E-utilities, OpenAlex, CrossRef, PMC. Set `NCBI_API_KEY` to lift the
  3 req/s rate limit. Optional env vars: `TIME_BUDGET`, `IF_THRESHOLD`, `PERCENTILE`,
  `RETRY_FAILED`, `GROBID_*`, … (see `step_1_publications.html`).
- **`grobid_xml.py` needs Docker + a GROBID server on `:8070`** (it can auto-launch Docker
  Desktop + the container). It is skippable when every article already has JATS full text.
- Output: the NER corpus `gpu_bundle/experimental_ner/PMC*.xml` — the input to stage 2.
- **Optional — `subtract.py`** (not one of the seven, not run by the orchestrator): reads two
  directory paths from **STDIN** and moves entries of `directory_1` whose names also appear in
  `directory_2` into `gpu_bundle/removed/` (relocated, not deleted; name collisions get a
  `_1`/`_2` suffix), writing `summaries/subtract_optional.html`. Handy for de-duplicating this
  project's `gpu_bundle/experimental_ner/` against another corpus. See `step_1_publications.html`.

### Stage 2 — triples / GPU bundle (Kaggle or local GPU)
`gpu_bundle/gpu.py` orchestrates a 16-step chain (RE-model training → BioBERT NER → GENETIC/
DISEASE/CHEMICAL normalization → rule triples → learned relation extraction → **zip**) in one
working directory.

**On Kaggle (recommended):**
1. Upload this bundle as a Kaggle Dataset (the `gpu_bundle/` scripts + `experimental_ner/`
   from stage 1 + the ontology DBs below).
2. *Settings → Accelerator → GPU* and enable *Internet*.
3. In a cell: `!pip install -q 'datasets<4'` then `!python gpu.py` (from the dataset dir).
4. Download the produced **`kaggle_working.zip`**.

**Locally:** `python run_pipeline.py --steps 2` (or `cd gpu_bundle && python gpu.py`). Needs
the GPU deps and the DB files present; preview with `python gpu_bundle/gpu.py --list`.

Output: `TRIPLES/` (incl. the scored + normalized triples) and `kaggle_working.zip`.

### Stage 3 — graph (local)
`high_confidence.py` filters the scored triples to the high-confidence gene–gene set in a
disease/chemical context and renders the interactive graph.

> **NB!** If stage 2 ran on Kaggle, **download `kaggle_working.zip` and unzip it here first** —
> stage 3 reads its inputs from that unzipped run directory.

```bash
# after unzipping kaggle_working.zip into ./kaggle_working
python run_pipeline.py --steps 3 --data-root kaggle_working
# or directly:
python high_confidence.py --data-root kaggle_working
```

Output: `<data-root>/summaries/high_confidence.html`, `<data-root>/TRIPLES/high_confidence_G_D_C.json`,
and a copy of the graph in the bundle root named after the current directory plus today's date
(e.g. `lung_large_2026_07_19.html`). The gene-only `high_confidence_g.py` writes the same set with a
`_G` suffix (e.g. `lung_large_2026_07_19_G.html`).

---

## Data you must provide

These are **git-ignored for size** (30 MB – 500 MB each) and are not in the bundle — place them
under `gpu_bundle/databases/` before running stage 2 (see
`gpu_bundle/databases/PLACE_DATABASES_HERE.md`):

- `hgnc_complete_set_2026-05-01.json` — HGNC gene symbols
- `mondo-clingen.json` — MONDO disease ontology
- `chebi.json` — ChEBI chemical ontology

`gpu_bundle/databases/pmc_years.json` is produced by stage 2 (or supply it to run the year
filter offline). The `experimental_ner/` corpus is produced by **stage 1** (or drop in your
own). Trained checkpoints (`ppi-biobert-re/`) and run outputs are generated, not committed.

---

## Repository layout

```
cancer_knowledge_graph/
├── run_pipeline.py            # end-to-end orchestrator (this bundle's entry point)
├── requirements.txt           # local deps (stages 1 & 3): requests
├── pubmed_query.py … pre_ner_xml_structure.py   # stage 1: the 7 publications scripts
├── subtract.py                # stage 1: optional dir-subtract utility (-> gpu_bundle/removed)
├── high_confidence.py         # stage 3: the graph
├── gpu_bundle/                # stage 2: the GPU pipeline
│   ├── gpu.py                 #   orchestrator (16 steps)
│   ├── requirements.txt       #   GPU deps: torch/transformers/datasets<4/numpy/lxml
│   ├── *.py                   #   the step scripts
│   └── databases/             #   ontology JSONs (provide these; git-ignored)
│       └── PLACE_DATABASES_HERE.md
├── step_1_publications.html   # rendered walk-throughs …
├── step_2_triples.html
├── step_3_graph.html
└── requirements.html
```

Generated corpora, model checkpoints, run trees (`kaggle_working/`), and `*.zip` bundles are
excluded via `.gitignore`.

---

## Notes

- All scripts resolve paths **relative to their own location**, so they run from any working
  directory and a copied tree runs in isolation.
- Stages 1 and 2 reach external services and (stage 1) need Docker; they are **not** meant to
  run unattended without those prerequisites.

## License

Released under the [MIT License](LICENSE) — © 2026 the cancer_knowledge_graph authors.
