# cancer_knowledge_graph

A biomedical relation-extraction pipeline: from a single **PubMed query** to an interactive
**gene–gene relationship graph** for a disease/chemical context, in three stages.

| Stage | What | Where it runs | Entry point |
|------:|------|---------------|-------------|
| **1** | Publications → full-text NER corpus | local (network + Docker/GROBID) | `step_1_orchestrator.py` (10 stages over 9 root scripts) |
| **2** | NER corpus → normalized, model-scored relation **triples** | GPU (Kaggle or local) | `gpu_bundle/gpu.py` (20-step chain) |
| **3** | Triples → high-confidence gene / disease / chemical **graph** | local | `high_confidence_g.py` |

Each stage hands off to the next **by files**. Rendered walk-throughs of every stage ship with
this bundle: [`step_1_publications.html`](step_1_publications.html),
[`step_2_triples.html`](step_2_triples.html), [`step_3_graph.html`](step_3_graph.html), and a
consolidated [`requirements.html`](requirements.html).
[`Step_2_updated_triples.html`](Step_2_updated_triples.html) covers what stage 2 gained when the
BioRED relation model was added alongside the PPI one — the two extra steps, typed/signed edges,
and the model comparison.

---

## Quick start (all three stages, locally)

```bash
# 1. install the light local deps (stage 1 + stage 3)
python -m pip install -r requirements.txt

# 2. (stage 2 only) GPU deps — normally you run stage 2 on Kaggle instead; see below
python -m pip install -r gpu_bundle/requirements.txt

# 3. stage 1 — query to NER corpus (prompts for the impact percentile)
python step_1_orchestrator.py "pancreatic cancer"

# 4. stage 2 — NER corpus to scored triples (needs a CUDA GPU; usually Kaggle)
python gpu_bundle/gpu.py

# 5. stage 3 — triples to the graph
python high_confidence_g.py --data-root kaggle_working
```

**There is no single all-stage runner** — each stage has its own entry point and they hand off
by files, because stage 2 almost always runs on different hardware than 1 and 3.
`step_1_orchestrator.py` chains stage 1's ten steps and aborts on the first non-zero exit;
run a sub-range with `--start` / `--stop` / `--only`. Its last step (`clean_up.py`) **deletes** the
intermediate XML/PDF directories once the corpus is built — `--stop 7` keeps them.

> **Realistically, stage 2 runs on Kaggle**, not your laptop — it needs a CUDA GPU and the
> large ontology databases. The typical flow is: **stage 1 locally → stage 2 on Kaggle →
> download `kaggle_working.zip`, unzip it here → stage 3 locally.** See below.

---

## Requirements at a glance

|          | Stage 1                                                                               | Stage 2                                                | Stage 3              |
|----------|---------------------------------------------------------------------------------------|--------------------------------------------------------|----------------------|
| Python   | 3.9+                                                                                  | 3.9+                                                   | 3.9+                 |
| Packages | `requests`                                                                            | `torch`, `transformers`, `datasets<4`, `numpy`, `lxml` | *stdlib only*        |
| Hardware | any                                                                                   | **CUDA GPU** (CPU = very slow)                         | any                  |
| Network  | NCBI / OpenAlex / CrossRef / PMC                                                      | HuggingFace (models + BigBIO), NCIt                    | optional (graph CDN) |
| Extra    | **Docker + GROBID** (for `grobid_xml.py`); *optional* local XML archive (steps 1b/6b) | ontology DB files (below)                              | —                    |

Full details per script are in the `step_*.html` docs and `requirements.html`.

---

## The three stages

### Stage 1 — publications (local)
Ten scripts in the bundle root, run as eleven steps by `step_1_orchestrator.py`:
`pubmed_query.py` **(1)** → `from_archive.py` **(1b)** → `high_impact_xml.py` **(2)** →
`xml_structure.py` **(3)** → `ncbi_pdf.py` **(4)** → `grobid_xml.py` **(5)** →
`named_entity_xml.py` **(6)** → `from_archive.py` **(6b)** → `pub_year_xml.py` **(6c)** →
`pre_ner_xml_structure.py` **(7)** → `clean_up.py` **(8)**.

Steps are addressed by **label**, not position, so the core stages keep the numbers they
always had (`--start 4` still resumes at `ncbi_pdf.py`) while `from_archive.py` slots in as 1b/6b
and `pub_year_xml.py` as 6c.

- Input: a PubMed query (read from **STDIN** by `pubmed_query.py`; the orchestrator's first
  bare argument pipes it in). The EXAMPLE-QUERY used for this project is:
  ```
  "non-small cell lung cancer"[Title/Abstract] NOT "small cell lung carcinoma"[Title/Abstract]
  ```
  > **Watch the hyphen.** Do *not* write the exclusion as `NOT "small cell lung cancer"`:
  > PubMed splits `non-small` into `non` + `small`, so the phrase `"small cell lung cancer"`
  > is a token-substring of every `"non-small cell lung cancer"` record and the `NOT`
  > excludes all of them → **0 hits**. Excluding `"small cell lung carcinoma"` (carcinoma,
  > not cancer) avoids the trap and returns the intended set. A query that matches 0 records
  > now aborts step 1 fast with an explanation instead of hanging.
- **Reuse a local archive — `from_archive.py` (steps 1b + 6b).** Matches the query result
  against a local XML archive (`ARCHIVE_DIR`, default `../xmls`) by **PMC id** and serves every
  publication already held there instead of re-fetching it, turning the pipeline's slowest part
  (one eFetch per article at ~3 req/s) into a file copy. It picks one representative per paper —
  **full-text JATS → GROBID TEI → abstract-only JATS** — **copies** it (the archive is opened
  read-only and never modified) into `archive_xmls/` and `gpu_bundle/experimental_ner/`, and
  writes `pmids/from_archive_pmcids.txt` so step 2 skips those ids. On the reference
  lung large-cell run, **309 of 331** PMC-bearing hits (93%) came from the archive and only 22
  were downloaded.
  > **Why it runs twice.** Step 6 clean-rebuilds `gpu_bundle/experimental_ner/` and would discard
  > the archive contribution, so step **6b** re-runs the same script to restore it from
  > `archive_xmls/`. It is idempotent — run it a third time and nothing changes.

  Steps 1b/6b are **skipped automatically, with a printed notice, when the archive directory does
  not exist**, so a checkout with no local corpus still runs the plain pipeline (the eight core
  stages plus 6c).
  Naming one explicitly (`--archive` / `ARCHIVE_DIR`) that is missing is a hard error instead —
  that is a typo, not an absence. `--no-archive` skips them outright.
- **Publication-year stamp — `pub_year_xml.py` (step 6c).** Runs once the corpus is final (after
  6/6b, which both rewrite it) and writes one processing instruction, `<?pub-year YYYY?>`, into the
  prolog of every `gpu_bundle/experimental_ner/*.xml` (right after the XML declaration — valid in
  both PMC JATS and GROBID TEI, ignored by parsers). The year comes from the `year` column of
  `pmids/pmid_pmc_ids.tsv` (PubMed; it wins), falling back to the date inside the XML (JATS
  `<pub-date><year>`, TEI `<publicationStmt><date when>`). Offline and idempotent (an existing
  stamp is replaced; files already correct are not rewritten); `--check` reports without writing.
  It **exits non-zero, listing the files, if any document stays undated**, aborting the pipeline —
  so stage 2 gets a fully dated corpus and never needs the PubMed table or NCBI. On the current
  corpus: 3,769 XML, all dated from the table (411 state a different year themselves; the table is
  kept), range 1949–2026.
- **Impact percentile prompt:** when `step_1_orchestrator.py` runs step 2 it prompts
  `Publication impact percentile (decimal: 0 <= impact <= 1):` on its own line right after
  the query, and hands the entered value to `high_impact_xml.py` on its **STDIN** and as the
  `PERCENTILE` env var (env wins there, so the validated value takes effect either way). It
  selects articles whose journal impact factor is at or above that percentile (e.g. `0.90` →
  top 10%; `0.01` → effectively everything). A blank line falls back to the built-in `0.90`
  default. Run standalone, the script reads the percentile from STDIN:
  For example: `echo 0.01 | python high_impact_xml.py`. Decimal = 0 includes all PubMed articles
  regardless of their impact.
- Reaches NCBI E-utilities, OpenAlex, CrossRef, PMC. Set `NCBI_API_KEY` to lift the
  3 req/s rate limit, and `OPENALEX_API_KEY` (free key) so the step-1 impact-factor
  lookup doesn't run out of OpenAlex's shared per-IP daily budget. Optional env vars: `TIME_BUDGET`, `IF_THRESHOLD`, `PERCENTILE`,
  `RETRY_FAILED`, `ARCHIVE_DIR`, `USE_ARCHIVE_SKIP`, `GROBID_*`, …
  (see `step_1_publications.html`).
- **`grobid_xml.py` needs Docker + a GROBID server on `:8070`** (it can auto-launch Docker
  Desktop + the container). It is skippable when every article already has JATS full text.
- Output: the NER corpus `gpu_bundle/experimental_ner/PMC*.xml` — the **union** of the downloaded
  and archive-served papers, de-duplicated by PMC id. This is the input to stage 2.
- **Clean-up — `clean_up.py` (step 8, the last one).** Once the corpus exists and step 7 has
  reported on it, the working directories are just bulk on disk (hundreds of MB of JATS XML, PDFs
  and TEI that no later stage reads), so step 8 **deletes them and their contents**:
  `archive_xmls/`, `grobid_xmls/`, `high_impact_xmls/`, `named_entity_xmls/`, `ncbi_pdfs_grobid/`.
  It **never touches** `gpu_bundle/experimental_ner/` (the corpus), `pmids/` (the query result and
  its resumable caches), `summaries/` (the HTML reports) or any script, and everything it removes
  is reproducible from `pmids/pmid_pmc_ids.tsv` by re-running steps 2–6. The removal list is a
  fixed literal of five names — nothing comes from arguments or the environment — and each must be
  a real directory directly inside the bundle root, so a symlink or an out-of-tree path is refused
  rather than followed. Absent directories are fine, so it is idempotent. `DRY_RUN=1` (or
  `--dry-run`) reports files and bytes per directory without deleting; `KEEP=archive_xmls` (comma
  separated) spares one; `--stop 7` skips the step entirely. Writes `summaries/clean_up.html`.
  > **Re-running after a clean-up** starts from step 2 — the inputs of steps 3–6b are gone. The
  > exception is `--start 6b --stop 6c`: `from_archive.py` re-copies the hit set out of the real
  > archive (`ARCHIVE_DIR`), so it still works, it just pays the file copy again — and 6c
  > re-stamps the years on the copies it replaced.
- **Optional — `subtract.py`** (not one of the eleven, not run by the orchestrator): reads two
  directory paths from **STDIN** and moves entries of `directory_1` whose names also appear in
  `directory_2` into `gpu_bundle/removed/` (relocated, not deleted; name collisions get a
  `_1`/`_2` suffix), writing `summaries/subtract_optional.html`. Handy for de-duplicating this
  project's `gpu_bundle/experimental_ner/` against another corpus. See `step_1_publications.html`.

### Stage 2 — triples / GPU bundle (Kaggle or local GPU)
`gpu_bundle/gpu.py` orchestrates a 20-step chain (RE-model training ×2 → **drug lexicon** →
BioBERT NER → GENETIC/DISEASE/CHEMICAL normalization → rule triples → learned relation extraction
→ model comparison → scores by publication year → **zip**) in one working directory. See
[`Step_2_updated_triples.html`](Step_2_updated_triples.html) for the two BioRED steps and
[`step_2_triples.html`](step_2_triples.html) for the original 16.

**On Kaggle (recommended):**
1. Upload this bundle as a Kaggle Dataset (the `gpu_bundle/` scripts + `experimental_ner/`
   from stage 1 + the ontology DBs below + — if you already have them — the four training dirs
   `ppi-biobert-re/`, `ppi_data/`, `biored-biobert-re/`, `biored_data/`, which make the gate
   skip steps 1–2; see below).
2. *Settings → Accelerator → GPU T4 x2* (both cards are used; see *Multi-GPU* below) and
   enable *Internet*.
3. In a cell: `!pip install -q 'datasets<4' bioc` then `!python gpu.py` (from the dataset dir).
   (`bioc` is what the BigBIO loading script for **BioRED** needs — that corpus is BioC XML;
   without it step 2 fails with `ModuleNotFoundError: No module named 'bioc'`.)
4. Download the produced **`kaggle_working.zip`**.

**Locally:** `cd gpu_bundle && python gpu.py`. Needs the GPU deps and the DB files present;
preview with `python gpu_bundle/gpu.py --list`.

Output: `TRIPLES/` (incl. the scored + normalized triples) and `kaggle_working.zip`.

**Multi-GPU (Kaggle *GPU T4 x2*).** Every GPU step uses all visible GPUs; there is nothing
to configure. With one GPU or none, each step behaves exactly as before.

| step                                | how it uses two GPUs                                                                                              | opt out                 |
|-------------------------------------|-------------------------------------------------------------------------------------------------------------------|-------------------------|
| 1 + 2 `re_pipeline(_biored)`        | the two trainings run **side by side**, PPI on GPU 0 and BioRED on GPU 1; log lines are prefixed with the step name | `NORM_PARALLEL_TRAIN=0` |
| 5 `sentences.py`                    | one copy of the scorer + NER models per GPU; 16-file chunks go to whichever GPU is free                            | `BIOBERT_GPUS=1`        |
| 17 `relation_extraction.py`         | each checkpoint is loaded once per GPU; batches go to whichever GPU is free                                         | `RE_GPUS=1`             |

Results are consumed in input order, so the output files are identical to a single-GPU run. Each
training run still uses one GPU: running both at once is faster than splitting one small BERT
across two cards (`--re-args "--gpus 0,1"` asks for that split, and also turns off the
side-by-side run). The two runs write separate reports, `summaries/calibration_ppi.html` and
`summaries/calibration_biored.html`. `gpu.py --list` shows the GPU count on its `accelerator`
line.

**Train the two checkpoints once, reuse them for every corpus.** Steps 1 and 2 never read
`experimental_ner/` or your PubMed query — they fine-tune BioBERT on fixed public BigBIO corpora
(**bioinfer** → `ppi-biobert-re/`, **biored** → `biored-biobert-re/`). The checkpoint is a
function of (dataset, seed, hyperparams) only, so a run on "<pubmed_query_1>" produces the same
model as one on "<pubmed_query_2>"; retraining per corpus is wasted GPU time. The same
holds for the normalization libraries (steps 6–14) — HGNC / ChEBI / MONDO are ontologies, not
corpus-derived.

**`gpu.py` therefore gates the two training steps** — they run **only** when the previous training
output is not already in the bundle. The gate tests the full content of four directories under
`gpu_bundle/`:

| directory                               | required content                                                                                                   |
|-----------------------------------------|--------------------------------------------------------------------------------------------------------------------|
| `ppi-biobert-re/`, `biored-biobert-re/` | `config.json`, `tokenizer.json`, `tokenizer_config.json`, `calibration.json`, weights (`*.safetensors` or `*.bin`) |
| `ppi_data/`, `biored_data/`             | `train.tsv`, `dev.tsv`, `test.tsv`                                                                                 |

All four complete → steps 1 **and** 2 are dropped from the plan and the checkpoints on disk are
used as-is. Anything absent or empty → both run, and the header names exactly which entries were
missing. It is **all four or none** by design: step 17 routes between the two checkpoints and step
18 compares them, so a run must never mix a reused PPI model with a freshly trained BioRED one. An
empty file counts as missing, so a half-copied bundle retrains instead of loading a truncated
checkpoint.

```bash
python gpu_bundle/gpu.py --list        # shows the verdict and its reason, runs nothing
python gpu_bundle/gpu.py --retrain     # train anyway (also: --steps re_pipeline / --re-args)
```

#### What the four directories hold

**The two checkpoints.** Both are BioBERT (`dmis-lab/biobert-v1.1`) fine-tuned as a
`BertForSequenceClassification` — same architecture, same 28,996-token WordPiece vocabulary, same
512-position limit. They differ only in the label head and the corpus behind it.

| file                                       | `ppi-biobert-re/`                                             | `biored-biobert-re/`                                                                               |
|--------------------------------------------|---------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| `model.safetensors`                        | 413 MB — fine-tuned encoder + classifier head                 | same                                                                                               |
| `config.json`                              | `id2label` = `interacts`, `false` (binary)                    | `associated`, `binds`, `downregulator/inhibitor`, `upregulator/activator`, `false` (5-way, signed) |
| `tokenizer.json` + `tokenizer_config.json` | BioBERT's WordPiece vocab, unchanged                          | same                                                                                               |
| `calibration.json`                         | `{"method": "platt", "a": 0.821, "b": -0.340, "clip": 0.001}` | `{"method": "platt", "a": 0.388, "b": -1.220, "clip": 0.001}`                                      |
| `training_args.bin`                        | HF `TrainingArguments` — provenance only                      | same                                                                                               |
| `checkpoint-N/`                            | `checkpoint-878` (epoch 2)                                    | `checkpoint-6357` (epoch 3)                                                                        |

The **first five** are what the gate requires and what `relation_extraction.py` loads:
weights + label map + tokenizer, and `calibration.json` mapping the raw softmax to the calibrated
`p_rel` (both Platt, so the two models' scores share one scale — see `step_2_triples.html` §6).

`training_args.bin` and `checkpoint-N/` are **neither gated nor read at inference**. `train_re.py`
runs with `load_best_model_at_end=True` and `save_total_limit=1`, then `trainer.save_model(out)`,
so the best epoch's weights are already at the directory root and the one surviving trainer
checkpoint is an audit/resume artifact. Deleting it takes each directory from **1.7 GB to ~414 MB**
— worth doing before a Kaggle upload. Its `trainer_state.json` is the only thing worth keeping: it
holds the per-epoch dev metrics that `EPOCH_DEFAULTS` in `train_re.py` was measured from (this
bundle's run: PPI best at epoch 2, dev F1 0.843; BioRED at epoch 3, dev F1 0.624).

**The two data dirs** are `bigbio_to_re.py`'s conversion output — the entity-blinded TSVs
`train_re.py` trains on. Three columns (`index`, `sentence`, `label`), one **entity pair** per row,
that pair's two mentions replaced by type markers:

```
0	The herpes simplex virus type 1 (HSV-1) @GENE$, an essential component of the viral DNA replication machinery, is a trimeric complex of the virus-coded UL5, UL8, and @GENE$ proteins.	1
0	Chloroacetaldehyde (CAA) is a metabolite of the alkylating agent @CHEMICAL$ (IFO) and putatively responsible for renal damage following anti-@DISEASE$ therapy with IFO.	Negative_Correlation
```

|                          | `ppi_data/` (BioInfer)                                     | `biored_data/` (BioRED)                                                                                          |
|--------------------------|------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| rows: train / dev / test | 7,018 / 779 / 1,604                                        | 33,889 / 9,891 / 9,202                                                                                           |
| markers                  | `@GENE$` only                                              | `@GENE$`, `@DISEASE$`, `@CHEMICAL$`, `@VARIANT$`                                                                 |
| train labels             | `1` interacts 1,897 · `0` no relation 5,121                | `Association` 10,908 · `Positive_Correlation` 5,186 · `Negative_Correlation` 4,089 · `Bind` 458 · `false` 13,248 |
| dev split                | carved from train (`--val-frac 0.1` — BioInfer ships none) | BioRED's own 400/100/100 abstracts (`--val-frac 0`)                                                              |
| on disk                  | 2.3 MB                                                     | 11 MB                                                                                                            |

BioRED's four rare chemical–chemical types (`Cotreatment`, `Comparison`, `Drug_Interaction`,
`Conversion`) fold into `Association` by default — `--biored-all-types` keeps all eight. Its row
count dwarfs "600 abstracts" because its relations are annotated at the **document** level: one
annotated pair emits a row for every sentence where both endpoints co-occur (`BioRED_task_mapping.md`
§5), which is also why its positives are noisier than BioInfer's sentence-scoped ones.

Neither data dir is needed to *score* triples. They are gated because they are what makes a reused
bundle re-calibratable in place: `run_re_pipeline.py --skip-convert --skip-train` re-runs the
checkpoint over `dev.tsv` (fit) and `test.tsv` (metrics + ECE) and rewrites `calibration.json` and
`summaries/calibration_<task>.html` — no BigBIO download, no retraining. `train.tsv` is read there only
for the row count in the summary.

All four dirs are git-ignored — they are inputs as often as they are outputs, so keep them in
`gpu_bundle/` between runs (and upload them with the Kaggle dataset) rather than in the
`kaggle_working/` run tree, which is regenerated each time.

> **What *is* corpus-sensitive is transfer quality, not checkpoint validity.** Both training
> corpora are abstracts, and `calibration.py` fits the calibrator on their own dev/test split — so
> `p_rel` is calibrated to *those* sentence distributions, not to your GROBID full text. The knobs
> worth re-tuning on a new corpus are therefore `RE_MIN_SCORE` and the stage-3 cutoffs, not the
> models. A corpus that is clinical or epidemiological rather than molecular is further out of
> domain, and BioRED's typed/signed labels degrade faster there than the PPI model's binary
> `interacts`. Retrain only to change the label space or add a source dataset — `bigbio_to_re.py`
> also supports `--task chemprot / gad / ddi`.

### Stage 3 — graph (local)
`high_confidence_g.py` filters the scored triples to the high-confidence set and renders the
interactive graph over gene, disease and chemical nodes. Edges are typed by the **relation**
the model predicted (activates / inhibits / binds / interacts / associated, dashed when
negated), and `--merge` collapses the two verdicts an additive stage-2 run writes per pair
(PPI + BioRED) into one — default `union`: every pair either model kept (BioRED-only edges
included), with the typed label preferred wherever BioRED fired and the higher of the two
scores. `--merge gate` is the stricter variant that keeps only pairs the binary PPI model
also claimed.

> **NB!** If stage 2 ran on Kaggle, **download `kaggle_working.zip` and unzip it here first** —
> stage 3 reads its inputs from that unzipped run directory.

```bash
# after unzipping kaggle_working.zip into ./kaggle_working
python high_confidence_g.py --data-root kaggle_working
python high_confidence_g.py --data-root kaggle_working --merge gate    # or typed / none
```

The drug-target colouring rests on a cross-link (corpus chemical → HGNC gene) built from
ChEBI roles and DGIdb. Three ways a gene could pick up a colour it had not earned are fixed
in the pipeline and guarded by a regression check — run it after any DGIdb or ChEBI update:

```bash
python audit_drug_targets.py --data-root kaggle_working     # exit 0 = clean
```

The graph spans gene, DISEASE and CHEMICAL nodes (diseases and chemicals identified by MONDO /
ChEBI label — or NCIt label, for the biologics ChEBI has no term for — shaped ◆ and ■), so the
gene–disease and chemical–gene edges BioRED contributes are drawn instead of discarded — on the
reference run, 129 nodes / 151 edges against 62 / 55 for a gene-only view. That gene-only view
and its `_G` outputs have been retired; the `_M` outputs below are the only ones written.

Output: `<data-root>/summaries/high_confidence_M.html`, `<data-root>/TRIPLES/high_confidence_M.json`,
and a copy of the graph in the bundle root named after the current directory plus today's date
(e.g. `lung_large_2026_07_19_M.html`).

**Reading the graph.** The page is self-contained — payload, library and all — and every edge
traces back to the sentences behind it, with clickable PMIDs. Beyond filtering, it ranks the
**nodes** carrying the current view (*Significance in view*: publications, partners or sentences,
with a z-score, bootstrap confidence intervals on count and rank, and Fisher enrichment against
the corpus with BH q-values), and a *Significance* slider hides what misses a q cutoff.
[`Graph_description.md`](Graph_description.md) is the full guide: what each control does, the
order the filters run in, and the limits worth carrying.

**One corpus for every `lung_*` project.** The enrichment denominator is the union of
publications across all sibling `lung_*` directories, so an odds ratio computed in one means the
same thing as one computed in another. Each project caches its slice in
`<data-root>/databases/corpus_contrib.json` and merges whatever the siblings have written; a
project that has never been run is reported as absent rather than silently assumed. Run stage 3
once per project to populate it.

> **`high_confidence.py` has been removed.** The older "G_D_C" script (gene–gene *in a
> disease/chemical context*, unsuffixed outputs) was deleted in `8f420f1`: it ignored
> `predicate.text`, so BioRED's signed labels collapsed into one edge category, and it had no
> multi-model merge, so its counts double-counted an additive run. Sections 1–7 of
> [`step_3_graph.html`](step_3_graph.html) still describe it, as history.

---

## Data you must provide

These are **git-ignored for size** (12 MB – 500 MB each) and are not in the bundle — place them
under `gpu_bundle/databases/` before running stage 2 (see
`gpu_bundle/databases/PLACE_DATABASES_HERE.md`):

- `hgnc_complete_set_2026-05-01.json` — HGNC gene symbols. **Required** by `roman.py`,
  `greek.py`, `controls.py`, `nonchemical.py`, `target_pharm.py`
- `mondo-clingen.json` — MONDO disease ontology. **Required** by `disease.py`
- `chebi.json` — ChEBI chemical ontology. **Required** by `chemical.py`, `target_pharm.py`
- `interactions.tsv` — DGIdb open drug–gene interactions. *Optional*, read by `chemical.py` for
  the drug–gene target layer. Absent, it degrades silently to an empty drug set — no error, just
  fewer CHEMICAL surfaces
- `ncit_drugs.json` — NCI Thesaurus drug names. *Optional but recommended*, and **generated**:
  `python gpu_bundle/drug_lexicon.py --build` (or stage-2 step 3) downloads NCIt and distils it.
  Two steps read it — `sentences.py` tags the drugs the BioBERT chemical model misses, and
  `chemical.py` normalizes the ones ChEBI cannot represent. Absent, both fall back to the
  INN-stem regex alone: antibody *names* are still recognised, but brand names, code names and
  synonyms are not, and none of them can be normalized

`gpu.py --list` preflights the three required ones and names any that are missing; the DGIdb
table is not preflighted, precisely because it is optional.
`gpu_bundle/databases/pmc_years.json` and `ncit_drugs.json` are produced by stage 2 (supply
`ncit_drugs.json` to run the drug lexicon offline; `pmc_years.json` is harvested offline from the
corpus's `<?pub-year?>` stamps, so supply it only for a run that reuses `sentences/` without
`experimental_ner/`). The `experimental_ner/` corpus is produced by **stage 1** (or drop in your
own). The trained checkpoints and their converted TSVs (`gpu_bundle/{ppi-biobert-re,ppi_data,
biored-biobert-re,biored_data}/`) are generated by stage-2 steps 1–2 on the first run and are
git-ignored — **keep them in `gpu_bundle/` afterwards**: their presence is what makes the training
gate skip those two steps on every later corpus. Run outputs are likewise generated, not committed.

**Optional — a local XML archive.** If you have a directory of previously fetched
`PMC*.xml` / `PMC*.grobid.tei.xml` files, point `ARCHIVE_DIR` (or `--archive`) at it and stage 1's
steps 1b/6b will serve matching publications from disk instead of re-downloading them. Nothing in
it is modified. Without one, stage 1 simply downloads everything.

---

## Repository layout

```
<parent_directory>/
├── step_1_orchestrator.py                       # stage 1 entry point (chains its 11 steps)
├── requirements.txt                             # local deps (stages 1 & 3): requests
├── pubmed_query.py … pre_ner_xml_structure.py   # stage 1: the 7 core publications scripts
├── from_archive.py                              # stage 1: steps 1b/6b — serve the query from a local XML archive
├── pub_year_xml.py                              # stage 1: step 6c — stamp <?pub-year YYYY?> into every corpus XML
├── clean_up.py                                  # stage 1: step 8 — deletes the intermediate XML/PDF dirs (last)
├── subtract.py                                  # stage 1: optional dir-subtract utility (-> gpu_bundle/removed)
├── high_confidence_g.py                         # stage 3: the graph (typed edges + --merge)
├── audit_drug_targets.py                        # stage 3: regression test for false drug targets
├── gpu_bundle/                                  # stage 2: the GPU pipeline
│   ├── gpu.py                                   #   orchestrator (20 steps)
│   ├── drug_lexicon.py                          #   NCIt drug lexicon: builder + matcher (step 3)
│   ├── requirements.txt                         #   GPU deps: torch/transformers/datasets<4/numpy/lxml
│   ├── *.py                                     #   the step scripts
│   ├── ppi-biobert-re/                          #   PPI checkpoint (BioInfer) — git-ignored, KEEP between runs
│   │   ├── model.safetensors                    # ‡   413 MB: fine-tuned BioBERT encoder + classifier head
│   │   ├── config.json                          # ‡   architecture + id2label: interacts / false
│   │   ├── tokenizer.json                       # ‡   BioBERT WordPiece vocab
│   │   ├── tokenizer_config.json                # ‡
│   │   ├── calibration.json                     # ‡   {method: platt, a, b, clip} — fit by calibration.py
│   │   ├── training_args.bin                    #     HF TrainingArguments — provenance only
│   │   └── checkpoint-878/                      #     best epoch (2); 1.3 GB, safe to delete before uploading
│   │       ├── config.json, model.safetensors, tokenizer{,_config}.json, training_args.bin
│   │       ├── optimizer.pt, scheduler.pt, rng_state.pth   #   resume state
│   │       └── trainer_state.json                          #   per-epoch dev metrics
│   ├── ppi_data/                                #   bigbio_to_re.py output: index / sentence / label TSVs
│   │   ├── train.tsv                            # ‡   7,018 rows  (@GENE$-blinded pairs, labels 1 / 0)
│   │   ├── dev.tsv                              # ‡     779 rows  — the calibrator is fit here
│   │   └── test.tsv                             # ‡   1,604 rows  — test metrics + ECE
│   ├── biored-biobert-re/                       #   BioRED checkpoint — same seven entries, 5-way signed head
│   │   ├── model.safetensors                    # ‡   413 MB
│   │   ├── config.json                          # ‡   id2label: associated / binds / downregulator-inhibitor /
│   │   │                                        #       upregulator-activator / false
│   │   ├── tokenizer.json                       # ‡
│   │   ├── tokenizer_config.json                # ‡
│   │   ├── calibration.json                     # ‡   its own Platt a/b (not interchangeable with the PPI one)
│   │   ├── training_args.bin                    #
│   │   └── checkpoint-6357/                     #     best epoch (3); same file set as checkpoint-878/
│   ├── biored_data/                             #   bigbio_to_re.py output, document-level → sentence-level
│   │   ├── train.tsv                            # ‡   33,889 rows (@GENE$/@DISEASE$/@CHEMICAL$/@VARIANT$)
│   │   ├── dev.tsv                              # ‡    9,891 rows — BioRED's own dev split, not carved
│   │   └── test.tsv                             # ‡    9,202 rows
│   └── databases/                               #   ontology + reference data (provide these; git-ignored)
│       ├── PLACE_DATABASES_HERE.md              # the only tracked file here: what to put in this dir
│       ├── hgnc_complete_set_2026-05-01.json    #  34 MB  HGNC gene symbols → roman.py,
│       │                                        #     greek.py, controls.py, nonchemical.py, target_pharm.py
│       ├── mondo-clingen.json                   #  83 MB  MONDO disease ontology → disease.py
│       ├── chebi.json                           # 507 MB  ChEBI chemical ontology → chemical.py, target_pharm.py
│       ├── interactions.tsv                     #  12 MB  DGIdb drug–gene interactions (open) → chemical.py
│                                                #     OPTIONAL: absent = empty drug set, no error
│       ├── ncit_drugs.json                      # 5.8 MB  NCIt drug names → sentences.py (NER lexicon pass) and
│                                                #     chemical.py (non-ChEBI fallback); generated by step 3
│       └── pmc_years.json                       #   generated by pub_years.py (step 4) from the corpus's
│                                                #     <?pub-year?> stamps (offline); read by relationships.py,
│                                                #     triples_per_year.py and stage 3
├── Graph_description.md                         # what the graph viewer draws and how to read it
├── monoclonal_antibody_NER.md                   # why biologics were missed, and what now recognises them
├── step_1_publications.html                     # rendered walk-throughs …
├── step_2_triples.html
├── Step_2_updated_triples.html                  #   stage 2 after the BioRED model was added
├── step_3_graph.html
└── requirements.html
```

**‡** = required by the stage-2 **training gate**. All five per checkpoint and all three TSVs per
data dir must be present and non-empty in *all four* directories, or `gpu.py` retrains both models
(see Stage 2). `training_args.bin` and `checkpoint-N/` are unmarked: nothing reads them at
inference.

Generated corpora, model checkpoints, run trees (`kaggle_working/`), and `*.zip` bundles are
excluded via `.gitignore` — but the four training directories above are the exception worth
remembering: git-ignored, yet they belong in `gpu_bundle/` permanently, because their presence is
what stops every later run from retraining.

---

## Notes

- All scripts resolve paths **relative to their own location**, so they run from any working
  directory and a copied tree runs in isolation.
- Stages 1 and 2 reach external services and (stage 1) need Docker; they are **not** meant to
  run unattended without those prerequisites.

## License

Released under the [MIT License](LICENSE) — © 2026 the cancer_knowledge_graph authors.
