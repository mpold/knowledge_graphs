# Draft: adding a `biored` task to the RE pipeline

Adds BioRED (Luo et al. 2022, 600 PubMed abstracts, NCBI) as a fifth `--task`
alongside `ppi` / `chemprot` / `gad` / `ddi`, so the relation extractor emits
**typed, signed** edges instead of a single undirected `interacts`.

Status: **untested draft.** The BigBIO type strings below are from BioRED's
published schema and have NOT been verified against the installed loader. Run the
verification step first — it writes no files.

---

## Source

Luo L, Lai P-T, Wei C-H, Arighi CN, Lu Z.
**BioRED: a rich biomedical relation extraction dataset.**
*Briefings in Bioinformatics*, 23(5):bbac282, September 2022.

- Journal: <https://academic.oup.com/bib/article/23/5/bbac282/6645993>
- DOI: `10.1093/bib/bbac282`
- PubMed: [PMID 35849818](https://pubmed.ncbi.nlm.nih.gov/35849818/)
- Preprint: [arXiv:2204.04263](https://arxiv.org/abs/2204.04263)
  ("BioRED: A Rich Biomedical Relation Extraction Dataset")
- Data + code: <https://github.com/ncbi/BioRED>

From the NCBI/NLM group behind PubTator, which is why it fits the tooling this
pipeline already uses.

Two claims in the abstract bear directly on this document:

- *"600 PubMed abstracts with multiple entity types and relation pairs ... **at
  the document level**"* — confirms the annotation granularity that section 5
  treats as the main risk. Read the paper's annotation section before deciding
  how aggressively to filter.
- *"each relation is labeled as describing either a novel finding or previously
  known background knowledge"* — confirms the novelty signal discussed in
  section 7. Whether BigBIO carries it through is still open; the GitHub
  distribution above has it intact regardless.

```bibtex
@article{luo2022biored,
  author  = {Luo, Ling and Lai, Po-Ting and Wei, Chih-Hsuan and
             Arighi, Cecilia N and Lu, Zhiyong},
  title   = {{BioRED}: a rich biomedical relation extraction dataset},
  journal = {Briefings in Bioinformatics},
  volume  = {23},
  number  = {5},
  pages   = {bbac282},
  year    = {2022},
  month   = sep,
  doi     = {10.1093/bib/bbac282},
  pmid    = {35849818},
  url     = {https://academic.oup.com/bib/article/23/5/bbac282/6645993},
}
```

---

## 0. Verify before writing any code

```
pip install -q 'datasets<4'
python bigbio_to_re.py --task ppi --dataset bigbio/biored --print-types
```

`--task ppi` is only a placeholder so the flag parses; `print_types()` dumps the
distinct entity and relation types actually present. Check the printed strings
against the tables below and adjust before converting. The default config name
will resolve to `biored_bigbio_kb`, matching the script's `<name>_bigbio_kb`
convention.

---

## 1. BioRED schema (expected)

**Entity types → marker**

| BioRED entity type            | marker      |
|-------------------------------|-------------|
| `GeneOrGeneProduct`           | `@GENE$`     |
| `DiseaseOrPhenotypicFeature`  | `@DISEASE$`  |
| `ChemicalEntity`              | `@CHEMICAL$` |
| `SequenceVariant`             | `@VARIANT$`  |
| `OrganismTaxon`               | `@SPECIES$`  |
| `CellLine`                    | `@CELLLINE$` |

**Relation types → predicate**

| BioRED relation        | predicate (readable)        | KG meaning              |
|------------------------|-----------------------------|-------------------------|
| `Positive_Correlation` | `upregulator/activator`     | signed, positive        |
| `Negative_Correlation` | `downregulator/inhibitor`   | signed, negative        |
| `Bind`                 | `binds`                     | physical interaction    |
| `Association`          | `associated`                | unsigned, generic       |
| `Cotreatment`          | `cotreatment`               | chemical-chemical       |
| `Comparison`           | `comparison`                | chemical-chemical       |
| `Drug_Interaction`     | `drug-interaction`          | chemical-chemical       |
| `Conversion`           | `conversion`                | chemical-chemical       |
| *(no annotated relation)* | `false` (NEG)            | negative                |

`Positive_Correlation` / `Negative_Correlation` / `Bind` are the three that
actually upgrade the graph; the rest keep coverage without adding sign.

---

## 2. `bigbio_to_re.py` — three edits

### 2a. Module constants (next to `CHEMPROT_EVAL` / `DDI_TYPES`)

```python
# BioRED relation types (Luo et al. 2022). Unknown/other annotated relations fall
# back to the generic unsigned bucket rather than being dropped.
BIORED_TYPES = {"positive_correlation": "Positive_Correlation",
                "negative_correlation": "Negative_Correlation",
                "association": "Association", "bind": "Bind",
                "cotreatment": "Cotreatment", "comparison": "Comparison",
                "drug_interaction": "Drug_Interaction", "conversion": "Conversion"}

# Entity-type pairs BioRED actually annotates. Same-type pairs collapse to a
# 1-element frozenset on both sides, so GENE-GENE is written as one member.
BIORED_PAIRS = {frozenset(("GENE",)), frozenset(("CHEMICAL",)), frozenset(("VARIANT",)),
                frozenset(("GENE", "DISEASE")), frozenset(("CHEMICAL", "GENE")),
                frozenset(("CHEMICAL", "DISEASE")), frozenset(("VARIANT", "DISEASE")),
                frozenset(("CHEMICAL", "VARIANT"))}
```

### 2b. `marker_for()` — add before the `ppi` fallback

Follows the existing substring-matching idiom, so minor naming variants in the
loader do not break it. Order matters: `variant` is tested first.

```python
    if task == "biored":
        if "variant" in t:                       return "VARIANT"
        if "chem" in t or "drug" in t:           return "CHEMICAL"
        if "dis" in t or "phenotyp" in t:        return "DISEASE"
        if "cell" in t:                          return "CELLLINE"
        if "taxon" in t or "organism" in t or "species" in t: return "SPECIES"
        return "GENE"
```

### 2c. `valid_pair()` — add before the trailing `return True`

```python
    if task == "biored":
        return frozenset((ma, mb)) in BIORED_PAIRS
```

`SPECIES` and `CELLLINE` are deliberately absent from `BIORED_PAIRS`: they are
annotated as entities but are not relation endpoints we want as graph edges. They
still get blinded correctly if they ever appear as a target.

### 2d. `label_for()` — add before the trailing `return "1"`

```python
    if task == "biored":
        if rel_type is None:
            return "false"
        return BIORED_TYPES.get(rel_type.strip().lower(), "Association")
```

Negative label is `"false"`, matching the `chemprot` / `ddi` convention rather
than the `ppi` `"0"`. This matters downstream — see section 4.

### 2e. Two one-line updates

- `--task` choices: `choices=["ppi", "chemprot", "gad", "ddi", "biored"]`
- The `TASKS` block in the module docstring:

```
    biored   : @GENE$/@DISEASE$/@CHEMICAL$/@VARIANT$ ; the pair types BioRED
               annotates ; 8 relation types + false
```

---

## 3. `train_re.py` — one registry entry

Add to `TASKS` (the `label_names` value must map the negative raw label to `NEG`,
which is `"false"`):

```python
    # BioRED (Luo et al. 2022): typed, signed relations over gene / disease /
    # chemical / variant entities. Unlike ppi this is multi-class, so the graph
    # gets a direction and a sign instead of a bare "interacts".
    "biored": {
        "label_names": {"Positive_Correlation": "upregulator/activator",
                        "Negative_Correlation": "downregulator/inhibitor",
                        "Bind": "binds", "Association": "associated",
                        "Cotreatment": "cotreatment", "Comparison": "comparison",
                        "Drug_Interaction": "drug-interaction",
                        "Conversion": "conversion", "false": NEG},
        "known_raw": {"Positive_Correlation", "Negative_Correlation", "Bind",
                      "Association", "Cotreatment", "Comparison",
                      "Drug_Interaction", "Conversion", "false", "true"},
        "markers": ["@GENE$", "@DISEASE$", "@CHEMICAL$", "@VARIANT$"],
    },
```

No other change is needed in this file: the label space is built from whatever
`train.tsv` contains (`labels = sorted({l for _, l in train if l != NEG})`, with
`NEG` appended last), so multi-class flows through unchanged — exactly as it
already does for ChemProt and DDI.

---

## 4. Downstream — what needs nothing

- **`relation_extraction.py`** reads the label set off the checkpoint and drops
  no-relation predictions by auto-detecting the negative class. Its default set
  is `{"0", "false", "none", "no_relation", "negative", "not_related"}`, so
  choosing `"false"` in 2d means `RE_NEG_LABELS` never has to be set.
- **`calibration.py`** / `run_re_pipeline.py` step 3 are task-agnostic.
- **`run_re_pipeline.py`** already forwards `--task` to both subprocesses:

```
python run_re_pipeline.py --task biored --dataset bigbio/biored \
    --data biored_data --model biored-biobert-re
```

BioRED ships its own train/dev/test split, so `--val-frac` should be **0** here
(unlike the BioInfer path, which carves 10% because BioInfer has no dev split —
see `kaggle_working/ppi_data/dev_test_train_explained.md`). Confirm with the split names printed by the
converter; `SPLIT_FILE` already maps `validation` / `valid` / `dev` → `dev.tsv`.

---

## 5. The one real risk: BioRED relations are document-level

**BioRED annotates relations at the abstract level, not the sentence level.**
`iter_instances()` forms pairs strictly *within* a sentence, so:

- **Recall loss (harmless):** relations whose two endpoints never share a
  sentence simply produce no row. Nothing wrong is learned.
- **Precision loss (the real problem):** if entities A and B have a
  document-level relation, then *every* sentence where A and B co-occur yields a
  **positive** row — including sentences that merely mention both without
  asserting anything. Those are false positives in the training labels.

This is the standard document-level-to-sentence-level distant-supervision
artefact. It does not make the corpus unusable — plenty of published work trains
this way — but expect noisier positives than BioInfer, whose relations are
sentence-scoped.

Mitigations, cheapest first:

1. **Measure it before trusting it.** Convert, then hand-read ~50 positive rows
   and count how many actually assert the relation. That number is the ceiling on
   label quality.
2. **Require a relational cue** for a document-level positive: reuse the
   `RELCUE` stem list that `relation_extraction.py` already applies in
   `cue_factor`, and demote a positive with no cue in the connecting text to
   negative. This is a new filter in `iter_instances()`, ~5 lines.
3. **Keep the composite score as the safety net.** `relation_extraction.py`
   already multiplies `p_rel` by `cue_factor`, `margin_factor` and
   `section_factor`, so a confidently-wrong `p_rel` still gets damped at scoring
   time.

A cleaner but larger alternative is to restrict positives to pairs whose relation
BioRED marks with sentence-level evidence, if the BigBIO loader exposes it.

---

## 6. Is 600 abstracts enough to train on?

The reasonable objection to this whole plan. Short answer: for the three edge
types we actually want, probably yes — but "600 abstracts" is the wrong unit to
judge it by, and corpus size is not the thing that will bite us.

### 600 abstracts is not the training set size

The unit that matters is labelled instances. BioRED carries roughly 20,000
entity mentions and ~6,500 annotated relations across its 400/100/100 abstract
split (figures from memory — confirm against the paper or `--print-types`
before quoting them). Set that against what we train on today: BioInfer converts
to ~9,600 pairs, 7,028 of which land in `train.tsv`. Same order of magnitude.

The sentence-level conversion in `bigbio_to_re.py` also *inflates* that count
rather than shrinking it, in one respect: BioRED annotates relations between
normalized **concepts**, so one document-level relation emits a row for every
sentence where a mention of each endpoint co-occurs, and every non-related
co-occurring pair becomes a negative. Expect a row count well above 6,500. Those
extra rows are redundant rather than independent — but the model does not starve.

### Fine-tuning is not training from scratch

Same reason the 3-epoch default works (see
`kaggle_working/ppi_data/dev_test_train_explained.md`): only the classification
head is fresh, and the encoder already knows biomedical language. Thousands of
labelled examples is the normal operating range for BERT fine-tuning, not the
edge of it — several GLUE tasks BERT handles fine (RTE ~2.5k, MRPC ~3.7k) are
smaller than what BioRED converts to. Low-thousands corpus size is a legitimate
worry for training a model from nothing; it is not disqualifying for adapting one.

### What will actually hurt: the class tail, not the total

Eight relation types over ~6,500 instances is not 810 per class — the
distribution is steeply skewed. `Association`, `Positive_Correlation` and
`Negative_Correlation` dominate; `Conversion`, `Comparison`, `Drug_Interaction`
and `Cotreatment` are rare, plausibly tens of instances each. Split three ways, a
class with ~30 total instances has a handful in test, which makes its F1 noise
and its learned boundary close to arbitrary.

This lands where we can afford it. Section 1 already identifies
`Positive_Correlation` / `Negative_Correlation` / `Bind` as the three that
upgrade the graph — and those are the well-populated ones. The starved classes
are the chemical-chemical bookkeeping types described there as "keep coverage
without adding sign."

So the mitigation is to stop asking the model to learn them. **Collapse the
8-way label set to 4** — signed-positive / signed-negative / binds / associated —
by folding the rare chemical-chemical types into `Association` in `BIORED_TYPES`
(section 2a) and in the `train_re.py` registry entry (section 3). Per-class
support rises, the metric stops being dominated by classes we never put in the
graph, and nothing of value is lost.

### Two things to settle before committing

**The document-level artefact is the bigger threat.** Section 5 has this right,
and it interacts badly with size: noisy labels hurt a small corpus far more than
a large one, because there is no volume to average the noise out. The
hand-read-50-positives check (5.1) matters more than any argument about counts.

**Switch the calibrator to Platt for this task.** `train_re.py` defaults to
`--calibration isotonic`, which is fine against BioInfer's 780-row dev split.
BioRED's dev split is 100 abstracts. Isotonic regression is nonparametric and
overfits small calibration sets badly — it will fit a step function to a few
hundred points and feed confidently wrong `p_rel` into the composite score.
Platt scaling has two parameters and degrades gracefully.

*Implemented:* the calibrator is now resolved **per task** rather than being one
global default, so nothing has to be remembered at the call site.

| entry point                                | task                          | calibrator   |
|--------------------------------------------|-------------------------------|--------------|
| `run_re_pipeline.py` (gpu.py step 1)        | `ppi` (also chemprot/gad/ddi) | **isotonic** |
| `run_re_pipeline.py --task biored` (step 2) | `biored`                      | **platt**    |

In both `run_re_pipeline.py` and `train_re.py` the `--calibration` default is `None`
and is resolved after parsing: `platt` if `--task biored`, else `isotonic`. An
explicit `--calibration ...` still wins, and the resolved value is printed in the
startup `[task]` line. Where the fit actually happens is unchanged: the orchestrated
path always invokes `train_re.py --calibration none` and fits the calibrator in
`run_re_pipeline.py`'s step 3, so the trainer's own default only applies when
`train_re.py` is run directly — both now follow the same rule either way.

Not a code change to the calibration itself: `--calibration platt` was already a
supported choice in `calibration.py`; only which one is picked by default is new.

Related: model selection gets noisier on a 100-abstract dev split, so the
`load_best_model_at_end` epoch pick is less reliable here than on BioInfer.
Not a reason to avoid BioRED — a reason not to read small dev-F1 differences as
real.

### The practical answer

Do not frame it as replace-or-not; see step 6 of section 8. Run BioRED alongside
the PPI checkpoint: BioInfer keeps supplying well-labelled sentence-scoped binary
interactions, BioRED adds sign and type where it is confident. If the typed edges
prove noisy in practice, nothing has been lost.

If the goal is specifically more data, the cheapest route is not a bigger single
corpus but concatenating the label-compatible ones already supported —
ChemProt/DrugProt covers chemical-gene with far more instances, GAD covers
gene-disease, and both overlap BioRED's pair types.

---

## 7. Bonus: BioRED's novelty annotation

BioRED labels each relation as **Novel** vs not — i.e. whether the relation is a
new finding of that paper or restated prior knowledge. That is the same
distinction `sentences.py` currently approximates with hand-written anchor
centroids and a cosine margin, with no training data at all.

If the BigBIO KB schema carries this attribute through, it is a supervised signal
for the "is this an original result" question — potentially a trained replacement
for, or a validation set against, the `RESULT_MARGIN` heuristic.

**Verify before planning around it.** BigBIO's KB schema normalizes relations to
`{id, type, arg1_id, arg2_id, normalized}`, and a novelty attribute may simply be
dropped in that normalization. Check a raw document:

```python
import datasets
ds = datasets.load_dataset("bigbio/biored", name="biored_bigbio_kb", trust_remote_code=True)
print(ds["train"][0]["relations"][:3])
```

If novelty is absent there, it is still in the original BioRED distribution from
NCBI and could be joined back on relation identity.

---

## 8. Suggested order of work

1. `--print-types` against `bigbio/biored`; correct the type strings above.
2. Apply the `bigbio_to_re.py` edits; convert with `--val-frac 0`.
3. Hand-read ~50 positives (section 5.1) and decide whether the cue filter is
   needed.
4. Apply the `train_re.py` registry entry; train and calibrate.
5. Compare against the current PPI checkpoint on the same documents — the
   question is not raw F1 but whether the extra edge types are *right*.
6. Only then decide whether BioRED replaces the PPI model or runs alongside it.
   `relation_extraction.py` supports multiple RE models, so both is an option.
