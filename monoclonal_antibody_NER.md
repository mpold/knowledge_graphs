# Recognising antibody therapies and combination regimens

What the pipeline does with them today, and what would recognise them.

## What happens today

**Antibodies are mostly never recognised.** Twelve named mAbs appear **680 times** in the
candidate sentences of the lung_adeno corpus, but only **181 (27%)** were tagged as entities at
all — and 18 of those were labelled `GENETIC`, not `CHEMICAL`:

| bevacizumab | nivolumab | pembrolizumab | cetuximab | trastuzumab | atezolizumab |
|---|---|---|---|---|---|
| 213 | 150 | 90 | 87 | 59 | 39 |

Then normalisation removes most of the survivors: across the whole corpus only **2** `-mab`
surfaces got a ChEBI id, while durvalumab, dupilumab, inetetamab (×42) and tinurilimab (×20) were
dropped. ChEBI is a small-molecule ontology; biologics are largely outside its scope. That is the
second stage of the same loss, and it sits inside the 19,460 "chemical: no ChEBI id" drops.

**Combinations are not represented at all.** There are **5,486** adjacent CHEMICAL pairs joined by
a combination cue (`and` 4166, `/` 662, `+` 278, `plus` 159, `combined with` 153,
`combination with` 68). Each becomes two independent nodes; the regimen as a unit does not exist
in the graph.

## What would recognise antibodies

**INN stem rules** are the cheapest high-precision win and need no resource: `-mab` plus the 2021
WHO replacements `-tug`, `-bart`, `-mig`, and `-cept` for fusion proteins. This is what catches
agents no dictionary has yet — inetetamab and tinurilimab are recent and will not be in ChEBI or
DrugBank for a while.

**Dictionary/ontology recognition**, because drug naming is closed-vocabulary in a way that gene
naming is not:

- **NCI Thesaurus** — best oncology coverage, and it carries regimens as well as agents
- **ChEMBL** and **DrugBank** — both include biologics with INN names
- **RxNorm** — clinical drugs including biologics
- **MeSH Supplementary Concept Records** — where new drugs live before they get descriptors
- **UMLS** as the umbrella, via MetaMapLite, QuickUMLS, or scispaCy's `EntityLinker`

**Better taggers**, to swap the model instead: **PubTator3 / TaggerOne**, **BERN2** (tags *and*
normalises against several ontologies), **HunFlair2**. The current chemical NER behaves like a
BC5CDR-trained model — CTD-chemical oriented, so small molecules are well covered and biologics
systematically are not, which is exactly the 27% observed.

## What would recognise combination therapies

This is not an NER problem, and treating it as one is why it is missing. Three layers, in
increasing effort:

1. **A coordination post-processor** over whatever drug NER is used: match `A plus B`, `A + B`,
   `A/B`, `A in combination with B`, `A-based chemotherapy` between two drug mentions and emit a
   combination as its own entity, or an edge typed `co-administered`. The 5,486 cues make this
   immediately productive, and it is a rule layer, not a model.
2. **A regimen dictionary** — **HemOnc** is the resource: oncology regimens, their components and
   lines of therapy, in BioPortal with an OMOP mapping. It is what recognises `FOLFOX`, `CHOP`,
   `EP`, `carboplatin–pemetrexed–pembrolizumab` as single concepts. NCIt covers many regimens too.
   Without this every combination named by acronym or by class (`platinum-based chemotherapy`,
   `chemoimmunotherapy`) is missed, and no component-level tagger can reconstruct them.
3. **A relation label.** CHEMICAL–CHEMICAL pairs already route to a DDI model, but its labels
   (`mechanism`, `effect`, `advise`, `int`) mean drug *interaction*, not co-administration, and
   `rel_cat` folds them all into `associated`. A `combination` label on that route would put
   regimens into the graph as first-class edges — the chemical–chemical channel already exists
   there (535 edges today).

## Where to start

The ordering that buys the most per unit of work:

1. **INN-stem regex + NCIt dictionary at the NER stage.** Recovers most of the 73% miss without
   retraining anything.
2. **Normalise CHEMICAL against ChEMBL or NCIt, falling back to ChEBI.** Keeps small molecules
   where they are and stops discarding every biologic. Watch the "no ChEBI id" bucket shrink as
   the measure.
3. **Coordination post-processor** for combinations, then HemOnc if regimen names are wanted.

One caution worth carrying: 18 antibody mentions came back as `GENETIC`. Antibodies are routinely
written by target — `anti-PD-1 antibody`, `anti-VEGF` — so a gene tagger grabs them, and whatever
is added should resolve `pembrolizumab` and `anti-PD-1` to the drug while keeping `PDCD1` the
gene. The `CHEMICAL_IGNORE`-style exclusions and the `GENETIC` false positives both need a look
after the change.

There is also a free evaluation metric: known-drug mentions in text versus mentions tagged is a
two-minute recall check that can be re-run against any tagger.
