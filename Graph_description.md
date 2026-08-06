# The graph viewer

What `high_confidence_g.py` produces, and how to read it.

The output is a self-contained HTML page — payload, vis-network library and all — showing
relations that the RE checkpoints extracted from the corpus. Every edge traces back to the
sentences that produced it, with clickable PMIDs. Nothing is fetched at view time; the file
works offline.

**One edge per unordered entity pair**, no matter how much evidence sits behind it. Each edge
carries every sentence that mentions that pair, and each sentence carries its own score, year,
PMID, relation label and the model that produced it.

## What you see

**Nodes** — genes (circle), diseases (diamond), chemicals (square), so type reads without
relying on colour. Size and label size grow with the node's unique sentences *in view*. Gene
fill also encodes drug-target status: deeper colour means more corpus chemicals target it,
green for DGIdb approved anti-neoplastic, amber for approved (other), pink for ChEBI-only.
Hover gives the name, kind, sentence count and target count.

**Edges** — colour is the relation carried by most of its sentences *in view*: `activates`
green, `inhibits` red, `binds` blue, `interacts` purple, `associated` grey, negations amber and
dashed. Thickness and arrowhead size scale with **independent publications**, square-rooted and
capped, so a paper repeating itself never outweighs several papers agreeing once each. An edge
takes its best-supported direction.

**Co-mention links** — optional dashed grey lines to a chosen disease, drawn where sentences
name it but no model claimed a relation. Never coloured like a prediction.

## Controls

**Left column — how the picture looks.** The PubMed query behind the corpus; the year range;
whether disease and drug names are drawn (gene symbols always are); *Shrink periphery* and
*Expand center*, two radial reshapers that pull in outlying clusters or thin the crowded core
without reordering anything; tissue stacking, which pools same-tissue disease nodes onto one
overlapping spot; the co-mention picker; zoom/fit; *Significance in view*; and the details box,
headed *Sentences of interest*, which lists an edge's sentences with linked PMIDs when you click
it. Each section's explanation is behind a small round **i** beside its heading.

**Right column — what is in the picture.** Support thresholds (min unique sentences, min unique
publications); label font size; structural pruning (min cluster size, min connections); gene and
drug search; focus-on-gene with 1 or 2 hops; drug filter; free-text sentence match (substring, or
`/regex/` when wrapped in slashes); node-type, relation-type and training-set filters; the score
slider — 0.5 to 0.95 in steps of 0.05, then 0.96 to 0.99 — and a live stats line.

## Significance in view

The graph shows which relations survive the filters; this ranks the **nodes** behind them, for
genes or for drugs. It is computed from the drawn edge list and recomputed on every redraw, so it
always describes the picture in front of you: open every filter and it reads as the corpus,
narrow them and it answers the same question of a slice. Co-mention links are excluded, since
this counts relations.

Three measures, chosen from a dropdown:

- **publications** — distinct papers behind the node's relations; the measure edge thickness
  uses, and the one a single talkative paper cannot inflate
- **partners** — distinct entities it is related to: breadth, not weight
- **sentences** — unique sentences supporting those relations

They disagree, which is the point of offering all three. Wide open, EGFR leads on publications
(882) while MALAT1 sits seventh — but MALAT1 is second on partners (674) and leads on sentences
(3253), and at score ≥0.99 it tops the publications list outright, because its evidence scores
unusually high.

Each row carries **z**: standard deviations above the mean on a log₁₀ scale, computed among the
node's own kind, since a gene is only remarkable among genes. The percentile is available too
(hover a row, or read the table) but it saturates — every member of a top six reads 99.9%, while
z still separates them (EGFR +6.3, CDH1 +5.6).

### How firm is the ranking?

**Bootstrap CIs**, a button in the table, resamples the view's **publications** with replacement
300 times, recomputes the measure and the whole ordering each time, and reports 95% percentile
intervals for every node's count and rank. Papers are the unit because they are the independent
one — sentences within an article are the same authors saying the same thing twice, and
resampling them would give intervals several times too tight.

It costs about 1.5 s over 5,288 genes and runs only when asked, since it is only worth having
once you have settled on a view. Any change to the graph discards it, so what you read always
belongs to what is on screen.

Read the rank interval before believing an ordering:

| rank | gene | count | count 95% CI | rank 95% CI |
|---|---|---|---|---|
| 1 | EGFR | 880 | 832–933 | 1–2 |
| 2 | AKT1 | 866 | 814–918 | 1–2 |
| 30 | HIF1A | 120 | 101–142 | 26–37 |
| 500 | UCA1 | 8 | 3–15 | 264–1447 |

The top of the list is firm and the tail is meaningless: a gene sitting 500th could belong
anywhere from 264th to 1447th. The count intervals track √n as they should — the width at
n = 880 is 108, against 2·1.96·√880 ≈ 116 for a Poisson count.

This measures **sampling variability of the corpus only**. It says nothing about whether a gene
is enriched relative to the literature at large; for that you need a background corpus and a
Fisher exact test, which this page does not do.

### Is it enriched, or just big?

The bootstrap says how firm the ranking is; it does not say whether a gene is *characteristic*
of the view. The table's **corpus papers**, **OR vs corpus** and **q (BH)** columns do, from a
2×2 over publications:

| | mentions it | doesn't |
|---|---|---|
| in view | a | view − a |
| rest of the corpus | c − a | N − view − (c − a) |

`c` and `N` come from a per-entity document index built over **every normalized triple in the
file** — 10,662 publications, score-unfiltered on purpose, since a corpus already filtered by the
cutoff you are testing is no denominator at all. It is embedded in the page (~104 KB) and keyed
by the same normalized ids the nodes use, so no string matching is involved.

Fisher's exact test, two-sided; the odds ratio carries a Haldane–Anscombe 0.5 correction and a
Woolf interval (an approximation, not the conditional MLE); **q** is Benjamini–Hochberg over the
entities with at least 5 corpus papers — roughly 1,000–1,300 genes depending on the view, the
rest being untestable at any sensible rate. It recomputes with the view, in well under a second.

What it finds is recognisable. Typing `immunotherapy` (159 papers):

| gene | view | corpus | OR | q |
|---|---|---|---|---|
| PDCD1 | 46 | 508 | 8.89 | 1.7e-22 |
| CD274 | 56 | 794 | 7.22 | 1.7e-22 |
| CTLA4 | 16 | 254 | 4.95 | 3.8e-05 |
| EGFR | 7 | 1224 | **0.38** | 0.06 |

The checkpoint axis is enriched and EGFR is *depleted* — which is the biology. Restricting to
2024–2026 surfaces GPX4 (OR 5.4) and SLC7A11 (3.6), the ferroptosis pair.

Two limits worth stating whenever you quote a q-value:

- **It cannot say "specific to lung adenocarcinoma."** Every paper in the denominator is already
  a lung paper, so the contrast is view-against-corpus. Answering the specificity question needs
  an outside corpus this page does not have.
- **The counts are model output.** At score ≥0.99 the enriched gene is MALAT1 while EGFR, KRAS,
  TP53 and ALK come out depleted — that is the extraction behaviour of the high-confidence
  checkpoint as much as the literature. Run the test at two cutoffs and trust what survives both.

The population is also "documents with an extracted candidate pair", not "documents mentioning
the entity": a gene named only in a methods section enters neither column.

**Full table view** replaces the canvas with the whole ranking: every node of that kind, all
three counts, percentile, z and a bar, sortable by any column, with its own year handles and
kind/measure dropdowns that mirror the panel's. Table rows do not navigate — they are there to be
read and sorted; the panel's top-six list is what centres a node in the graph, and if the filters
have removed it the details box says so rather than moving the view. When drugs are
selected, both views state that the list mixes clinically used drugs with lab chemicals used only
in experiments — ChEBI recognises both, so LY294002 ranks among the leaders and has never been
given to a patient.

## How filtering works

The order is deliberate, and it is what keeps the panel honest:

- **Support** — an edge's evidence is what survives score, year, relation type and training set.
  `Min unique sentences` and `Min unique publications` judge *that*, measured blind to the text
  query.
- **The text query is a lens applied late** — it selects among edges that already qualify and
  shows only the matching sentences. Typing can only remove edges from the view you had; it
  cannot collapse the graph by shrinking an edge's apparent evidence.
- **Structure comes last** — `Min connections` and `Min cluster size` describe the picture, so
  they run on what is actually drawn.
- **Relation type and training set cut at the sentence**, not the edge, and an edge is drawn and
  counted as the relation and provenance of its *visible* sentences. A count reading "0 in view"
  can therefore never sit above a visible tag, and *PPI-only* equals *All with just `interacts`
  ticked*.
- **Focus re-expands after the lens**, so a gene-focus view either contains its seed or is
  empty — never the seed's neighbours drawn around a missing centre.

## Limits worth carrying

- **No disease–disease relations exist.** The routing pairs a disease with a gene or a chemical
  only. Ticking disease alone therefore draws the disease landscape as unconnected nodes rather
  than a blank canvas.
- **Co-mentions are co-occurrence, not claims.**
- **Score >=0.95 and >=0.99 are effectively PPI-only** — the BioRED checkpoint's composite tops
  out around 0.926, as the build log warns.
- **`Min connections` is a single pass**, so nodes that lose links in it can finish below the bar.
- **Tissue grouping and co-mention aliases are read from names**, so a disease whose name does
  not say its tissue is not grouped.
- **Counts sum past the edge total** for relation types and training sets, since an edge with
  mixed evidence counts under each.
- On dense graphs **labels fade in as you zoom** (above 60 nodes), so names appear as you move in.

## Regenerating

The page and its payload are baked at generation time — editing the script changes nothing until
it is re-run:

```bash
python high_confidence_g.py --data-root kaggle_working --nodes all   # gene + disease + chemical
python high_confidence_g.py --data-root kaggle_working               # gene-only
```

Publication years come from `databases/pmc_years.json`, built by `pub_years.py`; run that first
if the corpus has grown, or the year slider will silently drop undated evidence.
