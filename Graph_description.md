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
overlapping spot; the co-mention picker; zoom/fit; and the details box, headed *Sentences of
interest*, which lists an edge's sentences with linked PMIDs when you click it.

**Right column — what is in the picture.** Support thresholds (min unique sentences, min unique
publications); label font size; structural pruning (min cluster size, min connections); gene and
drug search; focus-on-gene with 1 or 2 hops; drug filter; free-text sentence match (substring, or
`/regex/` when wrapped in slashes); node-type, relation-type and training-set filters; the score
slider — 0.5 to 0.95 in steps of 0.05, then 0.96 to 0.99 — and a live stats line.

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
