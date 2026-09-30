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
green for DGIdb approved anti-neoplastic, amber for approved (other), pink for no approved
drug at all -- DGIdb investigational and/or a ChEBI target role only, which is a statement about
approval status, not about which database the gene came from.
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
publications); label font size; structural pruning (min cluster size, min connections, both
floored at 2); gene and drug search; focus-on-gene with 1 or 2 hops; drug filter; free-text
sentence match (substring, or `/regex/` when wrapped in slashes); node-type, relation-type and
training-set filters; the score slider — 0.5 to 0.95 in steps of 0.05, then 0.96 to 0.99 — and a
live stats line.

The floor of 2 on the two structural controls matters beyond tidiness: at *Min connections* 1 the
view holds every partner every gene has, which leaves the enrichment test below with nothing
outside the view to contrast against.

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
| 1 | EGFR | 867 | 812–925 | 1–2 |
| 2 | AKT1 | 864 | 811–920 | 1–2 |
| 30 | SNAI2 | 119 | 99–145 | 26–36 |
| 500 | STAT6 | 8 | 3–13 | 292–1408 |

The top of the list is firm and the tail is meaningless: a gene sitting 500th could belong
anywhere from 292nd to 1408th, and EGFR and AKT1 are statistically indistinguishable, so "the top
gene" is not a claim this corpus supports. The count intervals track √n as they should — the
width at n = 867 is 113, against 2·1.96·√867 ≈ 115 for a Poisson count.

This measures **sampling variability only** — whether the ranking would survive another draw of
the literature. Whether a gene is *characteristic* of the view is the next section's question.

### Is it enriched, or just big?

The bootstrap says how firm the ranking is; it does not say whether a gene is *characteristic*
of the view. The table's **corpus papers**, **OR vs corpus** and **q (BH)** columns do, from a
2×2:

| | mentions it | doesn't |
|---|---|---|
| in view | a | view − a |
| rest of the corpus | c − a | N − view − (c − a) |

**The unit follows the ranking.** By publications the cells count papers against the corpus'
papers; by **partners** they count entities against the corpus' entities, drawn from a separate
partner background — testing breadth against a paper denominator would compare two different
things. Ranking by sentences falls back to the publication test on purpose: sentences are
pseudo-replicates, and a Fisher test on them would be anticonservative by a large factor.

### One corpus for every lung project

`c` and `N` come from a per-entity index built over **every normalized triple**, score-unfiltered
on purpose, since a corpus already filtered by the cutoff you are testing is no denominator at
all. The index is the **union across all `lung_*` projects** — currently **15,538 publications
and 9,685 entities** from `lung_adeno`, `lung_large`, `lung_neuroendocrine`, `lung_pancoast`,
`lung_small` and `lung_squamous` — so an odds ratio computed in one directory means the same
thing as one computed in another.

A union, not a sum: the corpora overlap (one paper today, between `lung_adeno` and
`lung_pancoast`), and adding counts would corrupt the denominator the moment two queries pull the
same article. Each project caches its own slice as sets under `databases/corpus_contrib.json`,
keyed to the source file's size and mtime; a project that has never been run is *absent* rather
than assumed, and both the console and the table name the corpora that contributed. The merged
index rides in the page at ~193 KB, keyed by the same normalized ids the nodes use, so no string
matching is involved.

### Reading the columns

Fisher's exact test, two-sided; the odds ratio carries a Haldane–Anscombe 0.5 correction and a
Woolf interval (an approximation, not the conditional MLE); **q** is Benjamini–Hochberg over the
entities with at least 5 in the corpus. It recomputes with the view, in well under a second.

Rows that cannot be tested show their corpus count with a blank OR and q — never a blank
denominator, which would imply the corpus holds fewer papers than the view. Two reasons a row
goes untested, both named in the subtitle: it falls below the 5-in-corpus floor, or it has
**nothing outside the view** to contrast against. If every row hits the second case the table
says so outright rather than hiding the columns, and narrowing the view restores the contrast.

What it finds is recognisable. Typing `immunotherapy` (86 papers in view):

| gene | view | corpus | OR | q |
|---|---|---|---|---|
| PDCD1 | 34 | 741 | 13.70 | 2.7e-21 |
| CD274 | 40 | 1198 | 10.75 | 6.9e-21 |
| CTLA4 | 16 | 327 | 11.38 | 4.4e-10 |
| EGFR | 5 | 1396 | **0.68** | 0.465 |

The checkpoint axis is enriched and EGFR is not — which is the biology.

### The Significance slider

A q cutoff, in the left column above *Significance in view* and in the table below its summary;
either handle moves the other and both views apply the same cut. Ten stops on the conventional
thresholds — off, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6 — because q is read on a log
scale. On `lung_adeno` at defaults, `q ≤ 0.05` takes the graph from 29,310 edges over 3,960 genes
to 1,249 over 62.

Genes and drugs are judged and untested counts as non-significant; **diseases stay**, since they
carry no ranking of their own, and survive only where they still connect to something
significant. The ranking is computed from the view *before* the cut and the cut only hides what
it judged — otherwise each notch would move the very numbers it filters on, and the q beside a
node would not be the q it was judged by.

Two limits worth stating whenever you quote a q-value:

- **It cannot say "specific to lung cancer."** Every paper in the denominator is a lung paper, so
  the contrast is view-against-corpus. What it *can* now say is that an entity is characteristic
  of one lung subtype against the rest, since the denominator spans all six — running the same
  query in `lung_small` and reading it against the shared corpus is a subtype question, not a
  cancer-versus-everything one.
- **The counts are model output.** At score ≥0.99 the enriched gene is MALAT1 while EGFR, KRAS,
  TP53 and ALK come out depleted — that is the extraction behaviour of the high-confidence
  checkpoint as much as the literature. Run the test at two cutoffs and trust what survives both.

The population is also "documents with an extracted candidate pair", not "documents mentioning
the entity": a gene named only in a methods section enters neither column.

### Full table view

Replaces the canvas with the whole ranking: every node of that kind, all three counts,
percentile, z, the enrichment trio and a bar, sortable by any column, with its own year handles,
significance slider and kind/measure dropdowns mirroring the panel's. `Graph view` or `Esc`
returns you.

It stops short of the **right panel, which stays live** — every threshold, filter and the score
slider keep working while you read, and the ranking under them recomputes. Tighten to score ≥0.99
without leaving the table and it goes from 5,288 rows led by EGFR to 2,519 led by MALAT1.

Columns appear only once they carry something: the bootstrap pair after you press the button, the
enrichment trio when there is a contrast to test. Twelve columns otherwise run off the right edge,
which is how a computed enrichment can look like a missing one. The name column is pinned while
the rest scrolls sideways.

Table rows do not navigate — they are there to be read and sorted; the panel's top-six list is
what centres a node in the graph, and if the filters have removed it the details box says so
rather than moving the view. When drugs are selected, both views state that the list mixes
clinically used drugs with lab chemicals used only in experiments — ChEBI recognises both, so
LY294002 ranks among the leaders and has never been given to a patient.

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
- **Partner enrichment degenerates more readily than publication enrichment.** A gene's whole
  neighbourhood tends to survive into the view while its papers do not, so at a loose view many
  genes have no partner outside it and go untested. Prune first — the floor of 2 on
  `Min connections` exists for this.
- **A few chemicals are blocked from being nodes at all** (`CHEMICAL_IGNORE`: tyrosine, glucose).
  ChEBI resolves them correctly, but their mentions are residue names and culture conditions, not
  compounds under study. The match is exact, so `2-deoxy-D-glucose` and
  `tyrosine kinase inhibitor` survive.
- **Tissue grouping and co-mention aliases are read from names**, so a disease whose name does
  not say its tissue is not grouped.
- **Counts sum past the edge total** for relation types and training sets, since an edge with
  mixed evidence counts under each.
- On dense graphs **labels fade in as you zoom** (above 60 nodes), so names appear as you move in.

## Regenerating

The page and its payload are baked at generation time — editing the script changes nothing until
it is re-run:

```bash
python high_confidence_g.py --data-root kaggle_working   # gene + disease + chemical
```

Publication years come from `databases/pmc_years.json`, built by `pub_years.py` (stage-2
step 4); re-run it if the corpus has grown, or the year slider will silently drop undated
evidence. It only harvests the `<?pub-year YYYY?>` stamp that stage 1 (`pub_year_xml.py`,
step 6c) writes into every `experimental_ner/` XML, so it needs no network, no NCBI and no
`pmids/pmid_pmc_ids.tsv`.

The enrichment denominator is shared across projects, so a run also writes its own slice to
`databases/corpus_contrib.json` and merges whatever slices the sibling `lung_*` projects have
written. **Run the script once in each project** to populate them; until then the console names
the corpora still missing, and the denominator is honestly smaller. Re-running a project after
its pipeline changes refreshes its slice automatically — the cache is keyed to the source file's
size and mtime.
