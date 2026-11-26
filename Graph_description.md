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
Hover gives the name, kind, sentence count and target count. For a gene it also gives how the
text wrote it -- every surface normalized to that HGNC symbol, with its sentence count (*CDKN1A*:
*p21* 223, *CDKN1A* 41, *Cdkn1a* 15, ...; *AR*: *AR* 124, *AR-V7* 9, *androgen receptor* 3, ...),
also shown when the node is clicked -- and in every sentence listed for an edge, the words that
stand for its gene ends are set in bold blue, so you can see what was read as the gene.

**Variants** — a gene's sequence variants in HGVS notation (*p.* protein, *c.* coding DNA, *g.*
genomic DNA, *r.* RNA), from tmVar 3 as NCBI PubTator3 annotated the corpus documents. Hover shows
the most-reported ones with their paper counts; clicking the gene lists them all, grouped by level,
each with its rsID (linked to dbSNP) and the papers that report it. tmVar writes a bare change in
HGVS form (*T790M* → *p.T790M*) and names the gene the variant belongs to, which is the gene it is
listed under. They are read over each document's whole text (the abstract only, where PubTator3 has
no full text), so a variant is listed whether or not a sentence naming it became a triple.
Each gene's variants are also summarized by **HGVS type**, in the tooltip, at the top of the
click panel and in *All variants*: per level, how many are *missense*, *nonsense*, *synonymous*,
*unspecified (X)*, *deletion*, *insertion*, *duplication*, *delins*, *frameshift* or *extension*
(p.), and *substitution*, *deletion*, *insertion*, *duplication*, *delins* or *inversion* (c., g.,
r.); the list below is grouped the same way. *X* is left open because it means "any residue" in
current HGVS but a stop in older papers (*p.R216X*). *other* collects what tmVar wrote in a form
that is not one change (*p.L19,858R*, a patient carrying two). The type is set in
`gpu_bundle/variants.py` (`hgvs_type`), and `summaries/gene_variants.md` / `.tsv` carry it too.

The *Variants* button (next to *Disease names* / *Drug names*) writes the variants into the gene
nodes themselves, under the symbol in a small monospace face, most-reported first. *Variants per
gene* caps how many are drawn (3, 5, 10, 20 or all) and *in ≥ N papers* hides those fewer papers
report; a gene with more shows "+N more". *All variants* lists every gene's full list, grouped by
level, in the details box.
*Search gene or variant* takes a variant too (*T790M*, *p.T790M*, *c.2573T>G*): it selects every
gene carrying it and lists them with their paper counts.

**Edges** — colour is the relation carried by most of its sentences *in view*: `activates`
green, `inhibits` red, `binds` blue, `interacts` purple, `associated` grey, negations amber and
dashed. *activates* ends in an arrowhead, *inhibits* in an inhibition bar (T). Thickness and arrowhead / bar size scale with **independent publications**, square-rooted and
capped, so a paper repeating itself never outweighs several papers agreeing once each. An edge
takes its best-supported direction.

**Distance** — each edge is a spring whose rest length shrinks with its sentence count *in
view*, on a log scale: 1 sentence 200, 2 → 165, 4 → 130, 8 → 95, 16 → 60, floored at 45. Node
repulsion still spreads things out, so it is a tendency, not a ruler: in the default view the
median drawn length runs 224 (1 sentence), 183 (2), 171 (3–4), 134 (5–8), 82 (9–16). Nodes
already on screen keep their place when filters change (only newcomers are placed), so the
distances reflect the counts of the view in which each node first appeared.

**Co-mention links** — optional dashed grey lines to a chosen disease, drawn where sentences
name it but no model claimed a relation. Never coloured like a prediction.

## Controls

**Left column — how the picture looks.** The PubMed query behind the corpus; the year range;
whether disease and drug names are drawn (gene symbols always are); *Shrink periphery* and
*Expand center*, two radial reshapers that pull in outlying clusters or thin the crowded core
without reordering anything; tissue stacking, which pools same-tissue disease nodes onto one
overlapping spot; the *Keep diseases* list, with same-tissue diseases grouped under a tissue
heading that ticks them all at once; the co-mention picker; two sentence lenses -- a **Candidate driver**
toggle that keeps only edges with a sentence naming an amino-acid change, case-sensitive, one-letter
(*V600E*, *G12C*) or three-letter (*p.Gly12Cys*). Any punctuation may stand around it -- *KRAS
G12C,*, *V600E-mutant*, *KRAS-G12V*, *L858R/T790M*, *(G12C)*, HGVS *p.L747P* -- and it may be glued
to a whole gene symbol of the graph (*KRASG12C*, *TP53R273H*) or a Greek letter (*p110αH1047R*),
but not to the front of a longer symbol (*CDKN2A*, *PPARGC1A*); no letter or digit may follow, nor
*.X* (*H2A.X*). Figure / table panels are not counted -- any label right after the word, as in
*Fig. S6K*, *Fig.S2F*, *Fig. E3B* or *Figure A1C*, and supplementary panels in a list or PLOS-style
(*S8B Fig*) -- and neither is a stoplist of look-alikes: *S6K*, *S6R*, *E2F*, *E2A*, *G9A*, *G2M*,
*S1P*, *S100P*, *B2M*, the histones *H2A* / *H2B*, the drug code *H3B*, and the cell lines *T47D*,
*T98G*, *M059J*, *NCI-H295R*. (Until 2026-11-24 the change had to stand between spaces or in
parentheses, which caught 1,076 sentences; it now catches 2,110.) Beyond that list, general rules apply to one-letter changes: the two letters must
be different amino acids (not *U266R*, *LCAR-B38M*, *R10R*); cell-line derivatives (*H929R*: its stem
*H929* is used as a cell line; *ANBL-6/V10R*) and compound codes (*S3I-201*), both learned from the
corpus at build time, never count; a gene, chemical or cell line of that very name (HGNC, ChEBI, NCIt,
Cellosaurus: *C3G*, *H2S*) and any position below 10 (*Q2W* dosing) count only in mutation context -- *p.*,
a gene symbol right before, or *mutant* / *mutation* / *variant* / *allele* / *substitution* / *WT* and
the like within 30 characters. Phosphosite shorthand (*S15P-p53*, *S536P-p65*: the P is phospho, not
proline) is dropped for S / T / Y sites when *phospho* stands within 60 characters or a protein is
hyphened on behind, unless in mutation context. On the myeloma graph of 2026-11-24 this took it from 199 to 113 sentences; an **Enrich for fusions** toggle that keeps only
edges with a sentence saying *fusion* (as a word: not *infusion*, nor the organelle event: *mitochondrial fusion*, *fission and fusion*, nor a lab construct: *GST-CENP-N fusion protein*, *Myc-ER fusion*, *luciferase fusion reporter*) or *translocation* in the
chromosomal sense -- after a fusion pair or fusion gene (*EML4-ALK translocation*,
*ALK-translocated*), with a karyotype (*t(11;22)*), after *chromosomal* / *balanced* /
*reciprocal* / *truncating*, or as *translocation involving MLL*; or karyotype notation on its own, with no other word
(*t(4;14)*, *t(11;14)*, *t(9;22)(q34;q11)*, *ins(3;3)*, *inv(16)*, *der(22)*: myeloma and leukemia
papers write their translocations this way; a *t(28) = 3.1* t-statistic has no semicolon and never
matches); or a fusion named by its two genes in HGNC notation, *NUP98::NSD1*, *AGK::BRAF*,
*EWSR1::WT1* (both sides genes of the graph; 49 sentences say it without the word *fusion*), the
older single colon (*TMPRSS2:ERG*, *KIAA1549:BRAF*) only with no space around it and a fusion cue in
the sentence -- it also writes proximity pairs (*EGFR:GRB2* PLA), kinase:substrate (*MAPK8:JUN*),
interactions (*KRAS:CDKN2A*) and transgenes (*rag2:Myc*) -- never a protein moving
(*nuclear translocation of p65*, *translocation of YAP into the nucleus*); an **Enrich for ITD** toggle for internal tandem duplications (*ITD* in capitals as a word,
as in *FLT3-ITD*, or *internal tandem duplication* spelled out); an **Enrich for
overexpression** toggle (*overexpres* / *over-expres*, any case); an **Enrich for
amplification** toggle (*amplif*, any case, minus the laboratory sense: *PCR-amplified*, *amplified
using primers*, *amplification of the cDNA* -- but *MET amplification detected by ddPCR* counts); an **Enrich for deletion / loss** toggle (*deletion* / *deleted* / *co-deletion*, a gene's *del* or
*del(17p)*, *LOH* / *loss of heterozygosity* / *loss of function*, and *loss* tied to a gene or
chromosome such as *PTEN loss* or *9p loss*; not *EGFR exon 19 deletion*, engineered deletion
mutants, generic loss like *weight loss*, or a lost signal like *loss of ERK phosphorylation*); an **Enrich for truncation** toggle (*truncating*, *frameshift*, *nonsense*,
*stop-gain*, splice-site variants, and *Q55\** / *R213Ter* / *P191fs* notation; not engineered
truncated proteins, splicing biology or splicing-factor mutations, or *G719X*-style X notation -- X as a stop, *R363X*, counts only with a
truncation word within 60 characters -- or *nonsense-mediated decay*, an RNA process); an **Enrich for exon variants** toggle (*exon N* with a variant word close
after it -- *EGFR exon 19 deletion*, *HER2 exon 20 insertion*, *MET exon 14 skipping*, *exon 21
L858R* -- and compact forms such as *METex14*, *ex19del*, *EGFRex20ins*, *ERBB2ΔEx16*; not an exon
named without a variant, fusion breakpoints such as *Ex14-16*, the drug *EX527* or the gene
*Deltex1*); an **Enrich for copy number** toggle (*copy number* in any phrasing, *CNV* / *CNA* /
*SCNA* / *CN gain*, chromosome-arm and whole-chromosome gains and losses such as *1q gain* or
*loss of chromosome 10*, and *copy-neutral* LOH; not mtDNA, plasmid or viral copy number); an **Enrich for methylation** toggle (promoter / DNA / CpG
(hyper/hypo)methylation, methylomes, CpG islands, epigenetic silencing, and a bare *methylation*
only in a sentence that speaks of DNA; never RNA (*m6A*), histone or protein methylation) -- these
sit boxed under **Enrich for alteration:**; a second box, **Enrich for therapy:**, holds
**Resistance** (*resist\**, *chemo-* / *radioresistance*, *refractory*; not apoptosis, anoikis,
ferroptosis, oxidative-stress or insulin resistance, nor antibiotic selection, unless the sentence
also names a treatment), **Acquired resistance** (acquired / secondary / intrinsic / primary /
adaptive / de novo resistance, development or mechanisms of resistance, resistance mutations,
gatekeeper mutations, bypass signalling beside an inhibitor), **Synthetic lethality** (*synthetic(ally)
lethal*, *synthetic sick*, *collateral* / *selectively lethal*, *SL interaction / pair*, *MYC-SL*) and
**Biomarker** (*biomarker*, *predictive marker*, *companion diagnostic*, a feature that *predicts
response / sensitivity / benefit*; prognostic biomarkers included); a third box, **Enrich for
clinical relevance:**, holds **Clinical evidence** (clinical trials, a trial phase beside trial or
study wording, randomised / open-label / prospective / retrospective studies, NCT numbers, response
rates, *patients treated with / receiving*; plain *patients* is not enough), **Prognosis**
(*prognos\**, named survival endpoints and PFS / DFS / RFS / EFS, poor / shorter / median survival or
outcome, hazard ratios, Kaplan-Meier; *OS* only beside a survival cue and *HR* only as a number; cell
survival excluded) and **In vivo** (*in vivo*, xenografts / PDX / allografts, orthotopic or syngeneic
models, mice, rats, animals, mouse / murine models and tumours; not mouse cells, fibroblasts or
reagents); a fourth box, **Enrich for gene role and biology:**, holds **Tumor suppressor**
and **Oncogene**, which -- unlike every other Enrich button -- judge the GENE, not the sentence: an
edge counts when either end is a gene the **OncoKB Cancer Gene List** classes as a tumor suppressor
(*TSG*) / an oncogene (*ONCOGENE*), genes classed as both (*NOTCH1*, *TP63*, *EZH2*, ...) counting for
each; all its sentences pass and nothing is highlighted. The list (1,247 genes; 518 oncogenes and
408 tumor suppressors, 68 of them both) is read from `databases/oncokb_cancerGeneList.tsv`, fetched
2026-11-24 from `https://www.oncokb.org/api/v1/utils/cancerGeneList.txt`; 403 oncogenes and 296
tumor suppressors are genes of this graph. (Until then the two matched the words *oncogene* /
*tumor suppressor*, which few sentences use -- *Candidate driver* + *Oncogene* drew 8 nodes; it now draws
334.) The box's other two buttons still read the sentence: **Metastasis / EMT** (*metasta\**,
*EMT*, epithelial-mesenchymal transition and EMT markers; not the mesenchymal GBM subtype) and
**Immune** (*immune*, *immunity*, *immunotherapy*, checkpoint inhibitors, PD-1 / PD-L1 / CTLA-4 and
the other checkpoints, T / NK cells, CD8+, TILs, macrophages / TAMs, MDSCs, Tregs; not the lab
methods named *immuno-* such as immunoblot or immunohistochemistry, cell-cycle checkpoints, T-cell
lymphomas or a singular *TAM*, which is also tamoxifen). Within a box the toggles are ORed (a sentence matching any of those on counts); the boxes
are ANDed with each other (*Candidate driver* + *Resistance* keeps sentences such as *T790M confers resistance
to gefitinib*; *Candidate driver* + *Oncogene* only sentences that name a point mutation and speak of an
oncogene, so the view is small by design), and with the text match; zoom/fit (filters keep each node where it was and keep the camera unless something new
falls outside it, and Fit never zooms closer than 100%, so a small view is not blown up); and the details box,
headed *Sentences of interest*, which lists an edge's sentences with linked PMIDs when you click
it. Each section's explanation is behind a small round **i** beside its heading.

**Right column — what is in the picture.** Support thresholds (min unique sentences, min unique
publications); label font size; structural pruning (
*Min connections* at 2 when the page opens and *Min cluster size* off); gene and drug search; focus-on-gene with 1 or 2 hops (only the edges those hops travel); drug filter; free-text
sentence match (substring, or `/regex/` when wrapped in slashes); node-type, relation-type and
training-set filters; the score slider, which opens at 0.5 — 0.5 to 0.95 in steps of 0.05, then 0.96 to 0.99 — and a
live stats line.

## How filtering works

The order is deliberate, and it is what keeps the panel honest:

- **Support** — an edge's evidence is what survives score, year, relation type and training set.
  `Min unique sentences` and `Min unique publications` judge *that*.
- **The text query and the Enrich buttons are applied late.** They keep the edges with a
  matching sentence and show only those sentences. For the **Enrich buttons** `Min unique
  sentences` / `Min unique publications` are judged again on the lens-matching sentences alone,
  so every edge drawn shows at least the minimum of that kind of evidence. The **text query** is
  a search, not a lens: an edge that passes the thresholds on its own support is found if any one
  of its sentences contains the text (CTNNB1–SNAI1, five sentences with one naming OXCT1, is found
  by *OXCT1* at min 2 sentences). Typing can only remove edges from the view you had.
- **Min connections counts links of every relation type**, ticked or not (since 2026-11-25), so
  unticking types never hides a well-connected node; the bullet below otherwise holds. The stats
  line says how many edges it hides.
- **Structure comes last** — `Min connections` and `Min cluster size` describe the picture, so
  they run on what is actually drawn. `Min connections` is peeled until it holds, so every node
  left has at least that many links on screen, and a node that is not drawn (say, a chemical
  with a single link) cannot keep its neighbour in view.
- **Relation type and training set cut at the sentence**, not the edge, and an edge is drawn and
  counted as the relation and provenance of its *visible* sentences. A count reading "0 in view"
  can therefore never sit above a visible tag. Each training-set button shows *every* sentence
  that model kept, labelled with that model's own reading: *PPI model* shows only `interacts`
  (its one label), so `binds` and every other typed relation read 0 there. When both models keep
  a sentence they agree a relation is there, not on its type: PPI's `interacts` is not BioRED's
  `binds`. Since the two models share no relation type, there is no
  "kept by both" view: such a sentence appears under each model's button with that model's
  reading, and the tooltip gives both. Under *All* the same sentences carry both readings and
  match either relation box, so ticking just `interacts` equals *PPI model*, and a typed box
  equals *BioRED model*.
  Each button counts the edges holding at least one such sentence, so the counts sum past *All*.
- **Focus re-expands after the lens**, so a gene-focus view either contains its seed or is
  empty — never the seed's neighbours drawn around a missing centre.
- **A kept disease never vanishes silently.** Diseases are linked by `associated`, almost never
  by `binds`, so with `associated` unticked a disease ticked under *Keep diseases* has no link.
  It is then drawn on its own and the stats line says so (and to tick `associated`), instead of
  *Only their direct neighbours* blanking the canvas. The stats line also says when no relation
  type is ticked at all.
- **Direct neighbours are strictly one drawn step away, and on by default.** Ticking the first
  disease under *Keep diseases* switches *Only their direct neighbours* on; unticking it is then
  respected until the keep-list is emptied. The option keeps
  the kept diseases, every node with a *drawn* edge to one of them, and those edges only (since
  2026-11-24 not the neighbours' links with each other: glioblastoma has 180 direct edges, the old
  view drew hundreds more). It is judged after every other setting (relation boxes, Enrich buttons, min
  sentences / publications / cluster size; *Min connections* is not applied), so a gene whose own disease link is
  hidden is not shown, however it connects to the neighbours that are. The **text box** is the
  exception: it *searches* the neighbourhood rather than redrawing it. The neighbours are chosen
  from the view without the text, and the text then keeps the edges among them that match —
  otherwise typing *OXCT1* under *skin* emptied the canvas, because CTNNB1–SNAI1 names OXCT1
  but CTNNB1's and SNAI1's own links to melanoma do not. The cost: with
  `associated` unticked, gene-disease links are mostly gone, so the view shrinks to a few drugs
  (linked by `inhibits`) or to the kept diseases drawn on their own.
- **Kept diseases also scope the other edges' sentences.** An edge that does not touch a kept
  disease (gene–gene, gene–drug) keeps only sentences from papers about a kept disease — a
  paper's diseases being those its own sentences relate to anything — and is re-judged against
  `Min unique sentences / publications` on what is left. Under *brain*, PIK3CA–alpelisib goes
  from 49 sentences (breast cancer trials) to the one glioma paper, and the edges not touching a
  brain disease from 1,994 to 569, mentioning melanoma 5 times instead of 58 (those 5 are papers
  on melanoma brain metastases). A paper whose sentences relate no disease at all — about a third
  of them — is placed nowhere and does not appear on these edges in any disease view.

## Limits worth carrying

- **No disease–disease relations exist.** The routing pairs a disease with a gene or a chemical
  only. Ticking disease alone therefore leaves the canvas empty, and the stats line says to
  tick genes or chemicals. (Diseases you kept in the disease menu are still drawn, on their
  own.) Unticking a node type only ever removes nodes; it never brings back ones the other
  filters dropped.
- **Co-mentions are co-occurrence, not claims.**
- **`-i` shorthand is a drug class, not the gene.** `BRAFi`, `MEKi`, `CDK4 / 6i` name
  inhibitors, so `triples.py` turns them into chemical nodes (*BRAF inhibitor*, *MEK inhibitor*,
  *CDK4/6 inhibitor*). "p-AKT was upregulated by BRAFi" is therefore *BRAF inhibitor activates
  AKT1*, never *BRAF activates AKT1*. Genes that really end in `i` (`let-7i`, `Gtf2i`) are
  untouched, and so are knock-down constructs (`RNAi`, `CRISPRi`) and the CaSki cell line.
- **Score >=0.95 is effectively PPI-model-only, and >=0.99 selects nothing** — the BioRED
  checkpoint's composite tops out around 0.92 (PPI's around 0.95), as the build log warns.
- **A disease is a node only when its name resolves to exactly one MONDO term.** Abbreviations
  the normalizer cannot place (*PDAC*, *LUAD*) and plain names it resolves to several terms
  (*pancreatic cancer*, *CRC*) used to be dropped -- 95% of pancreatic mentions, which kept
  KRAS G12 nearly off pancreatic cancer. `DISEASE_SURFACE` now maps ~35 such surfaces that mean
  one thing in a cancer corpus (*PDAC*, *CRC*, *LUAD*, *pancreatic cancer*, *colon cancer*,
  *bladder cancer*, *RCC*, *HNSC*, ...; 7,802 mentions). Still dropped, as ambiguous: *ALL*,
  *BC*, *GC*, *PCa*, *PD*, *EC*, *OC*, *UC*, *PDA*. *RRMM* and *NDMM* (relapsed/refractory, newly
  diagnosed) go to *multiple myeloma*: they are stages of its treatment, not other diseases.
- **Drug abbreviations** *CFZ*, *Dex*, *BTZ* / *Btz* (any case) are carfilzomib, dexamethasone and
  bortezomib (`CHEMICAL_SURFACE`), also where the NER tagged them as genes: *BTZ* is an old alias
  of CASC3, which used to collect 506 bortezomib mentions.
- **A gene is a node only when its surface resolves to one HGNC symbol** -- and the normalizer
  often gives several: plain *AR* resolves to AR + FDXR + AKR1B1 + AREG (the approved symbol is an
  old alias of three other genes), and was dropped, so the AR node once had 7 sentences and no
  prostate cancer. `gene_symbol()` now takes the symbol the surface spells exactly (*AR*, *HGF*,
  *RAC1*, *BAP1*, *ATR*, *CDH1*; minus words such as *TF*, *cAMP*, *mOS*, *SAG*, *Pgp*), and
  `GENE_SURFACE` maps one-gene short names: *p21* → CDKN1A, *p16* → CDKN2A, *p27* → CDKN1B,
  *NRF2* → NFE2L2, *p38* → MAPK14, *p62* → SQSTM1, *ER* → ESR1, *TAZ* → WWTR1, *S6 / pS6* → RPS6,
  *TTF-1* → NKX2-1, *AR-V7* → AR. Families stay dropped (*ERK*, *NF-κB*, *FGFR*), and *AIS*
  (an AR alias, here adenocarcinoma in situ) is ignored. AR now has 79 edges and 184 sentences.
- **A few chemicals are blocked from being nodes at all** (`CHEMICAL_IGNORE`: tyrosine, glucose).
  ChEBI resolves them correctly, but their mentions are residue names and culture conditions, not
  compounds under study. The match is exact, so `2-deoxy-D-glucose` and
  `tyrosine kinase inhibitor` survive.
- **Tissue grouping and co-mention aliases are read from names**, so a disease whose name does
  not say its tissue is not grouped.
- **Counts sum past the edge total** for relation types and training sets, since an edge with
  mixed evidence counts under each.
- **Labels are sized in screen pixels, whatever the zoom or the view.** A label's size depends on
  its node's unique sentences in view only, on a log scale: at the default 50% a gene reads 11 px
  with one sentence, 12.5 with 2, 15.5 with 8, 17 with 16 and 20 px from 64 on (diseases and drugs
  at 80% of that). Zooming re-sizes the fonts so they stay that size on screen. The *Font size*
  slider stretches everything above an 8 px floor: about 8.6–10.4 px at 10%, 14–32 px at 100%.
  Labels are never faded out. (Sizes used to be in graph units, so on screen they were size × zoom:
  a big view zoomed out drew everything at the floor, a small one drew busy genes at 96 px, and
  4 → 5 sentences doubled a label from 13 to 26 — the same 50% looked different with every filter.)

## Regenerating

The page and its payload are baked at generation time — editing the script changes nothing until
it is re-run:

```bash
python high_confidence_g.py --data-root kaggle_working   # gene + disease + chemical
```

The Oncogene / Tumor suppressor roles come from `databases/oncokb_cancerGeneList.tsv` (the data
root is not in git, so a fresh checkout needs it fetched again:
`curl -o kaggle_working/databases/oncokb_cancerGeneList.tsv https://www.oncokb.org/api/v1/utils/cancerGeneList.txt`).
Without it the build warns and the two buttons keep nothing.

Cell-line names for the Candidate driver lens come from Cellosaurus (CC BY 4.0, SIB; optional,
122 MB, warned about when absent), looked for in `databases/` of the data root, then `gpu_bundle`:
`curl -o gpu_bundle/databases/cellosaurus.txt https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt`.

Publication years come from `databases/pmc_years.json`, built by `pub_years.py` (stage-2
step 4); re-run it if the corpus has grown, or the year slider will silently drop undated
evidence. It only harvests the `<?pub-year YYYY?>` stamp that stage 1 (`pub_year_xml.py`,
step 6c) writes into every `experimental_ner/` XML, so it needs no network, no NCBI and no
`pmids/pmid_pmc_ids.tsv`.

The variants come from `VARIANT/variants.json`, built by `gpu_bundle/variants.py` (stage-2
step 4b, internet). Locally, against the downloaded run:
`python gpu_bundle/variants.py --root kaggle_working --xml-dir gpu_bundle/experimental_ner`.
It caches every fetched document in `VARIANT/pubtator3_variants.jsonl`, so a second run only
fetches what is new (`--refresh` refetches everything). Without the file, gene nodes list no variants.

## Re-layout

*Re-layout this view* places the drawn nodes afresh from the drawn edges alone, with *Shrink
periphery* and *Expand center* at 0%. The original positions and slider values are kept: *Restore
layout*, or any change that draws a node the re-laid view did not have, brings them back.

## siRNA and shRNA constructs

*siMLL1* and *MLL1 shRNA* are drawn as their target gene, so a signed relation of theirs is
re-read from the sentence (`gpu_bundle/constructs.py`): *siMLL1 decreased AURKA* is drawn as KMT2A
`activates` AURKA. Where no sign can be read reliably (group labels such as *shX + DOX group*,
comparisons, no clear verb) the relation is drawn as `associated`. A gene named with its loss
(*knockdown of X*, *X-deficient*) is read the same way when its partner is a gene or a disease and
the verb is plain, also as *decreased in cells transfected with X siRNA* (but not when the siRNA
targets a fusion); with a drug the relation is `associated`; otherwise the model's label stays.

## Finding drugs and hidden nodes

*Search drug* and *Filter to drug/chemical* list the drugs of the curated target table (choosing one shows
its target genes) and, below them, every chemical node (choosing one shows its own edges). A name in
both is found in the target table first. When a searched node is not on the canvas, the panel names
the settings that hide it, such as *Relation type* or *Min connections*.

## Antibody and CAR-T treatments

*anti-PD-1*, *CD19 CAR T cells* and *MSLN-targeted ADCs* are drawn as the treatment (*anti-PD-1
antibody*, *CD19 CAR-T*, *MSLN ADC*, chemical nodes) and not as the target gene
(`gpu_bundle/antibodies.py`). Lab reagents, *anti-X siRNA* and *X blockade* keep their reading.

## Drugs as agents

The models predict no direction, and the first-mentioned entity is taken as the agent unless the
text says otherwise. Between a gene or disease and a drug the drug is made the agent
(`gpu_bundle/drug_direction.py`): *the BCL-2 inhibitor venetoclax* is drawn as venetoclax `inhibits`
BCL2. Sentences on drug response (*G expression and resistance to D*) labelled `activates` become
`associated`; those of a gene labelled `inhibits` (*FTO inhibits ATRA-induced differentiation*),
natural metabolites and negated relations keep their reading. Between two drugs, synergy or
combination wording (*the IC50 of A fell when combined with B*, *A+B group*) is `associated`;
antagonism (*A attenuated B-induced apoptosis*) keeps its sign. *phospho-*, methylcellulose and
*BetA* (read as *electron*) are not drawn, *IDA* is idarubicin, and *NOXA* is the gene PMAIP1.

## Stand-alone nodes

*Hide stand-alone nodes* takes off the canvas every node with no drawn edge to another node: a kept
disease whose links are all filtered out, or a node whose only edge is to itself. The button shows
how many the current filters leave, and it is judged afresh on every change, so a node returns once
it has a link again.

## List edges

When one gene of a list has a qualifying relation (score >= 0.5) with an entity outside the
list in the same sentence (before or after it: *X inhibits A, B and C* counts too), every other gene of that list that maps to an HGNC symbol gets the same
relation: *PDH1, IDH2, and SDHB gene expression is significantly under expressed in GBM* links
PDHA1 to glioblastoma although its own pair scored 0.387. Such a sentence is marked *from list*
with the genes it was taken from, carries the lowest donor score, and never changes the direction
or label of an edge the models read themselves. Not applied to *respectively* sentences, *between A
and B*, two names joined by a bare comma, slash-written pathways (*Rheb/mTOR axis*), +/− cell markers
(*CD4+PD-1+*), negative relations (*not associated*), or
donors read with conflicting signs. The exported JSON
holds model relations only; `--no-list-edges` turns the rule off. `audit_lists.py` lists what is
still missing.
