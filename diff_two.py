"""Build a difference graph over two knowledge graphs.

Reads two vis-network graphs, keeps the pairs that are prominent in either, and
renders them as one network in which node colour is a continuous balance between
the two, and edge colour is the pair's provenance.

    python diff_two.py

Everything file-related is asked for at the prompt — the two inputs, the output
and the graph's name. Nothing is hardcoded, so the script has no idea which
diseases it is diffing until it is run.

    input file name 1:   the FIRST graph  -> the red pole of the scale
    input file name 2:   the SECOND graph -> the blue pole, and the source of
                         the inlined vis-network runtime
    output file name:
    graph name:          the page's heading

Only the standard library is used. The vis-network runtime is lifted out of the
second input and inlined, so the output opens offline.

The minimum relationship score is a slider (SCORE_STEPS), so the whole model —
which pairs qualify, each gene's connectivity, the balance index — is recomputed
in the browser as it moves. Python therefore ships the *raw* evidence rather than
one baked snapshot: every sentence at or above the lowest step, plus per-gene
degree tables for each step.

Shipping only the lowest step's selection is sufficient for every step above it:
sentences passing 0.99 are a subset of those passing 0.95, so a pair drawn at any
step is always drawn at the floor.

Prominence is measured in *independent publications*, not sentences. In these
graphs the overwhelming majority of pairs rest on a single paper (96.5% of them
in the SCLC/LUAD pair this was first built for), so a sentence count mostly ranks
how often one author repeated themselves rather than how well attested a relation
is.
"""

import io
import json
import math
import os
from collections import defaultdict

# Bare names resolve here; an absolute path is taken as given.
HERE = os.path.dirname(os.path.abspath(__file__))

# Selectable minimum relationship scores, following the source graphs' own
# scheme: coarse 0.05 steps through the low band, then 0.01 steps from 0.95 up,
# where a single point still changes the graph materially.
# The floor drives the payload — every sentence at or above it ships.
SCORE_STEPS = ([round(0.75 + 0.05 * i, 2) for i in range(4)]      # 0.75 .. 0.90
               + [round(0.95 + 0.01 * i, 2) for i in range(5)])   # 0.95 .. 0.99
FLOOR = min(SCORE_STEPS)
# Where the slider starts. Resolved to the nearest step rather than a hardcoded
# index, so it survives edits to SCORE_STEPS.
DEFAULT_SCORE = 0.95
DEFAULT_STEP = min(range(len(SCORE_STEPS)),
                   key=lambda i: abs(SCORE_STEPS[i] - DEFAULT_SCORE))

# A pair is drawn if either graph backs it with this many distinct publications.
# Shared pairs are always drawn, whatever their support.
MIN_PUBS = 2


# --------------------------------------------------------------------------
# Reading the source graphs
# --------------------------------------------------------------------------

def extract_data(path):
    """Return the `const DATA={...}` object embedded in a source graph."""
    with io.open(path, encoding="utf-8") as fh:
        src = fh.read()
    start = src.index("{", src.index("const DATA="))
    return json.JSONDecoder().raw_decode(src[start:])[0]


def extract_vis_runtime(path):
    """Return the vis-network bundle — the first <script> of a source graph."""
    with io.open(path, encoding="utf-8") as fh:
        src = fh.read()
    start = src.index("<script>") + len("<script>")
    return src[start:src.index("</script>", start)]


def collapse_to_pairs(data):
    """Fold directed, category-split edges into undirected gene pairs.

    Keeps every sentence at or above FLOOR, tagged with its own score, so the
    page can re-filter at any step. Returns {(gene_a, gene_b): [sentence, ...]}
    with gene_a < gene_b, so A->B and B->A land on one key.
    """
    pairs = defaultdict(list)
    for edge in data["edges"]:
        key = tuple(sorted((edge["from"], edge["to"])))
        for s in edge["sents"]:
            if s["sc"] >= FLOOR:
                pairs[key].append({
                    "p": s["pmid"],
                    "t": s["text"],
                    "c": round(s["sc"], 4),
                    "y": s["yr"],
                    "g": edge["cat"],
                })
    return dict(pairs)


def present_at(sents, score):
    return any(s["c"] >= score for s in sents)


def pubs_at(sents, score):
    return len({s["p"] for s in sents if s["c"] >= score})


# --------------------------------------------------------------------------
# Colour — a diverging red<->blue ramp
# --------------------------------------------------------------------------
#
# Red pole = pure input 1, blue pole = pure input 2, neutral grey midpoint = a
# gene equally central to both. Blue<->red is the documented diverging pair: the
# poles read as opposites and stay apart under every simulated colour-vision
# deficiency (worst case protanopia, dE 15.1 — the >=8 target); blue<->green
# would not.
#
# The arms mirror the blue ramp's OKLCH lightness, and each mode's direction is
# deliberate: on a light surface the poles are dark and the midpoint recedes
# toward white; on a dark surface that inverts, poles bright and midpoint sinking
# into the background. Both directions are monotonic in lightness, which is the
# check a diverging ramp answers to (a categorical validator fails it by design).

RAMP_STEPS = 41

RAMP_MODES = {
    # anchors run red pole -> midpoint -> blue pole
    "light": {
        "anchors": ["#671215", "#a6272a", "#e14c4a", "#f1968e", "#f0efec",
                    "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"],
        # Outline shifts away from the surface. Without it the near-white
        # midpoint (1.12:1) would vanish; with it every node clears 2.67:1.
        "border_dl": -0.26,
    },
    "dark": {
        "anchors": ["#ed7f78", "#e14c4a", "#ba3334", "#902123", "#383835",
                    "#184f95", "#256abf", "#3987e5", "#6da7ec"],
        "border_dl": +0.24,
    },
}


def _srgb_to_linear(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _linear_to_srgb(c):
    c = min(1.0, max(0.0, c))
    return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


def _cbrt(x):
    return math.copysign(abs(x) ** (1.0 / 3.0), x)


def _hex_to_linear(h):
    h = h.lstrip("#")
    return [_srgb_to_linear(int(h[i:i + 2], 16) / 255.0) for i in (0, 2, 4)]


def _linear_to_oklab(rgb):
    r, g, b = rgb
    l = _cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b)
    m = _cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b)
    s = _cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b)
    return [
        0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
        1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
        0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
    ]


def _oklab_to_linear(lab):
    L, a, b = lab
    l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3
    m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3
    s = (L - 0.0894841775 * a - 1.2914855480 * b) ** 3
    return [
        4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
        -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
        -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s,
    ]


def _oklab_to_hex(lab):
    # floor(x + 0.5) to match the JS Math.round the ramp was authored against;
    # Python's round() breaks .5 ties to even and would shift some steps.
    chans = (
        int(math.floor(_linear_to_srgb(c) * 255 + 0.5))
        for c in _oklab_to_linear(lab)
    )
    return "#" + "".join("%02x" % min(255, max(0, c)) for c in chans)


def build_ramp():
    """Interpolate each mode's anchors in OKLab into fill/outline pairs."""
    ramp = {}
    for mode, cfg in RAMP_MODES.items():
        labs = [_linear_to_oklab(_hex_to_linear(h)) for h in cfg["anchors"]]
        fills, borders = [], []
        for i in range(RAMP_STEPS):
            pos = i / (RAMP_STEPS - 1.0) * (len(labs) - 1)
            k = min(int(math.floor(pos)), len(labs) - 2)
            t = pos - k
            lab = [a + (b - a) * t for a, b in zip(labs[k], labs[k + 1])]
            fills.append(_oklab_to_hex(lab))
            lightness = min(0.96, max(0.06, lab[0] + cfg["border_dl"]))
            borders.append(_oklab_to_hex([lightness, lab[1], lab[2]]))
        ramp[mode] = {"fill": fills, "border": borders}
    return ramp


# --------------------------------------------------------------------------
# The difference model
# --------------------------------------------------------------------------
#
# The balance index the page computes from these tables is
#
#     index = (degree_1/total_1 - degree_2/total_2)
#           / (degree_1/total_1 + degree_2/total_2)
#
# i.e. +1 pure input 1 .. 0 balanced .. -1 pure input 2.
#
# Expressing each degree as a share of its *own* graph's connectivity is the
# whole point, and it is not cosmetic. The two graphs are never the same size —
# one literature has simply been written about more — so on raw counts a gene
# central to BOTH is dragged toward whichever graph is larger and reads as if it
# belonged there. Normalising asks the same question of each graph ("what
# fraction of this graph's wiring runs through this gene?"), which is answerable
# on equal terms however lopsided the two literatures are.
#
# Dividing by the sum, rather than just subtracting, makes the result a
# proportion instead of a magnitude: 3-vs-1 partners and 30-vs-10 read as the
# same lean. Magnitude is not lost — node size carries it.
#
# The page's "How balance works" panel demonstrates the failure live: it shows
# the same contrast computed on unnormalised counts beside this one, and counts
# the genes whose labels disagree. It derives that from whatever two graphs were
# loaded, so it stays honest for any pair.

def select_pairs(pairs1, pairs2):
    """Pairs worth drawing, judged at FLOOR (a superset of every higher step)."""
    shared = set(pairs1) & set(pairs2)
    return shared | {
        key
        for pairs in (pairs1, pairs2)
        for key, sents in pairs.items()
        if pubs_at(sents, FLOOR) >= MIN_PUBS
    }


def degree_tables(pairs, genes):
    """Per-gene degree and whole-graph total degree, one entry per score step."""
    table = {g: [0] * len(SCORE_STEPS) for g in genes}
    totals = [0] * len(SCORE_STEPS)
    for key, sents in pairs.items():
        for i, score in enumerate(SCORE_STEPS):
            if not present_at(sents, score):
                continue
            totals[i] += 2                      # a pair contributes two endpoints
            for gene in key:
                if gene in table:
                    table[gene][i] += 1
    return table, totals


def build_model(pairs1, pairs2):
    selected = select_pairs(pairs1, pairs2)
    edges = [
        {"f": key[0], "t": key[1],
         "s": pairs1.get(key, []), "l": pairs2.get(key, [])}
        for key in sorted(selected)
    ]
    drawn = sorted({g for key in selected for g in key})
    table1, totals1 = degree_tables(pairs1, set(drawn))
    table2, totals2 = degree_tables(pairs2, set(drawn))
    nodes = [
        {"id": g, "sd": table1[g], "ld": table2[g]} for g in drawn
    ]
    return nodes, edges, {"s": totals1, "l": totals2}


# --------------------------------------------------------------------------
# Page
# --------------------------------------------------------------------------

TEMPLATE = r"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<script>__VIS__</script>
<style>
:root{
  color-scheme: light;
  --surface-1:#fcfcfb; --surface-2:#f2f2ef; --border:#dcdcd6;
  --text-primary:#0b0b0b; --text-secondary:#52514e; --text-muted:#77766f;
  --edge-sclc:#d13b3c; --edge-luad:#2a78d6; --edge-both:#a5a49d;
  --grad:__GRADL__;
  --net-bg:#fcfcfb;
  /* evidence tooltip — the source graphs' own palette */
  --tip-bg:#ffffff; --tip-fg:#1a1a1a; --tip-bd:#999999; --tip-rule:#e3e3e3;
  --tip-mut:#5b6677; --tip-more:#888888; --tip-head:#333333;
  --tip-pm-bg:#eef3fb; --tip-pm-fg:#2b6cb0; --tip-pm-bg-hover:#d6e6fb;
  --tip-mark:#ffe680;
}
@media (prefers-color-scheme: dark){
  :root:where(:not([data-theme="light"])){
    color-scheme: dark;
    --surface-1:#1a1a19; --surface-2:#232322; --border:#3a3a37;
    --text-primary:#ffffff; --text-secondary:#c3c2b7; --text-muted:#8f8e85;
    --edge-sclc:#e66767; --edge-luad:#3987e5; --edge-both:#7a7972;
    --grad:__GRADD__;
    --net-bg:#1a1a19;
    --tip-bg:#232322; --tip-fg:#e9e8e1; --tip-bd:#55554f; --tip-rule:#3a3a37;
    --tip-mut:#a9a89f; --tip-more:#8f8e85; --tip-head:#d8d7cf;
    --tip-pm-bg:#1e3555; --tip-pm-fg:#8fbdf0; --tip-pm-bg-hover:#27456e;
    --tip-mark:#ffe680;
  }
}
:root[data-theme="dark"]{
  color-scheme: dark;
  --surface-1:#1a1a19; --surface-2:#232322; --border:#3a3a37;
  --text-primary:#ffffff; --text-secondary:#c3c2b7; --text-muted:#8f8e85;
  --edge-sclc:#e66767; --edge-luad:#3987e5; --edge-both:#7a7972;
  --grad:__GRADD__;
  --net-bg:#1a1a19;
  --tip-bg:#232322; --tip-fg:#e9e8e1; --tip-bd:#55554f; --tip-rule:#3a3a37;
  --tip-mut:#a9a89f; --tip-more:#8f8e85; --tip-head:#d8d7cf;
  --tip-pm-bg:#1e3555; --tip-pm-fg:#8fbdf0; --tip-pm-bg-hover:#27456e;
  --tip-mark:#ffe680;
}
*{box-sizing:border-box}
body{margin:0;background:var(--surface-1);color:var(--text-primary);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
  /* a hard viewport box, so the rail scrolls internally instead of the page
     growing and pushing the pinned selection block out of sight */
  display:flex;flex-direction:column;height:100vh;overflow:hidden;}

/* One rail on the right carries title, controls and legend; the graph gets
   everything else. The table swaps in over the canvas rather than stacking
   below it — parked under a full-height canvas it opened off-screen, which
   read as a dead button. */
.stage{flex:1;display:flex;align-items:stretch;min-height:0;}
.canvas{flex:1;min-width:0;position:relative;background:var(--net-bg);}
#net{position:absolute;inset:0;}
/* Only the upper rail scrolls; the selection block is pinned to the rail's foot
   so a click's detail can never land below the fold. */
.side{flex:none;width:286px;overflow:hidden;
  border-left:1px solid var(--border);background:var(--surface-1);
  display:flex;flex-direction:column;}
.railtop{flex:1;min-height:0;overflow-y:auto;padding:16px 18px;
  display:flex;flex-direction:column;gap:12px;}
h1{margin:0;font-size:20px;font-weight:700;letter-spacing:-0.01em;line-height:1.25;}
.sub{margin:9px 0 0;color:var(--text-secondary);font-size:12.5px;}
.sub + .sub{margin-top:3px;}
.controls{display:flex;flex-wrap:wrap;gap:9px;align-items:center;}
select,button{font:inherit;font-size:13px;padding:5px 9px;border-radius:6px;
  border:1px solid var(--border);background:var(--surface-1);color:var(--text-primary);}
button{cursor:pointer} button:hover{background:var(--surface-2)}
/* both thresholds on one line each: label · track · value */
.sliders{display:flex;flex-direction:column;gap:7px;}
.srow{display:flex;align-items:center;gap:8px;}
.slab{flex:none;width:66px;font-size:12.5px;color:var(--text-secondary);}
/* these carry a prompt-supplied label, so they must survive any length */
.balrow .slab{width:80px;font-size:11.5px;line-height:1.25;white-space:nowrap;
  overflow:hidden;text-overflow:ellipsis;}
.srow input[type=range]{flex:1;min-width:0;accent-color:var(--edge-luad);cursor:pointer;}
.sval{flex:none;font-size:13px;font-weight:600;color:var(--text-primary);
  font-variant-numeric:tabular-nums;min-width:30px;text-align:right;}
.slabel{font-size:13px;font-weight:600;color:var(--text-secondary);display:block;
  margin-bottom:5px;}
input[type=search]{width:100%;font:inherit;font-size:13px;padding:5px 8px;border-radius:6px;
  border:1px solid var(--border);background:var(--surface-1);color:var(--text-primary);}
input[type=search]:focus{outline:2px solid var(--edge-luad);outline-offset:1px;}
/* Link stays in text ink with an underline as its affordance — blue is a pole of
   the scale here, and a blue link beside the graph would read as an encoding. */
/* context in secondary ink, the actionable instruction bold in primary */
.hgnc{margin:14px 0 0;font-size:12px;color:var(--text-secondary);line-height:1.45;}
.hgnc strong{font-weight:600;color:var(--text-primary);}
.hgnc a{color:inherit;text-decoration:underline;text-underline-offset:2px;}
.hgnc a:hover{text-decoration-thickness:2px;}
.hgnc a:focus-visible{outline:2px solid var(--edge-luad);outline-offset:2px;}
.hint{font-size:11.5px;color:var(--text-muted);margin:5px 0 0;}
.hint code{font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:11px;}
.bad{color:var(--edge-sclc);}
.pair{display:flex;gap:7px;align-items:center;}
.pair input[type=search]{flex:1;min-width:0;}
.pair select{flex:none;padding:5px 4px;}
.count{margin:10px 0 0;font-size:12.5px;color:var(--text-secondary);}
.count b{color:var(--text-primary);font-weight:600;}
.resetbtn{margin-top:10px;width:100%;font-size:12.5px;}
.sect{font-size:11px;font-weight:600;letter-spacing:0.04em;text-transform:uppercase;
  color:var(--text-muted);margin:0 0 9px;}
/* Collapsed legends keep a colour preview in the summary, so the encoding is
   still readable at a glance without spending 405px of rail on it. */
.leg summary{list-style:none;cursor:pointer;display:flex;align-items:center;gap:9px;
  padding:1px 0;}
.leg summary::-webkit-details-marker{display:none;}
.leg summary .sect{margin:0;}
.leg summary::after{content:'';margin-left:auto;width:6px;height:6px;flex:none;
  border-right:1.5px solid var(--text-muted);border-bottom:1.5px solid var(--text-muted);
  transform:rotate(45deg) translate(-2px,-2px);transition:transform .15s;}
.leg[open] summary::after{transform:rotate(-135deg);}
.leg summary:hover .sect{color:var(--text-secondary);}
.leg .body{padding-top:12px;}
.mgrad{width:56px;height:8px;border-radius:4px;flex:none;border:1px solid var(--border);
  background:linear-gradient(to right,var(--grad));}
.mkeys{display:flex;gap:4px;flex:none;}
.mkeys i{width:13px;height:3px;border-radius:2px;}
.leg[open] .mgrad,.leg[open] .mkeys{display:none;}
/* horizontal colourbar: a vertical ramp cost ~110px of rail for no more
   information, and this matches the mini gradient in the collapsed summary */
.hramp{position:relative;display:block;width:100%;height:12px;border-radius:6px;
  border:1px solid var(--border);background:linear-gradient(to right,var(--grad));}
/* the excluded ends of the scale, greyed back over the ramp itself */
.cut{position:absolute;top:0;bottom:0;background:var(--surface-1);opacity:.76;
  pointer-events:none;border-radius:6px;}
.mgrad{position:relative;}
.mcut{position:absolute;top:0;bottom:0;background:var(--surface-1);opacity:.76;
  pointer-events:none;border-radius:4px;}
.hticks{display:flex;justify-content:space-between;gap:6px;margin-top:5px;}
.hticks span{font-size:11px;color:var(--text-secondary);line-height:1.2;white-space:nowrap;}
.hticks b{color:var(--text-primary);font-weight:600;}
.skey{display:flex;flex-direction:column;gap:9px;}
.lg{display:flex;align-items:center;gap:9px;font-size:13px;color:var(--text-primary);}
.ln{width:22px;height:3px;border-radius:2px;flex:none;}
.ln.dash{background:repeating-linear-gradient(to right,
  var(--text-secondary) 0 5px,transparent 5px 9px);}
.dot{width:11px;height:11px;border-radius:50%;flex:none;}
.note{font-size:12px;color:var(--text-muted);margin:9px 0 0;}
/* selection detail: pinned to the foot of the rail so it can grow without
   shoving the legend around. Laid out as key/value rows — the old one-line
   form wrapped into nonsense at this width. */
.selbox{flex:none;padding:13px 18px 16px;border-top:1px solid var(--border);
  background:var(--surface-1);max-height:44%;overflow-y:auto;}
#info{font-size:12.5px;color:var(--text-secondary);}
#info .nm{display:block;font-size:14px;font-weight:600;color:var(--text-primary);
  line-height:1.3;word-break:break-word;}
#info .lb{display:block;margin:1px 0 9px;color:var(--text-secondary);}
#info .kv{display:flex;justify-content:space-between;gap:12px;padding:2px 0;}
#info .kv span:last-child{color:var(--text-primary);font-variant-numeric:tabular-nums;
  text-align:right;white-space:nowrap;}
/* narrow viewports: the rail moves above the graph and reads as a header */
@media (max-width:820px){
  body{height:auto;overflow:visible;}
  .stage{flex-direction:column;}
  .canvas{height:60vh;min-height:340px;flex:none;}
  .side{order:-1;width:auto;overflow:visible;
    border-left:none;border-bottom:1px solid var(--border);}
  .railtop{overflow:visible;}
  .selbox{max-height:none;overflow:visible;}
}
/* covers the canvas exactly, so the toggle is always visible */
.tblwrap,.docwrap{position:absolute;inset:0;display:none;overflow:auto;
  padding:14px 18px 24px;background:var(--surface-1);}
.tblwrap.on,.docwrap.on{display:block}
.docwrap{padding:24px 28px 48px;}
.doc{max-width:74ch;}
.doc h2{margin:0 0 4px;font-size:19px;font-weight:700;letter-spacing:-0.01em;}
.doc .std{margin:0 0 22px;color:var(--text-secondary);font-size:13px;}
.doc h3{margin:24px 0 6px;font-size:14px;font-weight:600;}
.doc h3 .n{display:inline-block;min-width:20px;color:var(--text-muted);font-weight:700;}
.doc p{margin:0 0 8px;font-size:13.5px;color:var(--text-secondary);max-width:70ch;}
.doc b{color:var(--text-primary);}
.doc .fml{display:block;margin:9px 0;padding:9px 12px;border-radius:7px;
  background:var(--surface-2);border:1px solid var(--border);
  font-family:ui-monospace,SFMono-Regular,Consolas,monospace;font-size:12.5px;
  color:var(--text-primary);white-space:pre;overflow-x:auto;}
.doc .why{margin:22px 0 0;padding:12px 14px;border-radius:8px;
  background:var(--surface-2);border:1px solid var(--border);}
.doc .why p{margin:0;}
.doc table{margin:10px 0 0;}
.doc td.g{font-weight:600;color:var(--text-primary);}
/* the discredited column: muted, never red — red already means a pole of the
   scale here, and a red "+1.00" would read as an encoding, not "ignore this" */
.doc .raw{color:var(--text-muted);text-decoration:line-through;
  text-decoration-color:var(--text-muted);text-decoration-thickness:1px;}
.doc th.raw{text-decoration:none;}
table{border-collapse:collapse;font-size:13px;min-width:640px;}
th,td{text-align:left;padding:5px 12px 5px 0;border-bottom:1px solid var(--border);white-space:nowrap;}
th{color:var(--text-secondary);font-weight:600;}
td.num{text-align:right;font-variant-numeric:tabular-nums;padding-right:18px;}
.tag{display:inline-flex;align-items:center;gap:6px;}

/* Evidence tooltip — same anatomy and styling as the source graphs (PMID chip,
   score/year in muted brackets, rule between sentences), themed for dark. */
div.vis-tooltip{max-width:480px!important;white-space:normal!important;
  background:var(--tip-bg)!important;color:var(--tip-fg)!important;
  border:1px solid var(--tip-bd)!important;border-radius:8px!important;
  padding:8px 10px!important;box-shadow:0 4px 16px rgba(0,0,0,.35)!important;
  font:12px/1.45 Segoe UI,Arial,sans-serif!important}
@media (max-width:600px){div.vis-tooltip{max-width:88vw!important}}
div.vis-tooltip .eth{font-size:13px;margin-bottom:6px}
div.vis-tooltip .esrc{margin-top:8px;font-size:12px;font-weight:600;color:var(--tip-head)}
div.vis-tooltip .stip{padding:3px 0;border-top:1px solid var(--tip-rule)}
div.vis-tooltip .mut{color:var(--tip-mut);font-size:12px}
div.vis-tooltip .more{margin-top:5px;color:var(--tip-more);font-style:italic}
div.vis-tooltip .pm{display:inline-block;background:var(--tip-pm-bg);color:var(--tip-pm-fg);
  border-radius:4px;padding:0 5px;margin-right:5px;font-weight:600;font-size:11px;
  text-decoration:none}
div.vis-tooltip .pm:hover{background:var(--tip-pm-bg-hover);text-decoration:underline}
/* forced dark ink: the source leaves this at `inherit`, which is fine on its
   always-white tooltip but would put pale text on yellow in dark mode */
div.vis-tooltip mark{background:var(--tip-mark);color:#1a1a1a;border-radius:2px;padding:0 1px}
</style></head><body>
<div class="stage">
  <div class="canvas">
    <div id="net"></div>
    <div class="tblwrap" id="tblwrap"><table id="tbl"></table></div>
    <div class="docwrap" id="docwrap"><div class="doc" id="doc"></div></div>
  </div>
  <aside class="side">
   <div class="railtop">
    <div>
      <h1>__TITLE__</h1>
      <p class="sub">Pairs need ≥ __MINPUBS__ publications.</p>
      <p class="sub">Hover a connector for its sentences.</p>
    </div>
    <div>
      <div class="sliders">
        <div class="srow">
          <label class="slab" for="conf">Min score</label>
          <input id="conf" type="range" min="0" max="__MAXSTEP__" step="1" value="__DEFSTEP__"
                 aria-label="Minimum relationship score">
          <span class="sval" id="confval">__DEFCONF__</span>
        </div>
        <div class="srow">
          <label class="slab" for="mincluster">Min cluster</label>
          <input id="mincluster" type="range" min="2" max="12" step="1" value="2"
                 aria-label="Minimum connected cluster size">
          <span class="sval" id="mcval">2</span>
        </div>
        <p class="hint" id="mchint" hidden></p>
        <button class="resetbtn" id="rbtn"
                title="Clears the saved score, cluster size, gene focus, text search and legend state">Reset all (score ≥ __DEFCONF__)</button>
      </div>
      <p class="hgnc">Gene symbols and names in the sentences are normalized to HGNC symbols displayed
      as nodes. <strong>Use the <a href="https://www.genenames.org/" target="_blank" rel="noopener">HGNC
      website</a> to look up HGNC symbols when the nodes and sentences do not use the same
      terminology.</strong></p>
      <div style="margin-top:12px">
        <label class="slabel" for="genefilter">Focus on gene</label>
        <div class="pair">
          <input id="genefilter" type="search" autocomplete="off" spellcheck="false"
                 list="genelist" placeholder="gene symbol">
          <select id="hops" aria-label="Neighbourhood depth">
            <option value="1">1 hop</option>
            <option value="2">2 hops</option>
            <option value="3">3 hops</option>
          </select>
        </div>
        <datalist id="genelist"></datalist>
        <p class="hint" id="ghint">Exact name, else prefix match.</p>
      </div>
      <div style="margin-top:12px">
        <label class="slabel" for="textfilter">Match text in sentence</label>
        <input id="textfilter" type="search" autocomplete="off" spellcheck="false"
               placeholder="e.g. phosphorylat or /inhibit(s|ed)?/">
        <p class="hint" id="thint">Substring, or <code>/regex/</code>.</p>
      </div>
      <div class="controls" style="margin-top:12px">
        <button id="tbtn" aria-expanded="false" aria-controls="tblwrap">Show table</button>
        <button id="dbtn" aria-expanded="false" aria-controls="docwrap">How balance works</button>
      </div>
      <p class="count" id="count"></p>
    </div>
    <details class="leg" id="legbal">
      <summary><span class="sect">Gene balance</span>
        <span class="mgrad" role="img" aria-label="red for __LAB1__ through grey to blue for __LAB2__"><i
          class="mcut" id="mcutL"></i><i class="mcut" id="mcutR"></i></span>
      </summary>
      <div class="body">
        <span class="hramp" role="img" aria-label="Diverging scale: pure __LAB1__ (red) at the left, balanced (grey) in the middle, pure __LAB2__ (blue) at the right"><i
          class="cut" id="cutL"></i><i class="cut" id="cutR"></i></span>
        <div class="hticks">
          <span><b>pure __LAB1__</b></span>
          <span>balanced</span>
          <span><b>pure __LAB2__</b></span>
        </div>
        <div class="sliders" style="margin-top:11px">
          <div class="srow balrow">
            <label class="slab" for="balhi" title="__LAB1__ end">__LAB1__ end</label>
            <input id="balhi" type="range" min="-1" max="1" step="0.05" value="1"
                   aria-label="Show genes up to this balance">
            <span class="sval" id="balhiv">+1.00</span>
          </div>
          <div class="srow balrow">
            <label class="slab" for="ballo" title="__LAB2__ end">__LAB2__ end</label>
            <input id="ballo" type="range" min="-1" max="1" step="0.05" value="-1"
                   aria-label="Show genes down to this balance">
            <span class="sval" id="ballov">−1.00</span>
          </div>
        </div>
        <p class="hint" id="balhint" hidden></p>
        <p class="note">Share of each gene's own graph — the bigger graph doesn't skew it.</p>
      </div>
    </details>
    <details class="leg" id="legrel">
      <summary><span class="sect">Relationship</span>
        <span class="mkeys" role="img" aria-label="red __LAB1__, blue __LAB2__, grey shared">
          <i style="background:var(--edge-sclc)"></i><i style="background:var(--edge-luad)"></i><i style="background:var(--edge-both)"></i>
        </span>
      </summary>
      <div class="body">
        <div class="skey">
          <span class="lg"><span class="ln" style="background:var(--edge-sclc)"></span>__LAB1__ pair</span>
          <span class="lg"><span class="ln" style="background:var(--edge-luad)"></span>__LAB2__ pair</span>
          <span class="lg"><span class="ln" style="background:var(--edge-both)"></span>shared pair</span>
          <span class="lg"><span class="ln dash"></span>negated</span>
        </div>
        <p class="note">Width = publications. Node size = partners.</p>
      </div>
    </details>
   </div>
   <div class="selbox">
     <p class="sect">Selection</p>
     <div id="info" aria-live="polite">Click a gene or a relationship for detail.</div>
   </div>
  </aside>
</div>
<script>
const DATA=__PAYLOAD__, RAMP=__RAMPJS__, STEPS=__STEPSJS__, TOT=__TOTALSJS__,
      MINPUBS=__MINPUBS__, DEFSTEP=__DEFSTEP__;
// the two sides, named at the prompt; every label below is built from these
const LAB1=__LAB1JS__, LAB2=__LAB2JS__;
const ENAME={sclc:LAB1+' only',luad:LAB2+' only',both:'Shared'};
// formula-safe form of a label: 'some disease' -> 'some_disease'
function slug(s){return s.trim().replace(/[^A-Za-z0-9]+/g,'_').replace(/^_|_$/g,'').toLowerCase()||'x';}
const NODE_BY_ID={}; DATA.nodes.forEach(n=>{NODE_BY_ID[n.id]=n;});
function cvar(n){return getComputedStyle(document.documentElement).getPropertyValue(n).trim();}
function isDark(){const t=document.documentElement.getAttribute('data-theme');
  if(t)return t==='dark';return matchMedia('(prefers-color-scheme: dark)').matches;}
function ramp(){return RAMP[isDark()?'dark':'light'];}
function stepIdx(){const v=parseInt(document.getElementById('conf').value);
  return isNaN(v)?DEFSTEP:Math.max(0,Math.min(STEPS.length-1,v));}
// balance at the active step: +1 pure LAB1 .. 0 balanced .. -1 pure LAB2
function balance(n,i){const ss=n.sd[i]/TOT.s[i], ls=n.ld[i]/TOT.l[i];
  return (ss+ls)?(ss-ls)/(ss+ls):0;}
function stepOf(i){const n=ramp().fill.length;
  return Math.max(0,Math.min(n-1,Math.round((1-i)/2*(n-1))));}
function label(i){return i>=0.999?'pure '+LAB1:i>=0.6?'strongly '+LAB1+'-leaning':i>=0.2?LAB1+'-leaning':
  i>-0.2?'balanced':i>-0.6?LAB2+'-leaning':i>-0.999?'strongly '+LAB2+'-leaning':'pure '+LAB2;}
function esc(s){return String(s==null?'':s).replace(/[&<>"]/g,
  c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));}
function pmA(p){return '<a class=pm target=_blank rel=noopener href="https://www.ncbi.nlm.nih.gov/pmc/articles/'+
  encodeURIComponent(p)+'/">'+esc(p)+'</a>';}
// --- text search, same contract as the source graphs -----------------------
// plain query = case-insensitive substring; /.../flags = regex (bad regex falls
// back to a literal match, as in the source)
function reEsc(s){return s.replace(/[.*+?^${}()|[\]\\]/g,'\\$&');}
function activeText(){return (document.getElementById('textfilter').value||'').trim();}
function textMatcher(q){
  if(!q)return null;
  const m=/^\/(.*)\/([a-z]*)$/.exec(q);
  let src=null,flags='i';
  if(m){try{new RegExp(m[1],m[2]);src=m[1];
    flags=(m[2].indexOf('i')>=0?m[2]:m[2]+'i').replace(/g/g,'');}catch(err){src=null;}}
  if(src===null)src=reEsc(q);
  try{return {test:t=>new RegExp(src,flags).test(t||''),hlre:new RegExp(src,flags+'g')};}
  catch(err){return null;}
}
let TM=null;   // matcher in force for the current view
// --- palette range ---------------------------------------------------------
// A brush over the colour scale: keep only the genes whose balance falls inside
// it. The handles are independent sliders, so they can cross — read them as a
// set rather than trusting which is which.
function balRange(){
  let hi=parseFloat(document.getElementById('balhi').value);
  let lo=parseFloat(document.getElementById('ballo').value);
  if(isNaN(hi))hi=1;
  if(isNaN(lo))lo=-1;
  return [Math.min(lo,hi), Math.max(lo,hi)];
}
function fmtBal(v){return (v>=0?'+':'\u2212')+Math.abs(v).toFixed(2);}
// the ramp runs +1 at the left to -1 at the right, so index v sits at (1-v)/2
function paintRange(lo,hi){
  const l=(1-hi)/2*100, r=(1-lo)/2*100;
  [['cutL','cutR'],['mcutL','mcutR']].forEach(ids=>{
    const A=document.getElementById(ids[0]), B=document.getElementById(ids[1]);
    if(!A||!B)return;
    A.style.left='0%';   A.style.width=Math.max(0,l)+'%';
    B.style.left=r+'%';  B.style.width=Math.max(0,100-r)+'%';
  });
}

// --- gene focus ------------------------------------------------------------
function activeGene(){return (document.getElementById('genefilter').value||'').trim();}
// resolve like the source graphs: exact name first, else the first prefix match
function findGene(q){
  q=q.toLowerCase();
  const ids=DATA.nodes.map(n=>n.id);
  return ids.find(id=>id.toLowerCase()===q) ||
         ids.find(id=>id.toLowerCase().indexOf(q)===0) || null;
}
// --- cluster size ----------------------------------------------------------
function minCluster(){
  const v=parseInt(document.getElementById('mincluster').value);
  return isNaN(v)?2:v;
}
// Drop connected groups smaller than `mc`. Components are found over the edges
// still standing, so the sizes are the ones on screen. At mc=2 this is a no-op:
// every component built from edges already has two nodes.
function dropSmallClusters(view,mc){
  const adj={};
  view.forEach(o=>{(adj[o.f]=adj[o.f]||[]).push(o.t);(adj[o.t]=adj[o.t]||[]).push(o.f);});
  const comp={};let cid=0;
  for(const n in adj){
    if(comp[n]!==undefined)continue;
    const stack=[n];comp[n]=cid;
    while(stack.length){const x=stack.pop();
      (adj[x]||[]).forEach(y=>{if(comp[y]===undefined){comp[y]=cid;stack.push(y);}});}
    cid++;
  }
  const size={};
  for(const n in comp)size[comp[n]]=(size[comp[n]]||0)+1;
  return view.filter(o=>size[comp[o.f]]>=mc);
}
// breadth-first over the edges that survived every other filter, so the
// neighbourhood is the one actually on screen — not the whole graph's
function focusOn(view,seed,hops){
  const adj={};
  view.forEach(o=>{(adj[o.f]=adj[o.f]||[]).push(o.t);(adj[o.t]=adj[o.t]||[]).push(o.f);});
  const seen=new Set([seed]);
  let front=[seed];
  for(let h=0;h<hops;h++){
    const next=[];
    front.forEach(x=>(adj[x]||[]).forEach(y=>{if(!seen.has(y)){seen.add(y);next.push(y);}}));
    front=next;
  }
  return view.filter(o=>seen.has(o.f)&&seen.has(o.t));
}
// mark hits on the RAW text, escaping each piece, so a query containing < or &
// still highlights without injecting markup
function hl(t){
  t=t||'';
  if(!TM)return esc(t);
  const re=TM.hlre;re.lastIndex=0;
  let out='',last=0,m;
  while((m=re.exec(t))!==null){
    if(!m[0].length){re.lastIndex++;continue;}   // zero-length match: never loop
    out+=esc(t.slice(last,m.index))+'<mark>'+esc(m[0])+'</mark>';
    last=m.index+m[0].length;
  }
  return out+esc(t.slice(last));
}
// evidence block for one source graph, mirroring the source files' tooltip.
// `all` = sentences at the current score; `shown` = those matching the search.
function srcBlock(name,all,shown){
  if(!all.length)return '';
  const np=new Set(shown.map(s=>s.p)).size, lim=10;
  const cats=[...new Set(all.map(s=>s.g))].join(', ');
  const cnt=TM?(shown.length+' of '+all.length+' sentences match')
             :(all.length+' sentence'+(all.length===1?'':'s'));
  let h='<div class=esrc>'+name+' &middot; '+cnt+' &middot; '+np+' PMID'+(np===1?'':'s')+
    ' &middot; '+esc(cats)+'</div>';
  shown.slice(0,lim).forEach(s=>{h+='<div class=stip>'+pmA(s.p)+' <span class=mut>['+
    s.c.toFixed(3)+(s.y?(' &middot; '+s.y):'')+']</span> '+hl(s.t)+'</div>';});
  if(shown.length>lim)h+='<div class=more>+'+(shown.length-lim)+' more</div>';
  return h;
}
function edgeTip(o){
  const d=document.createElement('div');
  d.innerHTML='<div class=eth><b>'+esc(o.f)+' &ndash; '+esc(o.t)+'</b> ('+ENAME[o.c]+')</div>'+
    srcBlock(LAB1,o.sv,o.svm)+srcBlock(LAB2,o.lv,o.lvm);
  return d;
}
let network=null, VIEW=[];
// Everything below the score slider is derived, so it all recomputes here.
function build(){
  const R=ramp(), si=stepIdx(), T=STEPS[si];
  document.getElementById('confval').textContent=T.toFixed(2);
  TM=textMatcher(activeText());
  const EC={sclc:cvar('--edge-sclc'),luad:cvar('--edge-luad'),both:cvar('--edge-both')};
  // balance is per-gene and per-step, so resolve it once rather than per edge
  const BAL={};
  DATA.nodes.forEach(n=>{BAL[n.id]=balance(n,si);});
  const [blo,bhi]=balRange();
  const full=(blo<=-1&&bhi>=1);
  document.getElementById('balhiv').textContent=fmtBal(bhi);
  document.getElementById('ballov').textContent=fmtBal(blo);
  paintRange(blo,bhi);
  const inBal=id=>{const b=BAL[id];return b>=blo-1e-9&&b<=bhi+1e-9;};
  VIEW=[];
  DATA.edges.forEach(e=>{
    // both endpoints must be in the brushed range — an edge with one end
    // outside would draw a node the palette says is not being shown
    if(!full&&(!inBal(e.f)||!inBal(e.t)))return;
    const sv=e.s.filter(x=>x.c>=T), lv=e.l.filter(x=>x.c>=T);
    if(!sv.length&&!lv.length)return;
    const sp=new Set(sv.map(x=>x.p)).size, lp=new Set(lv.map(x=>x.p)).size;
    // Provenance, width and balance stay measured on the FULL evidence at this
    // score. The text search picks which pairs to show; it must not restate what
    // they are — restating a shared pair as one-sided because the query happened
    // to hit only one side's sentences would be a false claim about the biology.
    const c=(sv.length&&lv.length)?'both':(sv.length?'sclc':'luad');
    if(c!=='both'&&sp<MINPUBS&&lp<MINPUBS)return;   // shared pairs always survive
    const svm=TM?sv.filter(x=>TM.test(x.t)):sv, lvm=TM?lv.filter(x=>TM.test(x.t)):lv;
    if(TM&&!svm.length&&!lvm.length)return;
    VIEW.push({f:e.f,t:e.t,c:c,sv:sv,lv:lv,svm:svm,lvm:lvm,sp:sp,lp:lp,
      ss:sv.length,ls:lv.length,neg:sv.concat(lv).some(x=>x.g==='negated')});
  });
  // gene focus runs last, over whatever the other filters left standing
  const gq=activeGene(), gh=document.getElementById('ghint');
  if(gq){
    const seed=findGene(gq), hops=parseInt(document.getElementById('hops').value)||1;
    if(!seed){VIEW=[];gh.className='hint bad';gh.textContent='No gene named or starting with “'+gq+'”.';}
    else{
      VIEW=focusOn(VIEW,seed,hops);
      gh.className='hint';
      gh.textContent=VIEW.length?('Showing '+seed+' + '+hops+' hop'+(hops>1?'s':'')+'.')
        :(seed+' has no relationships in the current view.');
    }
  }else{gh.className='hint';
    gh.textContent='Exact name, else prefix match.';}
  // Cluster pruning is skipped while a gene is focused, as in the source graphs:
  // the focus already states which neighbourhood you want, and a size rule laid
  // over it could silently delete the very gene you asked for.
  const bh=document.getElementById('balhint');
  if(full)bh.hidden=true;
  else{bh.hidden=false;
    bh.textContent='Showing '+fmtBal(blo)+' to '+fmtBal(bhi)+' — '+label(blo)+' to '+label(bhi)+'.';}
  const mc=minCluster(), mch=document.getElementById('mchint');
  document.getElementById('mcval').textContent=mc;
  // the hint only earns its line when the slider is doing something
  if(gq){mch.hidden=false;mch.textContent='Cluster size not applied while a gene is focused.';}
  else if(mc>2){VIEW=dropSmallClusters(VIEW,mc);
    mch.hidden=false;mch.textContent='Hiding groups smaller than '+mc+' genes.';}
  else mch.hidden=true;
  const keep=new Set(); VIEW.forEach(o=>{keep.add(o.f);keep.add(o.t);});
  const nodes=[...keep].sort().map(id=>{
    const n=NODE_BY_ID[id], b=balance(n,si), s=stepOf(b), deg=Math.max(n.sd[si],n.ld[si]);
    return {id:id, label:id, size:8+Math.sqrt(Math.max(deg,1))*3.2,
      color:{background:R.fill[s],border:R.border[s],
             highlight:{background:R.fill[s],border:cvar('--text-primary')}},
      borderWidth:2,
      font:{color:cvar('--text-primary'),size:Math.max(14,Math.min(14+deg*0.2,26)),
            strokeWidth:5,strokeColor:cvar('--net-bg'),vadjust:-1},
      title:id+' — '+label(b)+' (balance '+b.toFixed(2)+')\npartners: '+LAB1+' '+
        n.sd[si]+', '+LAB2+' '+n.ld[si]};
  });
  const edges=VIEW.map((o,i)=>({id:i,from:o.f,to:o.t,
    width:Math.min(1.4+Math.max(o.sp,o.lp)*0.75,9),
    color:{color:EC[o.c],opacity:o.c==='both'?0.9:0.62},
    dashes:o.neg?[5,4]:false,
    title:edgeTip(o)}));
  const data={nodes:new vis.DataSet(nodes),edges:new vis.DataSet(edges)};
  const options={layout:{improvedLayout:false},
    physics:{stabilization:{iterations:400},
      barnesHut:{gravitationalConstant:-7000,centralGravity:0.62,springLength:95,
                 springConstant:0.035,damping:0.4,avoidOverlap:0.28}},
    interaction:{hover:true,tooltipDelay:120},
    nodes:{shape:'dot'},
    edges:{smooth:false,arrows:{to:{enabled:false}},hoverWidth:0,selectionWidth:0}};
  if(network)network.destroy();
  network=new vis.Network(document.getElementById('net'),data,options);
  network.on('stabilizationIterationsDone',()=>{
    network.setOptions({physics:false});network.fit({animation:false});network.redraw();});
  network.on('click',p=>{
    const info=document.getElementById('info');
    if(p.nodes.length){const n=NODE_BY_ID[p.nodes[0]], b=balance(n,si);
      info.innerHTML=nm(n.id)+lb(label(b))+kv('balance',b.toFixed(2))+
        kv(LAB1+' partners',n.sd[si])+kv(LAB2+' partners',n.ld[si]);}
    else if(p.edges.length){const o=VIEW[p.edges[0]];
      info.innerHTML=nm(o.f+'–'+o.t)+lb(ENAME[o.c]+(o.neg?' · has a negated relationship':''))+
        kv(LAB1,o.sp+' pubs / '+o.ss+' sent.')+kv(LAB2,o.lp+' pubs / '+o.ls+' sent.');}
    else info.textContent='Click a gene or a relationship for detail.';});
  document.getElementById('count').innerHTML='<b>'+nodes.length+'</b> genes · <b>'+
    edges.length+'</b> relationships'+(TM?' matching':'');
  const th=document.getElementById('thint');
  if(TM&&!VIEW.length){th.className='hint bad';th.textContent='No sentences match that query.';}
  else{th.className='hint';
    th.innerHTML='Substring, or <code>/regex/</code>.';}
  table(si);
  renderDoc(si);
}
// The explainer is generated from the live model, not written out as prose, so
// its worked numbers are always the ones behind the graph on screen.
// Exemplars are chosen from the data, not named here: a hardcoded gene list
// would be wrong for any other pair of graphs. Takes the best-connected genes
// and walks them across the scale, so the table always spans pole to pole.
function docGenes(si){
  const deg=n=>Math.max(n.sd[si],n.ld[si]);
  const pool=DATA.nodes.filter(n=>deg(n)>0).sort((a,b)=>deg(b)-deg(a)).slice(0,40);
  pool.sort((a,b)=>balance(b,si)-balance(a,si));
  if(pool.length<=7)return pool;
  const out=[];
  for(let i=0;i<7;i++)out.push(pool[Math.round(i*(pool.length-1)/6)]);
  return out;
}
// first value is a representative point for the swatch, not the bucket edge —
// keyed to the edge, "balanced" would draw as pale blue instead of neutral.
// The name comes from label(), so the buckets follow the prompt's labels.
function buckets(){return [
  [1,'index = +1 — no '+esc(LAB2)+' partners at all'],
  [0.8,'+0.6 to +1.0'],[0.4,'+0.2 to +0.6'],
  [0,'−0.2 to +0.2 — about equally central to both'],
  [-0.4,'−0.6 to −0.2'],[-0.8,'−1.0 to −0.6'],
  [-1,'index = −1 — no '+esc(LAB1)+' partners at all']];}
function fmtInt(n){return n.toString().replace(/\B(?=(\d{3})+(?!\d))/g,',');}
function renderDoc(si){
  const R=ramp(), T=STEPS[si], st=TOT.s[si], lt=TOT.l[si];
  const ex=docGenes(si).map(n=>{
    const sd=n.sd[si], ld=n.ld[si];
    return {id:n.id,sd:sd,ld:ld,ss:sd/st,ls:ld/lt,
            raw:(sd+ld)?(sd-ld)/(sd+ld):0, idx:balance(n,si)};
  });
  // whichever graph is larger does the dragging — never assume it is either one
  const bigL=lt>=st?LAB2:LAB1, smallL=lt>=st?LAB1:LAB2;
  const ratio=(Math.max(st,lt)/Math.max(1,Math.min(st,lt))).toFixed(2);
  const s1=slug(LAB1), s2=slug(LAB2);
  // illustrate step 1 with a gene that actually has both sides
  const eg=ex.filter(e=>e.sd>0&&e.ld>0).sort((a,b)=>Math.max(b.sd,b.ld)-Math.max(a.sd,a.ld))[0]||ex[0];
  let h='<h2>How the balance scale works</h2>'+
    '<p class=std>Node colour is a continuous '+esc(LAB1)+'↔'+esc(LAB2)+' axis. Every number below is '+
    'computed at the score you currently have selected (<b>≥ '+T.toFixed(2)+'</b>), so it describes the '+
    'graph on screen.</p>';

  h+='<h3><span class=n>1.</span>Count partners</h3><p>For each gene, count how many distinct '+
     'genes it is connected to, in each graph separately. At ≥ '+T.toFixed(2)+', '+
     '<b>'+esc(eg.id)+'</b> has <b>'+eg.sd+'</b> partners in '+esc(LAB1)+' and <b>'+eg.ld+'</b> in '+
     esc(LAB2)+'.</p>';

  h+='<h3><span class=n>2.</span>Divide by each graph’s own size <span class=mut>— the crux</span></h3>'+
     '<p>The two graphs are not the same size: '+esc(bigL)+'’s is <b>'+ratio+'×</b> larger at this '+
     'score ('+fmtInt(Math.max(st,lt)/2)+' pairs vs '+fmtInt(Math.min(st,lt)/2)+'). Raw counts are '+
     'therefore not comparable — 20 partners means far more in '+esc(smallL)+' than in '+esc(bigL)+
     '. So each degree becomes a <b>share of its own graph’s total connectivity</b>:</p>'+
     '<span class=fml>'+s1+'_share = '+s1+'_degree / '+fmtInt(st)+'\n'+
     s2+'_share = '+s2+'_degree / '+fmtInt(lt)+
     '\n\n(total = 2 × pairs — every pair has two endpoints)</span>'+
     '<p>This asks <i>“what fraction of '+esc(LAB1)+'’s wiring runs through this gene?”</i> against the '+
     'same question for '+esc(LAB2)+' — which is fair to both graphs regardless of how much has been '+
     'published.</p>';

  h+='<h3><span class=n>3.</span>Contrast the two shares</h3>'+
     '<span class=fml>index = ('+s1+'_share − '+s2+'_share) / ('+s1+'_share + '+s2+'_share)</span>'+
     '<p>Bounded to <b>−1 … +1</b>. Dividing by the sum is what makes this a <b>proportion</b> rather '+
     'than a magnitude: a gene with 3-vs-1 partners and one with 30-vs-10 read as the same lean. '+
     'Magnitude is not lost — it is carried separately, by node size.</p>';

  h+='<h3><span class=n>4.</span>Name the number</h3><p>Colour uses the continuous value; the words are '+
     'just bucketed for readability.</p><table><tbody>';
  buckets().forEach(b=>{const s=stepOf(b[0]);
    h+='<tr><td><span class=tag><span class=dot style="background:'+R.fill[s]+';border:1px solid '+
      R.border[s]+'"></span><b>'+esc(label(b[0]))+'</b></span></td>'+
      '<td style="padding-left:14px">'+b[1]+'</td></tr>';});
  h+='</tbody></table>';

  h+='<h3>Worked examples at ≥ '+T.toFixed(2)+'</h3><table><thead><tr><th>Gene</th>'+
     '<th class=num>'+esc(LAB1)+'</th><th class=num>'+esc(LAB2)+'</th>'+
     '<th class=num>'+s1+'_share</th><th class=num>'+s2+'_share</th>'+
     '<th class="num raw">raw</th><th class=num>index</th><th>reads as</th>'+
     '</tr></thead><tbody>';
  ex.forEach(e=>{const s=stepOf(e.idx);
    h+='<tr><td class=g>'+e.id+'</td><td class=num>'+e.sd+'</td><td class=num>'+e.ld+'</td>'+
      '<td class=num>'+e.ss.toFixed(5)+'</td><td class=num>'+e.ls.toFixed(5)+'</td>'+
      '<td class="num raw">'+(e.raw>=0?'+':'')+e.raw.toFixed(2)+'</td>'+
      '<td class=num><b>'+(e.idx>=0?'+':'')+e.idx.toFixed(2)+'</b></td>'+
      '<td><span class=tag><span class=dot style="background:'+R.fill[s]+';border:1px solid '+
      R.border[s]+'"></span>'+label(e.idx)+'</span></td></tr>';});
  h+='</tbody></table>';

  // "wrong" is measured, not asserted: where the two columns disagree on the label
  const misread=ex.filter(e=>label(e.raw)!==label(e.idx));
  h+='<div class=why><p><b>Why bother?</b> The struck-through <b>raw</b> column is the same contrast '+
     'computed on unnormalised counts. It disagrees with the normalised index for <b>'+misread.length+
     ' of these '+ex.length+' genes</b>'+
     (misread.length?' — '+misread.map(e=>esc(e.id)).join(', '):'')+
     '. A raw difference is biased toward whichever graph has published more — here <b>'+esc(bigL)+
     '</b>, whose graph is '+ratio+'× larger — so a gene central to both gets dragged toward it and '+
     'reads as though it belonged there. Normalising asks the same question of both graphs regardless '+
     'of how much has been written about either.</p></div>';

  h+='<h3>Two things to know</h3>'+
     '<p>• The index is recomputed at <b>every score step</b> — degrees and totals both move, so a '+
     'gene’s label can shift as you drag the slider.</p>'+
     '<p>• It is always measured over the <b>full</b> graphs, never the subgraph on screen. Focusing a '+
     'gene or matching sentence text changes what is drawn, but never what a gene <i>is</i>.</p>';
  document.getElementById('doc').innerHTML=h;
}
function nm(t){return '<span class=nm>'+esc(t)+'</span>';}
function lb(t){return '<span class=lb>'+esc(t)+'</span>';}
function kv(k,v){return '<div class=kv><span>'+esc(k)+'</span><span>'+esc(v)+'</span></div>';}
function table(si){
  const R=ramp(), shown=new Set();
  VIEW.forEach(o=>{shown.add(o.f);shown.add(o.t);});
  const rows=[...shown].map(id=>{const n=NODE_BY_ID[id];return {n:n,b:balance(n,si)};})
    .sort((a,b)=>b.b-a.b||a.n.id.localeCompare(b.n.id));
  let h='<thead><tr><th>Gene</th><th>Balance</th><th class=num>index</th>'+
        '<th class=num>'+esc(LAB1)+' partners</th><th class=num>'+esc(LAB2)+
        ' partners</th></tr></thead><tbody>';
  rows.forEach(r=>{const s=stepOf(r.b);
    h+='<tr><td>'+esc(r.n.id)+'</td><td><span class=tag><span class=dot style="background:'+R.fill[s]+
      ';border:1px solid '+R.border[s]+'"></span>'+label(r.b)+'</span></td><td class=num>'+
      r.b.toFixed(2)+'</td><td class=num>'+r.n.sd[si]+'</td><td class=num>'+r.n.ld[si]+'</td></tr>';});
  document.getElementById('tbl').innerHTML=h+'</tbody>';
}
let tdeb=null;   // rebuilding restarts physics, so don't do it on every keystroke
document.getElementById('textfilter').addEventListener('input',()=>{
  clearTimeout(tdeb);tdeb=setTimeout(()=>{saveFilters();build();},280);});
let gdeb=null;
document.getElementById('genefilter').addEventListener('input',()=>{
  clearTimeout(gdeb);gdeb=setTimeout(()=>{saveFilters();build();},280);});
document.getElementById('hops').addEventListener('change',()=>{saveFilters();build();});
(()=>{const dl=document.getElementById('genelist');
  DATA.nodes.forEach(n=>{const o=document.createElement('option');o.value=n.id;dl.appendChild(o);});
  // name a gene that is actually in these graphs rather than a hardcoded example
  const d=n=>Math.max(n.sd[DEFSTEP],n.ld[DEFSTEP]);
  const top=DATA.nodes.slice().sort((a,b)=>d(b)-d(a))[0];
  if(top)document.getElementById('genefilter').placeholder='e.g. '+top.id;})();
document.getElementById('mincluster').addEventListener('input',()=>{
  document.getElementById('mcval').textContent=minCluster();});
document.getElementById('mincluster').addEventListener('change',()=>{saveFilters();build();});
['balhi','ballo'].forEach(id=>{
  const el=document.getElementById(id);
  el.addEventListener('input',()=>{const r=balRange();
    document.getElementById('balhiv').textContent=fmtBal(r[1]);
    document.getElementById('ballov').textContent=fmtBal(r[0]);
    paintRange(r[0],r[1]);});
  el.addEventListener('change',()=>{saveFilters();build();});
});
document.getElementById('conf').addEventListener('input',()=>{
  document.getElementById('confval').textContent=STEPS[stepIdx()].toFixed(2);});
document.getElementById('conf').addEventListener('change',()=>{saveScore();build();});
// the two overlays share the canvas, so opening one closes the other
let PANEL=null;
function showPanel(which){
  PANEL=(PANEL===which)?null:which;
  const t=document.getElementById('tbtn'), d=document.getElementById('dbtn');
  document.getElementById('tblwrap').classList.toggle('on',PANEL==='tbl');
  document.getElementById('docwrap').classList.toggle('on',PANEL==='doc');
  t.textContent=PANEL==='tbl'?'Hide table':'Show table';
  d.textContent=PANEL==='doc'?'Hide explainer':'How balance works';
  t.setAttribute('aria-expanded',PANEL==='tbl');
  d.setAttribute('aria-expanded',PANEL==='doc');
}
// Remember which legends the reader left open. Wrapped in try/catch throughout:
// this page is opened over file://, where localStorage can be unavailable or
// throw outright — a legend preference must never take the graph down with it.
const LEGENDS=['legbal','legrel'], LS_KEY='diff-two.legends',
      LS_SCORE='diff-two.score', LS_FILTERS='diff-two.filters';
// The gene focus and text query come back too. A restored filter is never
// silent — the query sits in its own input and the hint below it says what is
// being shown, so a narrowed graph always explains itself.
function saveFilters(){
  if(QUIET)return;
  try{localStorage.setItem(LS_FILTERS,JSON.stringify({
    gene:document.getElementById('genefilter').value,
    hops:document.getElementById('hops').value,
    text:document.getElementById('textfilter').value,
    cluster:document.getElementById('mincluster').value,
    balhi:document.getElementById('balhi').value,
    ballo:document.getElementById('ballo').value}));}catch(err){}
}
function restoreFilters(){
  try{const raw=localStorage.getItem(LS_FILTERS);
    if(!raw)return;                              // first visit: empty filters
    const st=JSON.parse(raw), hops=document.getElementById('hops');
    if(typeof st.gene==='string')document.getElementById('genefilter').value=st.gene;
    if(typeof st.text==='string')document.getElementById('textfilter').value=st.text;
    // only accept a hops value the select actually offers
    if(typeof st.hops==='string'&&[...hops.options].some(o=>o.value===st.hops))
      hops.value=st.hops;
    // clamp: the slider's range may have changed since this was written
    if(typeof st.cluster==='string'){
      const mc=document.getElementById('mincluster'), v=parseInt(st.cluster);
      if(!isNaN(v)){mc.value=String(Math.max(+mc.min,Math.min(+mc.max,v)));
        document.getElementById('mcval').textContent=mc.value;}
    }
    ['balhi','ballo'].forEach(id=>{
      if(typeof st[id]!=='string')return;
      const el=document.getElementById(id), v=parseFloat(st[id]);
      if(!isNaN(v))el.value=String(Math.max(-1,Math.min(1,v)));
    });
  }catch(err){}
}
// The score is stored as the value (0.97), never the slider index: an index
// would quietly resolve to a different threshold if SCORE_STEPS is ever edited.
// Restored to the nearest step, the same way DEFAULT_SCORE is.
function saveScore(){
  if(QUIET)return;
  try{localStorage.setItem(LS_SCORE,String(STEPS[stepIdx()]));}catch(err){}
}
function restoreScore(){
  try{const raw=localStorage.getItem(LS_SCORE);
    if(raw===null)return;                        // first visit: keep DEFAULT_SCORE
    const v=parseFloat(raw);
    if(!isFinite(v))return;
    let best=0;
    for(let i=1;i<STEPS.length;i++){
      if(Math.abs(STEPS[i]-v)<Math.abs(STEPS[best]-v))best=i;}
    document.getElementById('conf').value=best;
    document.getElementById('confval').textContent=STEPS[best].toFixed(2);
  }catch(err){}
}
// Reset writes the defaults back into the DOM, which fires the same events the
// save handlers listen on — QUIET holds them off so the wipe isn't immediately
// undone. It outlives the current task because <details> fires `toggle`
// asynchronously; a same-task flag would already be cleared by the time the
// legend's save handler ran.
let QUIET=0;
// Snapshot the markup's defaults before restoreLegends() can overwrite them, so
// Reset always agrees with a first visit. Hardcoding them here once drifted out
// of step with the markup when the default flipped.
const LEGEND_DEFAULTS={};
LEGENDS.forEach(id=>{LEGEND_DEFAULTS[id]=document.getElementById(id).open;});
function saveLegends(){
  if(QUIET)return;
  try{const st={};
    LEGENDS.forEach(id=>{st[id]=document.getElementById(id).open;});
    localStorage.setItem(LS_KEY,JSON.stringify(st));}catch(err){}
}
function resetAll(){
  QUIET++;
  document.getElementById('conf').value=DEFSTEP;
  document.getElementById('confval').textContent=STEPS[DEFSTEP].toFixed(2);
  document.getElementById('genefilter').value='';
  document.getElementById('textfilter').value='';
  document.getElementById('hops').value=document.getElementById('hops').options[0].value;
  document.getElementById('mincluster').value=document.getElementById('mincluster').min;
  document.getElementById('mcval').textContent=document.getElementById('mincluster').min;
  document.getElementById('balhi').value='1';
  document.getElementById('ballo').value='-1';
  LEGENDS.forEach(id=>{document.getElementById(id).open=LEGEND_DEFAULTS[id];});
  build();
  // queued after the toggle tasks, so nothing re-saves before the keys go
  setTimeout(()=>{
    QUIET--;
    try{[LS_KEY,LS_SCORE,LS_FILTERS].forEach(k=>localStorage.removeItem(k));}catch(err){}
  },0);
}
function restoreLegends(){
  try{const raw=localStorage.getItem(LS_KEY);
    if(!raw)return;                              // first visit: keep the markup's defaults
    const st=JSON.parse(raw);
    LEGENDS.forEach(id=>{if(typeof st[id]==='boolean')document.getElementById(id).open=st[id];});
  }catch(err){}
}
// before first paint: this script blocks rendering, so no flash of the defaults
restoreLegends();
restoreScore();    // must precede build() — build() reads the slider
restoreFilters();  // and the filter inputs
LEGENDS.forEach(id=>document.getElementById(id).addEventListener('toggle',saveLegends));
document.getElementById('rbtn').addEventListener('click',resetAll);
document.getElementById('tbtn').addEventListener('click',()=>showPanel('tbl'));
document.getElementById('dbtn').addEventListener('click',()=>showPanel('doc'));
document.addEventListener('keydown',e=>{if(e.key==='Escape'&&PANEL)showPanel(PANEL);});
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',build);
build();
</script></body></html>"""


def _esc(text):
    """Escape user text destined for HTML."""
    return (text.replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;").replace('"', "&quot;"))


def render(nodes, edges, totals, ramp, vis_runtime, title, lab1, lab2):
    payload = json.dumps({"nodes": nodes, "edges": edges}, separators=(",", ":"))
    gradient = {
        mode: ",".join(ramp[mode]["fill"][i] for i in range(0, RAMP_STEPS, 5))
        for mode in ("light", "dark")
    }
    html = TEMPLATE
    for token, value in (
        ("__VIS__", vis_runtime),
        ("__PAYLOAD__", payload),
        ("__RAMPJS__", json.dumps(ramp, separators=(",", ":"))),
        ("__STEPSJS__", json.dumps(SCORE_STEPS, separators=(",", ":"))),
        ("__TOTALSJS__", json.dumps(totals, separators=(",", ":"))),
        ("__GRADL__", gradient["light"]),
        ("__GRADD__", gradient["dark"]),
        ("__MINPUBS__", str(MIN_PUBS)),
        ("__MAXSTEP__", str(len(SCORE_STEPS) - 1)),
        ("__DEFSTEP__", str(DEFAULT_STEP)),
        ("__DEFCONF__", "%.2f" % SCORE_STEPS[DEFAULT_STEP]),
        # last: these are user text, so substituting them first would let a
        # value containing another token corrupt the following replacements
        ("__LAB1JS__", json.dumps(lab1)),
        ("__LAB2JS__", json.dumps(lab2)),
        ("__LAB1__", _esc(lab1)),
        ("__LAB2__", _esc(lab2)),
        ("__TITLE__", _esc(title)),
    ):
        html = html.replace(token, value)
    return html


# --------------------------------------------------------------------------
# Prompts
# --------------------------------------------------------------------------
#
# Every name comes from stdin; there are no defaults to fall back on. All of
# this lives behind main() rather than at module level so `import diff_two`
# stays silent — importing must never block on a prompt.

def _ask(label):
    """One line from stdin. EOF is an abort, not an empty answer."""
    try:
        return input(label).strip().strip('"\'')
    except EOFError:
        raise SystemExit("\naborted: %s needs a value" % label.strip().rstrip(":"))


def _resolve(name):
    return name if os.path.isabs(name) else os.path.join(HERE, name)


def ask_input_graph(label):
    """Prompt until given the name of a graph that actually exists."""
    while True:
        name = _ask(label)
        if not name:
            print("  a file name is required")
            continue
        path = _resolve(name)
        if not os.path.isfile(path):
            print("  no such file: %s" % path)
            continue
        return path


def ask_output_path(inputs):
    """Prompt until given a writable name that is not one of the inputs."""
    while True:
        name = _ask("output file name: ")
        if not name:
            print("  a file name is required")
            continue
        if not os.path.splitext(name)[1]:
            name += ".html"
        path = _resolve(name)
        # The inputs are read, not written. A typo here would destroy megabytes
        # of upstream pipeline output with no way back.
        if any(os.path.abspath(path) == os.path.abspath(p) for p in inputs):
            print("  that is an input graph — pick another name")
            continue
        return path


def ask_graph_name():
    while True:
        name = _ask("graph name: ")
        if name:
            return name
        print("  a graph name is required")


def ask_label(prompt):
    """What to call one side of the diff, everywhere in the finished page."""
    while True:
        name = _ask(prompt)
        if name:
            return name
        print("  a label is required")


def main():
    src1 = ask_input_graph("input file name 1: ")
    src2 = ask_input_graph("input file name 2: ")
    out = ask_output_path((src1, src2))
    title = ask_graph_name()
    lab1 = ask_label("label 1: ")
    lab2 = ask_label("label 2: ")

    pairs1 = collapse_to_pairs(extract_data(src1))
    pairs2 = collapse_to_pairs(extract_data(src2))
    nodes, edges, totals = build_model(pairs1, pairs2)

    html = render(nodes, edges, totals, build_ramp(),
                  extract_vis_runtime(src2), title, lab1, lab2)
    with io.open(out, "w", encoding="utf-8") as fh:
        fh.write(html)
    print("wrote %s (%.1f MB)" % (out, len(html) / 1048576.0))


if __name__ == "__main__":
    main()
