#!/usr/bin/env python3
"""
high_confidence_g.py -- extract the high-confidence Gene-context relation triples from
the scored RE output, summarize the threshold statistics, and draw the brain-cancer
gene-gene relationship graph.

This is the "G" (gene-only) variant of high_confidence.py, run as a SEPARATE STEP: the
same step-3 graph pipeline, differing in exactly one place -- the qualifying filter. Where
high_confidence.py applies the "G_D_C" filter (gene-gene IN a disease/chemical context),
this step applies the "G" filter (gene-gene, context-agnostic). It is a standalone step you
run in addition to -- or instead of -- the G_D_C step; its "_G" output names never clobber
the other's, so both can target the same data root.

  The "G" step -- the ONE logic difference from high_confidence.py:
  the DISEASE-or-CHEMICAL sentence-context requirement (G_D_C rule 5) is DROPPED.
  A triple QUALIFIES (the "G" filter) when ALL of these hold:
    1. score >= SCORE                                   (export default 0.8)
    2. annotated control:no   -- >=1 GENETIC endpoint has control == "no"
    3. NOT annotated control:yes -- no endpoint has control == "yes"
    4. NOT hgnc_symbol == "MKI67" on either endpoint
  (The G_D_C filter's 5th condition -- the sentence must also carry a DISEASE or
  CHEMICAL entity -- is intentionally NOT applied in this step, so the "G" universe is
  strictly larger: every high-confidence gene-gene relation, disease/chemical context or not.)

All inputs are read from the pipeline output tree (the writable run dir gpu.py produced),
which defaults to ``kaggle_working/`` next to this script; override with ``--data-root``.

Inputs  : <data-root>/TRIPLES/triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json  (scored + normalized)
          <data-root>/TRIPLES/triples.json                                          (disease facets for the graph)
          <data-root>/{DISEASE,CHEMICAL,databases,sentences}/...                     (normalization + year + corpus size)
Outputs : <data-root>/TRIPLES/high_confidence_G.json     qualifying triples at --score (default 0.8)
          <data-root>/summaries/high_confidence_G.html    gene-gene graph with a 0.5..0.99 in-browser score slider
          <root>/<pubmed_query>_G.html                    a copy of that graph, named after the PubMed query
                                                           (whitespace -> underscore); e.g. pancreatic_cancer_G.html

The output filenames all differ from those written by high_confidence.py (which uses the
"_G_D_C" JSON, "high_confidence.html" graph, and "<query>.html" copy) so the two scripts
can be run against the same data root without clobbering each other's outputs.

Run::  python high_confidence_g.py [--data-root kaggle_working] [--score 0.8] [--thresholds 0.8,0.95,0.99] [--no-graph]
"""
import argparse
import collections
import html
import json
import re
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# The pipeline outputs (gpu.py's writable run dir) live under kaggle_working/ by default;
# every input below is resolved beneath this data root. Override with --data-root.
DATA_ROOT = ROOT / "kaggle_working"

# module-level paths; (re)bound to DATA_ROOT by set_data_root() so --data-root can retarget them
OUT_DIR = XML_DIR = SENT_DIR = None
RE_FILE = BASE_FILE = PMC_YEARS = TARGET_FILE = DISEASE_FILE = DISEASE_AMBIG = None
JSON_OUT = GRAPH_OUT = None


def set_data_root(data_root):
    """Point every input/output path at `data_root` (the pipeline's output tree)."""
    global DATA_ROOT, OUT_DIR, XML_DIR, SENT_DIR, RE_FILE, BASE_FILE, PMC_YEARS
    global TARGET_FILE, DISEASE_FILE, DISEASE_AMBIG, JSON_OUT, GRAPH_OUT
    DATA_ROOT = Path(data_root).resolve()
    OUT_DIR = DATA_ROOT / "TRIPLES"
    XML_DIR = DATA_ROOT / "experimental_ner"   # input XML corpus (may be empty in the bundle)
    SENT_DIR = DATA_ROOT / "sentences"         # one JSON per source document (corpus-size fallback)
    RE_FILE = OUT_DIR / "triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json"
    BASE_FILE = OUT_DIR / "triples.json"
    PMC_YEARS = DATA_ROOT / "databases" / "pmc_years.json"
    TARGET_FILE = DATA_ROOT / "CHEMICAL" / "chemical_to_target.json"   # gene -> corpus chemicals (in_corpus_GENETIC flag)
    DISEASE_FILE = DATA_ROOT / "DISEASE" / "disease.json"            # surface -> MONDO label (single)
    DISEASE_AMBIG = DATA_ROOT / "DISEASE" / "disease_ambiguous.json"  # surface -> MONDO labels (list)
    # "_G" (gene-only) output names, distinct from high_confidence.py's "_G_D_C"/"high_confidence.html".
    JSON_OUT = OUT_DIR / "high_confidence_G.json"
    GRAPH_OUT = DATA_ROOT / "summaries" / "high_confidence_G.html"


set_data_root(DATA_ROOT)

GRAPH_BASE = 0.5           # graph universe = qualifying triples at this score (lowest in-browser slider stop)
VIS_URL = "https://unpkg.com/vis-network@9.1.9/standalone/umd/vis-network.min.js"
PCOLOR = {"positive": "#2e9e5b", "negated": "#e0533d", "speculated": "#d59a2e"}

# Disease-label synonyms: map surface/label variants (matched case-insensitively) onto a
# single canonical disease name so the same disease is not split across several filter
# entries. Add new aliases here. Keys must be lowercase.
DISEASE_ALIASES = {
    "luad": "lung adenocarcinoma",
    "lung adenocarcinoma": "lung adenocarcinoma",
    "primary lung adenocarcinoma": "lung adenocarcinoma",
    "or": "osimertinib resistance",   # NER false-positive DISEASE: "OR" = osimertinib resistance
}

# Non-disease surface tokens the NER mislabels as DISEASE entities: dropped from every
# disease-filtered view. Entries are lowercase; matched case-insensitively.
DISEASE_IGNORE = {
    "os",   # "overall survival" (a survival metric), not a disease
}

# Ambiguous surface abbreviations resolved to their single intended disease in THIS corpus,
# overriding the multi-label mapping in disease_ambiguous.json. Keys lowercase; the value
# replaces the ambiguous MONDO label list for that surface token.
DISEASE_SURFACE = {
    "gc": ["gastric cancer"],   # "GC" is gastric cancer here, not gonorrhea
}


def canon_disease(label):
    """Canonicalize a disease label via DISEASE_ALIASES (case-insensitive), e.g. both
    'LUAD' and 'lung adenocarcinoma' -> 'lung adenocarcinoma'. Labels in DISEASE_IGNORE
    return None (caller drops them). Unknown labels pass through unchanged (only surrounding
    whitespace stripped)."""
    if not isinstance(label, str):
        return label
    key = label.strip().lower()
    if key in DISEASE_IGNORE:
        return None
    return DISEASE_ALIASES.get(key, label.strip())


# ----- qualifying filter -----------------------------------------------------
def syms(e):
    v = e.get("hgnc_symbol")
    return (set(v) if isinstance(v, list) else {v}) if v is not None else set()


def has_mki67(t):
    return "MKI67" in (syms(t["subject"]) | syms(t["object"]))


def ctrl_no(t):
    return t["subject"].get("control") == "no" or t["object"].get("control") == "no"


def ctrl_yes(t):
    return t["subject"].get("control") == "yes" or t["object"].get("control") == "yes"


def qualifies(t, T):
    # gene-only filter: no DISEASE/CHEMICAL sentence-context requirement (that was the "_G_D_C" filter)
    sc = t.get("score")
    return (isinstance(sc, (int, float)) and sc >= T
            and ctrl_no(t) and not ctrl_yes(t) and not has_mki67(t))


def stats(d, T):
    sub = [t for t in d if isinstance(t.get("score"), (int, float)) and t["score"] >= T]
    f = [t for t in sub if qualifies(t, T)]
    return {"T": T, "triples": len(sub), "sentences": len({t.get("sentence") for t in sub}),
            "control_no": sum(1 for t in sub if ctrl_no(t)),
            "control_yes": sum(1 for t in sub if ctrl_yes(t)),
            "mki67": sum(1 for t in sub if has_mki67(t)),
            "filt": len(f), "filt_sentences": len({t.get("sentence") for t in f})}


# ----- brain-cancer gene-gene graph ------------------------------------------
def single(v):
    return v.strip() if isinstance(v, str) and v.strip() else None


def graph_payload(triples, base):
    """Gene-gene (single HGNC symbol, non-self) graph; per-sentence scores so the
    in-browser confidence toggle can re-filter to >=0.99."""
    try:
        years = json.loads(PMC_YEARS.read_text(encoding="utf-8"))
    except Exception:
        years = {}
    try:
        c2t = json.loads(TARGET_FILE.read_text(encoding="utf-8"))
    except Exception:
        c2t = {}
    # restrict to the corpus: genes flagged in_corpus_GENETIC, scaled by their # of corpus chemicals
    tgt, chems_by_gene, tsrc, tcat = {}, {}, {}, {}
    for g, v in c2t.items():
        if v.get("in_corpus_GENETIC"):
            chems = v.get("chemicals") or []
            tgt[g] = v.get("n_chemicals") or len(chems)
            chems_by_gene[g] = sorted({c.get("chebi_label") for c in chems if c.get("chebi_label")})
            vrs = [str(r) for c in chems for r in (c.get("via_roles") or [])]
            has_db = any(r.startswith("DGIdb") for r in vrs)
            has_ch = any(not r.startswith("DGIdb") for r in vrs)
            tsrc[g] = "dgidb" if (has_db and not has_ch) else ("chebi+dgidb" if has_db else "chebi")
            # colour class: green = approved anti-neoplastic, amber = approved (other), pink = other target
            if any(r.startswith("DGIdb-antineoplastic") for r in vrs):
                tcat[g] = "green"
            elif any(r.startswith("DGIdb-approved") for r in vrs):
                tcat[g] = "amber"
            else:
                tcat[g] = "other"
    dmap = {}
    for df in (DISEASE_FILE, DISEASE_AMBIG):
        try:
            for k, v in json.loads(df.read_text(encoding="utf-8")).items():
                ml = v.get("mondo_label")
                dmap[k] = ml if isinstance(ml, list) else [ml]
        except Exception:
            pass
    sent_dis = collections.defaultdict(set)   # sentence -> disease labels (MONDO where known, else the surface text)
    for bt in base:
        for r in ("subject", "object"):
            be = bt[r]
            if be.get("type") != "DISEASE":
                continue
            txt = (be.get("text") or "").strip()
            if txt.lower() in DISEASE_SURFACE:
                labels = list(DISEASE_SURFACE[txt.lower()])   # forced disambiguation for ambiguous surface tokens
            else:
                labels = [ml for ml in dmap.get(txt, []) if ml]   # MONDO label(s) for this surface, nulls dropped
            if not labels and txt:
                # no MONDO match (mondo_label null / surface unmatched, e.g. "PDAC"): fall back to the
                # DISEASE entity text itself so the sentence still carries a disease facet instead of
                # being dropped from every disease-filtered view.
                labels = [txt]
            for ml in labels:
                cl = canon_disease(ml)
                if cl:   # None -> suppressed (DISEASE_IGNORE); skip
                    sent_dis[bt.get("sentence")].add(cl)
    dir_sent = collections.defaultdict(set)           # (s,o,pol) -> sentences (for direction)
    pair_sent = collections.defaultdict(dict)          # pair -> {sentence: [maxscore, pmid]}
    node_sent = collections.defaultdict(dict)          # node -> {sentence: maxscore}
    for t in triples:
        s, o = single(t["subject"].get("hgnc_symbol")), single(t["object"].get("hgnc_symbol"))
        if not (s and o) or s == o:
            continue
        sc = float(t.get("score") or 0.0)
        pol = t.get("polarity", "positive")
        sent = t.get("sentence", "")
        pm = (t.get("pmid") or "?").replace(".grobid.tei", "")
        dir_sent[(s, o, pol)].add(sent)
        cur = pair_sent[frozenset((s, o))].get(sent)
        if cur is None or sc > cur[0]:
            pair_sent[frozenset((s, o))][sent] = [sc, pm]
        for nd in (s, o):
            if node_sent[nd].get(sent, -1) < sc:
                node_sent[nd][sent] = sc
    edges = []
    for pr, sd in pair_sent.items():
        cands = [(len(ss), f, to, pol) for (f, to, pol), ss in dir_sent.items() if frozenset((f, to)) == pr]
        cands.sort(reverse=True)
        _, ff, ft, fpol = cands[0]
        sents = [{"pmid": pm, "text": sent[:300], "sc": round(sc, 4), "yr": years.get(pm), "dis": sorted(sent_dis.get(sent, []))} for sent, (sc, pm) in sd.items()]
        sents.sort(key=lambda z: (-z["sc"], z["pmid"]))
        edges.append({"from": ff, "to": ft, "cat": fpol, "color": PCOLOR.get(fpol, "#888"), "sents": sents})
    nodes = [{"id": nd, "label": nd, "bg": "#cfe3ff", "border": "#2b6cb0",
              "sent95": len(sd), "sent99": sum(1 for v in sd.values() if v >= 0.99),
              "target": tgt.get(nd, 0), "chems": chems_by_gene.get(nd, []), "tsource": tsrc.get(nd, ""), "tcat": tcat.get(nd, "other")}
             for nd, sd in node_sent.items()]
    return {"nodes": nodes, "edges": edges}


def read_pubmed_query():
    """Pull the PubMed query out of summaries/pubmed_query.html (the step-1 publications
    summary), i.e. the value after 'Query (read strictly from <STDIN>):'. Looks under the
    data root first, then next to this script; returns '' if the file/line is absent."""
    for cand in (DATA_ROOT / "summaries" / "pubmed_query.html",
                 ROOT / "summaries" / "pubmed_query.html"):
        try:
            txt = cand.read_text(encoding="utf-8")
        except Exception:
            continue
        m = re.search(r"Query \(read strictly from &lt;STDIN&gt;\):\s*<code>(.*?)</code>", txt, re.S)
        if m:
            return html.unescape(m.group(1)).strip()
    return ""


def get_vis_lib():
    try:
        with urllib.request.urlopen(VIS_URL, timeout=30) as r:
            return r.read().decode("utf-8")
    except Exception as e:
        print(f"  [graph] could not fetch vis-network ({type(e).__name__}); HTML will use the CDN (needs internet to view)")
        return None


GRAPH_TEMPLATE = r"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
__LIBTAG__
<style>
 html,body{margin:0;height:100%;background:#ffffff;color:#1c2330;font:14px/1.5 Segoe UI,Arial,sans-serif}
 #net{position:absolute;top:0;left:0;right:0;bottom:0;background:#ffffff}
 #panel{position:absolute;top:12px;right:12px;z-index:5;background:rgba(255,255,255,.97);border:1px solid #cdd5e0;border-radius:10px;padding:14px 16px;max-width:320px;box-shadow:0 2px 12px rgba(0,0,0,.18);color:#1c2330}
 #panel h1{font-size:13px;margin:0 0 8px;color:#1c2330;font-variant:small-caps;letter-spacing:.4px}
 .row{margin:8px 0}
 input[type=range]{width:150px;max-width:100%;vertical-align:middle}
 .legend{display:flex;flex-wrap:wrap;align-items:center;gap:4px 12px}
 .legend b{display:inline-block;width:11px;height:11px;border-radius:2px;margin-right:5px;vertical-align:-1px;border:1px solid #999}
 .sw{display:inline-block;width:10px;height:10px;border-radius:2px;vertical-align:-1px}
 .mut{color:#5b6677;font-size:12px} b{color:#2b6cb0}
 #conf{width:190px;cursor:pointer}
 select,#search,#genefilter,#drugsearch{max-width:100%;background:#fff;border:1px solid #cdd5e0;color:#1c2330;border-radius:5px;padding:3px 6px;font-size:13px}
 #search,#genefilter,#drugsearch{width:200px}
 #catfilters label{display:block;cursor:pointer;white-space:nowrap;font-size:12px;margin:1px 0}
 #catfilters{border:1px solid #cdd5e0;border-radius:6px;padding:4px 6px;max-height:140px;overflow:auto}
 #zoom button{background:#eef2f7;color:#1c2330;border:1px solid #cdd5e0;border-radius:6px;padding:4px 10px;cursor:pointer;margin-right:6px;font-size:13px}
 #zoom button:hover{background:#dde4ee}
 .vis-tooltip{max-width:480px!important;white-space:normal!important;background:#fff!important;color:#1a1a1a!important;border:1px solid #999!important;border-radius:8px!important;padding:8px 10px!important;box-shadow:0 4px 16px rgba(0,0,0,.35)!important;font:12px/1.45 Segoe UI,Arial,sans-serif!important}
 .eth{font-size:13px;margin-bottom:6px} .stip{padding:3px 0;border-top:1px solid #e3e3e3}
 .pm{display:inline-block;background:#eef3fb;color:#2b6cb0;border-radius:4px;padding:0 5px;margin-right:5px;font-weight:600;font-size:11px;text-decoration:none}
 a.pm:hover{background:#d6e6fb;text-decoration:underline} .more{margin-top:5px;color:#888;font-style:italic}
 #info{max-height:240px;overflow:auto} #info .stip{border-top:1px solid #e3e3e3}
 #toggle{position:absolute;top:12px;left:12px;z-index:6;background:#fff;color:#1c2330;border:1px solid #cdd5e0;border-radius:8px;padding:5px 11px;font-size:18px;line-height:1;cursor:pointer;box-shadow:0 2px 8px rgba(0,0,0,.18)}
 #panel.collapsed{display:none}
 @media (max-width:700px){
  #panel{left:12px;right:12px;max-width:none;max-height:62vh;overflow:auto;top:56px}
  .vis-tooltip{max-width:88vw!important}
 }
</style></head><body>
<button id="toggle" aria-label="Toggle controls">&#9776;</button>
<div id="panel">
 __PUBMED_QUERY__
 <h1>Gene&ndash;gene interactions</h1>
 <div class="row mut" id="pubinfo"></div>
 <div class="row legend"><b style="background:#cfe3ff;border-color:#2b6cb0"></b>gene <b style="background:#1b7837;border-color:#145a28"></b>approved anti-neoplastic <b style="background:#e08600;border-color:#9a6700"></b>approved (other) <b style="background:#c2185b;border-color:#7a0f3a"></b>ChEBI</div>
 <div class="row mut">Drug-target genes (corpus chemicals): deeper colour = more chemicals. <b style="color:#1b7837">Green</b> = DGIdb approved anti-neoplastic, <b style="color:#e08600">amber</b> = DGIdb approved (non-anti-neoplastic), <b style="color:#c2185b">pink</b> = ChEBI, absent in DGIdb.</div>
 <div class="row">Min unique sentences/edge: <b id="thv">1</b><br><input id="thr" type="range" min="1" max="10" value="1"></div>
 <div class="row">Min cluster size: <select id="mincluster"><option>3</option><option>4</option><option>5</option><option>6</option><option>7</option><option>8</option><option>9</option><option>10</option><option>11</option><option>12</option></select></div>
 <div class="row">Year: <b id="yrlab"></b><br><input id="yrlo" type="range" style="width:74px"> <input id="yrhi" type="range" style="width:74px"></div>
 <div class="row">Search gene: <input id="search" placeholder="e.g. EGFR" autocomplete="off"></div>
 <div class="row">Filter to gene:<br><input id="genefilter" placeholder="e.g. EGFR (+neighbors)" autocomplete="off"> <select id="hops"><option value="1">1 hop</option><option value="2">2 hops</option></select></div>
 <div class="row">Search drug: <input id="drugsearch" placeholder="e.g. nivolumab" autocomplete="off"></div>
 <div class="row">Filter to drug:<br><select id="chemfilter"><option value="">(all drugs)</option></select></div>
 <div class="row">Filter to disease:<br><select id="disfilter"><option value="">(no disease)</option></select></div>
 <div class="row mut">Polarity:</div><div id="catfilters"></div>
 <div class="row" id="zoom"><button id="zin">+ Zoom in</button><button id="zout">&minus; Zoom out</button><button id="zfit">Fit</button></div>
 <div class="row">Relationship score: <b id="scval">&ge;0.99</b><br><input id="conf" type="range" min="0" max="13" step="1" value="13" aria-label="Minimum relationship score"></div>
 <div class="row mut" id="stats"></div>
 <div class="row mut" id="info">Click a node or edge for details.</div>
</div>
<div id="net"></div>
<script>
const DATA=__PAYLOAD__;
const CCOLOR=__CCOLOR__;
const MINY=__MINY__, MAXY=__MAXY__;
const MAXTGT=Math.max(1,...DATA.nodes.map(n=>n.target||0));
function nodeColor(n){if(!n.target)return {background:'#cfe3ff',border:'#2b6cb0'};const t=n.target/MAXTGT,L=(a,b)=>Math.round(a+(b-a)*t);if(n.tcat==='green')return {background:'rgb('+L(200,27)+','+L(230,120)+','+L(201,55)+')',border:'#145a28'};if(n.tcat==='amber')return {background:'rgb('+L(255,224)+','+L(231,134)+','+L(179,0)+')',border:'#9a6700'};return {background:'rgb('+L(255,194)+','+L(217,24)+','+L(232,91)+')',border:'#7a0f3a'};}
const net=document.getElementById('net'); let network=null;
const labelById={};DATA.nodes.forEach(n=>{labelById[n.id]=n.label;});
function esc(s){return (s||'').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');}
function pmA(p){return '<a class=pm target=_blank rel=noopener href="https://www.ncbi.nlm.nih.gov/pmc/articles/'+p+'/">'+p+'</a>';}
const SCORE_STEPS=[0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99]; // 0.05 up to 0.95, then 0.01
function activeConf(){const s=document.getElementById('conf');let i=s?parseInt(s.value):SCORE_STEPS.length-1;if(isNaN(i))i=SCORE_STEPS.length-1;return SCORE_STEPS[Math.max(0,Math.min(SCORE_STEPS.length-1,i))];}
function activeCats(){return new Set(Array.from(document.querySelectorAll(".catf:checked")).map(c=>c.value));}
function activeYears(){const a=parseInt(document.getElementById('yrlo').value),b=parseInt(document.getElementById('yrhi').value);return [Math.min(a,b),Math.max(a,b)];}
function passYear(yr,lo,hi){return (yr!=null&&yr>=lo&&yr<=hi)||(yr==null&&lo<=MINY&&hi>=MAXY);}
function activeDisease(){return (document.getElementById('disfilter').value||'').trim();}
function isSCC(d){d=(d||'').toLowerCase();return d.indexOf('squamous cell carcinoma')>=0||/^[a-z]*sccs?$/.test(d)||d==='lusc'||d==='ecss';}
function disMatch(s,dis){if(!dis)return true;if(!s.dis)return false;return dis==='__ALL_SCC__'?s.dis.some(isSCC):s.dis.indexOf(dis)>=0;}
function disLabel(dis){return dis==='__ALL_SCC__'?'all squamous cell carcinomas':dis;}
function visSents(e,conf,lo,hi,dis){return e.sents.filter(s=>s.sc>=conf&&passYear(s.yr,lo,hi)&&disMatch(s,dis));}
function edgeHead(e,vis){const np=new Set(vis.map(s=>s.pmid)).size;return '<div class=eth><b>'+esc(labelById[e.from])+' &rarr; '+esc(labelById[e.to])+'</b> ('+vis.length+' sentences &middot; '+np+' PMIDs &middot; '+e.cat+')</div>';}
function edgeTip(e,vis){const d=document.createElement('div');let h=edgeHead(e,vis);const lim=20;vis.slice(0,lim).forEach(s=>{h+='<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span> '+esc(s.text)+'</div>';});if(vis.length>lim)h+='<div class=more>+'+(vis.length-lim)+' more</div>';d.innerHTML=h;return d;}
function scaleNode(s){return 6+Math.sqrt(s)*3.4;}
function fontSize(c){return c<5?13:2*Math.max(13,Math.min(Math.round(c*2.2),48));}
function activeMinCluster(){const v=parseInt((document.getElementById('mincluster')||{}).value);return isNaN(v)?3:v;}
function build(thr){
 const conf=activeConf(), cats=activeCats(); const [ylo,yhi]=activeYears(); const dis=activeDisease(); const mc=activeMinCluster();
 let edges=[];
 DATA.edges.forEach(e=>{ if(!cats.has(e.cat))return; const vis=visSents(e,conf,ylo,yhi,dis); if(vis.length>=thr) edges.push({e:e,vis:vis,w:vis.length}); });
 const gf=(document.getElementById('genefilter').value||'').trim().toLowerCase();
 const chemSel=document.getElementById('chemfilter').value;
 let focusActive=false, focusLabel='';
 if(gf||chemSel){
   focusActive=true;
   const seeds=new Set(), labs=[];
   if(gf){const fn=DATA.nodes.find(n=>n.label.toLowerCase()===gf)||DATA.nodes.find(n=>n.label.toLowerCase().indexOf(gf)===0);if(fn){seeds.add(fn.id);labs.push(fn.label);}else labs.push('(no gene: '+gf+')');}
   if(chemSel){DATA.nodes.forEach(n=>{if((n.chems||[]).indexOf(chemSel)>=0)seeds.add(n.id);});labs.push('chem: '+chemSel);}
   const hops=parseInt(document.getElementById('hops').value)||1;
   const ag={};edges.forEach(o=>{(ag[o.e.from]=ag[o.e.from]||[]).push(o.e.to);(ag[o.e.to]=ag[o.e.to]||[]).push(o.e.from);});
   const seen=new Set(seeds);let fr=[...seeds];
   for(let h=0;h<hops;h++){const nf=[];fr.forEach(x=>{(ag[x]||[]).forEach(y=>{if(!seen.has(y)){seen.add(y);nf.push(y);}});});fr=nf;}
   edges=edges.filter(o=>seen.has(o.e.from)&&seen.has(o.e.to));
   focusLabel=labs.join(', ');
 }
 if(!focusActive){
   const adj={};
   edges.forEach(o=>{(adj[o.e.from]=adj[o.e.from]||[]).push(o.e.to);(adj[o.e.to]=adj[o.e.to]||[]).push(o.e.from);});
   const comp={};let cid=0;
   for(const n in adj){if(comp[n]!==undefined)continue;const sk=[n];comp[n]=cid;while(sk.length){const x=sk.pop();(adj[x]||[]).forEach(y=>{if(comp[y]===undefined){comp[y]=cid;sk.push(y);}});}cid++;}
   const csz={};for(const n in comp)csz[comp[n]]=(csz[comp[n]]||0)+1;
   edges=edges.filter(o=>csz[comp[o.e.from]]>=mc);
 }
 const keep=new Set();edges.forEach(o=>{keep.add(o.e.from);keep.add(o.e.to);});
 const nss={};edges.forEach(o=>{o.vis.forEach(s=>{(nss[o.e.from]=nss[o.e.from]||new Set()).add(s.text);(nss[o.e.to]=nss[o.e.to]||new Set()).add(s.text);});});
 const nsz=id=>(nss[id]?nss[id].size:0);
 const allCatsSel=[...new Set(DATA.edges.map(e=>e.cat))].every(c=>cats.has(c));
 const dfltView=!focusActive&&!dis&&ylo<=MINY&&yhi>=MAXY&&thr<=1&&allCatsSel&&mc<=3; // pristine view at either confidence
 const nodes=DATA.nodes.filter(n=>keep.has(n.id)).map(n=>({id:n.id,label:(dfltView?undefined:n.label),value:nsz(n.id),size:scaleNode(nsz(n.id)),title:n.label+' — '+nsz(n.id)+' unique sentences (in view)'+(n.target?' · drug target: '+n.target+' chemicals'+(n.tcat==='green'?' (approved anti-neoplastic)':(n.tcat==='amber'?' (approved)':' (ChEBI)')):''),color:nodeColor(n),font:{color:'#1a1a1a',size:(dfltView?13:fontSize(nsz(n.id)))}}));
 const eds=edges.map((o,i)=>({id:i,from:o.e.from,to:o.e.to,value:o.w,width:Math.min(1+o.w*0.7,10),color:{color:o.e.color,opacity:0.6},title:edgeTip(o.e,o.vis)}));
 const vpub=new Set();edges.forEach(o=>o.vis.forEach(s=>vpub.add(s.pmid)));
 document.getElementById('stats').innerHTML='Showing <b>'+nodes.length+'</b> genes, <b>'+eds.length+'</b> edges, <b>'+vpub.size+'</b> publications (&ge;'+conf+')'+(dis?' &middot; disease: <b>'+esc(disLabel(dis))+'</b>':'')+(focusActive?' &middot; focus: <b>'+esc(focusLabel)+'</b>':'');
 const data={nodes:new vis.DataSet(nodes),edges:new vis.DataSet(eds)};
 const options={layout:{improvedLayout:false},physics:{stabilization:{iterations:200},barnesHut:{gravitationalConstant:-14000,springLength:130,springConstant:0.02,avoidOverlap:0.3}},interaction:{hover:true,tooltipDelay:120},nodes:{shape:'dot',scaling:{min:6,max:60}},edges:{smooth:false,arrowStrikethrough:false,hoverWidth:0,selectionWidth:0,arrows:{to:{enabled:true,scaleFactor:0.6}}}};
 if(network)network.destroy();
 network=new vis.Network(net,data,options);
 network.on('stabilizationIterationsDone',()=>{network.setOptions({physics:false});network.fit({animation:false});network.redraw();});
 const _e=edges;
 network.on('click',p=>{const info=document.getElementById('info');
   if(p.nodes.length){const n=DATA.nodes.find(x=>x.id===p.nodes[0]);info.innerHTML='<b>'+n.label+'</b>: '+nsz(n.id)+' unique sentences (in view)';}
   else if(p.edges.length){const o=_e[p.edges[0]];info.innerHTML=edgeHead(o.e,o.vis)+o.vis.map(s=>'<div class=stip>'+pmA(s.pmid)+' <span class=mut>['+s.sc.toFixed(3)+(s.yr?(' · '+s.yr):'')+']</span> '+esc(s.text)+'</div>').join('');}});
}
const thr=document.getElementById('thr');
function buildCatFilters(){const counts={};DATA.edges.forEach(e=>counts[e.cat]=(counts[e.cat]||0)+1);const cats=Object.keys(counts).sort((a,b)=>counts[b]-counts[a]);document.getElementById('catfilters').innerHTML=cats.map(c=>'<label><input type=checkbox class=catf value="'+c+'" checked> <span class=sw style="background:'+(CCOLOR[c]||'#888')+'"></span> '+c+' ('+counts[c]+')</label>').join('');document.querySelectorAll('.catf').forEach(c=>c.addEventListener('change',()=>build(+thr.value)));}
thr.addEventListener('input',()=>{document.getElementById('thv').textContent=thr.value;build(+thr.value);});
const mcl=document.getElementById('mincluster');
mcl.addEventListener('change',()=>build(+thr.value));
const confEl=document.getElementById('conf');function updScore(){document.getElementById('scval').innerHTML='&ge;'+activeConf().toFixed(2);}updScore();confEl.addEventListener('input',()=>{updScore();build(+thr.value);});
document.getElementById('genefilter').addEventListener('change',()=>build(+thr.value));
document.getElementById('genefilter').addEventListener('keydown',ev=>{if(ev.key==='Enter')build(+thr.value);});
document.getElementById('hops').addEventListener('change',()=>build(+thr.value));
const searchBox=document.getElementById('search');
function doSearch(q){q=(q||'').trim();const info=document.getElementById('info');if(!q)return;const hit=DATA.nodes.find(n=>n.label.toLowerCase()===q.toLowerCase())||DATA.nodes.find(n=>n.label.toLowerCase().indexOf(q.toLowerCase())===0);if(!hit){info.innerHTML='No gene matching "'+q+'"';return;}try{network.selectNodes([hit.id]);network.focus(hit.id,{scale:1.3,animation:true});info.innerHTML='<b>'+hit.label+'</b>';}catch(e){info.innerHTML='<b>'+hit.label+'</b> not in current view';}}
searchBox.addEventListener('keydown',ev=>{if(ev.key==='Enter')doSearch(searchBox.value);});
searchBox.addEventListener('change',()=>doSearch(searchBox.value));
function zoomBy(f){if(!network)return;const s=network.getScale();network.moveTo({scale:s*f,animation:{duration:200}});}
document.getElementById('zin').addEventListener('click',()=>zoomBy(1.25));
document.getElementById('zout').addEventListener('click',()=>zoomBy(0.8));
document.getElementById('zfit').addEventListener('click',()=>{if(network)network.fit({animation:true});});
const yl=document.getElementById('yrlo'),yh=document.getElementById('yrhi');
[yl,yh].forEach(el=>{el.min=MINY;el.max=MAXY;});yl.value=MINY;yh.value=MAXY;
function updYr(){const a=activeYears();document.getElementById('yrlab').textContent=a[0]+'–'+a[1];}
updYr();
yl.addEventListener('input',()=>{updYr();build(+thr.value);});
yh.addEventListener('input',()=>{updYr();build(+thr.value);});
const chemGenes={};DATA.nodes.forEach(n=>(n.chems||[]).forEach(c=>{(chemGenes[c]=chemGenes[c]||[]).push(n.label);}));const csel=document.getElementById('chemfilter');Object.keys(chemGenes).sort().forEach(c=>{const g=chemGenes[c].slice().sort();const o=document.createElement('option');o.value=c;o.textContent=c+' → '+g.join(', ');csel.appendChild(o);});csel.addEventListener('change',()=>build(+thr.value));
const drugBox=document.getElementById('drugsearch');function findDrug(q){q=(q||'').trim().toLowerCase();if(!q)return;const info=document.getElementById('info');const opts=[...csel.options].filter(o=>o.value);const m=opts.find(o=>o.value.toLowerCase()===q)||opts.find(o=>o.value.toLowerCase().indexOf(q)===0)||opts.find(o=>o.value.toLowerCase().indexOf(q)>=0);if(m){csel.value=m.value;build(+thr.value);info.innerHTML='Drug filter: <b>'+esc(m.value)+'</b>';}else{info.innerHTML='No drug matching "'+esc(q)+'"';}}drugBox.addEventListener('keydown',ev=>{if(ev.key==='Enter')findDrug(drugBox.value);});drugBox.addEventListener('change',()=>findDrug(drugBox.value));
const disCounts={},sccMembers={};let sccCount=0;DATA.edges.forEach(e=>e.sents.forEach(s=>{let sHasSCC=false;(s.dis||[]).forEach(d=>{disCounts[d]=(disCounts[d]||0)+1;if(isSCC(d)){sccMembers[d]=(sccMembers[d]||0)+1;sHasSCC=true;}});if(sHasSCC)sccCount++;}));const dsel=document.getElementById('disfilter');const sccLabels=Object.keys(sccMembers).sort((a,b)=>sccMembers[b]-sccMembers[a]);const disOpts=[];if(sccLabels.length)disOpts.push({value:'__ALL_SCC__',text:'all squamous cell carcinomas — '+sccLabels.join(', '),count:sccCount});Object.keys(disCounts).filter(d=>disCounts[d]>1&&!isSCC(d)).forEach(d=>disOpts.push({value:d,text:d,count:disCounts[d]}));disOpts.sort((a,b)=>b.count-a.count).forEach(x=>{const o=document.createElement('option');o.value=x.value;o.textContent=x.text+' ('+x.count+')';dsel.appendChild(o);});if(disOpts.length)dsel.value=disOpts[0].value;dsel.addEventListener('change',()=>build(+thr.value));
document.getElementById('toggle').addEventListener('click',()=>document.getElementById('panel').classList.toggle('collapsed'));
if(window.innerWidth<=700)document.getElementById('panel').classList.add('collapsed');
window.addEventListener('resize',()=>{if(network)network.redraw();});
(function(){const pm=new Set();DATA.edges.forEach(e=>e.sents.forEach(s=>pm.add(s.pmid)));document.getElementById('pubinfo').innerHTML='<b>'+pm.size+'</b> out of <b>'+__NXML__+'</b> produced high-score gene-gene interactions';})();
buildCatFilters();build(1);
</script></body></html>"""


def render_graph(payload, lib, miny, maxy, nxml, pubmed_query=""):
    if lib:
        libtag = "<script>\n" + lib.replace("</script>", "<\\/script>") + "\n</script>"
    else:
        libtag = f'<script src="{VIS_URL}"></script>'
    # legend line above the title: the PubMed query this corpus came from (empty -> nothing shown)
    qrow = (f'<div class="row" id="pubmedq" style="font-size:18px;font-weight:600;color:#2b6cb0;'
            f'margin-bottom:6px">PubMed query = {html.escape(pubmed_query)}</div>'
            if pubmed_query else "")
    title = html.escape(pubmed_query) if pubmed_query else "High-confidence brain-cancer gene-gene interactions (gene-only filter)"
    return (GRAPH_TEMPLATE.replace("__LIBTAG__", libtag)
            .replace("__TITLE__", title)
            .replace("__PUBMED_QUERY__", qrow)
            .replace("__PAYLOAD__", json.dumps(payload, ensure_ascii=False))
            .replace("__CCOLOR__", json.dumps(PCOLOR))
            .replace("__NXML__", str(nxml))
            .replace("__MINY__", str(miny)).replace("__MAXY__", str(maxy)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=str(DATA_ROOT),
                    help="pipeline output tree to read inputs from (default: kaggle_working/ next to this script)")
    ap.add_argument("--score", type=float, default=0.8, help="threshold for the exported JSON (default 0.8)")
    ap.add_argument("--thresholds", default="0.8,0.95,0.99", help="thresholds summarized in the console stats")
    ap.add_argument("--no-graph", action="store_true", help="skip the brain-cancer graph HTML")
    args = ap.parse_args()

    set_data_root(args.data_root)
    if not RE_FILE.exists():
        raise SystemExit(f"ERROR: {RE_FILE} not found under data root {DATA_ROOT} "
                         f"(run the gpu_bundle pipeline first, or pass --data-root).")

    d = json.loads(RE_FILE.read_text(encoding="utf-8"))
    base = json.loads(BASE_FILE.read_text(encoding="utf-8"))   # DISEASE facets for the graph's disease filter
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) JSON export at --score
    kept = [t for t in d if qualifies(t, args.score)]
    JSON_OUT.write_text(json.dumps(kept, ensure_ascii=False), encoding="utf-8")
    ksents = len({t.get("sentence") for t in kept})
    print(f"exported {len(kept):,} triples ({ksents:,} sentences) at score>={args.score} -> {JSON_OUT}")

    # 2) threshold stats (for the console summary below)
    rows = [stats(d, float(x)) for x in args.thresholds.split(",")]

    # 3) brain-cancer gene-gene graph (universe = qualifying at GRAPH_BASE; in-browser score slider 0.5..0.99)
    if not args.no_graph:
        universe = kept if args.score <= GRAPH_BASE else [t for t in d if qualifies(t, GRAPH_BASE)]
        payload = graph_payload(universe, base)
        yrs = [s["yr"] for e in payload["edges"] for s in e["sents"] if s.get("yr")]
        miny, maxy = (min(yrs), max(yrs)) if yrs else (2000, 2026)
        lib = get_vis_lib()
        # corpus size = # source documents; experimental_ner/ is empty in the bundle
        # (it was a runtime symlink), so fall back to the per-document sentences/ files.
        nxml = len(list(XML_DIR.glob("*.xml"))) or len(list(SENT_DIR.glob("*.json")))
        GRAPH_OUT.parent.mkdir(parents=True, exist_ok=True)
        pubmed_query = read_pubmed_query()
        GRAPH_OUT.write_text(render_graph(payload, lib, miny, maxy, nxml, pubmed_query), encoding="utf-8")
        n99 = sum(1 for e in payload["edges"] if any(s["sc"] >= 0.99 for s in e["sents"]))
        npubs = len({s["pmid"] for e in payload["edges"] for s in e["sents"]})
        print(f"wrote graph -> {GRAPH_OUT}  ({len(payload['nodes']):,} genes, {len(payload['edges']):,} edges; "
              f"{n99:,} edges have a >=0.99 sentence; {npubs:,} of {nxml:,} input XMLs produced gene-gene triples)")

        # also drop a copy in the root dir, named after the PubMed query (whitespace -> underscore;
        # other filesystem-illegal characters likewise underscored so the name is always valid).
        # "_G" suffix keeps this distinct from high_confidence.py's "<query>.html" copy.
        if pubmed_query:
            fname = re.sub(r'[\\/:*?"<>|\s]+', "_", pubmed_query.strip())
            dest = ROOT / f"{fname}_G.html"
            shutil.copy2(GRAPH_OUT, dest)
            print(f"copied graph -> {dest}")

    for r in rows:
        print(f"  score>={r['T']}: {r['triples']:,} triples / {r['sentences']:,} sent; "
              f"G {r['filt']:,} triples / {r['filt_sentences']:,} sent")


if __name__ == "__main__":
    main()
