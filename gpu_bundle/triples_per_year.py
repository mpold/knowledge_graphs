#!/usr/bin/env python3
"""triples_per_year.py -- plot the quality-score distribution of the triples per publication year.

Reads the final relation-extraction triples (relation_extraction.py's normalized
multi-type file), dates each one by the publication year of its source article, and
writes a self-contained HTML page with three panels:

  1. a box plot of the composite `score` per year (box = Q1-Q3, line = median,
     whiskers = 5th-95th percentile; years with fewer than --min-n triples are drawn
     faded, because their spread is not reliable)
  2. the number of triples per year, i.e. how much evidence sits behind each box
  3. a histogram of every dated score, in bins of 0.025

plus a per-year table (n, publications, min/P5/Q1/median/Q3/P95/max/mean).

Step 19 of gpu.py (optional; after relation_extraction writes the triples, before zip_work
bundles the run). Standalone it works on any finished run tree via --root.

YEARS come from databases/pmc_years.json, the PMC -> year cache pub_years.py (step 16)
builds, topped up from the stage-1 table pmids/pmid_pmc_ids.tsv (pubmed_query.py) when one
is found near the run -- the table has a year for every document in the corpus, unlike the
XML itself, where half the GROBID TEI files carry no date. Nothing here goes to the
network. A triple's pmid is reduced to its bare PMC accession first, since GROBID-derived
documents carry ids like "PMC123.grobid.tei". Triples whose article has no year in either
source are counted as undated and left off the plots; the page reports how many.

Run::  python triples_per_year.py [--root kaggle_working] [--triples FILE] [--pmid-tsv FILE] [--out FILE] [--min-n 100]
"""
import argparse
import collections
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# relative to the run root (--root), the working dir gpu.py runs every step in
TRIPLES = "TRIPLES/triples_re_GENETIC_DISEASE_CHEMICAL_normalized.json"   # relation_extraction.py (step 17)
YEARS_CACHE = "databases/pmc_years.json"                                   # pub_years.py (step 16)
PMID_TSV = "pmids/pmid_pmc_ids.tsv"                                        # stage-1 table (pubmed_query.py)
OUT = "summaries/triples_per_year.html"
NBINS = 40   # histogram bins over [0, 1]


def find_tsv(explicit, root):
    """The stage-1 table, or None -- searched where pub_years.py looks for it."""
    if explicit:
        p = Path(explicit)
        return p if p.exists() else None
    for base in (root, root.parent, ROOT, ROOT.parent):
        if (base / PMID_TSV).exists():
            return base / PMID_TSV
    return None


def load_years(root, tsv):
    """PMC accession -> publication year: the step-16 cache first, the table for the rest."""
    years, cache = {}, root / YEARS_CACHE
    if cache.exists():
        years = {k: int(v) for k, v in json.loads(cache.read_text(encoding="utf-8")).items() if v}
        print(f"  {YEARS_CACHE}: {len(years):,} accessions with a year")
    if tsv:
        with open(tsv, encoding="utf-8", newline="") as fh:
            table = {r["pmc_id"]: int(r["year"]) for r in csv.DictReader(fh, delimiter="\t")
                     if r.get("pmc_id") and (r.get("year") or "").strip().isdigit()}
        added = {k: v for k, v in table.items() if k not in years}
        years.update(added)
        print(f"  {tsv.name}: {len(table):,} accessions with a year; {len(added):,} not in the cache")
    if not years:
        raise SystemExit(f"no publication years: neither {cache} nor {PMID_TSV} was found "
                         f"(run pub_years.py first, or pass --pmid-tsv)")
    return years


def quantile(sorted_vals, p):
    """Linear-interpolated quantile of an already sorted list."""
    i = (len(sorted_vals) - 1) * p
    lo = int(i)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (i - lo)


def summarise(vals):
    """[n, min, p5, q1, median, q3, p95, max, mean] -- the row shape the page's script expects."""
    v = sorted(vals)
    if not v:
        return [0] + [0.0] * 8
    return [len(v)] + [round(quantile(v, p), 4) for p in (0, .05, .25, .5, .75, .95, 1)] + [round(sum(v) / len(v), 4)]


def collect(triples, years):
    by_year = collections.defaultdict(list)
    pubs = collections.defaultdict(set)
    undated = []
    for t in triples:
        s = t.get("score")
        if s is None:
            continue
        acc = (t.get("pmid") or "").split(".")[0]   # "PMC123.grobid.tei" -> "PMC123"
        y = years.get(acc)
        if y is None:
            undated.append(s)
            continue
        by_year[y].append(s)
        pubs[y].add(acc)
    dated = [s for v in by_year.values() for s in v]
    hist = [0] * NBINS
    for s in dated:
        hist[min(int(s * NBINS), NBINS - 1)] += 1
    return {"rows": [[y, len(pubs[y])] + summarise(by_year[y]) for y in sorted(by_year)],
            "undated": summarise(undated), "all": summarise(dated), "hist": hist,
            "total": len(triples)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", "--data-root", default=str(ROOT),
                    help="run tree holding TRIPLES/ and databases/ (default: next to this script)")
    ap.add_argument("--triples", help=f"triples JSON to plot (default: <root>/{TRIPLES})")
    ap.add_argument("--pmid-tsv", default=None,
                    help=f"stage-1 table with pmc_id and year columns, topping up <root>/{YEARS_CACHE} "
                         f"(default: {PMID_TSV} found near the root or this script)")
    ap.add_argument("--out", help=f"HTML to write (default: <root>/{OUT})")
    ap.add_argument("--min-n", type=int, default=100, help="years with fewer triples are drawn faded (default 100)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    src = Path(args.triples) if args.triples else root / TRIPLES
    if not src.exists():
        raise SystemExit(f"no triples file at {src} (relation_extraction.py writes it)")
    with open(src, encoding="utf-8") as fh:
        triples = json.load(fh)
    data = collect(triples, load_years(root, find_tsv(args.pmid_tsv, root)))
    if not data["rows"]:
        raise SystemExit(f"no dated triples with a score in {src}")

    html = (TEMPLATE.replace("__SOURCE__", src.name)
                    .replace("__MIN_N__", str(args.min_n))
                    .replace("__DATA__", json.dumps(data)))
    out = Path(args.out) if args.out else root / OUT
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    a = data["all"]
    print(f"{data['total']:,} triples: {a[0]:,} dated, {data['undated'][0]:,} undated; "
          f"median score {a[4]:.3f} (IQR {a[3]:.3f}-{a[5]:.3f}), {data['rows'][0][0]}-{data['rows'][-1][0]}")
    print(f"wrote {out}")


TEMPLATE = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Triple Scores by Year</title>
<style>
:root{color-scheme:light;--surface:#fcfcfb;--text-primary:#0b0b0b;--text-secondary:#52514e;--muted:#8a8984;--grid:#e6e5e0;--series-1:#2a78d6;--series-1-soft:#a9c8ef;--tip-bg:#fff;--tip-border:#d9d8d2}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])){color-scheme:dark;--surface:#1a1a19;--text-primary:#fff;--text-secondary:#c3c2b7;--muted:#8f8e86;--grid:#2e2e2c;--series-1:#3987e5;--series-1-soft:#2a4a72;--tip-bg:#262625;--tip-border:#3a3a38}}
:root[data-theme="dark"]{color-scheme:dark;--surface:#1a1a19;--text-primary:#fff;--text-secondary:#c3c2b7;--muted:#8f8e86;--grid:#2e2e2c;--series-1:#3987e5;--series-1-soft:#2a4a72;--tip-bg:#262625;--tip-border:#3a3a38}
body{margin:0;background:var(--surface);color:var(--text-primary);font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1000px;margin:0 auto;padding:24px 16px}
h1{font-size:20px;margin:0 0 4px}h2{font-size:15px;margin:28px 0 2px}
.sub{color:var(--text-secondary);margin:0 0 20px}.cap{color:var(--text-secondary);font-size:12px;margin:0 0 8px}
.stats{display:flex;flex-wrap:wrap;gap:28px;margin-bottom:8px}
.stat b{display:block;font-size:22px;font-variant-numeric:tabular-nums}.stat span{color:var(--text-secondary);font-size:12px}
.chart{position:relative}svg{display:block;width:100%;height:auto}
.axis text{fill:var(--muted);font-size:11px;font-variant-numeric:tabular-nums}
.grid line{stroke:var(--grid)}
.box{fill:var(--series-1)}.whisk{stroke:var(--series-1);stroke-width:2}.med{stroke:var(--surface);stroke-width:2}
.thin .box{fill:var(--series-1-soft)}.thin .whisk{stroke:var(--series-1-soft)}.thin .med{stroke:var(--series-1)}
.bar{fill:var(--series-1)}.hit{fill:transparent}.on{opacity:.75}
.legend{display:flex;flex-wrap:wrap;gap:18px;font-size:12px;color:var(--text-secondary);margin:4px 0 6px}
.legend i{display:inline-block;width:12px;height:12px;border-radius:2px;vertical-align:-2px;margin-right:6px}
.tip{position:absolute;pointer-events:none;background:var(--tip-bg);border:1px solid var(--tip-border);border-radius:6px;padding:6px 10px;font-size:12px;box-shadow:0 2px 8px rgba(0,0,0,.12);display:none;white-space:nowrap;font-variant-numeric:tabular-nums}
.tip b{font-size:13px}.tip .m{color:var(--text-secondary)}
.note{color:var(--text-secondary);font-size:12px;margin-top:10px}
details{margin-top:24px}summary{cursor:pointer;color:var(--text-secondary)}
.tw{overflow-x:auto}table{border-collapse:collapse;margin-top:8px;font-variant-numeric:tabular-nums;font-size:13px}
th,td{padding:3px 14px 3px 0;text-align:right}th{color:var(--text-secondary);font-weight:600}td:first-child,th:first-child{text-align:left}
</style></head><body><main>
<h1>Triple quality scores by publication year</h1>
<p class="sub">Composite <code>score</code> of each triple in <code>__SOURCE__</code>, grouped by the publication year of the source article</p>
<div class="stats" id="stats"></div>

<h2>Score distribution per year</h2>
<div class="legend"><span><i style="background:var(--series-1)"></i>Box: 25th–75th percentile, light line = median</span><span>Whiskers: 5th–95th percentile</span><span><i style="background:var(--series-1-soft)"></i>Fewer than __MIN_N__ triples (unreliable)</span></div>
<div class="chart"><svg id="box" role="img" aria-label="Box plot of triple scores per publication year"></svg><div class="tip"></div></div>

<h2>Triples per year</h2>
<p class="cap">How many triples are behind each box</p>
<div class="chart"><svg id="cnt" role="img" aria-label="Bar chart of triple counts per publication year"></svg><div class="tip"></div></div>

<h2>Overall score distribution</h2>
<p class="cap">All dated triples, in bins of 0.025</p>
<div class="chart"><svg id="hist" role="img" aria-label="Histogram of all triple scores"></svg><div class="tip"></div></div>

<p class="note" id="note"></p>
<details><summary>Show data table</summary><div class="tw"><table id="t"><thead><tr><th>Year</th><th>Triples</th><th>Pubs</th><th>Min</th><th>P5</th><th>Q1</th><th>Median</th><th>Q3</th><th>P95</th><th>Max</th><th>Mean</th></tr></thead><tbody></tbody></table></div></details>
</main>
<script>
const D=__DATA__, MIN_N=__MIN_N__;
const fmt=n=>n.toLocaleString('en-US'), f3=v=>v.toFixed(3);
// row: [year,pubs,n,min,p5,q1,med,q3,p95,max,mean]
const byY=new Map(D.rows.map(r=>[r[0],r]));
const y0=D.rows[0][0],y1=D.rows[D.rows.length-1][0];
const rows=[];for(let y=y0;y<=y1;y++)rows.push(byY.get(y)||[y,0,0]);
const A=D.all;
document.getElementById('stats').innerHTML=[[fmt(A[0]),'dated triples'],[f3(A[4]),'median score'],[f3(A[3])+'–'+f3(A[5]),'interquartile range'],[f3(A[8]),'mean score'],[fmt(D.undated[0]),'undated (not plotted)']].map(([v,l])=>`<div class="stat"><b>${v}</b><span>${l}</span></div>`).join('');
document.getElementById('note').textContent=(D.undated[0]?`${fmt(D.undated[0])} triples from articles with no known publication year are left out (their median score is ${f3(D.undated[4])}, IQR ${f3(D.undated[3])}–${f3(D.undated[5])}). `:'Every triple has a publication year. ')+`Scores run from ${f3(A[1])} to ${f3(A[7])}.`;
document.querySelector('#t tbody').innerHTML=rows.filter(r=>r[2]).map(r=>`<tr><td>${r[0]}</td><td>${fmt(r[2])}</td><td>${fmt(r[1])}</td>${r.slice(3).map(v=>`<td>${f3(v)}</td>`).join('')}</tr>`).join('');

const W=960,m={t:10,r:8,b:28,l:56},iw=W-m.l-m.r;
function frame(svg,H,yv,ticks,tf,xs){svg.setAttribute('viewBox',`0 0 ${W} ${H}`);let s='<g class="grid axis">';
 ticks.forEach(v=>s+=`<line x1="${m.l}" x2="${W-m.r}" y1="${yv(v)}" y2="${yv(v)}"/><text x="${m.l-8}" y="${yv(v)+4}" text-anchor="end">${tf(v)}</text>`);
 s+='</g><g class="axis">';xs.forEach(([x,l])=>s+=`<text x="${x}" y="${H-8}" text-anchor="${x>W-m.r-20?"end":"middle"}">${l}</text>`);return s+'</g>'}
function hover(svg,html){const box=svg.parentNode,tip=box.querySelector('.tip');
 const clear=()=>{tip.style.display='none';svg.querySelectorAll('.on').forEach(b=>b.classList.remove('on'))};
 svg.addEventListener('mousemove',e=>{const t=e.target.closest('.hit');clear();if(!t)return;const h=html(+t.dataset.i);if(!h)return;
  svg.querySelectorAll(`[data-g="${t.dataset.i}"]`).forEach(b=>b.classList.add('on'));tip.innerHTML=h;tip.style.display='block';
  const bb=box.getBoundingClientRect();let lx=e.clientX-bb.left+12;if(lx+tip.offsetWidth>bb.width)lx=e.clientX-bb.left-tip.offsetWidth-12;
  tip.style.left=lx+'px';tip.style.top=Math.max(0,e.clientY-bb.top-tip.offsetHeight-8)+'px'});
 svg.addEventListener('mouseleave',clear)}
function barPath(x,w,base,top){const h=base-top,r=Math.min(4,w/2,h);if(h<=0)return'';
 return `M${x},${base}V${top+r}Q${x},${top} ${x+r},${top}H${x+w-r}Q${x+w},${top} ${x+w},${top+r}V${base}Z`}
const bw=iw/rows.length, xYears=rows.map((r,i)=>[m.l+bw*i+bw/2,r[0]]).filter(([,y])=>y%5===0);

// box plot
{const svg=document.getElementById('box'),H=400,ih=H-m.t-m.b,lo=0.2,hi=1.0,yv=v=>m.t+ih-ih*(v-lo)/(hi-lo);
 let s=frame(svg,H,yv,[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],v=>v.toFixed(2),xYears);
 rows.forEach((r,i)=>{if(!r[2])return;const cx=m.l+bw*i+bw/2,w=Math.max(4,bw*0.6),x=cx-w/2;
  s+=`<g class="${r[2]<MIN_N?'thin':''}" data-g="${i}"><line class="whisk" x1="${cx}" x2="${cx}" y1="${yv(r[4])}" y2="${yv(r[5])}"/><line class="whisk" x1="${cx}" x2="${cx}" y1="${yv(r[7])}" y2="${yv(r[8])}"/>
   <rect class="box" x="${x}" y="${yv(r[7])}" width="${w}" height="${Math.max(1,yv(r[5])-yv(r[7]))}" rx="3"/><line class="med" x1="${x}" x2="${x+w}" y1="${yv(r[6])}" y2="${yv(r[6])}"/></g>`});
 rows.forEach((r,i)=>s+=`<rect class="hit" data-i="${i}" x="${m.l+bw*i}" y="${m.t}" width="${bw}" height="${ih}"/>`);
 svg.innerHTML=s;
 hover(svg,i=>{const r=rows[i];if(!r[2])return`<b>${r[0]}</b><br><span class="m">no triples</span>`;
  return `<b>${r[0]}</b> <span class="m">· ${fmt(r[2])} triples${r[2]<MIN_N?' (few)':''}</span><br>Median ${f3(r[6])}<br>IQR ${f3(r[5])}–${f3(r[7])}<br>P5–P95 ${f3(r[4])}–${f3(r[8])}<br><span class="m">Mean ${f3(r[10])} · range ${f3(r[3])}–${f3(r[9])}</span>`})}

// counts
{const svg=document.getElementById('cnt'),H=200,ih=H-m.t-m.b,max=Math.max(...rows.map(r=>r[2])),tk=10000,yMax=Math.ceil(max/tk)*tk,yv=v=>m.t+ih-ih*v/yMax;
 const ticks=[];for(let v=0;v<=yMax;v+=tk)ticks.push(v);
 let s=frame(svg,H,yv,ticks,v=>v?(v/1000)+'k':'0',xYears);const gap=Math.min(2,bw*.2);
 rows.forEach((r,i)=>{s+=`<path class="bar" data-g="${i}" d="${barPath(m.l+bw*i+gap/2,bw-gap,yv(0),yv(r[2]))}"/><rect class="hit" data-i="${i}" x="${m.l+bw*i}" y="${m.t}" width="${bw}" height="${ih}"/>`});
 svg.innerHTML=s;hover(svg,i=>`<b>${rows[i][0]}</b><br>${fmt(rows[i][2])} triples<br><span class="m">${fmt(rows[i][1])} publications</span>`)}

// histogram
{const svg=document.getElementById('hist'),H=220,ih=H-m.t-m.b,h=D.hist,nb=h.length,lo=0.2,off=Math.round(lo*nb),bins=h.slice(off);
 const hb=iw/bins.length,max=Math.max(...bins),tk=10000,yMax=Math.ceil(max/tk)*tk,yv=v=>m.t+ih-ih*v/yMax;
 const ticks=[];for(let v=0;v<=yMax;v+=tk)ticks.push(v);
 const xs=[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0].map(v=>[m.l+iw*(v-lo)/(1-lo),v.toFixed(2)]);
 let s=frame(svg,H,yv,ticks,v=>v?(v/1000)+'k':'0',xs);
 bins.forEach((c,i)=>{s+=`<path class="bar" data-g="${i}" d="${barPath(m.l+hb*i+1,hb-2,yv(0),yv(c))}"/><rect class="hit" data-i="${i}" x="${m.l+hb*i}" y="${m.t}" width="${hb}" height="${ih}"/>`});
 svg.innerHTML=s;hover(svg,i=>{const a=(i+off)/nb,b=(i+off+1)/nb;return `<b>${a.toFixed(3)}–${b.toFixed(3)}</b><br>${fmt(bins[i])} triples<br><span class="m">${(100*bins[i]/A[0]).toFixed(1)}% of dated</span>`})}
</script></body></html>
'''


if __name__ == "__main__":
    main()
