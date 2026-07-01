#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""keys_values.py -- summary reports over the HGNC-linkage libraries.

Reads the seven JSON libraries produced by roman.py and greek.py (each a JSON
object keyed by the source entity, every entry carrying a nested `hgnc_symbol`
and `occurrences`, among other keys) and writes three self-contained HTML
reports under GENETIC/:

  * hgnc.html         unique `hgnc_symbol` values across ALL seven libraries and,
                      for each, the summed `occurrences` (symbol-attributed:
                      list-valued ambiguous/complex and '|'-joined cosine values
                      are split into member symbols, each credited the full
                      occurrence count) -- with a top-40 bar plot.
  * hgnc_1_to_1.html  the same, restricted to the two single-gene (1-to-1)
                      libraries roman.json + greek.json.
  * nested_keys.html  the nested keys of each entry object per file (coverage,
                      value type(s), example) + a cross-file presence matrix, and
                      the `match_mode` values + frequency per file.

Run from anywhere (paths resolve relative to this file)::

    python keys_values.py
"""

import html
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SUM = ROOT / "GENETIC"

ALL_FILES = ["roman.json", "roman_ambiguous.json", "roman_cosine.json",
             "greek.json", "greek_ambiguous.json", "greek_complex.json",
             "greek_cosine.json"]
ROMAN_FILES = ["roman.json", "roman_ambiguous.json", "roman_cosine.json"]
GREEK_FILES = ["greek.json", "greek_ambiguous.json", "greek_complex.json",
               "greek_cosine.json"]
ONE_TO_ONE_FILES = ["roman.json", "greek.json"]   # single-gene links only

TOPN = 40
esc = html.escape


# ---------------------------------------------------------------- helpers
def load(f):
    return json.loads((SUM / f).read_text(encoding="utf-8"))


def syms_of(v):
    """Yield individual HGNC symbols from an hgnc_symbol value: a str, a
    '|'-joined str (cosine), or a list of str (ambiguous / complex)."""
    for it in (v if isinstance(v, list) else [v]):
        if isinstance(it, str):
            for s in it.split("|"):
                s = s.strip()
                if s:
                    yield s


def to_int(o):
    if isinstance(o, int):
        return o
    if isinstance(o, str) and o.isdigit():
        return int(o)
    return 0


def typename(v):
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, int):
        return "int"
    if isinstance(v, float):
        return "float"
    if isinstance(v, str):
        return "str"
    if isinstance(v, list):
        inner = sorted({typename(x) for x in v}) or ["empty"]
        return f"list[{'|'.join(inner)}]"
    if isinstance(v, dict):
        return "object"
    if v is None:
        return "null"
    return type(v).__name__


def short(v, n=60):
    s = json.dumps(v, ensure_ascii=False)
    return s if len(s) <= n else s[:n - 1] + "…"


def css(headline_bg, headline_border):
    return """
  body { font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;
         max-width:1000px; margin:2rem auto; padding:0 1rem; color:#1a1a1a; }
  h1 { font-size:1.5rem; } h2 { font-size:1.15rem; margin-top:1.8rem;
       border-bottom:1px solid #ddd; padding-bottom:.3rem; }
  h3 { font-size:1rem; margin-top:1.3rem; }
  table { border-collapse:collapse; margin:.7rem 0; font-size:.9em; }
  th,td { border:1px solid #ccc; padding:.3rem .7rem; text-align:left; vertical-align:top; }
  th { background:#f7f7f7; }
  .num { text-align:right; font-variant-numeric:tabular-nums; }
  .big { font-size:2rem; font-weight:700; }
  .headline { background:%s; border:1px solid %s; border-radius:8px;
              padding:.8rem 1rem; margin:1rem 0; }
  code { background:#f3f3f3; padding:1px 5px; border-radius:3px; font-size:.9em; }
  .plot { border:1px solid #e3e3e3; border-radius:8px; padding:.5rem; margin:1rem 0; }
  .yes { text-align:center; color:#1a7f37; font-weight:700; }
  .no  { text-align:center; background:#fafafa; }
  details { margin:.6rem 0; } summary { cursor:pointer; color:#357; }
""" % (headline_bg, headline_border)


def svg_topn(ranked, bar_color):
    """An inline horizontal bar chart (no external deps) of the top-N
    (symbol, summed-occurrences) pairs."""
    top = ranked[:TOPN]
    maxv = top[0][1] if top else 1
    row_h, pad_l, pad_r, pad_t, bar_w = 22, 110, 70, 10, 640
    w_total = pad_l + bar_w + pad_r
    h_total = pad_t * 2 + row_h * len(top)
    bars = []
    for i, (s, n) in enumerate(top):
        y = pad_t + i * row_h
        w = max(1, round(bar_w * n / maxv))
        bars.append(
            f'<text x="{pad_l - 6}" y="{y + row_h / 2 + 4}" text-anchor="end" '
            f'font-size="11" fill="#222">{esc(s)}</text>'
            f'<rect x="{pad_l}" y="{y + 3}" width="{w}" height="{row_h - 7}" '
            f'fill="{bar_color}" rx="2"><title>{esc(s)}: {n:,}</title></rect>'
            f'<text x="{pad_l + w + 5}" y="{y + row_h / 2 + 4}" font-size="10.5" '
            f'fill="#333">{n:,}</text>')
    return (f'<svg viewBox="0 0 {w_total} {h_total}" width="100%" role="img" '
            f'aria-label="Top {TOPN} HGNC symbols by summed occurrences" '
            f'xmlns="http://www.w3.org/2000/svg" '
            f'font-family="Segoe UI, sans-serif">' + "".join(bars) + "</svg>")


def symbol_stats(files):
    """(per_file: f->(n_entries, set(symbols)), occ: symbol->summed occurrences,
    ranked: [(symbol, occ)] desc)."""
    per_file = {}
    occ = defaultdict(int)
    for f in files:
        data = load(f)
        s = set()
        for entry in data.values():
            hs = entry.get("hgnc_symbol")
            if hs is None:
                continue
            n = to_int(entry.get("occurrences"))
            for sym in syms_of(hs):
                s.add(sym)
                occ[sym] += n
        per_file[f] = (len(data), s)
    ranked = sorted(occ.items(), key=lambda kv: (-kv[1], kv[0].casefold()))
    return per_file, dict(occ), ranked


# ---------------------------------------------------------------- reports
def write_hgnc_report(files, out, h1, intro, headline_bg, headline_border,
                      bar_color, corpus_groups=None):
    """Unique-symbol + summed-occurrence report for `files`. `corpus_groups`, if
    given, is a list of (label, [files]) rendered as a by-group breakdown."""
    per_file, occ, ranked = symbol_stats(files)
    union = set(occ)
    total_occ = sum(occ.values())

    prows = "".join(
        f'<tr><td><code>{esc(f)}</code></td><td class="num">{per_file[f][0]:,}</td>'
        f'<td class="num">{len(per_file[f][1]):,}</td></tr>' for f in files)
    toprows = "".join(
        f'<tr><td class="num">{i}</td><td><code>{esc(s)}</code></td>'
        f'<td class="num">{n:,}</td></tr>' for i, (s, n) in enumerate(ranked[:TOPN], 1))
    fulllist = "".join(
        f'<tr><td><code>{esc(s)}</code></td><td class="num">{n:,}</td></tr>'
        for s, n in ranked)

    if corpus_groups is not None:
        unions = {label: set().union(*(per_file[f][1] for f in fs))
                  for label, fs in corpus_groups}
        (la, lb) = [u for u in unions.values()]
        grp_rows = "".join(
            f'<tr><td>{esc(label)} ({len(fs)} file{"s" if len(fs) > 1 else ""})</td>'
            f'<td class="num">{len(unions[label]):,}</td></tr>'
            for label, fs in corpus_groups)
        grp_rows += (f'<tr><td>shared by both groups</td>'
                     f'<td class="num">{len(la & lb):,}</td></tr>')
        grp_rows += (f'<tr><td><strong>grand union</strong></td>'
                     f'<td class="num"><strong>{len(union):,}</strong></td></tr>')
        breakdown = (f'<h2>By corpus</h2><table>'
                     f'<tr><th>group</th><th class="num">unique hgnc_symbol</th></tr>'
                     f'{grp_rows}</table>')
        per_file_footer = ""
    else:
        sets = [per_file[f][1] for f in files]
        shared = sets[0].intersection(*sets[1:]) if len(sets) > 1 else sets[0]
        per_file_footer = (
            f'<tr><td>shared by all</td><td class="num">&mdash;</td>'
            f'<td class="num">{len(shared):,}</td></tr>'
            f'<tr><td><strong>union</strong></td><td class="num">&mdash;</td>'
            f'<td class="num"><strong>{len(union):,}</strong></td></tr>')
        breakdown = ""

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{esc(h1)}</title>
<style>{css(headline_bg, headline_border)}</style></head><body>
<h1>{h1}</h1>
<p>{intro}</p>
<div class="headline"><span class="big">{len(union):,}</span> unique HGNC symbols
&nbsp;&middot;&nbsp; <span class="big">{total_occ:,}</span> summed occurrences
(symbol-attributed).</div>

<h2>Top {TOPN} symbols by summed occurrences</h2>
<div class="plot">{svg_topn(ranked, bar_color)}</div>
<table>
  <tr><th class="num">#</th><th>hgnc_symbol</th><th class="num">&Sigma; occurrences</th></tr>
  {toprows}
</table>

<h2>Per file</h2>
<table>
  <tr><th>file</th><th class="num">entries</th><th class="num">unique hgnc_symbol</th></tr>
  {prows}{per_file_footer}
</table>
{breakdown}

<h2>All {len(union):,} symbols (&Sigma; occurrences, descending)</h2>
<details><summary>show / hide the full ranked list</summary>
<table><tr><th>hgnc_symbol</th><th class="num">&Sigma; occurrences</th></tr>
{fulllist}</table></details>
</body></html>
"""
    (SUM / out).write_text(doc, encoding="utf-8")
    return len(union), total_occ


def write_nested_keys(files, out):
    """Per-file nested-key summary + cross-file presence matrix + per-file
    match_mode frequencies."""
    key_sections, mm_sections, overview = [], [], []
    for f in files:
        data = load(f)
        n = len(data)
        present, types, example, order, mm = Counter(), defaultdict(set), {}, [], Counter()
        for entry in data.values():
            if not isinstance(entry, dict):
                continue
            for k, v in entry.items():
                if k not in present:
                    order.append(k)
                present[k] += 1
                types[k].add(typename(v))
                example.setdefault(k, v)
            if "match_mode" in entry:
                mm[entry["match_mode"]] += 1
        overview.append((f, n, order))
        rows = "".join(
            f'<tr><td><code>{esc(k)}</code></td><td class="num">{present[k]:,}/{n:,}</td>'
            f'<td><code>{esc(" | ".join(sorted(types[k])))}</code></td>'
            f'<td><code>{esc(short(example[k]))}</code></td></tr>' for k in order)
        key_sections.append(
            f'<h3><code>{esc(f)}</code> &mdash; {n:,} entries, {len(order)} '
            f'nested key(s)</h3><table><tr><th>nested key</th><th class="num">present</th>'
            f'<th>type(s)</th><th>example value</th></tr>{rows}</table>')
        if mm:
            mrows = "".join(
                f'<tr><td><code>{esc(m)}</code></td><td class="num">{c:,}</td>'
                f'<td class="num">{100 * c / n:.1f}%</td></tr>' for m, c in mm.most_common())
            mm_sections.append(
                f'<h3><code>{esc(f)}</code> &mdash; {len(mm)} distinct '
                f'<code>match_mode</code> value(s) over {n:,} entries</h3>'
                f'<table><tr><th>match_mode</th><th class="num">count</th>'
                f'<th class="num">share</th></tr>{mrows}</table>')
        else:
            mm_sections.append(f'<h3><code>{esc(f)}</code></h3>'
                               f'<p>no <code>match_mode</code> key.</p>')

    allkeys = []
    for _, _, ks in overview:
        for k in ks:
            if k not in allkeys:
                allkeys.append(k)
    hdr = "".join(f'<th>{esc(f.replace(".json", ""))}</th>' for f, _, _ in overview)
    matrix_rows = "".join(
        f'<tr><td><code>{esc(k)}</code></td>' +
        "".join('<td class="yes">&#10003;</td>' if k in ks else '<td class="no"></td>'
                for _, _, ks in overview) + '</tr>' for k in allkeys)

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>nested-key &amp; match_mode summary of the HGNC-linkage libraries</title>
<style>{css("#f3f7ff", "#cdddf5")}</style></head><body>
<h1>Nested-key &amp; match_mode summary</h1>
<p>The seven HGNC-linkage libraries are JSON objects keyed by the source entity.
This summarises (1) the <strong>nested keys</strong> inside each entry object and
(2) the <strong><code>match_mode</code></strong> values and their frequency, per
file.</p>

<h2>1 &middot; Nested keys</h2>
<h3>Cross-file presence matrix</h3>
<table><tr><th>nested key</th>{hdr}</tr>{matrix_rows}</table>
{''.join(key_sections)}

<h2>2 &middot; match_mode values &amp; frequency (per file)</h2>
{''.join(mm_sections)}
</body></html>
"""
    (SUM / out).write_text(doc, encoding="utf-8")
    return overview


def main():
    n1, o1 = write_hgnc_report(
        ALL_FILES, "hgnc.html",
        "HGNC symbols &mdash; unique count &amp; summed occurrences",
        "Across all seven HGNC-linkage libraries (<code>roman.json</code>, "
        "<code>roman_ambiguous.json</code>, <code>roman_cosine.json</code>, "
        "<code>greek.json</code>, <code>greek_ambiguous.json</code>, "
        "<code>greek_complex.json</code>, <code>greek_cosine.json</code>): the "
        "distinct values under the nested <code>hgnc_symbol</code> key, and for "
        "each, the sum of the nested <code>occurrences</code> over every entry in "
        "which it appears. List values (ambiguous / complex) and "
        "<code>|</code>-joined cosine values are split into member symbols; an "
        "entry that resolves to several symbols contributes its full occurrence "
        "count to <em>each</em> of them.",
        "#f3f7ff", "#cdddf5", "#3b76d1",
        corpus_groups=[("roman_*.json", ROMAN_FILES), ("greek_*.json", GREEK_FILES)])

    n2, o2 = write_hgnc_report(
        ONE_TO_ONE_FILES, "hgnc_1_to_1.html",
        "HGNC symbols &mdash; 1-to-1 single-gene links",
        "Restricted to the two <strong>single-gene</strong> libraries "
        "<code>roman.json</code> and <code>greek.json</code> (each entry resolves "
        "to exactly one gene, so <code>hgnc_symbol</code> is a single value &mdash; "
        "a true 1-to-1 entity&rarr;gene link). Reports the distinct "
        "<code>hgnc_symbol</code> values and, for each, the sum of the nested "
        "<code>occurrences</code> over every entry mapping to it.",
        "#f1faf3", "#c7e7d2", "#2e8b57")

    ov = write_nested_keys(ALL_FILES, "nested_keys.html")

    print(f"hgnc.html        : {n1:,} unique symbols, {o1:,} summed occ")
    print(f"hgnc_1_to_1.html : {n2:,} unique symbols, {o2:,} summed occ")
    print(f"nested_keys.html : {len(ov)} files summarised")
    for f, n, ks in ov:
        print(f"   {f:24} {n:>6} entries, {len(ks)} nested keys")
    print(f"\nWrote: {SUM / 'hgnc.html'}\n       {SUM / 'hgnc_1_to_1.html'}"
          f"\n       {SUM / 'nested_keys.html'}")


if __name__ == "__main__":
    main()
