"""pre_ner_xml_structure.py

Summarize the XML *section types* across the ENTIRE corpus in
``experimental_ner/`` and write a self-contained HTML report to
``summaries/xml_section_types_experimental.html``.

The directory holds TWO XML formats, both of which are parsed:

  * GROBID **TEI** files (``*.grobid.tei.xml``) -- sections are ``<div>``s
    labelled by a ``<head>`` element; structural back-matter is tagged with
    ``<div type="...">``.
  * NLM/PMC **JATS** files (``PMC*.xml``) -- sections are ``<sec>`` elements
    labelled by a child ``<title>``; many carry a ``sec-type`` attribute.

Reporting two complementary notions of "section type":
  1. Section TITLES -- the human-readable label of each section
     (Introduction, Methods, Results, Discussion, Conclusions, ...), pooled
     across both formats. Figure/table caption titles (nested in
     ``<figure>`` / ``<fig>`` / ``<table-wrap>``) are detected and EXCLUDED.
  2. TYPED sections -- explicit machine type attributes: JATS ``sec-type``
     and TEI ``<div type>``.

Run from anywhere (paths resolve relative to this file)::

    python pre_ner_xml_structure.py
"""
import glob
import html
import os
import re
from collections import Counter


# --- paths (resolved relative to this script, not the CWD) ----------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(SCRIPT_DIR, "gpu_bundle", "experimental_ner")  # corpus written here by named_entity_xml.py
OUT_DIR = os.path.join(SCRIPT_DIR, "summaries")
OUT_PATH = os.path.join(OUT_DIR, "pre_ner_xml_structure.html")

# --- regexes --------------------------------------------------------------
# TEI
HEAD_RE = re.compile(r"<head\b[^>]*>(.*?)</head>", re.IGNORECASE | re.DOTALL)
TEI_FIGURE_RE = re.compile(r"<figure\b([^>]*)>(.*?)</figure>",
                           re.IGNORECASE | re.DOTALL)
DIVTYPE_RE = re.compile(r"<div\b[^>]*\btype=\"([^\"]*)\"", re.IGNORECASE)
# JATS
JATS_FIG_RE = re.compile(r"<(fig|table-wrap)\b[^>]*>(.*?)</\1>",
                         re.IGNORECASE | re.DOTALL)
# a <sec ...> opener, optional <label>, then its immediate <title>
SEC_RE = re.compile(
    r"<sec\b([^>]*)>\s*(?:<label[^>]*>.*?</label>\s*)?<title[^>]*>(.*?)</title>",
    re.IGNORECASE | re.DOTALL,
)
SECTYPE_RE = re.compile(r'sec-type="([^"]*)"', re.IGNORECASE)
TITLE_RE = re.compile(r"<title[^>]*>(.*?)</title>", re.IGNORECASE | re.DOTALL)
# generic
TAG_RE = re.compile(r"<[^>]+>")
NUM_PREFIX_RE = re.compile(r"^\s*[\divxlcIVXLC]+(\.[\divxlcIVXLC]+)*[.)]?\s+")
WS_RE = re.compile(r"\s+")


def clean_title(raw):
    """Strip inner tags, unescape entities, drop leading numbering, normalize ws."""
    txt = TAG_RE.sub(" ", raw)
    txt = html.unescape(txt)
    txt = WS_RE.sub(" ", txt).strip()
    txt = NUM_PREFIX_RE.sub("", txt)  # drop "1." / "2.3" / "IV." prefixes
    return txt.strip()


def norm_key(txt):
    """Normalization key for grouping equivalent headings."""
    return txt.lower().strip(" .:;-")


def is_tei(data):
    head = data[:800]
    return "www.tei-c.org" in head or "<TEI" in head


def collect(files):
    """Scan files (both TEI and JATS) and return aggregated stats."""
    title_counter = Counter()   # normalized title -> document-frequency count
    title_display = {}          # normalized title -> representative display form
    typed_counter = Counter()   # "sec-type" / "div type" value -> occurrences
    stats = {
        "total_files": len(files),
        "tei_files": 0,
        "jats_files": 0,
        "files_with_section": 0,
        "total_section_instances": 0,    # genuine section labels (captions excl.)
        "fig_caption_count": 0,          # figure caption titles excluded
        "table_caption_count": 0,        # table caption titles excluded
        "parse_errors": 0,
    }

    for path in files:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                data = fh.read()
        except OSError:
            stats["parse_errors"] += 1
            continue

        seen_keys = set()
        has_section = False

        if is_tei(data):
            stats["tei_files"] += 1

            # Remove <figure> blocks: their <head> is a caption, not a section.
            def _strip_tei_fig(m):
                attrs, inner = m.group(1), m.group(2)
                n = len(HEAD_RE.findall(inner))
                if 'type="table"' in attrs.lower():
                    stats["table_caption_count"] += n
                else:
                    stats["fig_caption_count"] += n
                return " "

            body = TEI_FIGURE_RE.sub(_strip_tei_fig, data)
            labels = [(None, raw) for raw in HEAD_RE.findall(body)]
            for dt in DIVTYPE_RE.findall(data):
                typed_counter[dt.lower()] += 1
        else:
            stats["jats_files"] += 1

            # Remove <fig>/<table-wrap> blocks: their <title> is a caption.
            def _strip_jats_fig(m):
                kind, inner = m.group(1).lower(), m.group(2)
                n = len(TITLE_RE.findall(inner))
                if kind == "table-wrap":
                    stats["table_caption_count"] += n
                else:
                    stats["fig_caption_count"] += n
                return " "

            body = JATS_FIG_RE.sub(_strip_jats_fig, data)
            labels = []
            for attrs, raw in SEC_RE.findall(body):
                labels.append((attrs, raw))
                m = SECTYPE_RE.search(attrs)
                if m:
                    typed_counter[m.group(1).lower()] += 1

        for _attrs, raw in labels:
            cleaned = clean_title(raw)
            if not cleaned:
                continue
            stats["total_section_instances"] += 1
            key = norm_key(cleaned)
            if not key:
                continue
            has_section = True
            if key not in seen_keys:           # count once per document
                seen_keys.add(key)
                title_counter[key] += 1
                if key not in title_display or (
                    cleaned[:1].isupper() and not title_display[key][:1].isupper()
                ):
                    title_display[key] = cleaned

        if has_section:
            stats["files_with_section"] += 1

    stats["distinct_titles"] = len(title_counter)
    return stats, title_counter, title_display, typed_counter


def pct(n, d):
    return f"{(100.0 * n / d):.1f}%" if d else "0%"


def render_html(stats, title_counter, title_display, typed_counter):
    top_titles = title_counter.most_common(60)
    typed_rows = typed_counter.most_common()
    fws = stats["files_with_section"]
    total = stats["total_files"]

    title_rows = "\n".join(
        f"<tr><td class='label'>{html.escape(title_display.get(k, k))}</td>"
        f"<td class='num'>{c:,}</td>"
        f"<td class='num'>{pct(c, fws)}</td></tr>"
        for k, c in top_titles
    )
    typed_html = "\n".join(
        f"<tr><td class='label'>{html.escape(t)}</td>"
        f"<td class='num'>{c:,}</td>"
        f"<td class='num'>{pct(c, total)}</td></tr>"
        for t, c in typed_rows
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>XML Section Types &mdash; experimental_ner</title>
<style>
  :root {{ --bg:#ffffff; --card:#f6f8fa; --fg:#1f2328; --muted:#57606a;
           --accent:#0969da; --border:#d0d7de; }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--bg); color:var(--fg);
          font:15px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif; }}
  .wrap {{ max-width:980px; margin:0 auto; padding:32px 24px 64px; }}
  h1 {{ font-size:26px; margin:0 0 4px; }}
  h2 {{ font-size:19px; margin:36px 0 12px; border-bottom:1px solid var(--border);
        padding-bottom:6px; }}
  p.sub {{ color:var(--muted); margin:0 0 20px; }}
  .cards {{ display:flex; flex-wrap:wrap; gap:14px; margin:18px 0 8px; }}
  .card {{ background:var(--card); border:1px solid var(--border); border-radius:10px;
           padding:14px 18px; min-width:140px; flex:1; }}
  .card .k {{ font-size:26px; font-weight:700; color:var(--accent); }}
  .card .v {{ color:var(--muted); font-size:13px; margin-top:2px; }}
  table {{ width:100%; border-collapse:collapse; background:var(--card);
           border:1px solid var(--border); border-radius:10px; overflow:hidden; }}
  th,td {{ padding:8px 12px; text-align:left; border-bottom:1px solid var(--border); }}
  th {{ background:#eaeef2; color:var(--muted); font-weight:600; font-size:13px;
        text-transform:uppercase; letter-spacing:.03em; }}
  tr:last-child td {{ border-bottom:none; }}
  td.num {{ text-align:right; font-variant-numeric:tabular-nums; color:var(--muted);
            width:110px; }}
  td.label {{ font-weight:500; }}
  tbody tr:hover {{ background:#eaeef2; }}
  .note {{ color:var(--muted); font-size:13px; margin-top:10px; }}
  code {{ background:#eaeef2; padding:1px 5px; border-radius:4px; font-size:13px; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>XML Section Types &mdash; <code>experimental_ner</code></h1>
  <p class="sub">Section structure across the <strong>entire</strong> corpus,
  covering both XML formats present in the directory, produced by
  <code>pre_ner_xml_structure.py</code>.</p>

  <div class="cards">
    <div class="card"><div class="k">{total:,}</div>
      <div class="v">XML files scanned</div></div>
    <div class="card"><div class="k">{stats['jats_files']:,}</div>
      <div class="v">JATS files (&lt;sec&gt;)</div></div>
    <div class="card"><div class="k">{stats['tei_files']:,}</div>
      <div class="v">GROBID TEI files (&lt;head&gt;)</div></div>
    <div class="card"><div class="k">{fws:,}</div>
      <div class="v">files with sections ({pct(fws, total)})</div></div>
    <div class="card"><div class="k">{stats['distinct_titles']:,}</div>
      <div class="v">distinct section titles</div></div>
    <div class="card"><div class="k">{stats['fig_caption_count']:,}</div>
      <div class="v">figure captions (excluded)</div></div>
    <div class="card"><div class="k">{stats['table_caption_count']:,}</div>
      <div class="v">table captions (excluded)</div></div>
  </div>

  <h2>1. Section titles &mdash; top 60 across the whole corpus</h2>
  <p class="note">Section labels pooled from JATS <code>&lt;sec&gt;&lt;title&gt;</code>
  and TEI <code>&lt;div&gt;&lt;head&gt;</code>. Figure/table caption titles
  (nested in <code>&lt;fig&gt;</code>, <code>&lt;table-wrap&gt;</code> or
  <code>&lt;figure&gt;</code>) are <strong>excluded</strong> &mdash;
  {stats['fig_caption_count']:,} figure and {stats['table_caption_count']:,} table
  captions were filtered out. Titles are normalized (case-folded, leading numbering
  such as <code>1.</code> / <code>2.3</code> removed) and counted by
  <em>document frequency</em> &mdash; how many of the {fws:,} files containing any
  section have at least one section with that title.</p>
  <table>
    <thead><tr><th>Section title</th><th class="num">Files</th>
      <th class="num">% of files</th></tr></thead>
    <tbody>
{title_rows}
    </tbody>
  </table>

  <h2>2. Machine-typed sections</h2>
  <p class="note">Explicit type attributes: JATS <code>sec-type</code> and TEI
  <code>&lt;div type&gt;</code>. Counted as total occurrences across all files;
  percentage is of all {total:,} files.</p>
  <table>
    <thead><tr><th>type value</th><th class="num">Occurrences</th>
      <th class="num">% of files</th></tr></thead>
    <tbody>
{typed_html}
    </tbody>
  </table>

  <p class="note" style="margin-top:28px">Generated by
  <code>pre_ner_xml_structure.py</code>. The directory mixes two scholarly XML
  schemas; the script auto-detects each file (TEI vs JATS) and pools their section
  labels, so standard publication sections (Introduction, Methods, Results,
  Discussion, Conclusions, &hellip;) are captured across the full corpus rather than
  the TEI subset only.</p>
</div>
</body>
</html>
"""


def main():
    os.makedirs(INPUT_DIR, exist_ok=True)    # create input dir if missing
    files = glob.glob(os.path.join(INPUT_DIR, "*.xml"))
    if not files:
        raise SystemExit(f"No XML files found in {INPUT_DIR!r}")

    stats, title_counter, title_display, typed_counter = collect(files)
    report = render_html(stats, title_counter, title_display, typed_counter)

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        fh.write(report)

    print(f"files scanned          : {stats['total_files']:,}")
    print(f"  TEI (grobid)         : {stats['tei_files']:,}")
    print(f"  JATS (pmc)           : {stats['jats_files']:,}")
    print(f"files with sections    : {stats['files_with_section']:,} "
          f"({pct(stats['files_with_section'], stats['total_files'])})")
    print(f"section label instances: {stats['total_section_instances']:,}")
    print(f"distinct titles        : {stats['distinct_titles']:,}")
    print(f"figure captions excl.  : {stats['fig_caption_count']:,}")
    print(f"table captions excl.   : {stats['table_caption_count']:,}")
    print(f"typed-section values   : {len(typed_counter)}")
    print(f"read errors            : {stats['parse_errors']}")
    print(f"wrote                  : {OUT_PATH}")


if __name__ == "__main__":
    main()
