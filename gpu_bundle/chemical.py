#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""chemical.py -- normalize CHEMICAL named entities to the ChEBI ontology.

The CHEMICAL analog of the GENETIC (HGNC) and DISEASE (MONDO) pipelines. Reads the
BioBERT CHEMICAL spans from sentences/*.json and maps each to a ChEBI term by exact
equality of a transformed key against ChEBI labels + synonyms.

Inputs
------
  sentences/*.json                  BioBERT NER output (CHEMICAL entities)
  databases/chebi.json              ChEBI ontology (OBO-Graph JSON: graphs[].nodes[])

Outputs (under CHEMICAL/)
-------
  clean_chemical_ne.tsv             aggregated CHEMICAL surfaces + occurrences
  chemical.json                     surfaces linked to exactly ONE ChEBI term
  chemical_ambiguous.json           surfaces linked to >=2 ChEBI terms
  unmatched_chemical.json           surfaces linked to no ChEBI term
  chemical.html                     self-documenting strategy + results report,
                                    incl. a per-file nested-key summary (chemical.json
                                    / chemical_ambiguous.json) and the distinct ChEBI
                                    terms (chebi_id -> chebi_label + summed
                                    occurrences) for chemical.json alone AND
                                    single + ambiguous

================================ STRATEGY ===================================
(0) UPSTREAM. Collect every CHEMICAL entity from sentences/*.json, dash-normalize,
    and sum occurrences per normalized surface (the STAGE 1 roman.py performs).

(0b) ChEBI INDEXING. For each non-deprecated CHEBI_* node, index its label and
    every synonym (hasExact/Related/Broad/NarrowSynonym) under four key folds:
    literal, casefold, hyphen<->space, separator-deletion. Surface key -> the set
    of (ChEBI id, label); field priority label > exact > related > narrow > broad.

(1) MATCH CASCADE, first hit wins (match_mode):
      1. curated synonym       hand map for corpus abbreviations ChEBI lacks
      2. case-sensitive        exact equality against label/synonym surfaces
      3. case-insensitive      casefold equality
      4. hyphen/whitespace     casefold + '-'<->space
      5. normalized            plural-strip (acids->acid) + British->American
                               chemical spelling (sulphate->sulfate,
                               oestradiol->estradiol), then re-match c-i
      6. separator-deletion    casefold + drop hyphens/spaces

(2) ROUTING by distinct ChEBI-term count: 1 -> chemical.json,
    2+ -> chemical_ambiguous.json, 0 -> unmatched_chemical.json
=============================================================================

Run from anywhere (paths resolve relative to this file)::

    python chemical.py
"""

import argparse
import csv
import glob
import html
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(__file__).resolve().parent
SENT_DIR = ROOT / "sentences"
OUT_DIR = ROOT / "CHEMICAL"
CHEBI_PATH = ROOT / "databases" / "chebi.json"
LABEL = "CHEMICAL"

CLEAN_TSV = OUT_DIR / "clean_chemical_ne.tsv"
SINGLE_OUT = OUT_DIR / "chemical.json"
AMBIG_OUT = OUT_DIR / "chemical_ambiguous.json"
DGIDB_PATH = ROOT / "databases" / "interactions.tsv"   # DGIdb (open) drug-gene table
DGIDB_DRUGS_OUT = OUT_DIR / "dgidb_drugs.json"          # DGIdb drugs in corpus NOT in ChEBI
UNMATCHED_OUT = OUT_DIR / "unmatched_chemical.json"
HTML_OUT = OUT_DIR / "chemical.html"

_SYN_FIELD = {"hasExactSynonym": "exact", "hasRelatedSynonym": "related",
              "hasNarrowSynonym": "narrow", "hasBroadSynonym": "broad"}
FIELD_RANK = {"label": 0, "exact": 1, "related": 2, "narrow": 3, "broad": 4}

DASH_VARIANTS = "-‐‑‒–—―−⁃"
_TO_ASCII = str.maketrans({c: "-" for c in DASH_VARIANTS})
_DESPACE = re.compile(r"\s*-\s*")


def dash_normalize(text):
    return _DESPACE.sub("-", (text or "").translate(_TO_ASCII))


def hsfold(s):
    return s.casefold().replace("-", " ")


def delsep(s):
    return s.casefold().replace("-", "").replace(" ", "")


# British -> American chemical spelling folds (ChEBI labels use American)
_UK_US = [(re.compile(p), r) for p, r in [
    (r"sulphn", "sulfn"), (r"sulph", "sulf"), (r"aluminium", "aluminum"),
    (r"oestr", "estr"), (r"caesi", "cesi"), (r"haem", "hem"),
    (r"glutamin", "glutamin"), (r"oxidis", "oxidiz"), (r"ised\b", "ized")]]
_IRREGULAR = {"analyses": "analysis", "indices": "index"}


def britishize(s):
    for rx, rep in _UK_US:
        s = rx.sub(rep, s)
    return s


def singular_candidates(s):
    out = []
    if s in _IRREGULAR:
        out.append(_IRREGULAR[s])
    if s.endswith("ies") and len(s) > 4:
        out.append(s[:-3] + "y")
    if s.endswith("es") and len(s) > 4:
        out.append(s[:-2])
    if s.endswith("s") and not s.endswith("ss") and len(s) > 3:
        out.append(s[:-1])
    return out


# ---- corpus abbreviations ChEBI lacks as a synonym (surface -> ChEBI id
#      'CHEBI:xxxxx' pinned to one term, or a label/synonym string matched c-i),
#      read in this cancer/neuro-oncology corpus's sense. PTX is a fix: ChEBI lists
#      'PTX' as a synonym of palytoxin, but here PTX = paclitaxel. ----
CURATED = {
    "PTX":  "CHEBI:45863",   # paclitaxel  (NOT palytoxin -- corpus fix)
    "DOX":  "CHEBI:28748", "Dox": "CHEBI:28748", "dox": "CHEBI:28748",  # doxorubicin
    "CQ":   "CHEBI:3638",    # chloroquine
    "HCQ":  "CHEBI:5801",    # hydroxychloroquine
    "DTX":  "CHEBI:4672",    # docetaxel (anhydrous)
    "3-MA": "CHEBI:38635",   # 3-methyladenine (autophagy inhibitor)
    "CHX":  "CHEBI:27641",   # cycloheximide
    "BTZ":  "CHEBI:52717",   # bortezomib
    "ATRA": "CHEBI:15367",   # all-trans-retinoic acid (tretinoin)
    "CBD":  "CHEBI:69478",   # cannabidiol
    "ICG":  "CHEBI:31696",   # indocyanine green
    "SFN":  "CHEBI:47807",   # sulforaphane
    "SOR":  "CHEBI:50924",   # sorafenib
    "2-HG": "CHEBI:17084",   # 2-hydroxyglutaric acid (IDH oncometabolite)
    "5-aza-dC": "CHEBI:50131",  # 5-aza-2'-deoxycytidine (decitabine)
    "TFP":  "CHEBI:45951",   # trifluoperazine
    "Ce6":  "CHEBI:168630",  # chlorin e6
}


# ============================================================ (0a) ChEBI gazetteer
# WHO INN stems for drug names the BioBERT chemical NER tends to miss (see useful_lessons.txt).
_DRUG_SUFFIXES = ("mab", "nib", "tinib", "ciclib", "parib", "lisib", "degib", "afenib",
                  "metinib", "rafenib", "platin", "taxel", "rubicin", "tecan", "mustine",
                  "citabine", "arabine", "trexate", "fosfamide", "zomib", "limus", "sartan",
                  "prazole", "statin", "mycin", "floxacin", "conazole", "cycline", "cillin",
                  "dipine", "gliptin", "gliflozin", "zolomide", "semide")
_GAZ_TOKEN = re.compile(r"^[a-z][a-z']{6,}$")        # single token, length >= 7, letters only


def _norm_drug(s):
    """Normalize a DGIdb drug name to a comparison key: lowercase, drop combination
    products, strip radiolabel (' 111in') and biosimilar ('-awwb') suffixes."""
    s = (s or "").strip().lower()
    if not s or s == "null" or "+" in s:
        return ""
    s = re.sub(r"\s+\d+[a-z]+$", "", s)     # radiolabel, e.g. 'bevacizumab 111in'
    s = re.sub(r"-[a-z]{4}$", "", s)          # biosimilar 4-letter code, e.g. 'bevacizumab-awwb'
    return s.strip()


def build_dgidb_drugset(path):
    """Set of clean DGIdb drug names (single-token, alpha, len>=7) from drug_name/
    drug_claim_name -- the surfaces scanned for in corpus text to recover drugs the
    NER misses AND that ChEBI does not contain (e.g. bevacizumab)."""
    p = Path(path)
    if not p.exists():
        return set()
    names = set()
    with open(p, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="	")
        cols = {(c or "").lower(): c for c in (rd.fieldnames or [])}
        dn, dc = cols.get("drug_name"), cols.get("drug_claim_name")
        for row in rd:
            for col in (dn, dc):
                if not col:
                    continue
                k = _norm_drug(row.get(col))
                if re.fullmatch(r"[a-z]{7,}", k):
                    names.add(k)
    return names


def build_gazetteer(IX):
    """Drug-name dictionary seeded from ChEBI surfaces (labels + synonyms), restricted to
    single-token, INN-suffixed names (e.g. -mab, -nib, -platin). Backstops the BioBERT
    chemical NER, which silently misses many long drug names. Returns a set of casefolded
    surfaces; match_cascade() then maps each recovered surface to its ChEBI id, so a hit
    only enters the corpus if it normalizes -- drugs absent from ChEBI cannot be recovered."""
    gaz = set()
    for surf in IX["idx"]:                            # idx keys = ChEBI labels + synonyms
        cf = surf.casefold()
        if _GAZ_TOKEN.match(cf) and cf.endswith(_DRUG_SUFFIXES):
            gaz.add(cf)
    return gaz


# ============================================================ (0) upstream
def stage1(gaz=None, dgidb=None):
    counts = Counter()
    gaz_counts = Counter()                            # surfaces recovered by the ChEBI gazetteer
    dgidb_hits = {}                                   # normalized DGIdb drug -> {occurrences, surfaces}
    n_files = 0
    word_re = re.compile(r"[A-Za-z][A-Za-z']{6,}")    # single tokens, length >= 7
    for fp in sorted(glob.glob(str(SENT_DIR / "*.json"))):
        try:
            rec = json.load(open(fp, encoding="utf-8"))
        except Exception:
            continue
        n_files += 1
        for sent in rec.get("sentences", []):
            tagged = set()
            for ent in sent.get("entities", []):
                if ent.get("label") == LABEL:
                    t = ent.get("text", "")
                    counts[t] += 1
                    tagged.add(dash_normalize(t).casefold())
            if gaz or dgidb:
                seen_g, seen_d = set(), set()
                for m in word_re.finditer(sent.get("text") or ""):
                    w = m.group(0); cf = w.casefold()
                    if gaz and cf in gaz and cf not in seen_g:
                        seen_g.add(cf)
                        if dash_normalize(w).casefold() not in tagged:   # NER already caught it here
                            gaz_counts[w] += 1
                    if dgidb and cf in dgidb and cf not in seen_d:
                        seen_d.add(cf)
                        d = dgidb_hits.setdefault(cf, {"occurrences": 0, "surfaces": set()})
                        d["occurrences"] += 1
                        d["surfaces"].add(w)
    for w, n in gaz_counts.items():                   # fold gazetteer hits into the surface tally
        counts[w] += n
    agg = Counter()
    forms = defaultdict(set)
    for form, n in counts.items():
        key = dash_normalize(form)
        agg[key] += n
        forms[key].add(form)
    rows = sorted(agg.items(), key=lambda kv: (-kv[1], kv[0]))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"{k}\t{n}\t{len(forms[k])}" for k, n in rows)
    CLEAN_TSV.write_text("clean_chemical_ne\toccurrences\tn_source_forms\n"
                         + body + "\n", encoding="utf-8")
    info = {"n_files": n_files, "occ": sum(counts.values()),
            "unique_before": len(counts), "unique_after": len(agg),
            "gaz_surfaces": len(gaz_counts), "gaz_occ": sum(gaz_counts.values()),
            "dgidb_found": len(dgidb_hits)}
    return rows, info, dgidb_hits


# ============================================================ (0b) ChEBI index
def chebi_curie(node_id):
    return node_id.rsplit("/", 1)[-1].replace("_", ":")


def build_indexes(chebi_path):
    docs = json.loads(Path(chebi_path).read_text(encoding="utf-8"))
    nodes = docs["graphs"][0]["nodes"]
    idx, idx_ci, idx_hs, idx_del = ({}, {}, {}, {})
    label_of = {}
    n_terms = 0

    def add(index, key, cid, field):
        slot = index.setdefault(key, {})
        if cid not in slot or FIELD_RANK[field] < FIELD_RANK[slot[cid]]:
            slot[cid] = field

    for nd in nodes:
        nid = nd.get("id", "")
        if "CHEBI_" not in nid:
            continue
        meta = nd.get("meta", {})
        if meta.get("deprecated"):
            continue
        lbl = nd.get("lbl")
        if not lbl:
            continue
        cid = chebi_curie(nid)
        label_of[cid] = lbl
        n_terms += 1
        surfaces = [(lbl, "label")]
        for syn in meta.get("synonyms", []):
            val = syn.get("val")
            fld = _SYN_FIELD.get(syn.get("pred"))
            if val and fld:
                surfaces.append((val, fld))
        for s, field in surfaces:
            add(idx, s, cid, field)
            add(idx_ci, s.casefold(), cid, field)
            add(idx_hs, hsfold(s), cid, field)
            add(idx_del, delsep(s), cid, field)
    return {"idx": idx, "ci": idx_ci, "hs": idx_hs, "del": idx_del,
            "label_of": label_of, "n_terms": n_terms}


def best_field(slot):
    return min((f for f in slot.values()), key=lambda f: FIELD_RANK[f])


def resolve(slot, mode):
    """slot = {chebi_id: field}. Label-preference guard: a surface that hits >=2
    terms but matches exactly ONE via the ChEBI *label* (the rest only via a
    synonym) resolves to that base term -- collapsing ChEBI's base-vs-charge/
    tautomer duplicates (glucose vs D-glucopyranose, ATP vs ATP(4-))."""
    if len(slot) > 1:
        lab = [c for c, f in slot.items() if f == "label"]
        if len(lab) == 1:
            return {lab[0]}, "label", mode
    return set(slot), best_field(slot), mode


# ============================================================ (1) cascade
def match_cascade(value, IX):
    cur = CURATED.get(value)
    if cur is not None:
        if cur.startswith("CHEBI:"):
            if cur in IX["label_of"]:
                return {cur}, "curated", "curated synonym"
        else:
            slot = IX["ci"].get(cur.casefold())
            if slot:
                return resolve(slot, "curated synonym")
    slot = IX["idx"].get(value)
    if slot:
        return resolve(slot, "case-sensitive")
    slot = IX["ci"].get(value.casefold())
    if slot:
        return resolve(slot, "case-insensitive")
    slot = IX["hs"].get(hsfold(value))
    if slot:
        return resolve(slot, "hyphen/whitespace")
    base = britishize(value.casefold())
    cands = [base] + singular_candidates(base)
    if base != value.casefold():
        for c in [value.casefold()] + singular_candidates(value.casefold()):
            if c not in cands:
                cands.append(c)
    for c in cands:
        slot = IX["ci"].get(c)
        if slot:
            return resolve(slot, "normalized")
    slot = IX["del"].get(delsep(value))
    if slot:
        return resolve(slot, "separator-deletion")
    return set(), None, None


_MODES = ["curated synonym", "case-sensitive", "case-insensitive",
          "hyphen/whitespace", "normalized", "separator-deletion"]


def chebi_term_stats(single, ambiguous):
    """Aggregate the linked libraries by ChEBI TERM: distinct chebi_id -> label and
    summed occurrences. A single entry credits its one term; an ambiguous entry
    credits EACH of its candidate terms the entry's full occurrence count. Returns
    (ranked [(chebi_id, label, occ)] desc, total_occ)."""
    occ_by_id = defaultdict(int)
    label_of = {}
    for e in single.values():
        occ_by_id[e["chebi_id"]] += e["occurrences"]
        label_of[e["chebi_id"]] = e["chebi_label"]
    for e in ambiguous.values():
        for cid, lbl in zip(e["chebi_id"], e["chebi_label"]):
            occ_by_id[cid] += e["occurrences"]
            label_of[cid] = lbl
    ranked = sorted(((cid, label_of[cid], occ) for cid, occ in occ_by_id.items()),
                    key=lambda r: (-r[2], r[0]))
    return ranked, sum(occ_by_id.values())


# ============================================================ driver
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chebi", default=str(CHEBI_PATH))
    ap.parse_args()

    IX = build_indexes(CHEBI_PATH)
    print(f"ChEBI    {IX['n_terms']:,} non-deprecated terms indexed "
          f"(label + synonyms, 4 key folds)")
    gaz = build_gazetteer(IX)
    print(f"GAZETTEER {len(gaz):,} ChEBI drug-name surfaces (INN-suffix) for NER backstop")
    dgidb = build_dgidb_drugset(DGIDB_PATH)
    print(f"DGIDB    {len(dgidb):,} clean DGIdb drug names scanned for in corpus text")

    rows, info, dgidb_hits = stage1(gaz, dgidb)
    print(f"STAGE 1  files={info['n_files']:,}  CHEMICAL occ={info['occ']:,}  "
          f"unique {info['unique_before']:,}->{info['unique_after']:,}  "
          f"(+gazetteer: {info['gaz_surfaces']:,} surfaces, {info['gaz_occ']:,} occ) "
          f"-> {CLEAN_TSV}")

    single, ambiguous, unmatched = {}, {}, {}
    mode_counts = Counter()
    occ = {"single": 0, "ambiguous": 0, "unmatched": 0}
    for value, n in rows:
        ids, field, mode = match_cascade(value, IX)
        if not ids:
            unmatched[value] = {"occurrences": n}
            occ["unmatched"] += n
            continue
        cids = sorted(ids)
        labels = [IX["label_of"][c] for c in cids]
        base = {"occurrences": n, "match_field": field, "match_mode": mode}
        if len(cids) == 1:
            single[value] = {**base, "chebi_id": cids[0], "chebi_label": labels[0]}
            occ["single"] += n
        else:
            ambiguous[value] = {**base, "chebi_id": cids, "chebi_label": labels,
                                "n_terms": len(cids)}
            occ["ambiguous"] += n
        mode_counts[mode] += 1

    for path, lib in ((SINGLE_OUT, single), (AMBIG_OUT, ambiguous),
                      (UNMATCHED_OUT, unmatched)):
        path.write_text(json.dumps(lib, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")

    # DGIdb drugs seen in corpus text that are NOT representable in ChEBI (e.g. bevacizumab):
    # these become the DGIdb-only drug layer in target_pharm.py.
    dgidb_only = {}
    for name, d in dgidb_hits.items():
        ids, _, _ = match_cascade(name, IX)
        if ids:                                       # ChEBI represents it -> handled by ChEBI path
            continue
        dgidb_only[name] = {"occurrences": d["occurrences"], "surfaces": sorted(d["surfaces"])}
    dgidb_only = dict(sorted(dgidb_only.items(), key=lambda kv: -kv[1]["occurrences"]))
    DGIDB_DRUGS_OUT.write_text(json.dumps(dgidb_only, ensure_ascii=False, indent=2) + "\n",
                               encoding="utf-8")
    print(f"DGIDB-ONLY in corpus (non-ChEBI) -> {DGIDB_DRUGS_OUT.name}: {len(dgidb_only):,} "
          f"drugs (of {info['dgidb_found']:,} DGIdb names seen)")

    tot = len(single) + len(ambiguous) + len(unmatched)
    print(f"\nMATCH    single -> {SINGLE_OUT.name}: {len(single):,} "
          f"({occ['single']:,} occ)")
    print(f"         ambiguous -> {AMBIG_OUT.name}: {len(ambiguous):,} "
          f"({occ['ambiguous']:,} occ)")
    print(f"         unmatched -> {UNMATCHED_OUT.name}: {len(unmatched):,} "
          f"({occ['unmatched']:,} occ)")
    print(f"         partition: {len(single)} + {len(ambiguous)} + "
          f"{len(unmatched)} = {tot:,} (of {len(rows):,} surfaces)")
    print(f"         by mode: " + ", ".join(f"{m} {mode_counts[m]}" for m in _MODES))
    linked_occ = occ["single"] + occ["ambiguous"]
    print(f"         linked occ: {linked_occ:,} / {info['occ']:,} "
          f"({100*linked_occ/info['occ']:.1f}%)")
    sterm_rank, sterm_occ = chebi_term_stats(single, {})
    print(f"         distinct ChEBI terms (chemical.json, single): {len(sterm_rank):,}; "
          f"summed occurrences: {sterm_occ:,}")
    term_rank, term_occ = chebi_term_stats(single, ambiguous)
    print(f"         distinct ChEBI terms (single+ambiguous): {len(term_rank):,}; "
          f"summed term occurrences: {term_occ:,}")

    render_html(rows, info, IX, single, ambiguous, unmatched, occ, mode_counts)
    print(f"\nWrote {HTML_OUT}")


def _typename(v):
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, int):
        return "int"
    if isinstance(v, float):
        return "float"
    if isinstance(v, str):
        return "str"
    if isinstance(v, list):
        inner = sorted({_typename(x) for x in v}) or ["empty"]
        return f"list[{'|'.join(inner)}]"
    if isinstance(v, dict):
        return "object"
    if v is None:
        return "null"
    return type(v).__name__


def _short(v, n=60):
    s = json.dumps(v, ensure_ascii=False)
    return s if len(s) <= n else s[:n - 1] + "…"


def render_nested_keys(named_libs):
    """HTML summary of the nested keys per library: a cross-file presence matrix
    plus, for each file, every key with its coverage, value type(s) and an example.
    named_libs = [(filename, lib_dict), ...]."""
    esc = html.escape
    perfile = []
    allkeys = []
    for fname, d in named_libs:
        present, types, example, order = Counter(), defaultdict(set), {}, []
        for e in d.values():
            if not isinstance(e, dict):
                continue
            for k, v in e.items():
                if k not in present:
                    order.append(k)
                present[k] += 1
                types[k].add(_typename(v))
                example.setdefault(k, v)
        perfile.append((fname, len(d), order, present, types, example))
        for k in order:
            if k not in allkeys:
                allkeys.append(k)
    hdr = "".join(f'<th>{esc(f)}</th>' for f, _, _, _, _, _ in perfile)
    matrix = "".join(
        '<tr><td><code>' + esc(k) + '</code></td>' + "".join(
            '<td class="yes">&#10003;</td>' if k in pf[2] else '<td class="no"></td>'
            for pf in perfile) + '</tr>' for k in allkeys)
    parts = ["<h2>Nested keys (per file)</h2>",
             "<p>The nested keys inside each entry of <code>chemical.json</code> and "
             "<code>chemical_ambiguous.json</code> &mdash; coverage, value type(s) and "
             "an example. (<code>nonchemical.py</code> later adds a "
             "<code>\"non_chemical\"</code> key to the libraries.)</p>",
             f'<table><tr><th>nested key</th>{hdr}</tr>{matrix}</table>']
    for fname, n, order, present, types, example in perfile:
        rows = "".join(
            f'<tr><td><code>{esc(k)}</code></td>'
            f'<td class="num">{present[k]:,}/{n:,}</td>'
            f'<td><code>{esc(" | ".join(sorted(types[k])))}</code></td>'
            f'<td><code>{esc(_short(example[k]))}</code></td></tr>' for k in order)
        parts.append(
            f'<h3><code>{esc(fname)}</code> &mdash; {n:,} entries, {len(order)} '
            f'nested key(s)</h3><table><tr><th>nested key</th>'
            f'<th class="num">present</th><th>type(s)</th><th>example</th></tr>'
            f'{rows}</table>')
    return "\n".join(parts)


# ============================================================ report
def render_html(rows, info, IX, single, ambiguous, unmatched, occ, mode_counts):
    esc = html.escape
    n_surf = len(rows)
    linked = len(single) + len(ambiguous)
    linked_occ = occ["single"] + occ["ambiguous"]

    def occ_of(d):
        return sorted(d.items(), key=lambda kv: -kv[1]["occurrences"])

    srows = "".join(
        f'<tr><td><code>{esc(v)}</code></td><td class="num">{e["occurrences"]:,}</td>'
        f'<td><code>{esc(e["chebi_id"])}</code></td><td>{esc(e["chebi_label"])}</td>'
        f'<td><code>{esc(e["match_mode"])}</code>/{esc(e["match_field"])}</td></tr>'
        for v, e in occ_of(single)[:25])
    arows = "".join(
        f'<tr><td><code>{esc(v)}</code></td><td class="num">{e["occurrences"]:,}</td>'
        f'<td class="num">{e["n_terms"]}</td>'
        f'<td>{esc(", ".join(e["chebi_label"][:4]))}'
        f'{"&hellip;" if len(e["chebi_label"])>4 else ""}</td></tr>'
        for v, e in occ_of(ambiguous)[:20])
    urows = "".join(
        f'<tr><td><code>{esc(v)}</code></td><td class="num">{e["occurrences"]:,}</td></tr>'
        for v, e in occ_of(unmatched)[:30])
    mrows = "".join(
        f'<tr><td><code>{esc(m)}</code></td><td class="num">{mode_counts[m]:,}</td></tr>'
        for m in _MODES if mode_counts[m])

    # distinct ChEBI terms (single + ambiguous): id -> label, summed occurrences
    term_rank, term_occ = chebi_term_stats(single, ambiguous)
    trows = "".join(
        f'<tr><td><code>{esc(cid)}</code></td><td>{esc(lbl)}</td>'
        f'<td class="num">{o:,}</td></tr>' for cid, lbl, o in term_rank)
    # distinct ChEBI terms in chemical.json alone (single matches only)
    sterm_rank, sterm_occ = chebi_term_stats(single, {})
    strows = "".join(
        f'<tr><td><code>{esc(cid)}</code></td><td>{esc(lbl)}</td>'
        f'<td class="num">{o:,}</td></tr>' for cid, lbl, o in sterm_rank)
    nested_section = render_nested_keys([("chemical.json", single),
                                         ("chemical_ambiguous.json", ambiguous)])

    style = (
        " body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;"
        "margin:2rem auto;padding:0 1rem;color:#1a1a1a;}"
        " h1{font-size:1.45rem;} h2{font-size:1.15rem;margin-top:1.8rem;"
        "border-bottom:1px solid #ddd;padding-bottom:.3rem;}"
        " table{border-collapse:collapse;margin:.6rem 0;font-size:.9em;}"
        " th,td{border:1px solid #ccc;padding:.28rem .6rem;text-align:left;vertical-align:top;}"
        " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
        " .big{font-size:2rem;font-weight:700;}"
        " .headline{background:#f1fbf6;border:1px solid #cdeede;border-radius:8px;"
        "padding:.8rem 1rem;margin:1rem 0;}"
        " code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.9em;}"
        " .yes{text-align:center;color:#1a7f37;font-weight:700;} .no{text-align:center;background:#fafafa;}"
        " details{margin:.6rem 0;} summary{cursor:pointer;color:#357;}")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>chemical &mdash; CHEMICAL entity normalization to ChEBI</title>
<style>{style}</style></head><body>
<h1>chemical &mdash; CHEMICAL entity &rarr; ChEBI normalization</h1>
<p>The CHEMICAL analog of the GENETIC (HGNC) and DISEASE (MONDO) pipelines:
BioBERT CHEMICAL spans from <code>sentences/*.json</code> linked to the
<strong>ChEBI</strong> ontology by exact equality of a transformed key against
ChEBI labels and synonyms. Produced by <code>chemical.py</code>.</p>
<div class="headline"><span class="big">{linked:,}</span> of {n_surf:,} distinct
CHEMICAL surfaces link to ChEBI ({100*linked/n_surf:.1f}%);
<span class="big">{100*linked_occ/info['occ']:.1f}%</span> of the
{info['occ']:,} CHEMICAL occurrences are linked.
Indexed against {IX['n_terms']:,} ChEBI terms.</div>

<h2>Buckets</h2>
<table><tr><th>bucket</th><th class="num">surfaces</th><th class="num">occurrences</th></tr>
<tr><td>single ChEBI term &rarr; <code>chemical.json</code></td><td class="num">{len(single):,}</td><td class="num">{occ['single']:,}</td></tr>
<tr><td>ambiguous (&ge;2 terms) &rarr; <code>chemical_ambiguous.json</code></td><td class="num">{len(ambiguous):,}</td><td class="num">{occ['ambiguous']:,}</td></tr>
<tr><td>unmatched &rarr; <code>unmatched_chemical.json</code></td><td class="num">{len(unmatched):,}</td><td class="num">{occ['unmatched']:,}</td></tr>
<tr><td><strong>total</strong></td><td class="num"><strong>{n_surf:,}</strong></td><td class="num"><strong>{info['occ']:,}</strong></td></tr>
</table>

<h2>Yield by matching pass</h2>
<table><tr><th>match_mode</th><th class="num">surfaces</th></tr>{mrows}</table>

<h2>Top single-term matches</h2>
<table><tr><th>surface</th><th class="num">occ</th><th>ChEBI id</th><th>ChEBI label</th><th>via</th></tr>
{srows}</table>

<h2>Top ambiguous matches</h2>
<table><tr><th>surface</th><th class="num">occ</th><th class="num">#terms</th><th>candidate labels</th></tr>
{arows}</table>

<h2>Top unmatched surfaces</h2>
<p>High-frequency CHEMICAL surfaces with no ChEBI label/synonym match
(abbreviations ChEBI lacks, non-compound tokens, or descriptors) &mdash;
candidates for curation.</p>
<table><tr><th>surface</th><th class="num">occ</th></tr>{urows}</table>

{nested_section}

<h2>Distinct ChEBI terms in <code>chemical.json</code> (single matches)</h2>
<p>The distinct <code>chebi_id</code> values in <code>chemical.json</code> (the
single-term matches), each mapped to its <code>chebi_label</code> and the sum of the
<code>occurrences</code> of every surface that resolves to it (several surfaces can
share one term, e.g. <code>TMZ</code> and <code>temozolomide</code>).</p>
<div class="headline"><span class="big">{len(sterm_rank):,}</span> distinct ChEBI
terms &nbsp;&middot;&nbsp; <span class="big">{sterm_occ:,}</span> summed occurrences.</div>
<details><summary>show / hide all {len(sterm_rank):,} ChEBI terms (by &Sigma; occ)</summary>
<table><tr><th>chebi_id</th><th>chebi_label</th><th class="num">&Sigma; occ</th></tr>
{strows}</table></details>

<h2>Distinct ChEBI terms linked (single + ambiguous)</h2>
<p>The distinct <code>chebi_id</code> values across <code>chemical.json</code> and
<code>chemical_ambiguous.json</code>, each mapped to its <code>chebi_label</code>
and the sum of the <code>occurrences</code> of every surface that resolves to it
(an ambiguous surface credits its full count to <em>each</em> candidate term).</p>
<div class="headline"><span class="big">{len(term_rank):,}</span> distinct ChEBI
terms &nbsp;&middot;&nbsp; <span class="big">{term_occ:,}</span> summed occurrences
(term-attributed).</div>
<details><summary>show / hide all {len(term_rank):,} ChEBI terms (by &Sigma; occ)</summary>
<table><tr><th>chebi_id</th><th>chebi_label</th><th class="num">&Sigma; occ</th></tr>
{trows}</table></details>

<h2>Reproduce</h2><pre><code>python chemical.py</code></pre>
</body></html>
"""
    HTML_OUT.write_text(doc, encoding="utf-8")


if __name__ == "__main__":
    main()
