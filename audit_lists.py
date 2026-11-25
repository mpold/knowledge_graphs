#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""audit_lists.py -- how often does a listed entity lose an edge its list-mates got?

"PDH1, IDH2, and SDHB gene expression is under expressed in GBM" puts one claim on three
genes, but each pair is scored, normalized and filtered on its own, so some members can reach
the graph and others not. For every list of same-type entities and every other entity X in the
sentence, this finds the lists where at least one member-X pair is in the graph's universe
(high_confidence_g.qualifies_multi, the same merge and cutoff) and at least one is not, and
says where each missing one was lost. It only reads the pipeline output.

A hit is a CANDIDATE: some list members really do differ. "A, but not B" is not read as a list.

Reasons, in pipeline order (the first one that applies):
  not scored          the pair never reached a model (too far apart, no model for it, ...)
  no relation         every model that scored it said "no relation"
  below cutoff        a model named a relation but the score missed RE_MIN_SCORE
  not kept (no log)   the run has no TRIPLES/re_rejected.jsonl.gz: one of the three above
  score < cutoff      kept, but under --score after the merge
  ambiguous symbol    the member maps to several HGNC symbols
  no symbol / id      the member has no HGNC / MONDO / ChEBI / NCIt id
  graph filter        it has an id, but a graph rule drops it (control, generic disease, ...)

Run:  python audit_lists.py [--data-root kaggle_working] [--score 0.5] [--out audit_lists.tsv]
"""

import argparse
import collections
import csv
import glob
import gzip
import hashlib
import json
import sys
from pathlib import Path

import high_confidence_g as H

sys.path.insert(0, str(Path(__file__).resolve().parent / "gpu_bundle"))
import triples as T   # noqa: E402  (entity_text: the surfaces the RE step wrote)

lists_in = H.lists_in


def load_rejected(path):
    """pair_id -> [record] from relation_extraction.py's log, or None when the run has none."""
    if not path.exists():
        return None
    out = collections.defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for ln in fh:
            r = json.loads(ln)
            out[r["pair_id"]].append(r)
    return out


def endpoint_reason(e, flags):
    """Why a normalized endpoint is not a node (None = it is one)."""
    if H.node_of(e, flags):
        return None
    if e.get("type") == "GENETIC" and isinstance(e.get("hgnc_symbol"), list):
        return "ambiguous symbol"
    ids = [e.get(k) for k in ("hgnc_symbol", "mondo_label", "chebi_label", "ncit_label")]
    return "graph filter" if any(ids) else "no symbol / id"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=str(H.DATA_ROOT))
    ap.add_argument("--score", type=float, default=H.GRAPH_BASE)
    ap.add_argument("--merge", choices=H.MERGE_POLICIES, default="union")
    ap.add_argument("--no-list-edges", action="store_true", help="audit the models' own pairs only")
    ap.add_argument("--out", default="audit_lists.tsv")
    ap.add_argument("--top", type=int, default=15, help="surfaces shown per reason")
    a = ap.parse_args()
    H.set_data_root(a.data_root)
    if not H.RE_FILE.exists():
        raise SystemExit(f"ERROR: {H.RE_FILE} not found")
    flags = H.load_type_flags()
    raw = json.loads(H.RE_FILE.read_text(encoding="utf-8"))
    if not raw or not raw[0].get("pair_id"):
        raise SystemExit("ERROR: these triples carry no pair_id (relation_extraction.py is too old)")
    for t in raw:
        H.D.orient(t)
    H.recheck_negation(raw)
    H.K.apply(raw)
    H.orient_drugs(raw, flags)
    merged, note = H.merge_models(raw, a.merge)
    print(note)
    by_pair = {t["pair_id"]: t for t in merged}
    if not a.no_list_edges:
        by_pair.update((t["pair_id"], t) for t in H.list_edges(merged, flags, a.score))
    element = {}                                  # (type, text) -> a normalized endpoint
    for t in raw:
        for e in (t["subject"], t["object"]):
            element.setdefault((e["type"], e["text"]), e)
    rejected = load_rejected(H.OUT_DIR / "re_rejected.jsonl.gz")
    if rejected is None:
        print("NOTE: no re_rejected.jsonl.gz -- pairs the models dropped cannot be explained")

    def key(e):
        return e["label"], T.entity_text(e)

    def node(e):
        el = element.get(key(e))
        return H.node_of(el, flags) if el else None

    def why(pid, m):
        """(reason, score, label) for one pair that is not in the graph."""
        t = by_pair.get(pid)
        if t is None:
            if rejected is None:
                return "not kept (no log)", "", ""
            rs = rejected.get(pid)
            if not rs:
                return "not scored", "", ""
            low = [r for r in rs if r["reason"] == "below_cutoff"]
            if low:
                r = max(low, key=lambda r: r["score"])
                return "below cutoff", r["score"], r["label"]
            return "no relation", max(r["p_any"] for r in rs), ""
        sc, lab = t["score"], t["predicate"]["text"]
        if sc < a.score:
            return "score < cutoff", sc, lab
        e = next((e for e in (t["subject"], t["object"])
                  if (e["type"], e["text"]) == key(m)), None)
        return (endpoint_reason(e, flags) if e else None) or "graph filter", sc, lab

    cases = whole = split = none_on = missing = 0
    reasons, pair_types = collections.Counter(), collections.Counter()
    surfaces = collections.defaultdict(collections.Counter)
    rows, sents = [], set()
    for fp in sorted(glob.glob(str(H.SENT_DIR / "*.json"))):
        try:
            with open(fp, encoding="utf-8") as fh:
                rec = json.load(fh)
        except (OSError, ValueError):
            continue
        pmid = Path(rec.get("source_file", Path(fp).name)).stem
        for s in rec.get("sentences", []):
            text = s.get("text", "")
            ents = sorted((e for e in s.get("entities", []) if e.get("start") is not None),
                          key=lambda e: e["start"])
            if len(ents) < 3:
                continue
            sig = hashlib.md5(text.encode("utf-8", "replace")).hexdigest()[:8]
            for run in lists_in(text, ents):
                ids = {id(e) for e in run}
                for x in ents:
                    if id(x) in ids:
                        continue
                    cases += 1
                    on, off = [], []
                    for m in run:
                        p, q = sorted((m, x), key=lambda e: e["start"])
                        pid = (f'{pmid}:{sig}:{p["start"]}-{p["end"]}'
                               f':{q["start"]}-{q["end"]}')
                        t = by_pair.get(pid)
                        ok = t is not None and H.qualifies_multi(t, a.score, flags)
                        (on if ok else off).append((m, pid))
                    if not on:
                        none_on += 1
                        continue
                    # a second name of a member already on the graph, or of X, is no loss
                    have = {node(m) for m, _ in on} | {node(x)}
                    off = [(m, pid) for m, pid in off
                           if node(m) is None or node(m) not in have]
                    if not off:
                        whole += 1
                        continue
                    split += 1
                    sents.add((pmid, sig))
                    for m, pid in off:
                        r, sc, lab = why(pid, m)
                        missing += 1
                        reasons[r] += 1
                        surfaces[r][T.entity_text(m)] += 1
                        pair_types[f'{m["label"]}-{x["label"]}'] += 1
                        el = element.get(key(m)) or {}
                        sym = (el.get("hgnc_symbol") or el.get("mondo_label")
                               or el.get("chebi_label") or el.get("ncit_label") or "")
                        rows.append([r, T.entity_text(m), m["label"],
                                     sym if isinstance(sym, str) else "|".join(sym),
                                     T.entity_text(x), x["label"], lab, sc,
                                     "; ".join(T.entity_text(o) for o, _ in on), pmid, text])
    with open(a.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reason", "missing", "type", "mapped_to", "with", "with_type", "label", "score",
                    "list_mates_on_graph", "pmid", "sentence"])
        w.writerows(sorted(rows, key=lambda r: (-reasons[r[0]], r[0], r[1])))
    print(f"list x partner cases: {cases:,}  (no member on graph {none_on:,}, "
          f"all on graph {whole:,}, SPLIT {split:,} in {len(sents):,} sentences)")
    print(f"missing member pairs: {missing:,}")
    for r, n in reasons.most_common():
        print(f"  {r:18s} {n:6,}  {100 * n / missing:4.1f}%   "
              + ", ".join(f"{s} x{c}" for s, c in surfaces[r].most_common(a.top)))
    print("  by type pair: " + ", ".join(f"{k} {n:,}" for k, n in pair_types.most_common()))
    print(f"wrote {len(rows):,} row(s) -> {a.out}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
