#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""target_pharm.py -- ChEBI chemicals that target a specific protein / gene,
and the cross-link of the corpus chemicals to their gene targets / GENETIC entities.

ChEBI has no gene entities; a chemical's molecular target is expressed through its
ROLES (the `has role` relation, RO_0000087) -- e.g. "EC 2.7.10.1 (receptor
protein-tyrosine kinase) inhibitor", "glucagon receptor antagonist".

Two things are produced:
(A) the set of ChEBI molecular entities whose roles name a specific protein target;
(B) the cross-link of the CORPUS chemicals (those in CHEMICAL/chemical.json +
    chemical_ambiguous.json that are also protein-targeting) to HGNC gene targets
    and on to the corpus GENETIC entities, emitted as the INVERSE view
    gene -> chemicals that target it.

Chain (B): CHEMICAL surface -> chebi_id -> ChEBI target role -> protein-target name
           -> HGNC gene symbol -> corpus GENETIC entity surface(s).

Inputs
------
  databases/chebi.json                         ChEBI (OBO-Graph JSON)
  databases/hgnc_complete_set_2026-05-01.json  HGNC (protein-name -> gene)
  CHEMICAL/chemical.json, chemical_ambiguous.json   corpus chemical surfaces/ids
  GENETIC/roman.json, roman_ambiguous.json, greek.json,
  GENETIC/greek_ambiguous.json, greek_complex.json  corpus GENETIC genes/surfaces

Outputs (under CHEMICAL/)
-------
  target_pharm.json        chebi_id -> {label, target_roles, classes, actions,
                           hgnc_targets:[{hgnc_symbol, via_roles}]}  (the specific
                           HGNC gene(s) each chemical's roles map to)
  chemical_to_target.json  INVERSE cross-link: gene -> {in_corpus_GENETIC,
                           genetic_surfaces, n_chemicals, chemicals[...]}
  target_pharm.html        summary (target classes/actions/roles + the inverse view)

Definitions / caveats
---------------------
A role is a "specific-protein target" role when its label ends in a pharmacological
action word (inhibitor/agonist/antagonist/activator/modulator/blocker/...), names a
protein (an "EC <number> (<enzyme>) ..." role or a protein-class word: receptor,
kinase, channel, transporter, synthase, polymerase, protease, ...), and is NOT a
pathway/process role (X biosynthesis/synthesis/signalling/uptake/... inhibitor).
The target name (enzyme inside an EC role, else role minus action word) is mapped to
HGNC by exact casefold / separator-deletion over the six HGNC fields -- high
precision, partial recall (viral/bacterial targets, broad EC classes, and generic
names do not resolve to a human gene). The GENETIC-surface side reflects the corpus's
own gene normalization.

Run from anywhere::  python target_pharm.py
"""

import html
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CHEM = ROOT / "CHEMICAL"
GEN = ROOT / "GENETIC"
CHEBI_PATH = ROOT / "databases" / "chebi.json"
HGNC_PATH = ROOT / "databases" / "hgnc_complete_set_2026-05-01.json"
DGIDB_PATH = ROOT / "databases" / "interactions.tsv"   # optional: DGIdb interactions TSV (open drug-gene targets; typed interactions)
DGIDB_DRUGS_PATH = CHEM / "dgidb_drugs.json"           # DGIdb drugs in corpus but NOT in ChEBI (from chemical.py)
JSON_OUT = CHEM / "target_pharm.json"
INV_OUT = CHEM / "chemical_to_target.json"
HTML_OUT = CHEM / "target_pharm.html"

HAS_ROLE = "http://purl.obolibrary.org/obo/RO_0000087"
HGNC_FIELDS = ("symbol", "alias_symbol", "prev_symbol", "name", "alias_name", "prev_name")

ACTION = re.compile(
    r"(inhibitor|inverse agonist|partial agonist|agonist|antagonist|activator|"
    r"modulator|blocker|opener|potentiator|stabiliser|stabilizer)$", re.I)
ACTION_STRIP = re.compile(
    r"\s*(inhibitor|inverse agonist|partial agonist|agonist|antagonist|activator|"
    r"modulator|blocker|opener|potentiator|stabiliser|stabilizer)$", re.I)
PROTEIN = re.compile(
    r"(^EC |receptor|kinase|channel|transporter|synthase|synthetase|polymerase|"
    r"reductase|oxidase|transaminase|aminotransferase|transferase|phosphatase|"
    r"protease|peptidase|dehydrogenase|isomerase|ligase|hydrolase|\blyase\b|"
    r"deacetylase|demethylase|methyltransferase|acetyltransferase|cyclase|esterase|"
    r"gtpase|atpase|\bpump\b|exchanger|symporter|antiporter|dioxygenase|"
    r"monooxygenase|hydroxylase|carboxylase|decarboxylase|topoisomerase|gyrase|"
    r"integrase|transcriptase|convertase|aromatase|elastase|lipase|nuclease|"
    r"helicase|sirtuin|caspase|cathepsin|telomerase)", re.I)
PROCESS = re.compile(
    r"biosynthesis|synthesis|production|signal|pathway|uptake|secretion|release|"
    r"aggregation|formation|metabolism|\btransport\b|replication|\brepair\b|"
    r"assembly|\bfusion\b|biosynthetic", re.I)

_ENZYME_WORDS = ("synthase", "synthetase", "polymerase", "reductase", "oxidase",
                 "transaminase", "aminotransferase", "transferase", "phosphatase",
                 "protease", "peptidase", "dehydrogenase", "isomerase", "ligase",
                 "hydrolase", "lyase", "deacetylase", "demethylase",
                 "methyltransferase", "cyclase", "esterase", "dioxygenase",
                 "monooxygenase", "hydroxylase", "carboxylase", "decarboxylase",
                 "topoisomerase", "gyrase", "integrase", "transcriptase",
                 "convertase", "aromatase", "elastase", "lipase", "nuclease",
                 "helicase", "sirtuin", "caspase", "cathepsin", "telomerase")
_ACTION_WORDS = ["inverse agonist", "partial agonist", "antagonist", "agonist",
                 "inhibitor", "activator", "modulator", "blocker", "opener",
                 "potentiator", "stabiliser", "stabilizer"]


def is_target_role(role):
    return (bool(ACTION.search(role)) and not PROCESS.search(role)
            and (role.startswith("EC ") or bool(PROTEIN.search(role))))


def target_class(role):
    r = role.lower()
    if role.startswith("EC "):
        return "EC enzyme"
    if "receptor" in r:
        return "receptor"
    if "kinase" in r:
        return "kinase"
    if "channel" in r:
        return "ion channel"
    if any(w in r for w in ("transporter", "pump", "exchanger", "symporter", "antiporter")):
        return "transporter"
    if any(w in r for w in _ENZYME_WORDS):
        return "other enzyme"
    return "other protein"


def action_of(role):
    rl = role.lower()
    for a in _ACTION_WORDS:
        if rl.endswith(a):
            return a
    return "other"


def target_names(role):
    """The protein-target name(s) a role denotes (for HGNC mapping)."""
    if role.startswith("EC "):
        m = re.search(r"[\(\[](.+)[\)\]]", role)        # enzyme name in (...) / [...]
        return [m.group(1).strip()] if m else []
    return [ACTION_STRIP.sub("", role).strip()]         # role minus its action word


def delsep(s):
    return s.casefold().replace("-", "").replace(" ", "")


def curie(node_id):
    return node_id.rsplit("/", 1)[-1].replace("_", ":")


def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def _norm_drug(s):
    """Normalize a DGIdb drug name to a key: lowercase, drop combos, strip radiolabel
    (' 111in') and biosimilar ('-awwb') suffixes. Keeps biosimilars/labelled forms with
    their parent (BEVACIZUMAB-AWWB -> bevacizumab) so typed rows are reachable by INN."""
    s = (s or "").strip().lower()
    if not s or s == "null" or "+" in s:
        return ""
    s = re.sub(r"\s+\d+[a-z]+$", "", s)
    s = re.sub(r"-[a-z]{4}$", "", s)
    return s.strip()


def load_dgidb(path, map_gene):
    """DGIdb interactions TSV -> {drug_name_lower: {hgnc_symbol: set(interaction_types)}}.

    ChEBI `has role` only covers small molecules; DGIdb is an open, gene-centric
    drug-gene interaction resource that recovers antibody/biologic targets
    (e.g. nivolumab -> PDCD1). Column names vary by DGIdb release, so the header is
    matched by name. An interaction is kept iff its `interaction_type` is
    non-NULL (target specificity); each is tagged by drug class --
    "DGIdb-antineoplastic" (approved & anti_neoplastic), "DGIdb-approved"
    (approved only), else "DGIdb-investigational" -- so the graph can colour
    green / amber / other. Returns {} if the file is absent (provider no-ops)."""
    import csv
    path = Path(path)
    if not path.exists():
        print(f"(A2) DGIdb: {path.name} not found in databases/ -- skipping DGIdb targets.")
        return {}
    out = defaultdict(lambda: defaultdict(set))    # name_lower -> gene -> set(DGIdb class tags)
    nrows = nlinked = 0
    with open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = {(c or "").lower(): c for c in (rd.fieldnames or [])}

        def pick(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        dncol = pick("drug_name")
        dccol = pick("drug_claim_name", "drug_claim_primary_name")
        gcol = pick("gene_name", "gene_claim_name")
        tcol = pick("interaction_type", "interaction_types")
        acol = pick("approved")
        ncol = pick("anti_neoplastic", "antineoplastic")
        if not ((dncol or dccol) and gcol):
            print(f"(A2) DGIdb: {path.name} has no drug/gene columns -- skipping.")
            return {}
        miss = [nm for nm, c in (("approved", acol), ("anti_neoplastic", ncol),
                                 ("interaction_type", tcol)) if c is None]
        if miss:
            print(f"(A2) DGIdb: column(s) absent: {', '.join(miss)} -- affected classes fall back to 'investigational'.")
        ncat = Counter()
        for row in rd:
            nrows += 1
            gene = (row.get(gcol) or "").strip()
            keys = {_norm_drug(row.get(dncol)), _norm_drug(row.get(dccol))} - {""}
            if not keys or not gene:
                continue
            raw = (row.get(tcol) or "").strip()          # target specificity: a defined interaction_type
            if not raw or raw.upper() == "NULL":
                continue
            approved = (row.get(acol) or "").strip().lower() == "true"
            antineo = (row.get(ncol) or "").strip().lower() == "true"
            cat = "antineoplastic" if (approved and antineo) else ("approved" if approved else "investigational")
            itypes = {t.strip().lower() for t in raw.replace("|", ",").split(",")
                      if t.strip() and t.strip().lower() not in ("n/a", "na", "none", "null")}
            tag = "DGIdb-" + cat + ((": " + ", ".join(sorted(itypes))) if itypes else "")
            for sym in (map_gene(gene) or {gene}):
                for key in keys:
                    out[key][sym].add(tag)
            ncat[cat] += 1
            nlinked += 1
    print(f"(A2) DGIdb: read {nrows:,} rows; {nlinked:,} typed drug-gene links "
          f"(antineoplastic={ncat['antineoplastic']:,}, approved={ncat['approved']:,}, "
          f"investigational={ncat['investigational']:,}) -> {len(out):,} drug names indexed")
    return {k: {g: a for g, a in v.items()} for k, v in out.items()}




# ============================================================ driver
def main():
    g = load(CHEBI_PATH)["graphs"][0]
    label = {n["id"]: (n.get("lbl") or "") for n in g["nodes"]}
    chem_roles = defaultdict(set)
    for e in g["edges"]:
        if e.get("pred") == HAS_ROLE:
            chem_roles[e["sub"]].add(label.get(e["obj"], ""))

    # ---- (A) protein-targeting chemicals -> target_pharm.json --------------
    n_with_roles = len(chem_roles)
    n_action = sum(1 for rs in chem_roles.values() if any(ACTION.search(r) for r in rs))
    tlib = {}
    for cid, roles in chem_roles.items():
        troles = sorted(r for r in roles if is_target_role(r))
        if troles:
            tlib[curie(cid)] = {
                "chebi_label": label.get(cid, ""),
                "n_target_roles": len(troles),
                "target_roles": troles,
                "target_classes": sorted({target_class(r) for r in troles}),
                "actions": sorted({action_of(r) for r in troles}),
            }
    tlib = dict(sorted(tlib.items(), key=lambda kv: kv[1]["chebi_label"].casefold()))

    # ---- HGNC index (protein-target name -> approved gene symbol) ----------
    docs = load(HGNC_PATH)["response"]["docs"]
    idx_ci, idx_del = defaultdict(set), defaultdict(set)
    for d in docs:
        sym = d.get("symbol")
        if not sym:
            continue
        for f in HGNC_FIELDS:
            v = d.get(f)
            for s in ([v] if isinstance(v, str) else (v or [])):
                if isinstance(s, str) and s:
                    idx_ci[s.casefold()].add(sym)
                    idx_del[delsep(s)].add(sym)

    def map_gene(name):
        return idx_ci.get(name.casefold()) or idx_del.get(delsep(name)) or set()

    # ---- annotate EVERY protein-targeting chemical with its HGNC gene target(s),
    #      derived from its ChEBI target roles -> written into target_pharm.json
    chem_genes = {}                              # chebi_id -> {gene: set(roles)}
    for cid, e in tlib.items():
        gr = defaultdict(set)
        for role in e["target_roles"]:
            for nm in target_names(role):
                for gene in map_gene(nm):
                    gr[gene].add(role)
        chem_genes[cid] = gr
        e["hgnc_targets"] = [{"hgnc_symbol": g, "via_roles": sorted(rs), "source": "chebi"}
                             for g, rs in sorted(gr.items())]
        e["n_hgnc_targets"] = len(gr)

    # ---- corpus CHEMICAL surfaces per chebi_id (used by the DGIdb merge + cross-link)
    chebi_surf = defaultdict(set)              # chebi_id -> corpus CHEMICAL surfaces
    for surface, e in load(CHEM / "chemical.json").items():
        chebi_surf[e["chebi_id"]].add(surface)
    for surface, e in load(CHEM / "chemical_ambiguous.json").items():
        for cid in e["chebi_id"]:
            chebi_surf[cid].add(surface)

    # ---- (A2) DGIdb target provider: add targets for corpus chemicals that ChEBI
    #      has no `has role` annotation for (antibodies/biologics, etc.) ----
    clabel = {curie(nid): lbl for nid, lbl in label.items()}
    db = load_dgidb(DGIDB_PATH, map_gene)              # {} when the file is absent
    n_db_chem = 0
    if db:
        def src_of(rs):
            has_db = any(r.startswith("DGIdb") for r in rs)
            has_ch = any(not r.startswith("DGIdb") for r in rs)
            return "chebi+dgidb" if (has_db and has_ch) else ("dgidb" if has_db else "chebi")
        for cid in list(chebi_surf):
            names = {(clabel.get(cid) or "").casefold()} | {x.casefold() for x in chebi_surf[cid]}
            names.discard("")
            hits = defaultdict(set)                    # gene -> set(DGIdb class tags)
            for nm in names:
                for gene, tags in db.get(nm, {}).items():
                    hits[gene].update(tags)
            if not hits:
                continue
            n_db_chem += 1
            ent = tlib.get(cid)
            if ent is None:                            # DGIdb-only (no ChEBI target role)
                ent = tlib[cid] = {"chebi_label": clabel.get(cid, ""), "n_target_roles": 0,
                                   "target_roles": [], "target_classes": [], "actions": []}
                chem_genes[cid] = defaultdict(set)
            gr = chem_genes[cid]
            for gene, tags in hits.items():
                gr[gene].update(tags)
            allacts = sorted({t for tags in hits.values() for t in tags})
            ent["actions"] = sorted(set(ent.get("actions", [])) | set(allacts))
            ent["hgnc_targets"] = [{"hgnc_symbol": gx, "via_roles": sorted(rs), "source": src_of(rs)}
                                   for gx, rs in sorted(gr.items())]
            ent["n_hgnc_targets"] = len(gr)
        tlib = dict(sorted(tlib.items(), key=lambda kv: kv[1]["chebi_label"].casefold()))
        print(f"(A2) DGIdb: added targets for {n_db_chem:,} corpus chemicals")

    # ---- (A3) DGIdb-only drug layer: drugs found in the corpus text that ChEBI does
    #      NOT contain (e.g. bevacizumab) -- detected by chemical.py -> dgidb_drugs.json.
    #      Keyed by a synthetic "DGIDB:<name>" id so they flow through the (B) cross-link. ----
    dd = load(DGIDB_DRUGS_PATH)                         # {name: {occurrences, surfaces}} ({} if absent)
    n_dd = 0
    if dd and db:
        corpus_chem_names = {(clabel.get(c) or "").casefold() for c in chebi_surf}
        corpus_chem_names |= {x.casefold() for ss in chebi_surf.values() for x in ss}
        corpus_chem_names.discard("")
        for name, meta in dd.items():
            key = name.casefold()
            if key in corpus_chem_names:                # already a ChEBI corpus chemical -> A2 handled it
                continue
            genes = db.get(key)                         # {gene: set(class tags)} (typed only)
            if not genes:
                continue
            sid = "DGIDB:" + key
            chebi_surf[sid] = set(meta.get("surfaces") or [name])
            gr = chem_genes[sid] = defaultdict(set)
            for gene, tags in genes.items():
                gr[gene].update(tags)
            ent = tlib[sid] = {
                "chebi_label": name, "n_target_roles": 0, "target_roles": [],
                "target_classes": [], "actions": sorted({t for ts in genes.values() for t in ts}),
                "hgnc_targets": [{"hgnc_symbol": gx, "via_roles": sorted(rs), "source": "dgidb"}
                                 for gx, rs in sorted(gr.items())],
                "n_hgnc_targets": len(gr)}
            n_dd += 1
        tlib = dict(sorted(tlib.items(), key=lambda kv: kv[1]["chebi_label"].casefold()))
        print(f"(A3) DGIdb-only drugs (non-ChEBI, in corpus): added {n_dd:,}")

    CHEM.mkdir(parents=True, exist_ok=True)
    JSON_OUT.write_text(json.dumps(tlib, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")

    by_class, by_action, role_chem = Counter(), Counter(), Counter()
    for e in tlib.values():
        for c in e["target_classes"]:
            by_class[c] += 1
        for a in e["actions"]:
            by_action[a] += 1
        for r in e["target_roles"]:
            role_chem[r] += 1
    n_chem_gene = sum(1 for e in tlib.values() if e["hgnc_targets"])

    print(f"ChEBI: {n_with_roles:,} chemicals with a has_role link; "
          f"{n_action:,} with an action-type role.")
    print(f"(A) protein-targeting chemicals: {len(tlib):,} "
          f"({n_chem_gene:,} annotated with >=1 HGNC gene target) -> {JSON_OUT.name}")

    # ---- (B) cross-link corpus chemicals -> genes -> GENETIC entities ------
    corpus_ids = [c for c in tlib if c in chebi_surf]   # corpus & (ChEBI- or DGIdb-)targeting set

    gene_surf = defaultdict(set)                        # corpus GENETIC gene -> surfaces
    def add_genetic(path):
        for surface, e in load(path).items():
            hs = e.get("hgnc_symbol")
            for sym in ([hs] if isinstance(hs, str) else (hs or [])):
                gene_surf[sym].add(surface)
    for p in ("roman.json", "roman_ambiguous.json", "greek.json",
              "greek_ambiguous.json", "greek_complex.json"):
        add_genetic(GEN / p)
    corpus_genes = set(gene_surf)

    # gene -> corpus chemicals that target it (INVERSE view), reusing chem_genes
    inv = defaultdict(dict)                              # gene -> {chebi_id -> roles}
    for cid in corpus_ids:
        for gene, roles in chem_genes[cid].items():
            inv[gene][cid] = roles

    inverse = {}
    for gene, d in inv.items():
        chems = [{
            "chebi_id": cid,
            "chebi_label": tlib[cid]["chebi_label"],
            "chemical_surfaces": sorted(chebi_surf[cid]),
            "via_roles": sorted(roles),
        } for cid, roles in sorted(d.items(),
                                   key=lambda kv: tlib[kv[0]]["chebi_label"].casefold())]
        inverse[gene] = {
            "in_corpus_GENETIC": gene in corpus_genes,
            "genetic_surfaces": sorted(gene_surf.get(gene, ())),
            "n_chemicals": len(chems),
            "chemicals": chems,
        }
    # sort: corpus GENETIC genes first, then by #chemicals, then symbol
    inverse = dict(sorted(inverse.items(),
                          key=lambda kv: (not kv[1]["in_corpus_GENETIC"],
                                          -kv[1]["n_chemicals"], kv[0])))
    INV_OUT.write_text(json.dumps(inverse, ensure_ascii=False, indent=2) + "\n",
                       encoding="utf-8")

    n_genes = len(inverse)
    n_corpus_genes = sum(1 for v in inverse.values() if v["in_corpus_GENETIC"])
    n_chems_linked = len({c for v in inverse.values() for c in
                          (x["chebi_id"] for x in v["chemicals"])})
    print(f"(B) cross-link: {len(corpus_ids):,} corpus protein-targeting chemicals; "
          f"{n_chems_linked} map to a gene; {n_genes} target genes "
          f"({n_corpus_genes} are corpus GENETIC entities) -> {INV_OUT.name}")

    render_html(tlib, n_with_roles, n_action, by_class, by_action, role_chem,
                inverse, len(corpus_ids), n_chems_linked, n_genes, n_corpus_genes)
    print(f"Wrote {HTML_OUT}")


def render_html(tlib, n_with_roles, n_action, by_class, by_action, role_chem,
                inverse, n_corpus_chem, n_chems_linked, n_genes, n_corpus_genes):
    esc = html.escape
    n = len(tlib)

    def table(title, counter, head):
        rows = "".join(f'<tr><td>{esc(str(k))}</td><td class="num">{v:,}</td></tr>'
                       for k, v in counter.most_common())
        return (f'<h2>{title}</h2><table><tr><th>{head}</th>'
                f'<th class="num">chemicals</th></tr>{rows}</table>')

    toprole = "".join(f'<tr><td>{esc(r)}</td><td class="num">{c:,}</td></tr>'
                      for r, c in role_chem.most_common(30))
    full = "".join(
        f'<tr><td><code>{esc(cid)}</code></td><td>{esc(e["chebi_label"])}</td>'
        f'<td>{esc(" | ".join(e["target_roles"]))}</td></tr>'
        for cid, e in tlib.items())

    # inverse view: the corpus-GENETIC genes (the cross-link to GENETIC entities)
    corpus_rows = []
    for gene, v in inverse.items():
        if not v["in_corpus_GENETIC"]:
            continue
        chems = ", ".join(esc(c["chebi_label"]) for c in v["chemicals"][:6])
        if v["n_chemicals"] > 6:
            chems += "&hellip;"
        sf = esc(", ".join(v["genetic_surfaces"][:4]))
        corpus_rows.append(
            f'<tr><td><strong>{esc(gene)}</strong></td>'
            f'<td class="sf">{sf}{"&hellip;" if len(v["genetic_surfaces"])>4 else ""}</td>'
            f'<td class="num">{v["n_chemicals"]}</td><td>{chems}</td></tr>')
    inv_full = "".join(
        f'<tr><td>{esc(gene)}{"*" if v["in_corpus_GENETIC"] else ""}</td>'
        f'<td class="num">{v["n_chemicals"]}</td>'
        f'<td>{esc(", ".join(c["chebi_label"] for c in v["chemicals"][:8]))}'
        f'{"&hellip;" if v["n_chemicals"]>8 else ""}</td></tr>'
        for gene, v in inverse.items())

    style = (
        " body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1040px;"
        "margin:2rem auto;padding:0 1rem;color:#1a1a1a;}"
        " h1{font-size:1.45rem;} h2{font-size:1.15rem;margin-top:1.8rem;"
        "border-bottom:1px solid #ddd;padding-bottom:.3rem;}"
        " table{border-collapse:collapse;margin:.6rem 0;font-size:.88em;}"
        " th,td{border:1px solid #ccc;padding:.28rem .6rem;text-align:left;vertical-align:top;}"
        " th{background:#f7f7f7;} .num{text-align:right;font-variant-numeric:tabular-nums;}"
        " .big{font-size:2rem;font-weight:700;}"
        " .headline{background:#f1fbf6;border:1px solid #cdeede;border-radius:8px;"
        "padding:.8rem 1rem;margin:1rem 0;}"
        " .hl2{background:#f3f0fb;border:1px solid #d8cdf0;}"
        " code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.9em;}"
        " .sf{color:#666;font-size:.9em;} details{margin:.6rem 0;} summary{cursor:pointer;color:#357;}"
        " p.note{color:#444;font-size:.92em;}")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>ChEBI chemicals targeting specific proteins &amp; gene cross-link</title>
<style>{style}</style></head><body>
<h1>ChEBI chemicals that target a specific protein / gene</h1>
<p>Molecular entities in <code>databases/chebi.json</code> whose ChEBI roles name a
<strong>specific protein target</strong>; plus the cross-link of the corpus
chemicals to HGNC gene targets and the corpus GENETIC entities. Produced by
<code>target_pharm.py</code> &mdash; <code>CHEMICAL/target_pharm.json</code>
(chemicals&rarr;targets) and <code>CHEMICAL/chemical_to_target.json</code>
(inverse: gene&rarr;chemicals).</p>
<div class="headline"><span class="big">{n:,}</span> chemicals target a specific
protein &mdash; out of {n_action:,} with any action-type role and {n_with_roles:,}
with any role.</div>

{table("By target class", by_class, "target class")}
{table("By pharmacological action", by_action, "action")}

<h2>Top target roles (by # chemicals)</h2>
<table><tr><th>target role</th><th class="num">chemicals</th></tr>{toprole}</table>

<h2>Inverse cross-link &mdash; gene &rarr; chemicals that target it</h2>
<div class="headline hl2">Of the {n_corpus_chem:,} corpus protein-targeting
chemicals, <strong>{n_chems_linked}</strong> map to an HGNC gene target across
<strong>{n_genes}</strong> genes; <strong>{n_corpus_genes}</strong> of those genes
are also <strong>GENETIC entities</strong> in this corpus. Full data in
<code>CHEMICAL/chemical_to_target.json</code>.</div>
<p>The {n_corpus_genes} target genes that are corpus GENETIC entities, with the
chemicals that target them. <span class="sf">[bracketed]</span> = the gene's corpus
GENETIC surface forms.</p>
<table><tr><th>gene</th><th>[GENETIC surfaces]</th><th class="num">#chem</th>
<th>chemicals targeting it</th></tr>{"".join(corpus_rows)}</table>
<details><summary>show / hide all {n_genes} target genes (* = corpus GENETIC)</summary>
<table><tr><th>gene</th><th class="num">#chemicals</th><th>chemicals</th></tr>
{inv_full}</table></details>

<h2>Method &amp; caveats</h2>
<p class="note">Target-role = action-word suffix + protein term (EC enzyme or
receptor/kinase/&hellip;), excluding pathway/process roles. The cross-link maps the
target name (enzyme inside an EC role, else role minus action word) to HGNC by exact
casefold / separator-deletion over the six HGNC fields (high precision, partial
recall &mdash; viral/bacterial/broad-EC/generic targets do not resolve to a human
gene). EC roles can be enzyme-family level. The GENETIC-surface side reflects the
corpus's own gene normalization (so broad aliases there carry through).</p>

<h2>All {n:,} protein-targeting chemicals</h2>
<details><summary>show / hide the full list</summary>
<table><tr><th>chebi_id</th><th>chebi_label</th><th>target role(s)</th></tr>
{full}</table></details>
</body></html>
"""
    HTML_OUT.write_text(doc, encoding="utf-8")


if __name__ == "__main__":
    main()
