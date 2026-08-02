#!/usr/bin/env python3
"""
gpu.py -- run the full pipeline (PPI relation-model training -> BioBERT NER ->
entity normalization -> base triples -> learned relation extraction) as ONE step,
on Kaggle (or locally). GPU-enabled for the model-training, NER and
relation-extraction steps.

Orchestrates eighteen steps in dependency order:

    1  run_re_pipeline.py train + calibrate the PPI relation model
                          (BigBIO bioinfer -> ppi-biobert-re/)             GPU, internet
    2  run_re_pipeline.py --task biored   train + calibrate the BioRED relation model
                          (BigBIO biored -> biored-biobert-re/)   GPU, internet [optional]
    3  sentences.py      BioBERT result-sentence selection + NER (experimental_ner/ -> sentences/)  GPU
    4  roman.py          GENETIC -> HGNC (roman key)                 CPU
    5  greek.py          GENETIC -> HGNC (greek key)                 CPU
    6  keys_values.py    HGNC coverage reports                       CPU
    7  controls.py       control-gene flag                           CPU
    8  disease.py        DISEASE -> MONDO                            CPU
    9  phenotypes.py     phenotype flag                              CPU
    10 chemical.py       CHEMICAL -> ChEBI                           CPU
    11 nonchemical.py    non-chemical flag                           CPU
    12 target_pharm.py   chemical<->gene-target cross-links          CPU
    13 triples.py        base triples + normalized variants          CPU
    14 relationships.py  gene-gene slice (genetic_genetic.json)      CPU  [optional]
    15 pub_years.py      PMC->year via NCBI (pmc_years.json)         CPU  [optional, internet]
    16 relation_extraction.py --normalize --route-mode additive
                         BioBERT-scored triples, every applicable model  GPU
    17 compare_re.py     PPI vs BioRED on the same pairs             CPU  [optional]
    18 zip_work.py       bundle working dir -> kaggle_working.zip    CPU

Dependency / strategy notes:
  * run_re_pipeline.py (1, MOST UPSTREAM) trains the relation-extraction model the
    learned RE step (16) depends on. It is itself a three-script pipeline, run here
    as ONE subprocess, that fine-tunes BioBERT into ppi-biobert-re/:
        bigbio_to_re.py  BigBIO KB corpus (HF `datasets`) -> entity-blinded TSVs
        train_re.py      fine-tune BioBERT -> ppi-biobert-re/ checkpoint (imports calibration.py)
        calibration.py   evaluate on test + fit a probability calibrator
                         (-> ppi-biobert-re/calibration.json, summaries/calibration.html)
    GPU; needs internet + the `datasets` library to download BigBIO bioinfer on
    first use. Its output, ppi-biobert-re/, is exactly what step 16 loads via
    RE_MODEL_PPI -- so the model dir no longer has to be uploaded. (Upload it and
    skip step 1 via --steps to reuse an existing checkpoint, mirroring how the
    'sentences' step can be skipped when sentences/ is supplied.)
  * re_pipeline_biored (2) is the SAME script with --task biored: BioRED (Luo et al.
    2022) is a typed, SIGNED corpus, so its checkpoint predicts upregulator/activator,
    downregulator/inhibitor, binds and associated instead of a bare "interacts", and
    one checkpoint covers every entity-type pair this corpus produces. It runs
    ALONGSIDE the PPI model rather than replacing it (BioRED_task_mapping.md section
    6/8): step 16 scores each pair with both, step 17 compares them, and only then is
    the replace-or-not decision worth making. OPTIONAL -- if the download or training
    fails, the run continues with the PPI model alone and nothing is lost.
  * sentences.py (3) runs BioBERT over experimental_ner/ (XML) to select
    original-result sentences and tag DISEASE/GENE/CHEMICAL entities, writing
    sentences/*.json -- the input every later step reads. GPU; needs the HF BioBERT
    models (dmis-lab/biobert-v1.1 + alvaroalon2/biobert_{diseases,genetic,chemical}_ner),
    fetched on first use (internet) or cached.
  * Steps 4-12 build the GENETIC/DISEASE/CHEMICAL normalization libraries.
  * triples.py (13) reads those libraries; relationships.py (14) reads triples.py's
    output; pub_years.py (15) reads relationships.py's genetic_genetic.json.
  * relation_extraction.py (16) is a GPU step (BioBERT inference, auto CUDA) that loads
    the step-1 model via RE_MODEL_PPI and, when step 2 produced it, the step-2 model
    via RE_MODEL_BIORED; --route-mode additive makes BOTH score every pair they cover,
    each triple tagged with predicate.model and a shared pair_id. Its --normalize pass
    imports triples.py to reuse the normalization chain. NOTE: with two models the
    file holds up to two triples per entity pair -- filter by predicate.model before
    building a graph from it.
  * compare_re.py (17) joins those two verdicts on pair_id and writes
    summaries/compare_re.html: coverage, label agreement, how many edges gained a sign,
    and samples to hand-read. This is the evidence for the replace-or-not decision.
  * Steps 2, 14, 15 and 17 are OPTIONAL: if any fails (internet off, missing input,
    only one model trained) the orchestrator warns and CONTINUES.
  * zip_work.py (18, LAST) packs the whole working dir into kaggle_working.zip so the
    entire run is one download; it must run after every other step has written its output.

Why an orchestrator (not one merged file): the scripts are standalone but resolve
paths two ways -- most relative to their own file, greek.py relative to the
current directory -- so they are staged in one writable working dir and run there,
unmodified.

KAGGLE USAGE
  1. Upload the project as a Kaggle Dataset (read-only at /kaggle/input/<ds>/):
       the 18 step scripts above PLUS bigbio_to_re.py + train_re.py (run by
                                                  run_re_pipeline.py in steps 1 and 2)
                                  AND calibration.py (shared: steps 1, 2 + step 16)
       experimental_ner/                             (XML corpus; sentences.py input)
       databases/hgnc_complete_set_2026-05-01.json   (HGNC)
       databases/mondo-clingen.json                  (MONDO)
       databases/chebi.json                          (ChEBI)
     ppi-biobert-re/ and biored-biobert-re/ are GENERATED by steps 1 and 2 -- no need
     to upload them (or upload them and skip those steps via --steps if you already
     have the trained models). sentences/ is likewise generated by step 3.
     pmc_years.json is produced by step 15; upload it under databases/ only to run the
     year filter with internet OFF.
  2. Notebook Settings -> Accelerator -> GPU (steps 1, 2, 3 and 16 use CUDA); enable
     Internet so steps 1-2 can download the BigBIO corpora, step 3 the BioBERT NER
     models, and step 15 the publication years. Setup cell:
         !pip install -q 'datasets<4' bioc
  3. Cell:  !python /kaggle/input/<ds>/gpu.py
     Outputs go to /kaggle/working/{ppi-biobert-re,biored-biobert-re,sentences,summaries,GENETIC,DISEASE,CHEMICAL,TRIPLES}/.

  Read-only input is handled automatically: scripts/modules are copied,
  experimental_ner/ is symlinked, databases/ is copied (writable), and the
  generated dirs (the two model dirs from steps 1-2, sentences/ from step 3) plus the
  output dirs are created in /kaggle/working (the run root). RE_MODEL_PPI and
  RE_MODEL_BIORED are set automatically for step 16 -- each only if that model dir
  actually exists, so a skipped or failed training step never breaks the RE step.

LIBRARIES
  Steps 1 & 2 (run_re_pipeline.py): torch + transformers + `datasets<4` (BigBIO
                  download, needs internet). torch/transformers are preinstalled on
                  Kaggle GPU images, but BigBIO ships as loading scripts (dropped in
                  datasets 4.0), so the newer preinstalled datasets must be pinned in a
                  setup cell -- and step 2's corpus (BioRED) is BioC XML, so its loading
                  script also imports `bioc`, which the image does not carry:
                      !pip install -q 'datasets<4' bioc
                  (`bioc` is needed only for step 2; bioinfer does not use it.)
  Steps 4-15,17 : Python standard library only; step 15 (pub_years.py) needs
                  internet (NCBI E-utilities).
  Steps 3 & 16  : torch + transformers (+ lxml for step 3); preinstalled on Kaggle
                  GPU images. Step 3 also pulls BioBERT models from Hugging Face.

OPTIONS / ENV
  --list / --dry-run        show the plan (roots, inputs, accelerator, steps)
  --steps a,b,c             run only these step names (also env NORM_STEPS)
  --re-args "..."           extra args for step 1 (also env RE_PIPELINE_ARGS)
  --biored-args "..."       extra args for step 2, e.g. "--require-cue"
                            (also env BIORED_PIPELINE_ARGS)
  --input-root PATH         dataset root (also env NORM_INPUT_ROOT)
  --work-root PATH          writable run dir (also env NORM_WORK_ROOT)
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

# each step: name, script, extra args, required ontology DBs, GPU?, model dirs, description
# `models` = checkpoints the step READS (staged + exported via MODEL_ENV); optional keys:
# support=[modules staged with the step], produces_model="dir it generates", optional=True
STEPS = [
    dict(name="re_pipeline", script="run_re_pipeline.py", args=[], dbs=[], gpu=True, models=[],
         support=["bigbio_to_re.py", "train_re.py", "calibration.py"], produces_model="ppi-biobert-re",
         desc="[GPU][internet] train + calibrate the PPI relation model (BigBIO bioinfer -> ppi-biobert-re/) via bigbio_to_re.py + train_re.py + calibration.py"),
    dict(name="re_pipeline_biored", script="run_re_pipeline.py",
         args=["--task", "biored"], dbs=[], gpu=True, models=[], optional=True,
         support=["bigbio_to_re.py", "train_re.py", "calibration.py"],
         produces_model="biored-biobert-re",
         desc="[GPU][internet] train + calibrate the BioRED relation model (BigBIO biored -> biored-biobert-re/): typed + SIGNED edges over every entity-type pair [optional]"),
    dict(name="sentences", script="sentences.py", args=[], dbs=[], gpu=True, models=[],
         desc="[GPU] BioBERT result-sentence selection + NER over experimental_ner/ -> sentences/ (needs HF BioBERT models)"),
    dict(name="roman", script="roman.py", args=[], dbs=["hgnc"], gpu=False, models=[],
         desc="GENETIC surfaces -> HGNC (roman key); writes clean_genetic_ne.tsv + greek_clean_genetic_ne.tsv"),
    dict(name="greek", script="greek.py", args=[], dbs=["hgnc"], gpu=False, models=[],
         desc="Greek-symbol-expanded GENETIC surfaces -> HGNC"),
    dict(name="keys_values", script="keys_values.py", args=[], dbs=[], gpu=False, models=[],
         desc="HGNC coverage reports over the GENETIC libraries"),
    dict(name="controls", script="controls.py", args=[], dbs=["hgnc"], gpu=False, models=[],
         desc="flag experimental-control / tool GENETIC entities (control: yes/no)"),
    dict(name="disease", script="disease.py", args=[], dbs=["mondo"], gpu=False, models=[],
         desc="DISEASE surfaces -> MONDO"),
    dict(name="phenotypes", script="phenotypes.py", args=[], dbs=[], gpu=False, models=[],
         desc="flag phenotype/process DISEASE surfaces (phenotype: yes/no)"),
    dict(name="chemical", script="chemical.py", args=[], dbs=["chebi"], gpu=False, models=[],
         desc="CHEMICAL surfaces -> ChEBI"),
    dict(name="nonchemical", script="nonchemical.py", args=[], dbs=["hgnc"], gpu=False, models=[],
         desc="flag non-chemical CHEMICAL surfaces (non_chemical: yes/no)"),
    dict(name="target_pharm", script="target_pharm.py", args=[], dbs=["chebi", "hgnc"], gpu=False, models=[],
         desc="chemical<->gene-target cross-links (needs GENETIC libs + chemical.json)"),
    dict(name="triples", script="triples.py", args=[], dbs=[], gpu=False, models=[],
         desc="extract base triples + normalized variants (triples.json, triples_*_normalized.json, triples.html)"),
    dict(name="relationships", script="relationships.py", args=[], dbs=[], gpu=False, models=[], optional=True,
         desc="gene-gene slice in disease/chemical sentences (genetic_genetic.json) [optional]"),
    dict(name="pub_years", script="pub_years.py", args=[], dbs=[], gpu=False, models=[], optional=True,
         desc="[internet] PMC->publication year via NCBI; caches databases/pmc_years.json [optional]"),
    dict(name="relation_extraction", script="relation_extraction.py",
         args=["--normalize", "--route-mode", "additive"],
         dbs=[], gpu=True, models=["ppi-biobert-re", "biored-biobert-re"],
         support=["triples.py", "calibration.py"],
         desc="[GPU] BioBERT-scored relations + normalized variant; every applicable checkpoint scores each pair (RE_MODEL_PPI + RE_MODEL_BIORED)"),
    dict(name="compare_re", script="compare_re.py", args=[], dbs=[], gpu=False, models=[], optional=True,
         desc="PPI vs BioRED on identical pairs -> summaries/compare_re.html (needs both models) [optional]"),
    dict(name="zip", script="zip_work.py", args=[], dbs=[], gpu=False, models=[],
         desc="bundle the whole working dir into a downloadable kaggle_working.zip (final step)"),
]
DB_FILES = {
    "hgnc": "hgnc_complete_set_2026-05-01.json",
    "mondo": "mondo-clingen.json",
    "chebi": "chebi.json",
}
RAW_DIR = "experimental_ner"   # sentences.py input (XML); sentences/ is generated from it
OUT_DIRS = ["GENETIC", "DISEASE", "CHEMICAL", "TRIPLES"]
MODEL_ENV = {"ppi-biobert-re": "RE_MODEL_PPI",   # model dir -> env var the RE step reads
             "biored-biobert-re": "RE_MODEL_BIORED"}


def on_kaggle():
    return Path("/kaggle/input").exists() or bool(os.environ.get("KAGGLE_KERNEL_RUN_TYPE"))


def cuda_available():
    try:
        import torch
        return bool(torch.cuda.is_available())
    except Exception:
        return False


def report_accelerator():
    try:
        import torch
        if torch.cuda.is_available():
            return f"GPU: {torch.cuda.get_device_name(0)} (torch {torch.__version__}) -- used by the re_pipeline (training), sentences (NER) + relation_extraction steps"
        return f"no CUDA (torch {torch.__version__}, CPU build) -- the GPU steps (re_pipeline, sentences, relation_extraction) run on CPU (slow)"
    except Exception:
        if shutil.which("nvidia-smi"):
            return "GPU present (nvidia-smi) but torch not importable -- install torch for the GPU steps"
        return "none (CPU) -- GPU steps (re_pipeline, sentences, relation_extraction) would run on CPU"


def support_modules(steps):
    """Modules to stage alongside the selected steps (deduped, in stable order)."""
    seen, mods = set(), []
    for st in steps:
        for m in st.get("support", []):
            if m not in seen:
                seen.add(m)
                mods.append(m)
    return mods


def produced_models(steps):
    """Model dirs that a selected step GENERATES (so they need not be uploaded)."""
    return {st["produces_model"] for st in steps if st.get("produces_model")}


def _looks_like_root(p: Path):
    return (((p / "sentences").is_dir() or (p / "experimental_ner").is_dir())
            and (p / "databases").is_dir() and (p / "roman.py").exists())


def find_input_root(cli):
    cand = cli or os.environ.get("NORM_INPUT_ROOT")
    if cand:
        return Path(cand).resolve()
    here = Path(__file__).resolve().parent
    if _looks_like_root(here):
        return here
    base = Path("/kaggle/input")
    if base.is_dir():
        for c in sorted(base.glob("*")) + sorted(base.glob("*/*")):
            if c.is_dir() and _looks_like_root(c):
                return c.resolve()
    sys.exit("ERROR: could not locate an input root with sentences/, databases/ and the scripts.\n"
             "       Pass --input-root PATH or set NORM_INPUT_ROOT.")


def resolve_work_root(cli, input_root):
    cand = cli or os.environ.get("NORM_WORK_ROOT")
    if cand:
        return Path(cand).resolve()
    if on_kaggle():
        return Path("/kaggle/working").resolve()
    return input_root  # run in place locally


def link_or_copy(src: Path, dst: Path):
    if dst.exists() or dst.is_symlink():
        return
    try:
        os.symlink(src, dst, target_is_directory=src.is_dir())
    except (OSError, NotImplementedError):
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)


def model_dir_ok(mdir: Path):
    """True if `mdir` looks like a usable checkpoint (config + weights)."""
    return ((mdir / "config.json").exists()
            and bool(list(mdir.glob("*.safetensors")) + list(mdir.glob("*.bin"))))


def preflight(input_root: Path, steps):
    missing = []
    produced = produced_models(steps)
    for st in steps:
        if not (input_root / st["script"]).exists():
            missing.append(st["script"])
    for m in support_modules(steps):
        if not (input_root / m).exists():
            missing.append(m + "  (support module)")
    if any(st["name"] == "sentences" for st in steps):
        if not (input_root / RAW_DIR).is_dir():
            missing.append(RAW_DIR + "/  (XML input for sentences.py)")
    elif not (input_root / "sentences").is_dir():
        missing.append("sentences/  (or include the 'sentences' step to generate it)")
    for db in sorted({db for st in steps for db in st["dbs"]}):
        if not (input_root / "databases" / DB_FILES[db]).exists():
            missing.append(f"databases/{DB_FILES[db]}")
    # A step's models are OPTIONAL individually (the RE step reads two checkpoints and
    # runs with either), but it needs at least one: generated by a selected step, or
    # already on disk. Missing-but-not-generated ones are just not exported.
    for st in steps:
        if not st["models"]:
            continue
        have = [m for m in st["models"] if m in produced or model_dir_ok(input_root / m)]
        if not have:
            missing.append(f"{' or '.join(st['models'])}/  (a model dir with config.json + "
                           f"weights, or include the step that trains it)")
    if missing:
        sys.exit("ERROR: missing required inputs under %s:\n  %s" % (input_root, "\n  ".join(missing)))


def stage_work(input_root: Path, work_root: Path, steps):
    work_root.mkdir(parents=True, exist_ok=True)
    produced = produced_models(steps)
    if work_root.resolve() != input_root.resolve():
        for st in steps:
            shutil.copy2(input_root / st["script"], work_root / st["script"])
        for m in support_modules(steps):
            shutil.copy2(input_root / m, work_root / m)
        # databases: writable copy, NOT a symlink (pub_years.py writes pmc_years.json
        # here; a symlink would point back into read-only /kaggle/input).
        ddst = work_root / "databases"
        if ddst.is_symlink():
            ddst.unlink()
        ddst.mkdir(parents=True, exist_ok=True)
        for f in (input_root / "databases").glob("*"):
            if f.is_file():
                shutil.copy2(f, ddst / f.name)
        # sentence source: when the 'sentences' step runs it generates sentences/
        # from read-only experimental_ner/; otherwise symlink the prebuilt sentences/.
        if any(st["name"] == "sentences" for st in steps):
            link_or_copy(input_root / RAW_DIR, work_root / RAW_DIR)
        else:
            link_or_copy(input_root / "sentences", work_root / "sentences")
        # models: stage only those NOT generated this run (steps 1-2 produce theirs) and
        # actually present -- an absent optional checkpoint is simply not exported.
        for st in steps:
            for m in st["models"]:
                if m not in produced and model_dir_ok(input_root / m):
                    link_or_copy(input_root / m, work_root / m)
    outs = list(OUT_DIRS) + (["sentences", "summaries"] if any(st["name"] == "sentences" for st in steps) else [])
    for d in outs:
        (work_root / d).mkdir(parents=True, exist_ok=True)


def main():
    ap = argparse.ArgumentParser(description="RE-model training + normalization + relation-extraction pipeline (Kaggle/local, GPU-enabled)")
    ap.add_argument("--list", "--dry-run", dest="dry", action="store_true",
                    help="show the plan and exit")
    ap.add_argument("--steps", default=os.environ.get("NORM_STEPS"),
                    help="comma-separated subset of step names")
    ap.add_argument("--re-args", default=os.environ.get("RE_PIPELINE_ARGS"),
                    help="extra args passed through to step 1 (run_re_pipeline.py), e.g. "
                         "--re-args \"--neg-ratio 2 --add-marker-tokens\" "
                         "(also env RE_PIPELINE_ARGS)")
    ap.add_argument("--biored-args", default=os.environ.get("BIORED_PIPELINE_ARGS"),
                    help="extra args passed through to step 2 (run_re_pipeline.py --task biored), "
                         "e.g. --biored-args \"--require-cue\" (also env BIORED_PIPELINE_ARGS)")
    ap.add_argument("--input-root", default=None)
    ap.add_argument("--work-root", default=None)
    args = ap.parse_args()

    only = {s.strip() for s in args.steps.split(",")} if args.steps else None
    steps = [dict(st) for st in STEPS if (only is None or st["name"] in only)]
    if not steps:
        sys.exit(f"ERROR: no matching steps in {args.steps!r}. Valid: {', '.join(s['name'] for s in STEPS)}")

    for flag, value, target in (("--re-args", args.re_args, "re_pipeline"),
                                ("--biored-args", args.biored_args, "re_pipeline_biored")):
        if not value:
            continue
        st = next((s for s in steps if s["name"] == target), None)
        if st is None:
            sys.exit(f"ERROR: {flag} was given but the '{target}' step is not in the plan "
                     f"(it is excluded by --steps).")
        st["args"] = list(st["args"]) + shlex.split(value)

    input_root = find_input_root(args.input_root)
    preflight(input_root, steps)
    work_root = resolve_work_root(args.work_root, input_root)
    if any(st["name"] == "sentences" for st in steps):
        nfiles = len(list((input_root / RAW_DIR).glob("*.xml")))
        src_line = f" input      : {nfiles:,} XML in {RAW_DIR}/ (sentences/ generated by step 2)"
    else:
        nfiles = len(list((input_root / "sentences").glob("*.json")))
        src_line = f" sentences  : {nfiles:,} *.json files"

    print("=" * 66)
    print(" kaggle pipeline  (Kaggle)" if on_kaggle() else " kaggle pipeline  (local)")
    print("=" * 66)
    print(f" input_root : {input_root}")
    print(f" work_root  : {work_root}" + ("  (in place)" if work_root == input_root else "  (staged)"))
    print(src_line)
    print(f" accelerator: {report_accelerator()}")
    print(f" steps      : {', '.join(s['name'] for s in steps)}")
    if args.re_args:
        print(f" re_args    : {args.re_args}  (-> step 1 run_re_pipeline.py)")
    if args.biored_args:
        print(f" biored_args: {args.biored_args}  (-> step 2 run_re_pipeline.py --task biored)")
    print()
    for i, st in enumerate(steps, 1):
        tags = []
        if st["dbs"]:
            tags.append("needs: " + ", ".join(st["dbs"]))
        if st["gpu"]:
            tags.append("GPU")
        if st.get("produces_model"):
            tags.append("makes: " + st["produces_model"])
        if st["models"]:
            tags.append("models: " + ", ".join(st["models"]))
        tag = ("  [" + " | ".join(tags) + "]") if tags else ""
        print(f"  {i:>2}. {st['name']:<19} {st['script']:<22} {st['desc']}{tag}")
    if args.dry:
        print("\n(--list) plan only; nothing executed.")
        return

    stage_work(input_root, work_root, steps)

    results, t0 = [], time.time()
    for st in steps:
        if st["gpu"] and not cuda_available():
            print(f"\n[!] step '{st['name']}' is a GPU step but no CUDA is visible -- it will run on CPU "
                  f"(much slower). Enable the Kaggle GPU accelerator to speed it up.", flush=True)
        env = os.environ.copy()
        # export each checkpoint this step reads -- but only the ones that exist by now,
        # so a skipped/failed optional training step cannot point the RE step at nothing
        for m in st["models"]:
            if m in MODEL_ENV and model_dir_ok(work_root / m):
                env[MODEL_ENV[m]] = str(work_root / m)
            elif m in MODEL_ENV:
                print(f"[!] {m}/ not present -- {MODEL_ENV[m]} not set for step '{st['name']}'.",
                      flush=True)
        print(f"\n{'=' * 66}\n[{st['name']}] {st['desc']}\n{'=' * 66}", flush=True)
        t = time.time()
        rc = subprocess.run([sys.executable, st["script"], *st["args"]], cwd=str(work_root), env=env).returncode
        dt = time.time() - t
        results.append((st["name"], rc, dt))
        if rc != 0:
            if st.get("optional"):
                print(f"\n[!] optional step '{st['name']}' failed (exit {rc}) after {dt:.1f}s -- continuing.", flush=True)
                continue
            print(f"\n!! step '{st['name']}' FAILED (exit {rc}) after {dt:.1f}s -- stopping.", flush=True)
            break
        print(f"-- {st['name']} done in {dt:.1f}s", flush=True)

    print(f"\n{'=' * 66}\n SUMMARY\n{'=' * 66}")
    for name, rc, dt in results:
        print(f"  {'OK  ' if rc == 0 else 'FAIL'}  {name:<20} {dt:8.1f}s")
    print(f"  total {time.time() - t0:.1f}s")
    for d in OUT_DIRS:
        files = sorted((work_root / d).glob("*"))
        if files:
            print(f"\n  {d}/  ({len(files)} files)")
            for f in files:
                print(f"    {f.stat().st_size:>13,}  {f.name}")
    # the step-16 download bundle lives at the work-root, outside OUT_DIRS
    bundle = work_root / "kaggle_working.zip"
    if bundle.exists():
        print(f"\n  bundle: {bundle.stat().st_size:>13,}  {bundle.name}")
    sys.exit(1 if any(rc != 0 for _, rc, _ in results) else 0)


if __name__ == "__main__":
    main()
