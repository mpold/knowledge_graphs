#!/usr/bin/env python3
"""
gpu.py -- run the full pipeline (PPI relation-model training -> BioBERT NER ->
entity normalization -> base triples -> learned relation extraction) as ONE step,
on Kaggle (or locally). GPU-enabled for the model-training, NER and
relation-extraction steps.

Orchestrates twenty steps in dependency order:

    1  run_re_pipeline.py train + calibrate the PPI relation model
                          (BigBIO bioinfer -> ppi-biobert-re/)   GPU, internet [gated]
    2  run_re_pipeline.py --task biored   train + calibrate the BioRED relation model
                          (BigBIO biored -> biored-biobert-re/)  GPU, internet [gated]
    3  drug_lexicon.py --build   NCIt drug names -> databases/ncit_drugs.json
                         CPU  [optional, internet]
    4  sentences.py      BioBERT result-sentence selection + NER (experimental_ner/ -> sentences/)  GPU
    5  roman.py          GENETIC -> HGNC (roman key)                 CPU
    6  greek.py          GENETIC -> HGNC (greek key)                 CPU
    7  keys_values.py    HGNC coverage reports                       CPU
    8  controls.py       control-gene flag                           CPU
    9  disease.py        DISEASE -> MONDO                            CPU
    10 phenotypes.py     phenotype flag                              CPU
    11 chemical.py       CHEMICAL -> ChEBI, then NCIt for the rest   CPU
    12 nonchemical.py    non-chemical flag                           CPU
    13 target_pharm.py   chemical<->gene-target cross-links          CPU
    14 triples.py        base triples + normalized variants          CPU
    15 relationships.py  gene-gene slice (genetic_genetic.json)      CPU  [optional]
    16 pub_years.py      PMC->year, table first (pmc_years.json)     CPU  [optional, internet if no table]
    17 relation_extraction.py --normalize --route-mode additive
                         BioBERT-scored triples, every applicable model  GPU
    18 compare_re.py     PPI vs BioRED on the same pairs             CPU  [optional]
    19 triples_per_year.py  triple scores by publication year        CPU  [optional]
    20 zip_work.py       bundle working dir -> kaggle_working.zip    CPU

THE DRUG LEXICON (step 3 -- why it is optional, and what is lost without it)
  ChEBI is a small-molecule ontology, so biologics fall out of the pipeline twice: the
  BioBERT chemical model, trained on the same kind of data, tags only 27% of the named
  antibody mentions in a lung corpus, and ChEBI then has no id for the survivors. Step 3
  distils the NCI Thesaurus into databases/ncit_drugs.json, which both the NER stage
  (sentences.py, a lexicon pass after the models) and the normalizer (chemical.py, a
  fallback after the ChEBI cascade) read.
  It is NOT a hard dependency of either: without the file both fall back to the INN-stem
  regex alone ("-mab", "-tug", "-bart", "-mig", "-cept"), which still recognises antibody
  names but cannot resolve brand names, code names or synonyms, and chemical.py cannot
  normalize any of them. Run it once with internet on and upload the JSON with the dataset
  thereafter, exactly like pmc_years.json.

THE TRAINING GATE (steps 1-2 are not run unconditionally)
  Neither training step reads the corpus: they fine-tune BioBERT on fixed public
  BigBIO datasets, so the checkpoint is a function of (dataset, seed, hyperparams)
  only and is IDENTICAL for every corpus this pipeline is pointed at. Training them
  again per run is wasted GPU time, so steps 1 and 2 run ONLY when the training
  output is not already sitting in the input root (normally gpu_bundle/):

      ppi-biobert-re/      config.json, tokenizer.json, tokenizer_config.json,
      biored-biobert-re/   calibration.json + weights (*.safetensors or *.bin)
      ppi_data/            train.tsv, dev.tsv, test.tsv
      biored_data/

  All four complete    -> steps 1 AND 2 are dropped from the plan; the checkpoints
                          already on disk are staged and exported to step 16 instead.
  Anything missing/empty -> both training steps run, and the gate prints exactly which
                          entries were absent.

  It is all four or none, deliberately: step 16 routes between the two checkpoints and
  step 17 compares them, so a run must not mix a reused PPI model with a freshly trained
  BioRED one (or the reverse). An empty file counts as missing, so a half-copied bundle
  retrains rather than loading a truncated checkpoint. The *_data/ TSVs are gated too --
  they are what makes a reused bundle re-calibratable without a re-download.

  Overrides, any of which forces training back on: --retrain (env NORM_RETRAIN),
  naming a training step in --steps, or passing --re-args / --biored-args.

Dependency / strategy notes:
  * run_re_pipeline.py (1, MOST UPSTREAM) trains the relation-extraction model the
    learned RE step (16) depends on. It is itself a three-script pipeline, run here
    as ONE subprocess, that fine-tunes BioBERT into ppi-biobert-re/:
        bigbio_to_re.py  BigBIO KB corpus (HF `datasets`) -> entity-blinded TSVs
        train_re.py      fine-tune BioBERT -> ppi-biobert-re/ checkpoint (imports calibration.py)
        calibration.py   evaluate on test + fit a probability calibrator
                         (-> ppi-biobert-re/calibration.json, summaries/calibration_ppi.html)
    GPU; needs internet + the `datasets` library to download BigBIO bioinfer on
    first use. Its output, ppi-biobert-re/, is exactly what step 16 loads via
    RE_MODEL_PPI -- so the model dir no longer has to be uploaded. When it IS
    uploaded (with the three dirs above), the gate skips this step by itself; --steps
    still works as the manual equivalent, mirroring how the 'sentences' step can be
    skipped when sentences/ is supplied.
  * re_pipeline_biored (2) is the SAME script with --task biored: BioRED (Luo et al.
    2022) is a typed, SIGNED corpus, so its checkpoint predicts upregulator/activator,
    downregulator/inhibitor, binds and associated instead of a bare "interacts", and
    one checkpoint covers every entity-type pair this corpus produces. It runs
    ALONGSIDE the PPI model rather than replacing it (BioRED_task_mapping.md section
    6/8): step 16 scores each pair with both, step 17 compares them, and only then is
    the replace-or-not decision worth making. OPTIONAL -- if the download or training
    fails, the run continues with the PPI model alone and nothing is lost. Its two
    output dirs are half of what the training gate above tests for.
  * sentences.py (3) runs BioBERT over experimental_ner/ (XML) to select
    original-result sentences and tag DISEASE/GENE/CHEMICAL entities, writing
    sentences/*.json -- the input every later step reads. GPU; needs the HF BioBERT
    models (dmis-lab/biobert-v1.1 + alvaroalon2/biobert_{diseases,genetic,chemical}_ner),
    fetched on first use (internet) or cached.
  * Steps 4-12 build the GENETIC/DISEASE/CHEMICAL normalization libraries. Steps 5-12
    run as three lanes side by side -- GENETIC (roman -> greek -> keys_values ->
    controls), DISEASE (disease -> phenotypes), CHEMICAL (chemical -> nonchemical) --
    because each reads only sentences/ + databases/ and writes only its own dir. Within a
    lane the order is kept: controls, phenotypes and nonchemical rewrite their lane's
    libraries in place. target_pharm (13) waits for all
    three. A failed step stops only its lane; the others finish, then the run stops.
    NORM_PARALLEL_CPU=0 runs them in turn. Watch memory: chemical.py's ChEBI load overlaps
    the MONDO and HGNC loads.
  * triples.py (14) reads those libraries; relationships.py (15) reads triples.py's
    output; pub_years.py (16) reads relationships.py's genetic_genetic.json.
  * relation_extraction.py (17) is a GPU step (BioBERT inference, auto CUDA) that loads
    the step-1 model via RE_MODEL_PPI and, when step 2 produced it, the step-2 model
    via RE_MODEL_BIORED; --route-mode additive makes BOTH score every pair they cover,
    each triple tagged with predicate.model and a shared pair_id. Its --normalize pass
    imports triples.py to reuse the normalization chain. NOTE: with two models the
    file holds up to two triples per entity pair -- filter by predicate.model before
    building a graph from it.
  * compare_re.py (18) joins those two verdicts on pair_id and writes
    summaries/compare_re.html: coverage, label agreement, how many edges gained a sign,
    and samples to hand-read. This is the evidence for the replace-or-not decision.
  * triples_per_year.py (19) dates step 17's normalized triples by publication year
    (databases/pmc_years.json from step 16, topped up from pmids/pmid_pmc_ids.tsv when it
    is shipped) and writes summaries/triples_per_year.html: a per-year box plot of the
    composite score, triple counts per year, and a score histogram. Offline.
  * Steps 2, 3, 15, 16, 18 and 19 are OPTIONAL: if any fails (internet off, missing input,
    only one model trained) the orchestrator warns and CONTINUES.
  * zip_work.py (20, LAST) packs the whole working dir into kaggle_working.zip so the
    entire run is one download; it must run after every other step has written its output.

Why an orchestrator (not one merged file): the scripts are standalone but resolve
paths two ways -- most relative to their own file, greek.py relative to the
current directory -- so they are staged in one writable working dir and run there,
unmodified.

KAGGLE USAGE
  1. Upload the project as a Kaggle Dataset (read-only at /kaggle/input/<ds>/):
       the 20 step scripts above PLUS bigbio_to_re.py + train_re.py (run by
                                                  run_re_pipeline.py in steps 1 and 2)
                                  AND calibration.py (shared: steps 1, 2 + step 17)
       experimental_ner/                             (XML corpus; sentences.py input)
       databases/hgnc_complete_set_2026-05-01.json   (HGNC)
       databases/mondo-clingen.json                  (MONDO)
       databases/chebi.json                          (ChEBI)
       databases/ncit_drugs.json                     (NCIt drugs; step 3 builds it)
     ppi-biobert-re/ + ppi_data/ and biored-biobert-re/ + biored_data/ are GENERATED
     by steps 1 and 2 on a first run -- but once you have them, upload all four with
     the dataset and the training gate drops both steps automatically (no --steps
     needed), which is the normal case for every corpus after the first.
     sentences/ is likewise generated by step 4.
     pmc_years.json is produced by step 16, ncit_drugs.json by step 3; upload both under
     databases/ to run the year filter and the drug lexicon with internet OFF.
  2. Notebook Settings -> Accelerator -> GPU (steps 1, 2, 4 and 17 use CUDA); enable
     Internet so steps 1-2 can download the BigBIO corpora, step 3 the BioBERT NER
     models, and step 15 the publication years. Setup cell:
         !pip install -q 'datasets<4' bioc
  3. Cell:  !python /kaggle/input/<ds>/gpu.py
     Outputs go to /kaggle/working/{ppi-biobert-re,biored-biobert-re,sentences,summaries,GENETIC,DISEASE,CHEMICAL,TRIPLES}/.

  Read-only input is handled automatically: scripts/modules are copied,
  experimental_ner/ is symlinked, databases/ is copied (writable), model + *_data dirs
  the gate is reusing are symlinked/copied across, and the dirs a selected step
  generates (the two model dirs from steps 1-2, sentences/ from step 3) plus the
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
  Steps 4-15,17 : Python standard library only; pub_years.py needs internet
                  (NCBI E-utilities) ONLY for accessions pmids/pmid_pmc_ids.tsv
                  does not cover -- ship that table beside the run and it is offline.
  Steps 3 & 16  : torch + transformers (+ lxml for step 3); preinstalled on Kaggle
                  GPU images. Step 3 also pulls BioBERT models from Hugging Face.

OPTIONS / ENV
  --list / --dry-run        show the plan (roots, inputs, accelerator, steps) --
                            including the training gate's verdict and its reason, and
                            one OK/BUILD/ABSENT line per optional cache the plan reads
                            (ncit_drugs.json, pmc_years.json): ABSENT names what is lost
  --steps a,b,c             run only these step names (also env NORM_STEPS)
  --retrain                 train even when the gate finds the artifacts complete
                            (also env NORM_RETRAIN=1)
  NORM_PARALLEL_TRAIN=0     with >= 2 GPUs, steps 1 and 2 run side by side (one GPU
                            each, output lines prefixed with the step name); this env
                            runs them one after the other instead
  NORM_PARALLEL_CPU=0       steps 5-12 normally run as three lanes side by side (GENETIC:
                            roman -> greek -> keys_values -> controls |
                            DISEASE: disease -> phenotypes | CHEMICAL: chemical ->
                            nonchemical; each writes only its own dir); this env runs them
                            one after the other instead
  BIOBERT_GPUS / RE_GPUS    cap the GPUs steps 4 / 17 spread over (default: all visible)
  BIOBERT_PARSE_WORKERS     processes step 4 parses XML with, ahead of the GPU (default
                            min(4, cores - 1); 0 = inline)
  --re-args "..."           extra args for step 1 (also env RE_PIPELINE_ARGS);
                            implies --retrain for that step
  --biored-args "..."       extra args for step 2, e.g. "--require-cue"
                            (also env BIORED_PIPELINE_ARGS); implies --retrain
  --input-root PATH         dataset root (also env NORM_INPUT_ROOT)
  --work-root PATH          writable run dir (also env NORM_WORK_ROOT)
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

# each step: name, script, extra args, required ontology DBs, GPU?, model dirs, description
# `models` = checkpoints the step READS (staged + exported via MODEL_ENV); optional keys:
# support=[modules staged with the step], produces_model="dir it generates", optional=True,
# reads_sentences=False for a step that does not consume sentences/ (default True),
# lane="DIR" for a normalization step that reads sentences/ + databases/ and writes ONLY
# under DIR/ (see run_lanes). Steps sharing a lane run in plan order -- controls, phenotypes
# and nonchemical rewrite their lane's libraries in place -- while different lanes overlap.
STEPS = [
    dict(name="re_pipeline", script="run_re_pipeline.py", args=[], dbs=[], gpu=True, models=[],
         support=["bigbio_to_re.py", "train_re.py", "calibration.py"], produces_model="ppi-biobert-re",
         desc="[GPU][internet] train + calibrate the PPI relation model (BigBIO bioinfer -> ppi-biobert-re/) via bigbio_to_re.py + train_re.py + calibration.py"),
    dict(name="re_pipeline_biored", script="run_re_pipeline.py",
         args=["--task", "biored"], dbs=[], gpu=True, models=[], optional=True,
         support=["bigbio_to_re.py", "train_re.py", "calibration.py"],
         produces_model="biored-biobert-re",
         desc="[GPU][internet] train + calibrate the BioRED relation model (BigBIO biored -> biored-biobert-re/): typed + SIGNED edges over every entity-type pair [optional]"),
    dict(name="drug_lexicon", script="drug_lexicon.py", args=["--build"], dbs=[], gpu=False,
         models=[], optional=True, reads_sentences=False,
         desc="[internet] distil the NCI Thesaurus into databases/ncit_drugs.json: the drug names "
              "the BioBERT chemical model misses (biologics) and ChEBI cannot represent [optional]"),
    dict(name="sentences", script="sentences.py", args=[], dbs=[], gpu=True, models=[],
         support=["drug_lexicon.py"],
         desc="[GPU] BioBERT result-sentence selection + NER over experimental_ner/ -> sentences/ "
              "(needs HF BioBERT models; adds an INN-stem + NCIt lexicon pass when the lexicon is present)"),
    dict(name="roman", script="roman.py", args=[], lane="GENETIC", dbs=["hgnc"], gpu=False, models=[],
         desc="GENETIC surfaces -> HGNC (roman key); writes clean_genetic_ne.tsv + greek_clean_genetic_ne.tsv"),
    dict(name="greek", script="greek.py", args=[], lane="GENETIC", dbs=["hgnc"], gpu=False, models=[],
         desc="Greek-symbol-expanded GENETIC surfaces -> HGNC"),
    dict(name="keys_values", script="keys_values.py", args=[], lane="GENETIC", dbs=[], gpu=False, models=[],
         desc="HGNC coverage reports over the GENETIC libraries"),
    dict(name="controls", script="controls.py", args=[], lane="GENETIC", dbs=["hgnc"], gpu=False, models=[],
         desc="flag experimental-control / tool GENETIC entities (control: yes/no)"),
    dict(name="disease", script="disease.py", args=[], lane="DISEASE", dbs=["mondo"], gpu=False, models=[],
         desc="DISEASE surfaces -> MONDO"),
    dict(name="phenotypes", script="phenotypes.py", args=[], lane="DISEASE", dbs=[], gpu=False, models=[],
         desc="flag phenotype/process DISEASE surfaces (phenotype: yes/no)"),
    dict(name="chemical", script="chemical.py", args=[], lane="CHEMICAL", dbs=["chebi"], gpu=False, models=[],
         support=["drug_lexicon.py"],
         desc="CHEMICAL surfaces -> ChEBI, then NCIt for the biologics ChEBI cannot represent"),
    dict(name="nonchemical", script="nonchemical.py", args=[], lane="CHEMICAL", dbs=["hgnc"], gpu=False, models=[],
         desc="flag non-chemical CHEMICAL surfaces (non_chemical: yes/no)"),
    dict(name="target_pharm", script="target_pharm.py", args=[], dbs=["chebi", "hgnc"], gpu=False, models=[],
         desc="chemical<->gene-target cross-links (needs GENETIC libs + chemical.json)"),
    dict(name="triples", script="triples.py", args=[], dbs=[], gpu=False, models=[],
         desc="extract base triples + normalized variants (triples.json, triples_*_normalized.json, triples.html)"),
    dict(name="relationships", script="relationships.py", args=[], dbs=[], gpu=False, models=[], optional=True,
         desc="gene-gene slice in disease/chemical sentences (genetic_genetic.json) [optional]"),
    dict(name="pub_years", script="pub_years.py", args=[], dbs=[], gpu=False, models=[], optional=True,
         desc="PMC->publication year from pmids/pmid_pmc_ids.tsv, NCBI only for what it lacks; "
              "caches databases/pmc_years.json [optional]"),
    dict(name="relation_extraction", script="relation_extraction.py",
         args=["--normalize", "--route-mode", "additive"],
         dbs=[], gpu=True, models=["ppi-biobert-re", "biored-biobert-re"],
         support=["triples.py", "calibration.py"],
         desc="[GPU] BioBERT-scored relations + normalized variant; every applicable checkpoint scores each pair (RE_MODEL_PPI + RE_MODEL_BIORED)"),
    dict(name="compare_re", script="compare_re.py", args=[], dbs=[], gpu=False, models=[], optional=True,
         desc="PPI vs BioRED on identical pairs -> summaries/compare_re.html (needs both models) [optional]"),
    dict(name="triples_per_year", script="triples_per_year.py", args=[], dbs=[], gpu=False, models=[],
         optional=True, reads_sentences=False,
         desc="triple score distribution per publication year (pmc_years.json) -> "
              "summaries/triples_per_year.html [optional]"),
    dict(name="zip", script="zip_work.py", args=[], dbs=[], gpu=False, models=[],
         desc="bundle the whole working dir into a downloadable kaggle_working.zip (final step)"),
]
DB_FILES = {                    # hard requirements: a selected step aborts if one is missing
    "hgnc": "hgnc_complete_set_2026-05-01.json",
    "mondo": "mondo-clingen.json",
    "chebi": "chebi.json",
}
# Caches produced by a step here and read by later ones, so NOT in DB_FILES: gating on one
# would abort every run that has not built it yet, when the steps that read it degrade instead
# of failing. They are reported in the header (see optional_db_report) so a degraded run says
# so up front -- absent is a quieter result, not a broken one.
# read_by: steps that consume the file; built_by: the step that writes it; without: what is lost.
OPTIONAL_DBS = {
    "ncit_drugs": dict(file="ncit_drugs.json", built_by="drug_lexicon",
                       read_by=("sentences", "chemical"),
                       without="INN-stem regex only -- no brand/code names, and chemical.py "
                               "cannot normalize the biologics"),
    "pmc_years": dict(file="pmc_years.json", built_by="pub_years",
                      read_by=("pub_years", "triples_per_year"),
                      without="every PMC id is re-fetched from NCBI (needs internet)"),
}
RAW_DIR = "experimental_ner"   # sentences.py input (XML); sentences/ is generated from it
OUT_DIRS = ["GENETIC", "DISEASE", "CHEMICAL", "TRIPLES"]
MODEL_ENV = {"ppi-biobert-re": "RE_MODEL_PPI",   # model dir -> env var the RE step reads
             "biored-biobert-re": "RE_MODEL_BIORED"}
# converted-TSV dir each training step writes next to its checkpoint (run_re_pipeline
# --data); staged alongside the models so a reused bundle can be re-calibrated in place
MODEL_DATA = {"ppi-biobert-re": "ppi_data", "biored-biobert-re": "biored_data"}

# ---- the training gate ----------------------------------------------------------
TRAINING_STEPS = ("re_pipeline", "re_pipeline_biored")
WEIGHTS = ("*.safetensors", "*.bin")   # either satisfies "the weights are here"
CHECKPOINT_CONTENTS = ("config.json", "tokenizer.json", "tokenizer_config.json",
                       "calibration.json", WEIGHTS)
TSV_CONTENTS = ("train.tsv", "dev.tsv", "test.tsv")
# What "the training output is already here" means: BOTH steps' checkpoints AND both
# converted-TSV dirs, complete. Steps 1-2 run iff anything here is absent or empty --
# all four or none, so the two checkpoints the RE step routes between always come from
# the same generation (see THE TRAINING GATE above).
TRAINED_CONTENTS = {
    "ppi-biobert-re": CHECKPOINT_CONTENTS,
    "ppi_data": TSV_CONTENTS,
    "biored-biobert-re": CHECKPOINT_CONTENTS,
    "biored_data": TSV_CONTENTS,
}


def on_kaggle():
    return Path("/kaggle/input").exists() or bool(os.environ.get("KAGGLE_KERNEL_RUN_TYPE"))


def cuda_available():
    try:
        import torch
        return bool(torch.cuda.is_available())
    except Exception:
        return False


def gpu_count():
    try:
        import torch
        return torch.cuda.device_count() if torch.cuda.is_available() else 0
    except Exception:
        return 0


def parallel_training(steps):
    """Training steps to launch side by side, one per GPU -- or () to run them in turn.

    The two runs are independent (own BigBIO corpus, TSV dir, checkpoint dir and report),
    so with two GPUs they cost the wall time of the slower one instead of the sum. Needs
    both steps planned and >= 2 visible GPUs; an explicit --gpus in --re-args/--biored-args,
    or NORM_PARALLEL_TRAIN=0, keeps the sequential run."""
    planned = [st for st in steps if st["name"] in TRAINING_STEPS]
    if (len(planned) < 2 or gpu_count() < 2
            or os.environ.get("NORM_PARALLEL_TRAIN", "1").strip().lower() in ("0", "false", "no")
            or any("--gpus" in st["args"] for st in planned)):
        return ()
    return tuple(st["name"] for st in planned)


def _pump(tag, argv, env, cwd, lock):
    """Run one subprocess, echoing each output line as "[tag] line" under `lock` so
    concurrent steps' logs interleave line by line. Returns (returncode, seconds)."""
    t = time.time()
    env = dict(env, PYTHONUNBUFFERED="1",
               TQDM_MININTERVAL="30")   # a pipe turns every tqdm refresh into a line
    proc = subprocess.Popen(argv, cwd=cwd, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                            errors="replace")
    for line in proc.stdout:
        line = line.rstrip()
        if line:
            with lock:
                print(f"[{tag}] {line}", flush=True)
    return proc.wait(), time.time() - t


def _run_threads(targets):
    threads = [threading.Thread(target=fn, args=args) for fn, args in targets]
    for th in threads:
        th.start()
    for th in threads:
        th.join()


def run_side_by_side(jobs, cwd):
    """Run [(tag, argv, env), ...] concurrently; every output line is prefixed "[tag] " so
    the interleaved logs stay readable. Returns {tag: (returncode, seconds)}."""
    lock, out = threading.Lock(), {}

    def one(tag, argv, env):
        out[tag] = _pump(tag, argv, env, cwd, lock)

    _run_threads([(one, job) for job in jobs])
    return out


def lane_block(steps, i):
    """The run of consecutive lane-tagged steps starting at steps[i], grouped by lane in
    plan order -> {lane: [step, ...]}; {} when steps[i] has no lane."""
    lanes = {}
    for st in steps[i:]:
        if not st.get("lane"):
            break
        lanes.setdefault(st["lane"], []).append(st)
    return lanes


def parallel_lanes_enabled():
    return os.environ.get("NORM_PARALLEL_CPU", "1").strip().lower() not in ("0", "false", "no")


def run_lanes(lanes, cwd, step_env):
    """Run each lane's steps in order, the lanes side by side (one thread each).

    The lanes are independent -- each reads sentences/ + databases/ and writes only its own
    GENETIC/, DISEASE/ or CHEMICAL/ dir -- so they cost the wall time of the slowest lane
    instead of the sum. Within a lane a failed required step stops that lane (its later
    steps build on it); the other lanes run to completion. Returns {step name: (rc, seconds)}
    for every step that ran."""
    lock, out = threading.Lock(), {}

    def chain(lane_steps):
        for st in lane_steps:
            with lock:
                env = step_env(st)
            rc, dt = _pump(st["name"], [sys.executable, st["script"], *st["args"]], env, cwd, lock)
            out[st["name"]] = (rc, dt)
            with lock:
                print(f"-- {st['name']} {'done' if rc == 0 else f'FAILED (exit {rc})'} "
                      f"in {dt:.1f}s", flush=True)
            if rc != 0 and not st.get("optional"):
                return

    _run_threads([(chain, (ls,)) for ls in lanes.values()])
    return out


def report_accelerator():
    try:
        import torch
        if torch.cuda.is_available():
            n = torch.cuda.device_count()
            return (f"GPU: {n} x {torch.cuda.get_device_name(0)} (torch {torch.__version__}) -- "
                    f"sentences (NER) + relation_extraction use all {n} (BIOBERT_GPUS / RE_GPUS "
                    f"to cap); the two training runs take one GPU each when both run")
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


def missing_contents(d: Path, required):
    """Entries of `required` absent (or empty) under directory `d`; [] = full content.

    A pattern containing '*' is satisfied by any non-empty matching file; a tuple of
    patterns by any ONE of them (model weights are .safetensors OR .bin).
    """
    if not d.is_dir():
        return [f"{d.name}/  (directory absent)"]

    def present(pat):
        if "*" in pat:
            return any(f.is_file() and f.stat().st_size > 0 for f in d.glob(pat))
        f = d / pat
        return f.is_file() and f.stat().st_size > 0

    out = []
    for req in required:
        alts = req if isinstance(req, tuple) else (req,)
        if not any(present(p) for p in alts):
            out.append(f"{d.name}/" + " or ".join(alts))
    return out


def missing_training_artifacts(input_root: Path):
    """What TRAINED_CONTENTS is missing under `input_root` ([] = training can be skipped)."""
    return [m for name, req in TRAINED_CONTENTS.items()
            for m in missing_contents(input_root / name, req)]


def apply_training_gate(steps, input_root: Path, forced_by: str):
    """Drop steps 1-2 unless the already-trained artifacts are incomplete.

    Returns (steps, report_lines). `forced_by` is a non-empty reason string when an
    override (--retrain / --steps / --re-args / --biored-args) keeps training on.
    """
    planned = [st["name"] for st in steps if st["name"] in TRAINING_STEPS]
    if not planned:
        return steps, []
    gated = ", ".join(n + "/" for n in TRAINED_CONTENTS)
    if forced_by:
        return steps, [f" training   : RUN ({', '.join(planned)}) -- forced by {forced_by}"]

    gap = missing_training_artifacts(input_root)
    if gap:
        lines = [f" training   : RUN ({', '.join(planned)}) -- training output incomplete "
                 f"under {input_root}:"]
        lines += [f"              missing  {g}" for g in gap]
        return steps, lines

    kept = [st for st in steps if st["name"] not in TRAINING_STEPS]
    return kept, [f" training   : SKIP ({', '.join(planned)}) -- complete in {input_root}:",
                  f"              {gated}",
                  "              both checkpoints are reused as-is (--retrain to train anyway)"]


def optional_db_report(input_root: Path, steps):
    """Header lines for the optional caches -- advisory, never fatal.

    One line per OPTIONAL_DBS entry a selected step actually reads, so a plan that
    cannot use a cache does not report on it. Three verdicts: OK (uploaded with the
    dataset), BUILD (absent, but the step that writes it is in the plan -- needs
    internet), ABSENT (absent and nothing will build it, so say what the run loses).
    """
    names = {st["name"] for st in steps}
    lines = []
    for key, spec in OPTIONAL_DBS.items():
        if not names.intersection(spec["read_by"]):
            continue
        where = f"databases/{spec['file']}"
        if (input_root / "databases" / spec["file"]).exists():
            lines.append(f" {key:<11}: OK     {where}")
        elif spec["built_by"] in names:
            lines.append(f" {key:<11}: BUILD  {where} -- step '{spec['built_by']}' writes it "
                         f"(needs internet)")
        else:
            lines.append(f" {key:<11}: ABSENT {where} -- {spec['without']}")
    return lines


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
    # Only demand sentences/ for steps that actually read it. drug_lexicon downloads a
    # vocabulary and touches no corpus, so requiring the corpus to build it is a false gate --
    # and building it BEFORE a first run, when sentences/ does not exist yet, is the normal case.
    elif (any(st.get("reads_sentences", True) for st in steps)
          and not (input_root / "sentences").is_dir()):
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
        # actually present -- an absent optional checkpoint is simply not exported. The
        # matching *_data/ TSVs come along so a reused bundle stays re-calibratable.
        for st in steps:
            for m in st["models"]:
                if m not in produced and model_dir_ok(input_root / m):
                    link_or_copy(input_root / m, work_root / m)
                    data = input_root / MODEL_DATA.get(m, "")
                    if MODEL_DATA.get(m) and data.is_dir():
                        link_or_copy(data, work_root / data.name)
    outs = list(OUT_DIRS) + (["sentences", "summaries"] if any(st["name"] == "sentences" for st in steps) else [])
    for d in outs:
        (work_root / d).mkdir(parents=True, exist_ok=True)


def main():
    ap = argparse.ArgumentParser(description="RE-model training + normalization + relation-extraction pipeline (Kaggle/local, GPU-enabled)")
    ap.add_argument("--list", "--dry-run", dest="dry", action="store_true",
                    help="show the plan and exit")
    ap.add_argument("--steps", default=os.environ.get("NORM_STEPS"),
                    help="comma-separated subset of step names")
    ap.add_argument("--retrain", action="store_true", default=bool(os.environ.get("NORM_RETRAIN")),
                    help="run the training steps even when the four training dirs "
                         "({ppi,biored}-biobert-re/ + {ppi,biored}_data/) are already complete "
                         "in the input root (also env NORM_RETRAIN)")
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

    # The training gate needs the input root, so resolve it before the plan is final.
    input_root = find_input_root(args.input_root)
    forced_by = ", ".join(r for r in (
        "--retrain" if args.retrain else "",
        "--steps" if only and (only & set(TRAINING_STEPS)) else "",
        "--re-args" if args.re_args else "",
        "--biored-args" if args.biored_args else "") if r)
    steps, gate_report = apply_training_gate(steps, input_root, forced_by)

    for flag, value, target in (("--re-args", args.re_args, "re_pipeline"),
                                ("--biored-args", args.biored_args, "re_pipeline_biored")):
        if not value:
            continue
        st = next((s for s in steps if s["name"] == target), None)
        if st is None:
            sys.exit(f"ERROR: {flag} was given but the '{target}' step is not in the plan "
                     f"(it is excluded by --steps).")
        st["args"] = list(st["args"]) + shlex.split(value)

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
    for line in gate_report:
        print(line)
    for line in optional_db_report(input_root, steps):
        print(line)
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

    def step_env(st):
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
        return env

    together = parallel_training(steps)
    ran_early = {}      # step name -> (rc, seconds) for steps already run side by side
    laned = set()       # steps run by run_lanes (which prints their own done/failed line)
    results, t0 = [], time.time()
    for i, st in enumerate(steps):
        lanes = (lane_block(steps, i) if st.get("lane") and st["name"] not in ran_early
                 and parallel_lanes_enabled() else {})
        if st["name"] in ran_early:
            rc, dt = ran_early.pop(st["name"])
        elif len(lanes) >= 2:
            print(f"\n[parallel] lanes run side by side, each in order: "
                  + "  |  ".join(" -> ".join(s["name"] for s in ls) for ls in lanes.values())
                  + "  (NORM_PARALLEL_CPU=0 runs them in turn)", flush=True)
            ran_early = run_lanes(lanes, str(work_root), step_env)
            laned.update(ran_early)
            rc, dt = ran_early.pop(st["name"])
        elif st["name"] in together:
            group = [s for s in steps if s["name"] in together]
            jobs = [(s["name"], [sys.executable, s["script"], *s["args"], "--gpus", str(g)], step_env(s))
                    for g, s in enumerate(group)]
            print(f"\n[parallel] {' + '.join(together)} run side by side, one GPU each "
                  f"(NORM_PARALLEL_TRAIN=0 runs them in turn)", flush=True)
            ran_early = run_side_by_side(jobs, str(work_root))
            rc, dt = ran_early.pop(st["name"])
        else:
            env = step_env(st)
            t = time.time()
            rc = subprocess.run([sys.executable, st["script"], *st["args"]], cwd=str(work_root), env=env).returncode
            dt = time.time() - t
        results.append((st["name"], rc, dt))
        if rc != 0:
            if st.get("optional"):
                print(f"\n[!] optional step '{st['name']}' failed (exit {rc}) after {dt:.1f}s -- continuing.", flush=True)
                continue
            print(f"\n!! step '{st['name']}' FAILED (exit {rc}) after {dt:.1f}s -- stopping.", flush=True)
            # steps that already ran side by side (other lanes) still belong in the summary
            results += [(s["name"], *ran_early[s["name"]]) for s in steps if s["name"] in ran_early]
            break
        if st["name"] not in laned:
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
