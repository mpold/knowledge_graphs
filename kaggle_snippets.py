CELL 1:
!pip install -q "datasets<4" bioc
# needed only when the training gate lets steps 1-2 run (cell 2 tells you); harmless otherwise

CELL 2:
import glob, os
hits = glob.glob("/kaggle/input/**/gpu.py", recursive=True)
assert hits, "gpu.py not found under /kaggle/input -- is the dataset attached?"
INPUT_ROOT = os.path.dirname(hits[0])
print("INPUT_ROOT =", INPUT_ROOT)
for need in ("roman.py", "experimental_ner", "databases/hgnc_complete_set_2026-05-01.json",
             "databases/mondo-clingen.json", "databases/chebi.json"):
    print(("  OK   " if os.path.exists(f"{INPUT_ROOT}/{need}") else "  MISSING "), need)
# optional, not a hard requirement -- step 3 builds it (internet ON); without it sentences.py
# falls back to the INN-stem regex and chemical.py cannot normalize the biologics
for opt in ("databases/ncit_drugs.json", "databases/pmc_years.json"):
    print(("  OK   " if os.path.exists(f"{INPUT_ROOT}/{opt}") else "  ABSENT "), opt, "[optional]")
# the training gate: all four dirs complete -> steps 1-2 are skipped and the uploaded
# checkpoints are reused as-is. All four or none -- one TRAINS line retrains both models.
ckpt = ("config.json", "tokenizer.json", "tokenizer_config.json",
        "calibration.json", "model.safetensors")
tsv = ("train.tsv", "dev.tsv", "test.tsv")
gate = {"ppi-biobert-re": ckpt, "ppi_data": tsv,
        "biored-biobert-re": ckpt, "biored_data": tsv}
for d, files in gate.items():
    have = [f for f in files
            if os.path.exists(f"{INPUT_ROOT}/{d}/{f}") and os.path.getsize(f"{INPUT_ROOT}/{d}/{f}")]
    print(("  OK    " if len(have) == len(files) else "  TRAINS"),
          f"{d}/  ({len(have)}/{len(files)} files)")

CELL 3:
!nvidia-smi -L
!python {INPUT_ROOT}/gpu.py --input-root {INPUT_ROOT} --list
# the "training   :" line in the header prints SKIP or RUN, and its reason

CELL 4:
!python {INPUT_ROOT}/bigbio_to_re.py --task biored --dataset bigbio/biored --print-types

CELL 5:
!bash -c "set -o pipefail; PYTHONUNBUFFERED=1 TQDM_MININTERVAL=30 python {INPUT_ROOT}/gpu.py --input-root {INPUT_ROOT} 2>&1 | tee /tmp/run.log"
print("pipeline exit code:", _exit_code)   # 0 = every required step succeeded
import shutil, zipfile
shutil.copy("/tmp/run.log", "/kaggle/working/run.log")
with zipfile.ZipFile("/kaggle/working/kaggle_working.zip", "a", zipfile.ZIP_DEFLATED) as z:
    z.write("/kaggle/working/run.log", "run.log")
# --retrain trains the two RE models even when the gate finds them complete
# Output streams live AND is saved (stderr too: warnings, tracebacks). PYTHONUNBUFFERED keeps a
# piped run from printing in bursts; TQDM_MININTERVAL=30 stops progress bars flooding the log.
# pipefail makes the exit code the run's, not tee's -- '!' never stops a cell, so it is printed.
# The log is written to /tmp during the run: the pipeline's own zip step archives all of
# /kaggle/working, and would otherwise bundle a half-written copy. The finished log is then
# copied to /kaggle/working/run.log (Output panel) and added to kaggle_working.zip once.


==============================================================================
STEP 2 ONLY -- retrain the BioRED model with extra Bind rows from GE/PC events
(bigbio_to_re.py --extra-bind). Accelerator GPU + Internet ON. Afterwards copy
biored-biobert-re/ and biored_data/ from kaggle_working.zip into gpu_bundle/ and
re-upload; the training gate then reuses both models on a normal run.
==============================================================================

STEP 2 CELL 1:
!pip install -q "datasets<4" bioc
# datasets<4: BigBIO ships as loading scripts; bioc: BioRED's BioC XML parser

STEP 2 CELL 2:
import glob, os
hits = glob.glob("/kaggle/input/**/gpu.py", recursive=True)
assert hits, "gpu.py not found under /kaggle/input -- is the dataset attached?"
INPUT_ROOT = os.path.dirname(hits[0])
print("INPUT_ROOT =", INPUT_ROOT)
# step 2 only needs the training scripts -- no corpus, sentences/ or databases
for need in ("gpu.py", "run_re_pipeline.py", "bigbio_to_re.py", "train_re.py",
             "calibration.py", "zip_work.py"):
    print(("  OK   " if os.path.exists(f"{INPUT_ROOT}/{need}") else "  MISSING "), need)
# catch a stale upload: the old converter has no --extra-bind and would fail mid-run
for f, mark in (("bigbio_to_re.py", "def iter_bind_instances"),
                ("run_re_pipeline.py", "--extra-bind"),
                ("gpu.py", 'optional=True, reads_sentences=False,\n         support=["bigbio_to_re.py"')):
    ok = mark in open(f"{INPUT_ROOT}/{f}", encoding="utf-8").read()
    print(("  NEW  " if ok else "  OLD -- re-upload "), f)

STEP 2 CELL 3:
!nvidia-smi -L
!python {INPUT_ROOT}/gpu.py --input-root {INPUT_ROOT} --list --steps re_pipeline_biored,zip --biored-args "--extra-bind bigbio/bionlp_st_2013_ge,bigbio/bionlp_st_2013_pc"
# expect:  training   : RUN (re_pipeline_biored) -- forced by --steps, --biored-args
#          steps      : re_pipeline_biored, zip

STEP 2 CELL 4:
!bash -c "set -o pipefail; PYTHONUNBUFFERED=1 TQDM_MININTERVAL=30 python {INPUT_ROOT}/gpu.py --input-root {INPUT_ROOT} --steps re_pipeline_biored,zip --biored-args '--extra-bind bigbio/bionlp_st_2013_ge,bigbio/bionlp_st_2013_pc' 2>&1 | tee /tmp/run.log"
print("pipeline exit code:", _exit_code)   # 0 = step 2 and the zip succeeded
import shutil, zipfile
shutil.copy("/tmp/run.log", "/kaggle/working/run.log")
with zipfile.ZipFile("/kaggle/working/kaggle_working.zip", "a", zipfile.ZIP_DEFLATED) as z:
    z.write("/kaggle/working/run.log", "run.log")
# Look for two "[extra-bind] ... Bind + ... false -> train.tsv" lines (GE ~271, PC ~514 Bind)
# and "test: no events annotated -- skipped" for each -- that is expected.

STEP 2 CELL 5:
import collections, os
W = "/kaggle/working"
for split in ("train", "dev", "test"):
    labs = collections.Counter(l.rstrip("\n").split("\t")[-1]
                               for l in open(f"{W}/biored_data/{split}.tsv", encoding="utf-8").readlines()[1:])
    print(f"{split:5s}", dict(labs))
# train Bind should be ~458 + ~785 = ~1,240; dev/test Bind unchanged from the old run
for f in ("config.json", "model.safetensors", "calibration.json"):
    print(("  OK   " if os.path.exists(f"{W}/biored-biobert-re/{f}") else "  MISSING "), f"biored-biobert-re/{f}")
!ls {W}/summaries 2>/dev/null
# calibration_biored.html: compare Bind precision/recall and dev F1 against the old .6238
