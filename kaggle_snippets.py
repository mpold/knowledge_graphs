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
