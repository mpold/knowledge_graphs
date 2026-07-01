"""sentences.py

Find, in every XML file in ``experimental_ner/``, the sentences that **describe
original results** of the study -- scanning ANY section, plus **table headings
and figure legends** -- using **BioBERT** (the biomedical BERT encoder) as the
semantic engine, then run **BioBERT named-entity recognition** on each kept
sentence. Writes one JSON file per input document to ``sentences/`` and an HTML
summary to ``summaries/sentences.html`` with SEPARATE statistics for the JATS and
TEI corpora.

WHY BioBERT (and how it is used)
--------------------------------
BioBERT (``dmis-lab/biobert-v1.1``, a BERT-base encoder pre-trained on PubMed +
PMC) is a *bidirectional encoder*: unlike a generative model it cannot be prompted
to answer "Yes/No", and there is no off-the-shelf, fine-tuned "is-this-a-result"
head. It is therefore used in the natural unsupervised way for an encoder -- as a
sentence-embedding model:

  * Two small, hand-written ANCHOR sets are embedded once: POSITIVE anchors that
    exemplify original results ("We found that ... significantly increased ...",
    "Tumor volume was reduced by 42% ... (p < 0.01)", ...) and NEGATIVE anchors
    that exemplify background / methods / prior work ("Glioblastoma is the most
    common ...", "Cells were cultured in DMEM ...", "Previous studies have shown
    ...", ...). Their mean (L2-normalized) gives a positive and a negative
    centroid in BioBERT embedding space.
  * Each candidate sentence is embedded (mean-pooled, masked, L2-normalized) and
    scored by ``margin = cos(sentence, POS_centroid) - cos(sentence, NEG_centroid)``.
    The sentence is kept as an original result only when ``margin >= RESULT_MARGIN``
    -- i.e. it sits closer to the "result" region of BioBERT space than to the
    "background/method" region.

BioBERT is an encoder, so a single forward pass per sentence (batched) suffices --
no per-token generation -- which makes this CPU-feasible (BioBERT-base ~110M
parameters vs a multi-billion-parameter generative LM).

NAMED-ENTITY RECOGNITION (BioBERT)
----------------------------------
The base encoder above has no token-classification head, so NER uses BioBERT
models that were *fine-tuned* for token classification. Three BioBERT-based NER
models run over every KEPT result sentence and their entity spans are merged:

  * ``alvaroalon2/biobert_diseases_ner``  -> DISEASE entities
  * ``alvaroalon2/biobert_genetic_ner``   -> GENE / protein entities
  * ``alvaroalon2/biobert_chemical_ner``  -> CHEMICAL / drug entities

Each is a ``dmis-lab/biobert-v1.1`` fine-tune. They are run via the Transformers
``token-classification`` pipeline with ``aggregation_strategy="first"`` -- which
aggregates every WordPiece sub-token of a word together (one label per whole word,
taken from its first sub-token) BEFORE entity grouping, so fragments like
``"G" + "##lioblastoma"`` are stitched back into a single ``"Glioblastoma"`` span.
(``"simple"`` does NOT do this word-level merge: it groups only already-formed
tokens, so a noisy B-/I- prediction across sub-tokens leaks raw ``"##..."`` pieces
as separate entities.) Every entity records its surface text,
label, contributing model domain, confidence and character offsets within the
sentence. The model list is overridable via ``BIOBERT_NER_MODELS`` (comma list)
and NER can be switched off for a fast structural smoke test with
``BIOBERT_NER=0``.

Sentence *segmentation* (which BioBERT does not provide) uses a biomedical-
abbreviation-aware rule splitter. FRAGMENTS / poorly formed strings are excluded
deterministically by a strict surface filter (terminal punctuation, >= 4 words,
>= 25 chars, clean start, lowercase letters present, mostly alphabetic) before any
embedding -- BioBERT decides the *semantics* (result vs not), the surface filter
decides *well-formedness*.

PIPELINE
--------
1. SOURCE -- ANY section is scanned: every JATS <sec> / TEI <div> body section,
   every <abstract>, any loose <p> directly under <body>, PLUS -- per this step --
   every figure legend (JATS <fig><caption>, TEI <figure><figDesc>) and every
   table heading (JATS <table-wrap> label+caption and <thead>, TEI
   <figure type="table"> head/figDesc and header row). Reference lists, footnotes
   and title elements are never visited, so publication / section titles are
   structurally excluded.

2. SEGMENTATION -- each text block is split into sentences with an
   abbreviation-aware rule splitter.

3. PRE-SCREEN (cheap) -- a surface filter removes fragments and ALL-CAPS headings;
   prior-work attributions ("has been shown", "previous studies", ...) and
   references / URLs / DOIs are removed by regex.

4. ORIGINAL RESULT (cue gate + BioBERT) -- surviving candidates are embedded with
   BioBERT and kept only when they (a) match >= 1 result cue (finding / change /
   relation / statistic regex -- a hard precision gate) AND (b) reach the
   positive-vs-negative anchor margin >= RESULT_MARGIN. The cue gate supplies
   precision; BioBERT re-ranks within the cue-bearing candidates (raw mean-pooled
   embeddings are anisotropic, so the margin alone is too weak to use on its own).

5. NER (BioBERT) -- each kept result sentence is passed through the three BioBERT
   token-classification models above; merged DISEASE / GENE / CHEMICAL spans (with
   text, label, domain, score and char offsets) are attached to that sentence.

OUTPUT
------
  * ``sentences/<input-basename>.json`` -- one per input file (written even when
    zero sentences are extracted, so corpus coverage is auditable). Each kept
    sentence carries its result-scoring metadata AND an ``entities`` list from
    BioBERT NER.
  * ``summaries/sentences.html`` -- the generating prompt, the strategy, the
    BioBERT semantic + NER models / device / threshold used, SEPARATE JATS vs TEI
    summaries, the per-section distribution (incl. figure-legend / table-heading),
    the entity-label distribution per format, and the full list of zero-sentence
    output files.

ENVIRONMENT
-----------
Requires PyTorch + Transformers (already present in the project venv), the BioBERT
encoder weights (``dmis-lab/biobert-v1.1``) and the three BioBERT NER models
(``alvaroalon2/biobert_{diseases,genetic,chemical}_ner``) -- all cached locally on
first use. Set the optional ``HF_TOKEN`` env var for higher Hugging Face Hub rate
limits / faster downloads; ``hf_transfer`` is enabled automatically when installed.

Run from anywhere (paths resolve relative to this file)::

    python sentences.py

(Per the task, this script is NOT executed by the author.)
"""
import datetime
import glob
import html
import os
import re
import time

from lxml import etree

# ---------------------------------------------------------------------------
# Optional Hugging Face token (BioBERT is public, so a token is NOT required).
# When set, it raises Hub rate limits and speeds downloads. Read from the
# environment only -- no secret is embedded in source. ``hf_transfer`` (if
# installed) is enabled for accelerated downloads.
# ---------------------------------------------------------------------------
HF_TOKEN = os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN") or ""
try:
    import hf_transfer  # noqa: F401
    os.environ.setdefault("HF_HUB_ENABLE_HF_TRANSFER", "1")
except Exception:
    pass

# --- the generating prompt (embedded verbatim in the HTML summary) --------
PROMPT = """Using BioBERT, find sentences that describe original results;
Exclude titles;
Include table headings and figure legends;
Exclude fragments and poorly formed sentences;
Using BioBERT, carry out named entity recognition on each sentence;
write input-specific output .json files to 'sentences';
Include this prompt and stategy in 'summaries/sentences.html' that summarizes this step including output files with zero sentences;
Summarize the JATS and TEI xml files separately in 'summaries/sentences.html';
Set up BioBERT environment if needed;
Write the script called 'sentences.py';
Do not execute this script;"""

STRATEGY = [
    "Parse every *.xml in experimental_ner with lxml (recover=True, no DTD/network) "
    "and auto-detect the schema from the root element: TEI (GROBID) vs JATS/PMC.",
    "Scan ANY section as a source: every JATS <sec> / TEI <div> body section, every "
    "<abstract>, and any loose <p> directly under <body>. Nested sub-sections are "
    "collected with their enclosing top-level section (no double counting).",
    "INCLUDE table headings and figure legends: JATS figure legends (<fig><caption> "
    "+ label) and TEI figure legends (<figure><figDesc>/<head>); JATS table headings "
    "(<table-wrap> label+caption and <thead> cells) and TEI table headings "
    "(<figure type=\"table\"> head/figDesc and the header row). These are tagged with "
    "the section types 'figure-legend' and 'table-heading'.",
    "EXCLUDE titles: title elements (article-title, section <title>/<head>, "
    "teiHeader titles) are never read; any title-cased fragment that slips through "
    "lacks terminal punctuation and is rejected by the surface filter.",
    "Segment each text block into sentences with an abbreviation-aware rule splitter "
    "(decimals, 'e.g.', 'i.e.', 'et al.', 'vs.', single-letter initials, ... are "
    "protected) -- BioBERT is an encoder and provides no sentence boundaries.",
    "EXCLUDE fragments / poorly formed sentences deterministically with a surface "
    "filter (terminal punctuation, >= 4 words, >= 25 chars, clean start, lowercase "
    "letters present, mostly alphabetic) before any embedding; prior-work "
    "attributions and references / URLs / DOIs are removed by regex.",
    "Decide ORIGINAL RESULT with a cue hard-gate plus BioBERT embeddings. A "
    "candidate must FIRST match >= 1 result cue (finding / change / relation / "
    "statistic regex) -- this supplies precision. Positive anchor sentences "
    "(original results) and negative anchor sentences (background / methods / prior "
    "work) are embedded once into L2-normalized centroids; each cue-bearing "
    "candidate is mean-pooled, masked and L2-normalized, then scored by margin = "
    "cos(sent, POS) - cos(sent, NEG) and kept only when margin >= RESULT_MARGIN. "
    "BioBERT thus re-ranks within the cue-bearing set rather than classifying alone "
    "(raw mean-pooled embeddings are anisotropic, so the margin is a weak signal).",
    "Record, as metadata, the cosine similarities to each centroid, the matched "
    "result-cue groups and the source section's canonical type. Content -- not the "
    "section label -- decides whether a sentence is kept.",
    "De-duplicate identical sentence strings within each document.",
    "Run BioBERT NER on each kept result sentence. The base encoder has no "
    "token-classification head, so three BioBERT fine-tunes are used -- "
    "alvaroalon2/biobert_diseases_ner (DISEASE), biobert_genetic_ner (GENE/protein) "
    "and biobert_chemical_ner (CHEMICAL/drug) -- via the Transformers "
    "token-classification pipeline with aggregation_strategy='first' (every WordPiece "
    "sub-token of a word merged into one whole-word span before entity grouping, so "
    "fragments like 'G'+'##lioblastoma' rejoin into 'Glioblastoma'; 'simple' does not "
    "do this word-level merge and leaks raw '##' sub-tokens). Each sentence's merged "
    "entity spans (text, label, model domain, "
    "score, char offsets) are attached to it; the model list is overridable via "
    "BIOBERT_NER_MODELS and NER is skippable via BIOBERT_NER=0.",
    "Write one JSON per input file to sentences/ (even when empty), then write the "
    "HTML summary with SEPARATE JATS and TEI statistics, the per-section distribution "
    "of result sentences, the BioBERT-NER entity-label distribution per format, and "
    "the list of zero-sentence output files.",
]

# --- paths ----------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Paths are env-overridable so the script runs UNCHANGED on Kaggle, where inputs are
# read-only under /kaggle/input and outputs must go to /kaggle/working. Defaults resolve
# next to this file for a local run.
INPUT_DIR = os.environ.get("INPUT_DIR") or os.path.join(SCRIPT_DIR, "experimental_ner")
OUT_DIR = os.environ.get("OUT_DIR") or os.path.join(SCRIPT_DIR, "sentences")
SUMMARY_DIR = os.environ.get("SUMMARY_DIR") or os.path.join(SCRIPT_DIR, "summaries")
SUMMARY_PATH = os.environ.get("SUMMARY_PATH") or os.path.join(SUMMARY_DIR, "sentences.html")
# Optional: when OUTPUT_ZIP is set (e.g. on Kaggle), bundle the per-file JSONs in OUT_DIR plus
# the summary HTML into a single zip so the whole run downloads as ONE file from
# /kaggle/working. Empty (the default) = no archive, so local runs are unaffected.
OUTPUT_ZIP = os.environ.get("OUTPUT_ZIP", "").strip()

# --- BioBERT configuration -------------------------------------------------
MODEL_ID = os.environ.get("BIOBERT_MODEL", "dmis-lab/biobert-v1.1")
# Keep a candidate only when it (a) matches >= 1 result cue (hard precision gate)
# AND (b) cos(sent, POS) - cos(sent, NEG) >= this margin. Raw mean-pooled BioBERT
# embeddings are anisotropic (every sentence sits at cosine ~0.85-0.90 to both
# centroids), so the margin alone is a weak, near-zero signal -- the cue gate does
# the precision work and BioBERT re-ranks within the cue-bearing candidates.
RESULT_MARGIN = float(os.environ.get("BIOBERT_RESULT_MARGIN", "0.01"))
# Embedding micro-batch size and max tokens per sentence.
EMB_BATCH = int(os.environ.get("BIOBERT_BATCH", "32"))
EMB_MAXLEN = int(os.environ.get("BIOBERT_MAXLEN", "256"))
# How many input files to process per inference chunk. Candidates from this many
# files are gathered, then embedded in ONE pass and NER-tagged in ONE pass, so the
# BioBERT batches stay full and per-call overhead is amortized across files. This is
# purely a throughput knob -- output is identical to processing one file at a time
# (each sentence is scored independently and padding is masked). Tune for memory.
BATCH_FILES = int(os.environ.get("BIOBERT_FILE_BATCH", "16"))
# Optional cap on number of input files (handy for a smoke test); empty = all.
_MAX = os.environ.get("BIOBERT_MAX_FILES", "").strip()
MAX_FILES = int(_MAX) if _MAX else None
# Optional cap on candidates SCORED per file (smoke-test throttle); empty = all.
_MAXC = os.environ.get("BIOBERT_MAX_CANDS", "").strip()
MAX_CANDS = int(_MAXC) if _MAXC else None

_INSTALL_HINT = (
    "BioBERT requires PyTorch + Transformers. Install with:\n"
    "  pip install torch transformers huggingface_hub\n"
    "  pip install hf_transfer            # optional, faster downloads\n"
    "The weights (dmis-lab/biobert-v1.1 plus the alvaroalon2 BioBERT NER models)\n"
    "download automatically on first use."
)

# --- BioBERT NER configuration ---------------------------------------------
# The base encoder has no NER head, so NER uses BioBERT *fine-tunes* for token
# classification. Default = three BioBERT NER models (disease / gene / chemical);
# overridable via BIOBERT_NER_MODELS (comma-separated HF ids). NER is run over the
# kept result sentences only; set BIOBERT_NER=0 for a fast structural smoke test.
_NER_DEFAULT = (
    "alvaroalon2/biobert_diseases_ner,"
    "alvaroalon2/biobert_genetic_ner,"
    "alvaroalon2/biobert_chemical_ner"
)
NER_MODELS = [m.strip() for m in
              os.environ.get("BIOBERT_NER_MODELS", _NER_DEFAULT).split(",") if m.strip()]
NER_ENABLED = os.environ.get("BIOBERT_NER", "1").strip().lower() not in ("0", "false", "no")
# Drop low-confidence entity spans below this aggregated score.
NER_MIN_SCORE = float(os.environ.get("BIOBERT_NER_MIN_SCORE", "0.5"))
NER_BATCH = int(os.environ.get("BIOBERT_NER_BATCH", "16"))
# Outside / non-entity classes to drop. The pipeline only ignores "O" by default,
# but alvaroalon2/biobert_diseases_ner names its outside class "0" (the digit), so
# its non-entity spans (e.g. "in adults.") would otherwise leak in as entities.
_OUTSIDE_LABELS = {"O", "0", ""}

# Anchor sentences that define the "original result" vs "background/method/prior"
# regions of BioBERT embedding space. Hand-written, intentionally generic.
POS_ANCHORS = [
    "We found that knockdown of the gene significantly increased apoptosis in the cells.",
    "Tumor volume was reduced by 42% in treated mice compared with controls (p < 0.01).",
    "Expression of the marker was significantly higher in tumor tissue than in normal tissue.",
    "Overall survival was significantly longer in the treatment group than in the control group.",
    "We observed a significant correlation between marker expression and tumor grade.",
    "Knockout of the gene abolished cell proliferation in vitro.",
    "The combination treatment yielded a three-fold increase in response rate.",
    "Protein levels decreased markedly following drug exposure.",
]
NEG_ANCHORS = [
    "Glioblastoma is the most common malignant primary brain tumor in adults.",
    "Previous studies have shown that the receptor is frequently amplified in this cancer.",
    "Cells were cultured in DMEM supplemented with 10% fetal bovine serum.",
    "The aim of this study was to investigate the role of the gene in tumor growth.",
    "Total RNA was extracted according to the manufacturer's protocol.",
    "Brain tumors represent a major cause of cancer-related mortality worldwide.",
    "It is well known that the blood-brain barrier limits drug delivery to the brain.",
    "Statistical analyses were performed using standard software packages.",
]


def _pick_device(torch):
    """Return (device_str, torch_dtype, pipeline_device) for the best backend."""
    if torch.cuda.is_available():
        return "cuda", torch.float16, 0
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return "mps", torch.float32, "mps"
    return "cpu", torch.float32, -1


class BioBERTScorer:
    """Wrap BioBERT for original-result scoring via embedding-anchor margins."""

    def __init__(self, model_id=MODEL_ID, token=HF_TOKEN, batch_size=EMB_BATCH,
                 max_len=EMB_MAXLEN):
        try:
            import torch
            from transformers import AutoModel, AutoTokenizer
        except ImportError:
            raise SystemExit(_INSTALL_HINT)

        self.torch = torch
        self.batch_size = batch_size
        self.max_len = max_len
        self.model_id = model_id

        tok_kwargs = {"token": token} if token else {}
        self.tok = AutoTokenizer.from_pretrained(model_id, **tok_kwargs)

        self.device, self.dtype, _ = _pick_device(torch)

        self.model = AutoModel.from_pretrained(
            model_id, dtype=self.dtype, low_cpu_mem_usage=True, **tok_kwargs)
        self.model.to(self.device)
        self.model.eval()

        # Build the positive / negative centroids once.
        F = torch.nn.functional
        self.pos = F.normalize(self._embed(POS_ANCHORS).mean(0), p=2, dim=0)
        self.neg = F.normalize(self._embed(NEG_ANCHORS).mean(0), p=2, dim=0)

    def _embed(self, texts):
        """Mean-pooled, masked, L2-normalized BioBERT embeddings -> tensor [N, H]."""
        torch = self.torch
        F = torch.nn.functional
        out = []
        for i in range(0, len(texts), self.batch_size):
            chunk = texts[i:i + self.batch_size]
            enc = self.tok(chunk, padding=True, truncation=True,
                           max_length=self.max_len, return_tensors="pt")
            enc = {k: v.to(self.device) for k, v in enc.items()}
            with torch.no_grad():
                hidden = self.model(**enc).last_hidden_state          # [B, T, H]
            mask = enc["attention_mask"].unsqueeze(-1).to(hidden.dtype)
            summed = (hidden * mask).sum(dim=1)
            counts = mask.sum(dim=1).clamp(min=1e-9)
            mean = summed / counts
            out.append(F.normalize(mean.float(), p=2, dim=1).cpu())
        return torch.cat(out, dim=0) if out else torch.empty(0)

    def result_sims(self, sentences):
        """Return [(sim_pos, sim_neg), ...] cosine similarities to each centroid."""
        if not sentences:
            return []
        emb = self._embed(sentences)               # [N, H], L2-normalized
        sim_pos = (emb @ self.pos).tolist()
        sim_neg = (emb @ self.neg).tolist()
        return list(zip(sim_pos, sim_neg))


# --- domain inferred from each NER model id (for tagging merged entities) --
def _model_domain(model_id):
    mid = model_id.lower()
    if "disease" in mid:
        return "disease"
    if "genetic" in mid or "gene" in mid:
        return "gene"
    if "chemical" in mid or "drug" in mid:
        return "chemical"
    if "species" in mid or "organism" in mid:
        return "species"
    return mid.rsplit("/", 1)[-1]


class BioBERTNER:
    """Run BioBERT token-classification model(s) over sentences and merge spans.

    Base BioBERT has no NER head, so each configured model is a BioBERT *fine-tune*
    for token classification. Models are loaded with the Transformers
    ``token-classification`` pipeline (``aggregation_strategy="first"``) so every
    WordPiece sub-token of a word is merged into one whole-word span (label taken
    from the first sub-token) before entity grouping, with character offsets --
    avoiding the ``"simple"`` failure mode where a noisy B-/I- prediction across
    sub-tokens emits raw ``"##..."`` fragments as separate entities. The spans from
    every model are unioned per sentence and de-duplicated by (start, end, label).
    Best-effort: a model that fails to load is recorded in ``self.load_errors`` and
    skipped, so a single bad id never aborts the run.
    """

    def __init__(self, model_ids=NER_MODELS, token=HF_TOKEN, min_score=NER_MIN_SCORE,
                 batch_size=NER_BATCH):
        try:
            import torch
            from transformers import pipeline
        except ImportError:
            raise SystemExit(_INSTALL_HINT)

        self.min_score = min_score
        self.batch_size = batch_size
        self.pipes = []           # list of (domain, model_id, pipeline)
        self.load_errors = []     # list of (model_id, error string)
        self.device, _, pipe_dev = _pick_device(torch)
        tok_kwargs = {"token": token} if token else {}

        for mid in model_ids:
            try:
                pipe = pipeline(
                    "token-classification", model=mid, tokenizer=mid,
                    aggregation_strategy="first", device=pipe_dev, **tok_kwargs)
                self.pipes.append((_model_domain(mid), mid, pipe))
            except Exception as exc:                       # pragma: no cover
                self.load_errors.append((mid, f"{type(exc).__name__}: {exc}"))

    @property
    def model_ids(self):
        return [mid for _, mid, _ in self.pipes]

    def annotate(self, sentences):
        """Return a list (parallel to `sentences`) of merged entity-span lists."""
        per_sentence = [dict() for _ in sentences]         # key -> entity dict
        if not sentences or not self.pipes:
            return [[] for _ in sentences]

        for domain, _mid, pipe in self.pipes:
            try:
                results = pipe(sentences, batch_size=self.batch_size)
            except Exception:                              # pragma: no cover
                continue
            # pipeline returns a list-of-lists when given a list of inputs; a
            # single-input call would return one list -- normalize either way.
            if results and isinstance(results[0], dict):
                results = [results]
            for idx, ents in enumerate(results):
                bucket = per_sentence[idx]
                for e in ents or []:
                    score = float(e.get("score", 0.0))
                    if score < self.min_score:
                        continue
                    start = e.get("start")
                    end = e.get("end")
                    label = (e.get("entity_group") or e.get("entity") or "").upper()
                    if label in _OUTSIDE_LABELS:   # non-entity class -> not an entity
                        continue
                    text = e.get("word") or (
                        sentences[idx][start:end] if start is not None else "")
                    text = normalize_ws(text)
                    if not text:
                        continue
                    key = (start, end, label, domain)
                    prev = bucket.get(key)
                    if prev is None or score > prev["score"]:
                        bucket[key] = {
                            "text": text,
                            "label": label,
                            "domain": domain,
                            "score": round(score, 4),
                            "start": start,
                            "end": end,
                        }

        out = []
        for bucket in per_sentence:
            ents = sorted(bucket.values(),
                          key=lambda d: (d["start"] if d["start"] is not None else 0,
                                         d["end"] if d["end"] is not None else 0))
            out.append(ents)
        return out


# --- section-type canonicalization ----------------------------------------
_SEP_RE = re.compile(r"\s*(?:\band\b|&|/|\+|\||,)\s*")


def norm_type(value):
    """Canonicalize a section identifier (sec-type / div type / head / title)."""
    if not value:
        return ""
    v = value.strip().lower().replace("-", " ").replace("_", " ")
    v = _SEP_RE.sub("|", v)
    v = re.sub(r"\s+", " ", v).strip()
    if "|" in v:
        v = "|".join(sorted(p.strip() for p in v.split("|") if p.strip()))
    return v


ABSTRACT_LABEL = "abstract"
BODY_LABEL = "(body)"
UNTITLED_LABEL = "(untitled)"
FIGURE_LABEL = "figure-legend"
TABLE_LABEL = "table-heading"

# --- sentence segmentation (BioBERT has none) ------------------------------
WS_RE = re.compile(r"\s+")
WORD_RE = re.compile(r"[A-Za-z]{2,}")
_ABBREV = (
    "e.g|i.e|et al|vs|fig|figs|eq|eqs|ref|refs|approx|ca|cf|no|vol|nos|sp|ssp|"
    "gen|dr|prof|mr|mrs|ms|st|jr|sr|al|inc|ltd|co|dept|univ|p.o|i.v|i.p|s.c|"
    "i.c.v|q.d|b.i.d|t.i.d"
)
_ABBREV_RE = re.compile(r"\b(" + _ABBREV + r")\.", re.IGNORECASE)
_DECIMAL_RE = re.compile(r"(\d)\.(\d)")
_INITIAL_RE = re.compile(r"\b([A-Z])\.")
_SPLIT_RE = re.compile(r"(?<=[.!?])\s+(?=[\"“(\[]?[A-Z0-9])")
_DOT = "\x00DOT\x00"


def normalize_ws(text):
    return WS_RE.sub(" ", text or "").strip()


def split_sentences(text):
    """Abbreviation-aware sentence splitter for biomedical prose."""
    t = normalize_ws(text)
    if not t:
        return []
    t = _DECIMAL_RE.sub(lambda m: m.group(1) + _DOT + m.group(2), t)
    t = _ABBREV_RE.sub(lambda m: m.group(1) + _DOT, t)
    t = _INITIAL_RE.sub(lambda m: m.group(1) + _DOT, t)
    parts = _SPLIT_RE.split(t)
    return [p.replace(_DOT, ".").strip() for p in parts if p.strip()]


# --- surface filters -------------------------------------------------------
def _alpha_ratio(s):
    letters = sum(c.isalpha() for c in s)
    nonspace = sum(not c.isspace() for c in s)
    return (letters / nonspace) if nonspace else 0.0


def is_sentence(s):
    """Surface quality filter -- excludes short / malformed strings / titles."""
    s = s.strip()
    if len(s) < 25 or len(s.split()) < 4:
        return False
    if _alpha_ratio(s) < 0.5 or len(WORD_RE.findall(s)) < 2:
        return False
    if not re.search(r"[.!?][\"’”)\]]*$", s):     # terminal punctuation
        return False
    if not re.match(r'^["“(\[A-Z0-9]', s):         # clean start
        return False
    if not any(c.islower() for c in s):            # drop ALL-CAPS headings/titles
        return False
    return True


# --- result-cue metadata + prior-work / reference exclusion ----------------
CUE_GROUPS = [
    ("finding", re.compile(
        r"\b(found|observed|identified|detected|revealed|showed|shown|demonstrated|"
        r"indicated|exhibited|displayed|confirmed|noted|measured|quantified|"
        r"determined|yielded)\b|\bresulted in\b", re.IGNORECASE)),
    ("change", re.compile(
        r"\b(increased|decreased|reduced|elevated|enhanced|diminished|declined|"
        r"higher|lower|greater|upregulated|up-regulated|downregulated|down-regulated|"
        r"overexpressed|suppressed|induced|inhibited|enriched|depleted|"
        r"significantly)\b", re.IGNORECASE)),
    ("relation", re.compile(
        r"\b(associated with|correlated with|linked to|compared (with|to)|"
        r"relative to|versus|vs\.?)\b|\bno significant difference\b|\bsignificant\b",
        re.IGNORECASE)),
    ("stat", re.compile(
        r"\bp\s*[<=>]\s*0?\.\d+|\b\d+(\.\d+)?\s*%|\bn\s*=\s*\d+|"
        r"\b\d+(\.\d+)?[- ]?fold\b|±|\b95\s*%\s*ci\b|\b(or|hr|rr)\s*=\s*\d|"
        r"\br\s*=\s*-?0?\.\d+", re.IGNORECASE)),
]
PRIOR_WORK_RE = re.compile(
    r"\b(ha(s|ve) been (shown|reported|demonstrated|described|associated)|"
    r"it has been|previously (shown|reported|demonstrated|described|observed)|"
    r"previous(ly)? studies|prior studies|earlier studies|other studies|"
    r"studies have (shown|reported|demonstrated)|it is (well )?known|"
    r"are known to)\b", re.IGNORECASE)
REFERENCE_RE = re.compile(r"https?://|\bdoi:\s*10\.|doi\.org|www\.", re.IGNORECASE)


def is_prior_or_reference(s):
    return bool(REFERENCE_RE.search(s) or PRIOR_WORK_RE.search(s))


def matched_cues(s):
    return [name for name, rx in CUE_GROUPS if rx.search(s)]


# --- lxml helpers ----------------------------------------------------------
PARSER = etree.XMLParser(recover=True, resolve_entities=False, load_dtd=False,
                         no_network=True, huge_tree=True)


def lname(el):
    return etree.QName(el).localname if isinstance(el.tag, str) else None


def text_of(el):
    return normalize_ws("".join(el.itertext()))


def child_label(parent, tag_names):
    el = next((c for c in parent if lname(c) == tag_names[0]), None)
    return text_of(el) if el is not None else ""


def first_descendant(el, name):
    return next((d for d in el.iter() if lname(d) == name and d is not el), None)


def has_excluded_ancestor(el, stop, excluded):
    a = el.getparent()
    while a is not None and a is not stop:
        if lname(a) in excluded:
            return True
        a = a.getparent()
    return False


def in_ancestor(el, name):
    a = el.getparent()
    while a is not None:
        if lname(a) == name:
            return True
        a = a.getparent()
    return False


def top_level_matches(matched):
    matchset = {id(e) for e in matched}
    out = []
    for e in matched:
        a = e.getparent()
        nested = False
        while a is not None:
            if id(a) in matchset:
                nested = True
                break
            a = a.getparent()
        if not nested:
            out.append(e)
    return out


# --- candidate text-block collection ---------------------------------------
# Body prose excludes figures/tables/footnotes/captions: those are collected
# separately (as figure legends / table headings) so they are tagged correctly
# and never double-counted.
_JATS_BODY_EXCL = {"fig", "table-wrap", "table", "table-wrap-foot", "fn", "caption"}
_TEI_BODY_EXCL = {"figure"}


def _collect_paras(el, stype, stitle, body_excl, paras):
    """Append (section_type, section_title, paragraph_text) for each <p> in `el`."""
    for node in el.iter():
        if lname(node) != "p":
            continue
        if has_excluded_ancestor(node, el.getparent(), body_excl):
            continue
        txt = text_of(node)
        if txt:
            paras.append((stype, stitle, txt))


def _figures_tables_jats(root, paras):
    """Collect JATS figure legends and table headings (this step INCLUDES them)."""
    for fig in [e for e in root.iter() if lname(e) == "fig"]:
        label = child_label(fig, ["label"])
        cap = next((c for c in fig if lname(c) == "caption"), None)
        cap_txt = text_of(cap) if cap is not None else ""
        legend = " ".join(t for t in (label, cap_txt) if t).strip()
        if legend:
            paras.append((FIGURE_LABEL, label or "Figure", legend))

    for tw in [e for e in root.iter() if lname(e) == "table-wrap"]:
        label = child_label(tw, ["label"])
        cap = next((c for c in tw if lname(c) == "caption"), None)
        cap_txt = text_of(cap) if cap is not None else ""
        heading = " ".join(t for t in (label, cap_txt) if t).strip()
        if heading:
            paras.append((TABLE_LABEL, label or "Table", heading))
        for thead in [e for e in tw.iter() if lname(e) == "thead"]:
            htxt = text_of(thead)
            if htxt:
                paras.append((TABLE_LABEL, label or "Table", htxt))


def _figures_tables_tei(root, paras):
    """Collect TEI/GROBID figure legends and table headings."""
    for fg in [e for e in root.iter() if lname(e) == "figure"]:
        head = child_label(fg, ["head"])
        desc = next((c for c in fg if lname(c) == "figDesc"), None)
        desc_txt = text_of(desc) if desc is not None else ""
        is_table = (fg.get("type") or "").strip().lower() == "table"
        body = " ".join(t for t in (head, desc_txt) if t).strip()
        if is_table:
            if body:
                paras.append((TABLE_LABEL, head or "Table", body))
            tbl = next((c for c in fg if lname(c) == "table"), None)
            if tbl is not None:
                first_row = first_descendant(tbl, "row") or \
                    next((r for r in tbl.iter() if lname(r) == "row"), None)
                if first_row is not None:
                    htxt = text_of(first_row)
                    if htxt:
                        paras.append((TABLE_LABEL, head or "Table", htxt))
        else:
            if not body:
                body = text_of(fg)
            if body:
                paras.append((FIGURE_LABEL, head or "Figure", body))


def paragraphs_jats(root):
    paras = []
    for ab in [e for e in root.iter() if lname(e) == "abstract"]:
        atype = (ab.get("abstract-type") or "").strip()
        stitle = (child_label(ab, ["title"])
                  or (f"Abstract ({atype})" if atype else "Abstract"))
        _collect_paras(ab, ABSTRACT_LABEL, stitle, _JATS_BODY_EXCL, paras)

    secs = [e for e in root.iter()
            if lname(e) == "sec" and not in_ancestor(e, "abstract")]
    for sec in top_level_matches(secs):
        title = child_label(sec, ["title"])
        stype = norm_type(sec.get("sec-type")) or norm_type(title) or UNTITLED_LABEL
        _collect_paras(sec, stype, title or stype, _JATS_BODY_EXCL, paras)

    for p in root.iter():
        if (lname(p) == "p" and p.getparent() is not None
                and lname(p.getparent()) == "body"):
            txt = text_of(p)
            if txt:
                paras.append((BODY_LABEL, "Body", txt))

    _figures_tables_jats(root, paras)
    return paras


def paragraphs_tei(root):
    paras = []
    for ab in [e for e in root.iter() if lname(e) == "abstract"]:
        stitle = child_label(ab, ["head"]) or "Abstract"
        _collect_paras(ab, ABSTRACT_LABEL, stitle, _TEI_BODY_EXCL, paras)

    divs = [e for e in root.iter()
            if lname(e) == "div" and not in_ancestor(e, "abstract")]
    for div in top_level_matches(divs):
        head_text = child_label(div, ["head"])
        stype = norm_type(div.get("type")) or norm_type(head_text) or UNTITLED_LABEL
        _collect_paras(div, stype, head_text or stype, _TEI_BODY_EXCL, paras)

    for p in root.iter():
        if (lname(p) == "p" and p.getparent() is not None
                and lname(p.getparent()) == "body"):
            txt = text_of(p)
            if txt:
                paras.append((BODY_LABEL, "Body", txt))

    _figures_tables_tei(root, paras)
    return paras


def detect_and_collect(root):
    """Return (format, [(section_type, section_title, text_block), ...])."""
    rootname = lname(root)
    if rootname == "TEI":
        return "tei", paragraphs_tei(root)
    if rootname in ("pmc-articleset", "article"):
        return "jats", paragraphs_jats(root)
    names = {lname(e) for e in root.iter()}
    if "sec" in names:
        return "jats", paragraphs_jats(root)
    if "div" in names:
        return "tei", paragraphs_tei(root)
    return "unknown", []


# --- per-file processing ---------------------------------------------------
def out_name(basename):
    return re.sub(r"\.xml$", ".json", basename, flags=re.IGNORECASE)


def _base_record(basename, ner):
    return {
        "source_file": basename,
        "format": "unknown",
        "model": MODEL_ID,
        "ner_models": (ner.model_ids if ner is not None else []),
        "extraction": "original-results sentences from any section + figure legends "
                      "& table headings (BioBERT embedding-anchor margin), then "
                      "BioBERT NER on each kept sentence",
        "thresholds": {"result_margin": RESULT_MARGIN, "require_cue": True,
                       "ner_min_score": NER_MIN_SCORE},
        "section_types_matched": [],
        "n_candidates": 0,
        "n_scored": 0,
        "n_sentences": 0,
        "n_entities": 0,
        "sentences": [],
        "parse_error": None,
    }


def collect_candidates(path, ner=None):
    """Parse one file and build its candidate list (the pure-Python, no-BioBERT part).

    Returns ``(record, candidates)`` where ``candidates`` is a list of
    ``(text, section_type, section_title, cues)`` tuples, already de-duplicated,
    surface-filtered and MAX_CANDS-capped. No embedding or NER happens here -- that is
    done in batch by ``process_batch`` so the BioBERT passes span many files at once.
    """
    record = _base_record(os.path.basename(path), ner)
    try:
        root = etree.parse(path, PARSER).getroot()
    except Exception as exc:           # pragma: no cover - defensive
        record["parse_error"] = f"{type(exc).__name__}: {exc}"
        return record, []

    fmt, paras = detect_and_collect(root)
    record["format"] = fmt
    if not paras:
        return record, []

    # segment + cheap pre-screen --> candidate list (carry section metadata)
    seen = set()
    candidates = []          # (text, stype, stitle, cues)
    for stype, stitle, block in paras:
        for sent in split_sentences(block):
            text = normalize_ws(sent)
            if text in seen:
                continue
            if not is_sentence(text):              # fragments / titles -> drop
                continue
            if is_prior_or_reference(text):
                continue
            seen.add(text)
            candidates.append((text, stype, stitle, matched_cues(text)))

    record["n_candidates"] = len(candidates)
    if MAX_CANDS and len(candidates) > MAX_CANDS:
        candidates = candidates[:MAX_CANDS]          # smoke-test throttle
    record["n_scored"] = len(candidates)
    return record, candidates


def process_batch(paths, scorer, ner=None):
    """Process a chunk of files with ONE BioBERT embedding pass and ONE NER pass.

    Candidates from every file in ``paths`` are flattened into a single list, embedded
    together (``scorer.result_sims``), gated (cue + margin), and the kept sentences from
    the whole chunk are flattened again for a single ``ner.annotate`` call. Results are
    sliced back to their originating files. This is throughput-only: because each
    sentence is scored independently and padding is masked, the per-file records are
    byte-for-byte identical to calling ``process_file`` on each path in turn.
    """
    records, per_file_cands = [], []
    for path in paths:
        rec, cands = collect_candidates(path, ner)
        records.append(rec)
        per_file_cands.append(cands)

    # 1) ONE embedding pass over every candidate in the chunk
    flat_texts = [c[0] for cands in per_file_cands for c in cands]
    if flat_texts:
        print(f"    embedding {len(flat_texts)} candidate(s) across "
              f"{len(paths)} file(s) with BioBERT ...", flush=True)
    flat_sims = scorer.result_sims(flat_texts)

    # 2) slice sims back per file and apply the cue + margin gate
    pos = 0
    for rec, cands in zip(records, per_file_cands):
        sims = flat_sims[pos:pos + len(cands)]
        pos += len(cands)
        matched_types = set()
        for (text, stype, stitle, cues), (sim_pos, sim_neg) in zip(cands, sims):
            if not cues:                           # cue hard-gate (precision) -> drop
                continue
            margin = sim_pos - sim_neg
            if margin < RESULT_MARGIN:             # weak BioBERT re-rank -> drop
                continue
            matched_types.add(stype)
            rec["sentences"].append({
                "text": text,
                "section_type": stype,
                "section_title": stitle,
                "cues": cues,
                "sim_pos": round(sim_pos, 4),
                "sim_neg": round(sim_neg, 4),
                "result_margin": round(margin, 4),
                "entities": [],
            })
        rec["section_types_matched"] = sorted(matched_types)
        rec["n_sentences"] = len(rec["sentences"])

    # 3) ONE NER pass over every kept sentence in the chunk
    if ner is not None:
        flat_kept = [s["text"] for rec in records for s in rec["sentences"]]
        if flat_kept:
            print(f"    NER on {len(flat_kept)} kept sentence(s) across "
                  f"{len(paths)} file(s) with {len(ner.pipes)} model(s) ...",
                  flush=True)
            flat_ents = ner.annotate(flat_kept)
            k = 0
            for rec in records:
                n_ent = 0
                for s in rec["sentences"]:
                    s["entities"] = flat_ents[k]
                    k += 1
                    n_ent += len(s["entities"])
                rec["n_entities"] = n_ent

    return records


def process_file(path, scorer, ner=None):
    """Single-file convenience wrapper around :func:`process_batch`."""
    return process_batch([path], scorer, ner)[0]


# --- HTML summary ----------------------------------------------------------
def pct(n, d):
    return f"{(100.0 * n / d):.1f}%" if d else "0%"


def empty_fmt_stats():
    return {
        "files": 0,
        "files_with_sentences": 0,
        "zero_files": 0,
        "total_candidates": 0,
        "total_sentences": 0,
        "total_entities": 0,
        "sentences_with_entities": 0,
        "cue": {"finding": 0, "change": 0, "relation": 0, "stat": 0},
        "type_files": {},
        "type_sentences": {},
        "entity_labels": {},
    }


def render_html(agg, generated_at, model_name, device):
    total = agg["total_files"]
    jats = agg["by_format"]["jats"]
    tei = agg["by_format"]["tei"]
    with_sent = jats["files_with_sentences"] + tei["files_with_sentences"]
    total_sent = jats["total_sentences"] + tei["total_sentences"]
    total_cand = jats["total_candidates"] + tei["total_candidates"]
    total_ent = jats["total_entities"] + tei["total_entities"]
    zero_files = sorted(agg["zero_files"])

    strategy_html = "\n".join(f"<li>{html.escape(s)}</li>" for s in STRATEGY)

    def _fmt_row(name, fs):
        c = fs["cue"]
        return (
            f"<tr><td class='label'>{name}</td>"
            f"<td class='num'>{fs['files']:,}</td>"
            f"<td class='num'>{fs['files_with_sentences']:,}</td>"
            f"<td class='num'>{fs['zero_files']:,}</td>"
            f"<td class='num'>{fs['total_candidates']:,}</td>"
            f"<td class='num'>{fs['total_sentences']:,}</td>"
            f"<td class='num'>{c['finding']:,}</td>"
            f"<td class='num'>{c['change']:,}</td>"
            f"<td class='num'>{c['relation']:,}</td>"
            f"<td class='num'>{c['stat']:,}</td>"
            f"<td class='num'>{fs['total_entities']:,}</td></tr>"
        )
    combined = {
        "files": jats["files"] + tei["files"],
        "files_with_sentences": with_sent,
        "zero_files": jats["zero_files"] + tei["zero_files"],
        "total_candidates": total_cand,
        "total_sentences": total_sent,
        "total_entities": total_ent,
        "cue": {k: jats["cue"][k] + tei["cue"][k]
                for k in ("finding", "change", "relation", "stat")},
    }
    fmt_html = "\n".join([
        _fmt_row("JATS", jats),
        _fmt_row("TEI", tei),
        _fmt_row("<strong>All</strong>", combined),
    ])

    # --- entity-label distribution by format (BioBERT NER) ---
    all_labels = set(jats["entity_labels"]) | set(tei["entity_labels"])
    ranked_labels = sorted(all_labels,
                           key=lambda l: -(jats["entity_labels"].get(l, 0)
                                           + tei["entity_labels"].get(l, 0)))
    ent_rows = []
    for lbl in ranked_labels:
        je = jats["entity_labels"].get(lbl, 0)
        te = tei["entity_labels"].get(lbl, 0)
        ent_rows.append(
            f"<tr><td class='label'>{html.escape(lbl)}</td>"
            f"<td class='num'>{je:,}</td><td class='num'>{te:,}</td>"
            f"<td class='num'>{je + te:,}</td></tr>"
        )
    ent_html = "\n".join(ent_rows) or "<tr><td colspan='4'>(none)</td></tr>"
    jats_sent_ent = jats["sentences_with_entities"]
    tei_sent_ent = tei["sentences_with_entities"]
    sent_with_ent = jats_sent_ent + tei_sent_ent

    ner_models = agg.get("ner_models", [])
    ner_pills = "".join(
        f'<span class="pill">NER: {html.escape(m)}</span>\n' for m in ner_models)
    if not ner_models:
        ner_pills = ('<span class="pill">NER: disabled</span>'
                     if not agg.get("ner_enabled")
                     else '<span class="pill">NER: no models loaded</span>')
    ner_err_html = "".join(
        f"{html.escape(m)} -> {html.escape(e)}\n" for m, e in agg.get("ner_errors", []))

    all_types = set(jats["type_sentences"]) | set(tei["type_sentences"])
    ranked = sorted(all_types,
                    key=lambda t: -(jats["type_sentences"].get(t, 0)
                                    + tei["type_sentences"].get(t, 0)))
    type_rows = []
    for t in ranked:
        jf, js = jats["type_files"].get(t, 0), jats["type_sentences"].get(t, 0)
        tf, ts = tei["type_files"].get(t, 0), tei["type_sentences"].get(t, 0)
        type_rows.append(
            f"<tr><td class='label'>{html.escape(t)}</td>"
            f"<td class='num'>{jf:,}</td><td class='num'>{js:,}</td>"
            f"<td class='num'>{tf:,}</td><td class='num'>{ts:,}</td>"
            f"<td class='num'>{js + ts:,}</td></tr>"
        )
    type_html = "\n".join(type_rows) or "<tr><td colspan='6'>(none)</td></tr>"

    zero_list_html = "".join(f"{html.escape(f)}\n" for f in zero_files) or "(none)"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Original-result sentences (BioBERT) &mdash; experimental_ner</title>
<style>
  :root {{ --bg:#0f1117; --card:#1a1d27; --fg:#e6e8ee; --muted:#9aa3b2;
           --accent:#6ea8fe; --border:#2a2e3a; }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--bg); color:var(--fg);
          font:15px/1.55 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif; }}
  .wrap {{ max-width:980px; margin:0 auto; padding:32px 24px 64px; }}
  h1 {{ font-size:24px; margin:0 0 4px; }}
  h2 {{ font-size:19px; margin:36px 0 12px; border-bottom:1px solid var(--border);
        padding-bottom:6px; }}
  p.sub {{ color:var(--muted); margin:0 0 16px; }}
  .prompt {{ background:#15212e; border:1px solid #25425c; border-left:4px solid var(--accent);
             border-radius:8px; padding:12px 16px; margin:12px 0; color:#cfe0f5;
             font-size:13px; white-space:pre-wrap; word-break:break-word;
             font-family:ui-monospace,Consolas,monospace; }}
  .prompt .tag {{ display:block; color:var(--accent); font-weight:700;
                  text-transform:uppercase; letter-spacing:.05em; font-size:11px;
                  margin-bottom:6px; font-family:inherit; }}
  .pill {{ display:inline-block; background:#1d2b1f; border:1px solid #2f5a36;
           color:#bfe6c6; border-radius:999px; padding:2px 10px; font-size:12px; }}
  ol.strategy li {{ margin:6px 0; }}
  .cards {{ display:flex; flex-wrap:wrap; gap:14px; margin:18px 0 8px; }}
  .card {{ background:var(--card); border:1px solid var(--border); border-radius:10px;
           padding:14px 18px; min-width:140px; flex:1; }}
  .card .k {{ font-size:26px; font-weight:700; color:var(--accent); }}
  .card .v {{ color:var(--muted); font-size:13px; margin-top:2px; }}
  table {{ width:100%; border-collapse:collapse; background:var(--card);
           border:1px solid var(--border); border-radius:10px; overflow:hidden; }}
  th,td {{ padding:8px 12px; text-align:left; border-bottom:1px solid var(--border); }}
  th {{ background:#222633; color:var(--muted); font-weight:600; font-size:13px;
        text-transform:uppercase; letter-spacing:.03em; }}
  th.num {{ text-align:right; }}
  tr:last-child td {{ border-bottom:none; }}
  td.num {{ text-align:right; font-variant-numeric:tabular-nums; color:var(--muted);
            white-space:nowrap; }}
  td.label {{ font-weight:500; }}
  tbody tr:hover {{ background:#20242f; }}
  .note {{ color:var(--muted); font-size:13px; margin-top:10px; }}
  code {{ background:#222633; padding:1px 5px; border-radius:4px; font-size:13px; }}
  details {{ background:var(--card); border:1px solid var(--border); border-radius:10px;
             padding:10px 14px; margin-top:12px; }}
  details pre {{ max-height:340px; overflow:auto; background:#0c0e14; padding:10px;
                 border-radius:8px; font-size:12px; color:#b7c0d0; }}
</style>
</head>
<body>
<div class="wrap">
  <h1>Original-result sentences + NER (BioBERT) &mdash; <code>experimental_ner</code></h1>
  <p class="sub">Sentences describing original study results, drawn from any
  section plus figure legends &amp; table headings, judged by BioBERT embedding
  anchors, then annotated with BioBERT named-entity recognition.
  Generated {html.escape(generated_at)} by <code>sentences.py</code> &middot;
  <span class="pill">model: {html.escape(model_name)}</span>
  <span class="pill">device: {html.escape(device)}</span>
  <span class="pill">cue-gated</span>
  <span class="pill">margin&ge;{RESULT_MARGIN:g}</span>
  {ner_pills}</p>

  <div class="prompt"><span class="tag">Generating prompt</span>{html.escape(PROMPT)}</div>

  <h2>Strategy</h2>
  <ol class="strategy">
{strategy_html}
  </ol>

  <div class="cards">
    <div class="card"><div class="k">{total:,}</div>
      <div class="v">input XML files</div></div>
    <div class="card"><div class="k">{jats['files']:,}</div>
      <div class="v">JATS files</div></div>
    <div class="card"><div class="k">{tei['files']:,}</div>
      <div class="v">TEI files</div></div>
    <div class="card"><div class="k">{with_sent:,}</div>
      <div class="v">files with &ge;1 sentence ({pct(with_sent, total)})</div></div>
    <div class="card"><div class="k">{len(zero_files):,}</div>
      <div class="v">files with 0 sentences</div></div>
    <div class="card"><div class="k">{total_sent:,}</div>
      <div class="v">total result sentences</div></div>
    <div class="card"><div class="k">{total_ent:,}</div>
      <div class="v">BioBERT NER entities</div></div>
  </div>

  <h2>What <code>result_margin</code> (and the <code>{RESULT_MARGIN:g}</code> threshold) means</h2>
  <p class="note">Every candidate sentence is embedded with BioBERT (mean-pooled,
  masked, L2-normalized) and compared by cosine similarity to two centroids built
  once from the hand-written anchor sets: a <strong>positive</strong> centroid
  (the average of anchors that exemplify original results, e.g.
  <em>"Tumor volume was reduced by 42% &hellip; (p &lt; 0.01)"</em>) and a
  <strong>negative</strong> centroid (background / methods / prior-work anchors,
  e.g. <em>"Glioblastoma is the most common &hellip;"</em>). The per-sentence score
  recorded in each JSON is their difference:</p>
  <div class="prompt"><span class="tag">result_margin</span>result_margin = sim_pos &minus; sim_neg = cos(sentence, POSITIVE) &minus; cos(sentence, NEGATIVE)</div>
  <p class="note">A sentence is kept only when it BOTH (a) matches &ge;1 result cue
  (the hard precision gate) AND (b) reaches
  <code>result_margin &ge; {RESULT_MARGIN:g}</code> &mdash; i.e. it sits at least
  <code>{RESULT_MARGIN:g}</code> closer (in cosine terms) to the "result" region of
  BioBERT space than to the "background/method" region. The
  <code>thresholds.result_margin</code> field in every output file is just this
  <em>cutoff</em> that was in force for the run; each kept sentence additionally
  records its <em>own</em> <code>sim_pos</code>, <code>sim_neg</code> and
  <code>result_margin</code> so the decision is auditable
  (kept &hArr; <code>result_margin &ge; {RESULT_MARGIN:g}</code>).</p>
  <p class="note"><strong>Why the threshold is so small ({RESULT_MARGIN:g}, not e.g.
  0.5):</strong> raw mean-pooled BioBERT embeddings are <em>anisotropic</em> &mdash;
  almost every sentence sits at cosine ~0.85&ndash;0.90 to <em>both</em> centroids, so
  the margin is a weak, near-zero signal centred around 0. The result-cue regex does
  the real precision work; BioBERT only <em>re-ranks within</em> the cue-bearing
  candidates, nudging out the few that lean toward the background region. A higher
  per-sentence margin therefore is a weak tie-breaker, not a strong confidence score.
  The cutoff is tunable via the <code>BIOBERT_RESULT_MARGIN</code> environment
  variable.</p>

  <h2>Summary statistics by format (JATS vs TEI)</h2>
  <p class="note">Each input file is classified as JATS (NLM/PMC) or TEI (GROBID)
  and its original-result sentences are tallied separately. "Candidates" are the
  surface-clean, non-reference sentences embedded by BioBERT; "Sentences" are those
  that passed the result-cue hard gate AND reached the positive-vs-negative anchor
  margin. The next four columns count how many accepted sentences matched each
  result-cue group (every accepted sentence matches at least one); the last column
  is the total BioBERT-NER entities found in those sentences.</p>
  <table>
    <thead><tr><th>format</th><th class="num">Files</th>
      <th class="num">Files w/ sent.</th><th class="num">0-sentence files</th>
      <th class="num">Candidates</th><th class="num">Sentences</th>
      <th class="num">finding</th><th class="num">change</th>
      <th class="num">relation</th><th class="num">stat</th>
      <th class="num">entities</th></tr></thead>
    <tbody>
{fmt_html}
    </tbody>
  </table>

  <h2>Result sentences by source section (JATS vs TEI)</h2>
  <p class="note">Because any section is scanned -- plus figure legends
  (<code>figure-legend</code>) and table headings (<code>table-heading</code>) --
  this shows WHERE original-result sentences were found (canonical section type).
  "Files" = documents contributing at least one result sentence from that section
  type.</p>
  <table>
    <thead><tr><th>section type</th>
      <th class="num">JATS files</th><th class="num">JATS sent.</th>
      <th class="num">TEI files</th><th class="num">TEI sent.</th>
      <th class="num">Total sent.</th></tr></thead>
    <tbody>
{type_html}
    </tbody>
  </table>

  <h2>Named entities by BioBERT NER label (JATS vs TEI)</h2>
  <p class="note">Every kept result sentence is passed through the BioBERT NER
  models ({html.escape(', '.join(ner_models) or 'none loaded')}); merged spans are
  counted by label below. Of {total_sent:,} result sentences,
  {sent_with_ent:,} carry &ge;1 entity ({pct(sent_with_ent, total_sent)}) --
  {jats_sent_ent:,} in JATS and {tei_sent_ent:,} in TEI -- for {total_ent:,}
  entities in total (min score &ge; {NER_MIN_SCORE:g}). Each sentence's full entity
  list (text, label, domain, score, char offsets) is in its JSON file.</p>
  <table>
    <thead><tr><th>entity label</th>
      <th class="num">JATS ent.</th><th class="num">TEI ent.</th>
      <th class="num">Total ent.</th></tr></thead>
    <tbody>
{ent_html}
    </tbody>
  </table>
  {("<details><summary>BioBERT NER model load errors</summary><pre>"
    + html.escape(ner_err_html) + "</pre></details>") if ner_err_html else ""}

  <h2>Output files with zero sentences ({len(zero_files):,})</h2>
  <p class="note">A JSON file is still written for each of these to
  <code>sentences/</code> (with an empty <code>sentences</code> array), so corpus
  coverage is auditable. A file lands here when none of its candidate sentences,
  in any section, both match a result cue and reach the BioBERT margin.</p>

  <h3 style="margin:18px 0 8px">Why a file with real results can still yield zero</h3>
  <p class="note">This is a deliberately <strong>precision-over-recall</strong>
  pipeline: a sentence is kept only if it survives a chain of AND-ed hard gates
  (source collection &rarr; segmentation &rarr; surface filter &rarr;
  not-prior/reference &rarr; <em>&ge;1 result cue</em> &rarr;
  <em>BioBERT margin &ge; {RESULT_MARGIN:g}</em>). Any single gate can reject a
  genuine result, so a file that clearly contains result sentences can still produce
  zero. The two dominant causes:</p>
  <ul class="note" style="margin-top:0">
    <li><strong>The BioBERT margin gate.</strong> Because the embeddings are
    anisotropic, per-sentence margins are near-zero noise centred around 0, so the
    <code>{RESULT_MARGIN:g}</code> cutoff sits above the bulk of the distribution and
    rejects most cue-bearing candidates. It is compounded by <em>domain mismatch</em>:
    all positive/negative anchors are biomedical, so off-domain papers (e.g. a
    physics / imaging article) lean toward the negative centroid and score
    <em>negative</em> margins &mdash; every cue-bearing sentence can fail even when
    the paper is full of findings.</li>
    <li><strong>The result-cue hard gate.</strong> The cue regexes are narrow &mdash;
    they match mainly past-tense verbs (<code>increased</code>, <code>reduced</code>)
    and digit-based statistics. Results phrased in the present tense
    (<em>"the method increases speed"</em>), as nominalizations
    (<em>"a four-fold enhancement"</em>), or with spelled-out numbers
    (<em>"four-fold"</em> rather than <em>"4-fold"</em>) match no cue and are dropped
    before BioBERT ever sees them.</li>
  </ul>
  <p class="note" style="margin-top:0">Secondary contributors: the surface filter
  (requires terminal <code>.!?</code>, &ge;4 words and &ge;50% alphabetic characters,
  so number-dense or encoding-mangled lines are rejected &mdash; which is why TEI /
  GROBID text, often lacking final punctuation, has a higher zero-rate than JATS);
  and source-collection scope (only body sections, abstracts, loose paragraphs and
  figure/table captions are read, so results buried in atypical structures are never
  visited). To trade some precision for recall, lower
  <code>BIOBERT_RESULT_MARGIN</code> toward 0, broaden the cue regexes, and add
  domain-appropriate anchors.</p>

  <details><summary>Show {len(zero_files):,} file name(s)</summary>
<pre>{html.escape(zero_list_html)}</pre>
  </details>

  <p class="note" style="margin-top:28px">Parse errors: {agg['parse_errors']:,}.
  BioBERT semantic model: <code>{html.escape(model_name)}</code> on
  <code>{html.escape(device)}</code> (result-cue gate + margin &ge; {RESULT_MARGIN:g}).
  BioBERT NER models: <code>{html.escape(', '.join(ner_models) or 'none')}</code>.
  Output JSON per input file in <code>sentences/</code>.</p>
</div>
</body>
</html>
"""


# --- main ------------------------------------------------------------------
def main():
    files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.xml")))
    if not files:                                      # Kaggle datasets often nest the XMLs
        files = sorted(glob.glob(os.path.join(INPUT_DIR, "**", "*.xml"), recursive=True))
    if not files:
        raise SystemExit(f"No XML files found in {INPUT_DIR!r}")
    if MAX_FILES:
        files = files[:MAX_FILES]

    print(f"loading BioBERT ({MODEL_ID}) ...", flush=True)
    scorer = BioBERTScorer()               # exits with install hint if unavailable
    print(f"BioBERT model loaded  : {scorer.model_id}  (device={scorer.device})",
          flush=True)

    ner = None
    if NER_ENABLED:
        print(f"loading BioBERT NER ({', '.join(NER_MODELS)}) ...", flush=True)
        ner = BioBERTNER()                 # best-effort: bad ids recorded, not fatal
        if ner.pipes:
            print(f"BioBERT NER loaded    : {', '.join(ner.model_ids)} "
                  f"(device={ner.device})", flush=True)
        else:
            print("BioBERT NER           : NO models loaded "
                  "(entities will be empty) -- see load errors below", flush=True)
        for mid, err in ner.load_errors:
            print(f"    NER load error      : {mid} -> {err}", flush=True)
    else:
        print("BioBERT NER           : disabled (BIOBERT_NER=0)", flush=True)

    if MAX_FILES or MAX_CANDS:
        print(f"SMOKE TEST mode       : max_files={MAX_FILES} "
              f"max_cands_per_file={MAX_CANDS} (results are NOT a full run)",
              flush=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(SUMMARY_DIR, exist_ok=True)

    agg = {
        "total_files": len(files),
        "unknown_files": 0,
        "parse_errors": 0,
        "zero_files": [],
        "ner_enabled": NER_ENABLED,
        "ner_models": (ner.model_ids if ner is not None else []),
        "ner_device": (ner.device if ner is not None else "n/a"),
        "ner_errors": (ner.load_errors if ner is not None else []),
        "by_format": {
            "jats": empty_fmt_stats(),
            "tei": empty_fmt_stats(),
            "unknown": empty_fmt_stats(),
        },
    }

    import json
    n_total = len(files)
    done = 0
    for start in range(0, n_total, BATCH_FILES):
        chunk = files[start:start + BATCH_FILES]
        t0 = time.perf_counter()
        print(f"[{start + 1}-{start + len(chunk)}/{n_total}] processing "
              f"{len(chunk)} file(s) in one inference batch ...", flush=True)
        recs = process_batch(chunk, scorer, ner)
        dt = time.perf_counter() - t0

        for rec in recs:
            done += 1
            print(f"  [{done}/{n_total}] {rec['source_file']}: "
                  f"format={rec['format']} candidates={rec['n_candidates']} "
                  f"kept={rec['n_sentences']} entities={rec['n_entities']}"
                  + (f"  PARSE_ERROR={rec['parse_error']}"
                     if rec["parse_error"] else ""), flush=True)

            with open(os.path.join(OUT_DIR, out_name(rec["source_file"])), "w",
                      encoding="utf-8") as fh:
                json.dump(rec, fh, ensure_ascii=False, indent=2)

            fmt = rec["format"] if rec["format"] in ("jats", "tei") else "unknown"
            if fmt == "unknown":
                agg["unknown_files"] += 1
            if rec["parse_error"]:
                agg["parse_errors"] += 1

            fs = agg["by_format"][fmt]
            fs["files"] += 1
            fs["total_candidates"] += rec["n_candidates"]
            n = rec["n_sentences"]
            fs["total_sentences"] += n
            if n > 0:
                fs["files_with_sentences"] += 1
            else:
                fs["zero_files"] += 1
                agg["zero_files"].append(out_name(rec["source_file"]))

            fs["total_entities"] += rec["n_entities"]
            per_type = {}
            for s in rec["sentences"]:
                per_type[s["section_type"]] = per_type.get(s["section_type"], 0) + 1
                for g in s["cues"]:
                    if g in fs["cue"]:
                        fs["cue"][g] += 1
                if s["entities"]:
                    fs["sentences_with_entities"] += 1
                for e in s["entities"]:
                    lbl = e["label"] or "UNLABELED"
                    fs["entity_labels"][lbl] = fs["entity_labels"].get(lbl, 0) + 1
            for t, c in per_type.items():
                fs["type_sentences"][t] = fs["type_sentences"].get(t, 0) + c
                fs["type_files"][t] = fs["type_files"].get(t, 0) + 1

        print(f"    chunk done in {dt:.1f}s "
              f"({dt / max(len(chunk), 1):.1f}s/file)", flush=True)

    generated_at = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(SUMMARY_PATH, "w", encoding="utf-8") as fh:
        fh.write(render_html(agg, generated_at, scorer.model_id, scorer.device))

    jats, tei = agg["by_format"]["jats"], agg["by_format"]["tei"]
    total_sent = jats["total_sentences"] + tei["total_sentences"]
    with_sent = jats["files_with_sentences"] + tei["files_with_sentences"]

    total_ent = jats["total_entities"] + tei["total_entities"]

    def _report(name, fs):
        c = fs["cue"]
        print(f"[{name}] files={fs['files']:,}  with_sent={fs['files_with_sentences']:,}"
              f"  zero={fs['zero_files']:,}  cand={fs['total_candidates']:,}"
              f"  sentences={fs['total_sentences']:,}  entities={fs['total_entities']:,}"
              f"  (finding={c['finding']:,}, change={c['change']:,},"
              f" relation={c['relation']:,}, stat={c['stat']:,})")

    print(f"input files         : {agg['total_files']:,}")
    _report("JATS", jats)
    _report("TEI ", tei)
    if agg["by_format"]["unknown"]["files"]:
        _report("UNK ", agg["by_format"]["unknown"])
    print(f"files with sentences: {with_sent:,} ({pct(with_sent, agg['total_files'])})")
    print(f"files with 0 sent.  : {len(agg['zero_files']):,}")
    print(f"total result sent.  : {total_sent:,}")
    print(f"total NER entities  : {total_ent:,}")
    print(f"parse errors        : {agg['parse_errors']:,}")
    print(f"per-file JSON dir   : {OUT_DIR}")
    print(f"summary HTML        : {SUMMARY_PATH}")

    if OUTPUT_ZIP:                                      # bundle outputs for one-file download
        import zipfile
        os.makedirs(os.path.dirname(os.path.abspath(OUTPUT_ZIP)) or ".", exist_ok=True)
        nz = 0
        with zipfile.ZipFile(OUTPUT_ZIP, "w", zipfile.ZIP_DEFLATED) as zf:
            for fp in sorted(glob.glob(os.path.join(OUT_DIR, "*.json"))):
                zf.write(fp, os.path.join("sentences", os.path.basename(fp)))
                nz += 1
            if os.path.isfile(SUMMARY_PATH):
                zf.write(SUMMARY_PATH, os.path.join("summaries", os.path.basename(SUMMARY_PATH)))
                nz += 1
        print(f"output archive      : {OUTPUT_ZIP} ({nz:,} file(s))")


if __name__ == "__main__":
    main()
