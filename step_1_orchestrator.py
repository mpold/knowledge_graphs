#!/usr/bin/env python3
"""Step 1 orchestrator -- runs the publications pipeline in order.

Executes the scripts of the Step 1 publications pipeline in the exact
order given by section "1. Order of execution (data flow)" of
``step_1_publications.html``. Each stage consumes the files the previous one
wrote; the whole thing turns one PubMed query into a de-duplicated, full-text
corpus ready for named-entity recognition, then deletes the intermediate
directories that produced it (step 8).

The scripts are *not* run automatically by the docs because several reach the
network (NCBI E-utilities, OpenAlex, CrossRef) and step 5 needs Docker + a
running GROBID container. This orchestrator just chains them; the caller is
responsible for the environment (NCBI_API_KEY, GROBID on :8070, etc.).

Usage
-----
    # query on the command line
    python step_1_orchestrator.py "your pubmed query"

    # query piped in (fed to step 1's STDIN)
    echo "your pubmed query" | python step_1_orchestrator.py

    # query typed at the prompt (if neither of the above is given)
    python step_1_orchestrator.py

    # both inputs piped in (query first, percentile second)
    printf 'your pubmed query\n0.01\n' | python step_1_orchestrator.py

    # re-seed the archive contribution after a manual step 6 (and re-stamp its years:
    # 6b replaces stamped files with the unstamped archive copies)
    python step_1_orchestrator.py --start 6b --stop 6c

    # keep the intermediate directories (stop before the clean-up step)
    python step_1_orchestrator.py --stop 7 "your pubmed query"

Options
-------
    --start S     start from step S (1, 1b, 2 .. 6, 6b, 6c, 7, 8) instead of step 1
    --stop  S     stop after step S
    --only  S     run only step S
    --archive D   local XML archive for steps 1b/6b (default: ../../xmls)
    --no-archive  skip steps 1b/6b and download the whole selection
    --list        print the pipeline order and exit
    --dry-run     print what would run without executing anything

Step 1 (``pubmed_query.py``) reads a query from STDIN and step 2
(``high_impact_xml.py``) takes a publication-impact percentile (a decimal
between 0 and 1), prompted on its own line right after the query; the
remaining steps take no argument and discover their inputs from the previous
stage's output directory. Any non-zero exit code aborts the pipeline.

When step 2 is in the selected range this orchestrator prompts for the
percentile on a dedicated line with::

    Publication impact percentile (decimal between 0 and 1):

and hands the value to ``high_impact_xml.py`` on that script's STDIN, *and* as
its ``PERCENTILE`` env var (env wins over STDIN there, so the validated value
is what takes effect either way). A blank line is still fed as an empty STDIN
line -- so the child never blocks re-prompting -- and it falls back to its
built-in default of 0.90.

The archive steps (1b / 6b)
---------------------------
``from_archive.py`` appears **twice**. As step 1b it matches the query result
against a local XML archive and serves every publication already held there,
writing ``pmids/from_archive_pmcids.txt`` so step 2 does not re-download them.
Step 6 (``named_entity_xml.py``) clean-rebuilds ``gpu_bundle/experimental_ner/``
and would discard that contribution, so step 6b re-runs the same script to
restore it from ``archive_xmls/``.

Both are skipped automatically -- with a printed notice -- when the archive
directory does not exist, so a checkout with no local archive still runs the
plain eight-stage pipeline. Naming an archive explicitly (``--archive`` or the
``ARCHIVE_DIR`` env var) makes a missing directory a hard error instead, since
that is a typo rather than an absence. ``--no-archive`` skips them outright.

The publication-year stamp (6c)
-------------------------------
``pub_year_xml.py`` writes a ``<?pub-year YYYY?>`` processing instruction into
the prolog of every file in ``gpu_bundle/experimental_ner/``, from the year column
of ``pmids/pmid_pmc_ids.tsv`` (falling back to the date inside the XML). It runs
after 6/6b because both rewrite the corpus, and it fails the pipeline if any file
is left undated -- so stage 2 receives a corpus that carries its own dates and
never needs the PubMed table or NCBI to place a paper in time.

The clean-up step (8)
---------------------
``clean_up.py`` runs last and deletes the five intermediate directories --
``archive_xmls/``, ``grobid_xmls/``, ``high_impact_xmls/``,
``named_entity_xmls/``, ``ncbi_pdfs_grobid/`` -- once the corpus is assembled.
It leaves ``gpu_bundle/experimental_ner/``, ``pmids/`` and ``summaries/`` alone.
Stop before it with ``--stop 7`` when you want to keep the intermediates (to
re-run a middle stage by hand, say); ``DRY_RUN=1`` makes step 8 report what it
would delete without deleting anything.
"""

import os
import sys
import subprocess

# Scripts in execution order (matches step_1_publications.html section 1).
#
# Steps are keyed by LABEL, not by position: from_archive.py runs twice, as "1b"
# and "6b", so keeping the core stages numbered 1-8 means --start 4 still means
# what it always did and the labels match the documentation.
PIPELINE = [
    ("1",  "pubmed_query.py",          "PubMed query (STDIN) -> pmids/pmid_pmc_ids.tsv"),
    ("1b", "from_archive.py",          "query result vs local archive -> experimental_ner/ + skip-list"),
    ("2",  "high_impact_xml.py",       "pmid_pmc_ids.tsv (minus skip-list) -> high_impact_xmls/"),
    ("3",  "xml_structure.py",         "high_impact_xmls/ -> summaries/xml_structure.html"),
    ("4",  "ncbi_pdf.py",              "no-<body> XMLs -> ncbi_pdfs_grobid/PMC*.pdf"),
    ("5",  "grobid_xml.py",            "PDFs (Docker+GROBID) -> grobid_xmls/PMC*.grobid.tei.xml"),
    ("6",  "named_entity_xml.py",      "grobid+high_impact -> gpu_bundle/experimental_ner/ (REBUILD)"),
    ("6b", "from_archive.py",          "archive_xmls/ -> experimental_ner/ (re-seed after the rebuild)"),
    ("6c", "pub_year_xml.py",          "pmid_pmc_ids.tsv -> <?pub-year?> stamp in every experimental_ner/ XML"),
    ("7",  "pre_ner_xml_structure.py", "experimental_ner/ -> summaries/pre_ner_xml_structure.html"),
    ("8",  "clean_up.py",              "REMOVES the five intermediate XML/PDF dirs (runs last)"),
]

# Steps served by from_archive.py -- skipped together when there is no archive.
ARCHIVE_STEPS = {"1b", "6b"}

# label -> position in PIPELINE, for resolving --start/--stop/--only.
STEP_INDEX = {key: i for i, (key, _script, _desc) in enumerate(PIPELINE)}

# Directory this orchestrator lives in -- all paths are script-relative so the
# pipeline can be launched from any working directory.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Default archive for steps 1b/6b -- the same default from_archive.py itself uses.
# Two levels up: this script lives in <project>/relationship_graphs/lung_small/ and
# the shared archive is <project>/xmls/, a sibling of relationship_graphs/.
DEFAULT_ARCHIVE = os.path.join(BASE_DIR, os.pardir, os.pardir, "xmls")


def parse_step(value):
    """Resolve a step label ("1", "1b", "6b", "7") to its position in PIPELINE.

    Labels rather than bare indices, so that inserting from_archive.py as 1b/6b
    leaves the seven core stages numbered exactly as they always were (and as
    step_1_publications.html documents them).
    """
    key = str(value).strip().lower()
    if key not in STEP_INDEX:
        sys.exit("error: unknown step %r (valid: %s)"
                 % (value, ", ".join(k for k, _s, _d in PIPELINE)))
    return STEP_INDEX[key]


def print_list():
    print("Step 1 publications pipeline -- order of execution:")
    for key, script, desc in PIPELINE:
        print("  %-3s %-26s %s" % (key + ".", script, desc))
    print("\n  1b/6b are the same script (from_archive.py): 1b serves the query result")
    print("  from a local XML archive, 6b restores it after step 6 rebuilds the corpus.")
    print("  6c stamps every corpus XML with <?pub-year YYYY?> for stage 2.")
    print("  8 (clean_up.py) deletes archive_xmls/, grobid_xmls/, high_impact_xmls/,")
    print("  named_entity_xmls/ and ncbi_pdfs_grobid/ -- the corpus in")
    print("  gpu_bundle/experimental_ner/, pmids/ and summaries/ are kept. --stop 7 skips it.")


def get_query(cli_query):
    """Resolve the PubMed query for step 1.

    Priority: explicit CLI argument, else a single line read from STDIN (typed
    interactively or piped in). ``pubmed_query.py`` itself reads one line from
    STDIN via input(), so we hand it the query through the child's stdin.
    """
    if cli_query is not None:
        return cli_query
    # Always emit the prompt (flushed) *before* blocking on STDIN so it shows
    # even when STDIN is not a TTY -- e.g. PyCharm's run console, where
    # sys.stdin.isatty() is False, would otherwise wait with no visible prompt.
    sys.stdout.write("enter pubmed query: ")
    sys.stdout.flush()
    line = sys.stdin.readline()      # "" on EOF (closed/empty STDIN)
    return line.strip()


def get_percentile():
    """Read + validate the publication-impact percentile for step 2 from STDIN.

    Prompted on its own dedicated line *after* the query is entered; the value is
    handed to ``high_impact_xml.py`` on its STDIN and as its ``PERCENTILE`` env
    var (env wins there, so the validated value takes effect either way). A blank
    line (or EOF) is returned as "" and fed on as an empty STDIN line, so that
    script falls back to its 0.90 default without re-prompting.

    A non-blank entry must be a decimal in [0, 1]. The two common mis-entries --
    a whole-number percentile ("90" instead of "0.90") or a stray character --
    otherwise sail through as a valid ``PERCENTILE`` env value and make
    ``high_impact_xml.py`` silently select an unexpected subset, which reads as
    "my percentile was ignored". We reject them *here*, before the slow step 1
    runs, with a message that names the fix.
    """
    sys.stdout.write("Publication impact percentile (decimal between 0 and 1): ")
    sys.stdout.flush()
    line = sys.stdin.readline()      # "" on EOF (closed/empty STDIN)
    pctl = line.strip()
    if not pctl:
        return ""
    try:
        value = float(pctl)
    except ValueError:
        sys.exit("error: percentile must be a decimal between 0 and 1, got %r" % pctl)
    if not 0.0 <= value <= 1.0:
        sys.exit("error: percentile must be a fraction between 0 and 1 "
                 "(e.g. 0.90 for the 90th percentile), got %r" % pctl)
    return pctl


def resolve_archive(cli_archive, no_archive, wanted):
    """Decide the archive directory for steps 1b/6b. Returns (dir_or_None, skip).

    An archive named *explicitly* -- ``--archive`` or the ``ARCHIVE_DIR`` env var
    -- that does not exist is a hard error: the user meant a specific directory
    and got the path wrong, and silently downloading ~3 req/s instead is a poor
    way to find that out. An absent *default* archive is not an error at all: a
    fresh checkout simply has no local corpus, so 1b/6b are skipped with a notice
    and the plain eight-stage pipeline runs.
    """
    if not wanted:
        return None, True
    if no_archive:
        print("[orchestrator] --no-archive: skipping steps 1b/6b; the whole "
              "selection will be downloaded", flush=True)
        return None, True

    env_archive = os.environ.get("ARCHIVE_DIR", "").strip()
    explicit = cli_archive or env_archive
    path = os.path.abspath(explicit or DEFAULT_ARCHIVE)

    if os.path.isdir(path):
        print("[orchestrator] archive for steps 1b/6b: %s" % path, flush=True)
        return path, False
    if explicit:
        sys.exit("error: archive not found: %s\n"
                 "       (named via %s) -- fix the path or pass --no-archive"
                 % (path, "--archive" if cli_archive else "ARCHIVE_DIR"))
    print("[orchestrator] no local archive at %s -- skipping steps 1b/6b; every "
          "selected article will be downloaded" % path, flush=True)
    return None, True


def run_stage(step_key, script, stdin_text, dry_run, env_extra=None):
    """Run a single pipeline stage, returning its exit code.

    ``stdin_text`` is the text fed to the child's STDIN (with a trailing newline),
    or ``None`` for stages that take no STDIN input. ``env_extra`` is an optional
    dict of environment variables layered over the inherited environment for the
    child (e.g. ``PERCENTILE`` for ``high_impact_xml.py``, ``ARCHIVE_DIR`` for
    ``from_archive.py``).
    """
    path = os.path.join(BASE_DIR, script)
    label = "[step %s of %d] %s" % (step_key, len(PIPELINE), script)

    if not os.path.isfile(path):
        print("%s -- MISSING (%s)" % (label, path), file=sys.stderr)
        return 127

    # ``-u`` forces the child's stdout/stderr to be unbuffered. Without it, when
    # this orchestrator's stdout is a pipe rather than a TTY (e.g. PyCharm's run
    # console), Python block-buffers the child's output: its prompt line
    # ("enter command line argument: ") shows, then every progress line is held
    # in the 4 KB buffer until it fills or the process exits, making the stage
    # look stalled for minutes even though it is running. Unbuffered output
    # streams live so the console tracks the child's real progress.
    cmd = [sys.executable, "-u", path]
    feeds_stdin = stdin_text is not None
    child_env = {**os.environ, **env_extra} if env_extra else None

    print("=" * 70)
    tags = (" <- STDIN" if feeds_stdin else "") + (
        " <- " + ", ".join(env_extra) if env_extra else "")
    print(label + tags)
    print("=" * 70, flush=True)

    if dry_run:
        print("  (dry-run) would run: %s" % " ".join(cmd))
        return 0

    if script == "pubmed_query.py" and not (stdin_text or "").strip():
        print("%s -- no query provided for step 1" % label, file=sys.stderr)
        return 2

    if feeds_stdin:
        proc = subprocess.run(cmd, cwd=BASE_DIR, input=stdin_text + "\n",
                              text=True, env=child_env)
    else:
        proc = subprocess.run(cmd, cwd=BASE_DIR, env=child_env)
    return proc.returncode


def main(argv):
    args = argv[1:]

    start, stop, cli_query, dry_run = 0, len(PIPELINE) - 1, None, False
    archive, no_archive = None, False
    i = 0
    while i < len(args):
        a = args[i]
        if a in ("-h", "--help"):
            print(__doc__)
            return 0
        elif a == "--list":
            print_list()
            return 0
        elif a == "--dry-run":
            dry_run = True
        elif a == "--start":
            i += 1
            start = parse_step(args[i])
        elif a == "--stop":
            i += 1
            stop = parse_step(args[i])
        elif a == "--only":
            i += 1
            start = stop = parse_step(args[i])
        elif a == "--archive":
            i += 1
            archive = args[i]
        elif a == "--no-archive":
            no_archive = True
        elif a.startswith("--"):
            sys.exit("error: unknown option %r (try --help)" % a)
        else:
            # First bare argument is the PubMed query for step 1.
            if cli_query is not None:
                sys.exit("error: unexpected extra argument %r" % a)
            cli_query = a
        i += 1

    if start > stop:
        sys.exit("error: --start (%s) is after --stop (%s)"
                 % (PIPELINE[start][0], PIPELINE[stop][0]))

    selected = [(key, script) for key, script, _desc in PIPELINE[start:stop + 1]]
    archive_dir, skip_archive = resolve_archive(archive, no_archive,
                                                any(k in ARCHIVE_STEPS for k, _s in selected))
    if skip_archive:
        selected = [(k, s) for k, s in selected if k not in ARCHIVE_STEPS]

    # Collect each stage's input in the order the prompts are consumed: the query
    # (step 1, fed on STDIN) first, then the impact percentile (step 2).
    stdin_by_step = {}
    env_by_step = {}
    keys = {k for k, _s in selected}
    if "1" in keys:
        stdin_by_step["1"] = get_query(cli_query)
    if "2" in keys:
        pctl = get_percentile()
        # Fed on STDIN (the documented input) *and* as PERCENTILE, which wins in
        # high_impact_xml.py -- so the value validated here is what takes effect
        # either way. A blank entry is still fed as an empty STDIN line, so the
        # child falls through to its 0.90 default instead of re-prompting on a
        # TTY the orchestrator has already prompted on.
        stdin_by_step["2"] = pctl
        if pctl:
            env_by_step["2"] = {"PERCENTILE": pctl}
            # Echoing it makes "did my percentile get through?" answerable
            # without reading the child's [select] line.
            print("[orchestrator] step 2 (high_impact_xml.py) will run with "
                  "PERCENTILE=%s" % pctl, flush=True)
        else:
            print("[orchestrator] no percentile entered; step 2 (high_impact_xml.py) "
                  "will use its 0.90 default", flush=True)
    if archive_dir:
        for key in ARCHIVE_STEPS:
            env_by_step[key] = {"ARCHIVE_DIR": archive_dir}

    for key, script in selected:
        code = run_stage(key, script, stdin_by_step.get(key), dry_run,
                         env_by_step.get(key))
        if code != 0:
            print(
                "\nPIPELINE ABORTED at step %s (%s): exit code %d"
                % (key, script, code),
                file=sys.stderr,
            )
            return code

    print("\nPIPELINE COMPLETE: steps %s-%s finished successfully (%d stage%s run)."
          % (PIPELINE[start][0], PIPELINE[stop][0], len(selected),
             "" if len(selected) == 1 else "s"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
