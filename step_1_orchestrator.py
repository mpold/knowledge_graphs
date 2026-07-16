#!/usr/bin/env python3
"""Step 1 orchestrator -- runs the publications pipeline in order.

Executes the seven scripts of the Step 1 publications pipeline in the exact
order given by section "1. Order of execution (data flow)" of
``step_1_publications.html``. Each stage consumes the files the previous one
wrote; the whole thing turns one PubMed query into a de-duplicated, full-text
corpus ready for named-entity recognition.

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

Options
-------
    --start N     start from step N (1-7) instead of step 1
    --stop  N     stop after step N (1-7)
    --only  N     run only step N
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

and passes the value to ``high_impact_xml.py`` via its ``PERCENTILE`` env var
(which is the only place that script reads the percentile from). A blank line
lets that script fall back to any inherited ``PERCENTILE`` env var, or its
built-in default of 0.90.
"""

import os
import sys
import subprocess

# Scripts in execution order (matches step_1_publications.html section 1).
PIPELINE = [
    "pubmed_query.py",          # 1. PubMed query (STDIN) -> pmids/pmid_pmc_ids.tsv
    "high_impact_xml.py",       # 2. pmid_pmc_ids.tsv     -> high_impact_xmls/PMC*.xml
    "xml_structure.py",         # 3. high_impact_xmls/    -> summaries/xml_structure.html
    "ncbi_pdf.py",              # 4. no-<body> XMLs       -> ncbi_pdfs_grobid/PMC*.pdf
    "grobid_xml.py",            # 5. PDFs (Docker+GROBID) -> grobid_xmls/PMC*.grobid.tei.xml
    "named_entity_xml.py",      # 6. grobid+high_impact   -> gpu_bundle/experimental_ner/
    "pre_ner_xml_structure.py", # 7. experimental_ner/    -> summaries/pre_ner_xml_structure.html
]

# Directory this orchestrator lives in -- all paths are script-relative so the
# pipeline can be launched from any working directory.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_step(value):
    """Parse a 1-based step number, validating it against the pipeline length."""
    try:
        n = int(value)
    except ValueError:
        sys.exit("error: step must be an integer 1-%d, got %r" % (len(PIPELINE), value))
    if not 1 <= n <= len(PIPELINE):
        sys.exit("error: step out of range: %d (valid: 1-%d)" % (n, len(PIPELINE)))
    return n


def print_list():
    print("Step 1 publications pipeline -- order of execution:")
    for i, script in enumerate(PIPELINE, 1):
        print("  %d. %s" % (i, script))


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
    """Read the publication-impact percentile for step 2 from STDIN.

    Prompted on its own dedicated line *after* the query is entered; the value is
    passed to ``high_impact_xml.py`` via its ``PERCENTILE`` env var (the only place
    that script reads the percentile). A blank line (or EOF) is returned as "" so
    that script falls back to its inherited PERCENTILE env / 0.90 default.
    """
    sys.stdout.write("Publication impact percentile (decimal between 0 and 1): ")
    sys.stdout.flush()
    line = sys.stdin.readline()      # "" on EOF (closed/empty STDIN)
    return line.strip()


def run_stage(step_no, script, stdin_text, dry_run, env_extra=None):
    """Run a single pipeline stage, returning its exit code.

    ``stdin_text`` is the text fed to the child's STDIN (with a trailing newline),
    or ``None`` for stages that take no STDIN input. ``env_extra`` is an optional
    dict of environment variables layered over the inherited environment for the
    child (e.g. ``PERCENTILE`` for ``high_impact_xml.py``).
    """
    path = os.path.join(BASE_DIR, script)
    label = "[%d/%d] %s" % (step_no, len(PIPELINE), script)

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

    start, stop, cli_query, dry_run = 1, len(PIPELINE), None, False
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
        elif a.startswith("--"):
            sys.exit("error: unknown option %r (try --help)" % a)
        else:
            # First bare argument is the PubMed query for step 1.
            if cli_query is not None:
                sys.exit("error: unexpected extra argument %r" % a)
            cli_query = a
        i += 1

    if start > stop:
        sys.exit("error: --start (%d) is after --stop (%d)" % (start, stop))

    # Collect each stage's input in the order the prompts are consumed: the query
    # (step 1, fed on STDIN) first, then the impact percentile (step 2, passed as
    # the PERCENTILE env var since high_impact_xml.py reads it only from there).
    stdin_by_script = {}
    env_by_script = {}
    if start == 1:
        stdin_by_script["pubmed_query.py"] = get_query(cli_query)
    if start <= 2 <= stop:
        pctl = get_percentile()
        if pctl:
            env_by_script["high_impact_xml.py"] = {"PERCENTILE": pctl}

    for step_no in range(start, stop + 1):
        script = PIPELINE[step_no - 1]
        code = run_stage(step_no, script, stdin_by_script.get(script), dry_run,
                         env_by_script.get(script))
        if code != 0:
            print(
                "\nPIPELINE ABORTED at step %d (%s): exit code %d"
                % (step_no, script, code),
                file=sys.stderr,
            )
            return code

    print("\nPIPELINE COMPLETE: steps %d-%d finished successfully." % (start, stop))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
