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

Only step 1 (``pubmed_query.py``) reads a query from STDIN; the remaining
steps take no argument and discover their inputs from the previous stage's
output directory. Any non-zero exit code aborts the pipeline.
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

    Priority: explicit CLI argument, then piped STDIN, then an interactive
    prompt. ``pubmed_query.py`` itself reads a single line from STDIN via
    input(), so we hand it the query through the child's stdin.
    """
    if cli_query is not None:
        return cli_query
    if not sys.stdin.isatty():
        # Something was piped in -- forward the first line to step 1.
        line = sys.stdin.readline()
        return line.strip()
    try:
        return input("enter pubmed query: ").strip()
    except EOFError:
        return ""


def run_stage(step_no, script, query, dry_run):
    """Run a single pipeline stage, returning its exit code."""
    path = os.path.join(BASE_DIR, script)
    label = "[%d/%d] %s" % (step_no, len(PIPELINE), script)

    if not os.path.isfile(path):
        print("%s -- MISSING (%s)" % (label, path), file=sys.stderr)
        return 127

    cmd = [sys.executable, path]
    feeds_stdin = script == "pubmed_query.py"

    print("=" * 70)
    print(label + (" <- query on STDIN" if feeds_stdin else ""))
    print("=" * 70, flush=True)

    if dry_run:
        print("  (dry-run) would run: %s" % " ".join(cmd))
        return 0

    if feeds_stdin:
        if not query:
            print("%s -- no query provided for step 1" % label, file=sys.stderr)
            return 2
        proc = subprocess.run(cmd, cwd=BASE_DIR, input=query + "\n", text=True)
    else:
        proc = subprocess.run(cmd, cwd=BASE_DIR)
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

    # Only fetch a query if step 1 is actually in the selected range.
    query = get_query(cli_query) if start == 1 else None

    for step_no in range(start, stop + 1):
        script = PIPELINE[step_no - 1]
        code = run_stage(step_no, script, query, dry_run)
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
