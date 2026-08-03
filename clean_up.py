#!/usr/bin/env python3
"""clean_up.py -- remove the step 1 intermediate directories once the corpus is built.

Final stage (step 8) of the Step 1 publications pipeline. By the time step 7
(``pre_ner_xml_structure.py``) has run, the deliverable is
``gpu_bundle/experimental_ner/`` -- the de-duplicated NER corpus handed to
stage 2 -- and the working directories that produced it are just bulk on disk:
several GB of JATS XML, PDFs and TEI that no later stage reads.

Removed (directory *and* contents), each resolved relative to this script:

    archive_xmls/        staging copy of the archive hit set (step 1b/6b)
    grobid_xmls/         TEI produced from the PDFs (step 5)
    high_impact_xmls/    downloaded JATS full text (step 2)
    named_entity_xmls/   de-duplicated one-file-per-publication set (step 6)
    ncbi_pdfs_grobid/    PDFs for the no-<body> records (step 4)

Never touched: ``gpu_bundle/experimental_ner/`` (the corpus), ``pmids/`` (the
query result + resumable caches), ``summaries/`` (the HTML reports) and every
Python script. The removal list is a fixed literal -- nothing is derived from
arguments or the environment -- and each target must be a real directory
directly inside this script's own directory, so a stray symlink or an
out-of-tree path is refused rather than followed.

Usage
-----
    python clean_up.py              # remove them
    DRY_RUN=1 python clean_up.py    # report what would be removed, delete nothing
    python clean_up.py --dry-run    # same
    KEEP=archive_xmls python clean_up.py    # keep one (or a comma-separated list)

A directory that is already absent is not an error: it is reported as such and
the run still succeeds, so re-running the pipeline (or this script) is safe.

Outputs
-------
* writes ``summaries/clean_up.html`` -- what was removed, what was kept, how
  many files and bytes each directory held.

Re-running earlier steps afterwards
-----------------------------------
Step 8 deletes the inputs of steps 3-6b, so resuming mid-pipeline (e.g.
``--only 5``) needs those steps re-run from step 2. Step 6b is the exception
that matters in practice: ``from_archive.py`` re-copies the hit set out of the
real archive (``ARCHIVE_DIR``), so ``--only 6b`` still works after a clean-up --
it just pays the file copy again instead of reading the ``archive_xmls/``
staging copy. Keep it with ``KEEP=archive_xmls`` if you plan to re-seed often.
"""

import os
import sys
import html
import shutil
import datetime

# This script's directory -- every target is resolved against it (as in the rest
# of the bundle), so the pipeline can be launched from any working directory.
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# The directories this script removes. A fixed literal, in the order documented
# above: nothing here is built from CLI arguments or env vars, so the script can
# only ever delete these five names, and only inside BASE_DIR.
TARGETS = [
    "archive_xmls",
    "grobid_xmls",
    "high_impact_xmls",
    "named_entity_xmls",
    "ncbi_pdfs_grobid",
]

SUMMARY_DIR = os.path.join(BASE_DIR, "summaries")
SUMMARY_HTML = os.path.join(SUMMARY_DIR, "clean_up.html")


def measure(path):
    """Return (file_count, total_bytes) for the tree at ``path``.

    Unreadable entries are skipped rather than raising: this is a report of what
    is about to go, not an audit, and a locked file must not abort the clean-up.
    """
    files = 0
    total = 0
    for _root, _dirs, names in os.walk(path):
        for name in names:
            files += 1
            try:
                total += os.path.getsize(os.path.join(_root, name))
            except OSError:
                pass
    return files, total


def human(num_bytes):
    """Format a byte count as B / KB / MB / GB."""
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return "%.0f %s" % (size, unit) if unit == "B" else "%.1f %s" % (size, unit)
        size /= 1024


def classify(name, keep):
    """Decide what to do with target ``name``: returns (status, path, files, bytes).

    Status is one of ``remove`` (a real directory, in tree, not kept),
    ``kept`` (named in ``KEEP``), ``absent`` (nothing there) or ``refused``
    (present but not a plain directory inside BASE_DIR -- e.g. a symlink or a
    file). ``refused`` is a guard, not an expected outcome: it means the name
    points somewhere this script will not follow.
    """
    path = os.path.join(BASE_DIR, name)

    if name in keep:
        files, size = measure(path) if os.path.isdir(path) else (0, 0)
        return "kept", path, files, size
    if not os.path.exists(path) and not os.path.islink(path):
        return "absent", path, 0, 0
    # islink first: os.path.isdir() follows symlinks, and rmtree() on a link to a
    # directory outside the tree would delete somebody else's files.
    if os.path.islink(path) or not os.path.isdir(path):
        return "refused", path, 0, 0
    if os.path.dirname(os.path.realpath(path)) != os.path.realpath(BASE_DIR):
        return "refused", path, 0, 0
    files, size = measure(path)
    return "remove", path, files, size


CSS = """
 body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        margin: 2rem auto; max-width: 1100px; color: #222; line-height: 1.45; padding: 0 1rem; }
 h1 { margin-bottom: .25rem; } h2 { margin-top: 2.25rem; border-bottom: 1px solid #ddd; padding-bottom: .25rem; }
 .meta { color: #555; margin-bottom: 1rem; }
 .stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: .75rem; margin: 1rem 0 1.5rem; }
 .stat { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: .7rem 1rem; }
 .stat .v { font-size: 1.35rem; font-weight: 600; } .stat .k { color: #555; font-size: .82rem; }
 table { border-collapse: collapse; width: 100%; margin: .5rem 0 1rem; font-size: .9rem; }
 th, td { border: 1px solid #e1e4e8; padding: .3rem .55rem; text-align: left; vertical-align: middle; }
 th { background: #f6f8fa; }
 td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
 code { background: #f6f8fa; padding: 1px 4px; border-radius: 3px; font-size: .88em; }
 .dim { color: #888; font-size: .85em; }
 .key  { background: #ddf4ff; border-left: 4px solid #0969da; padding: .6rem .9rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 .warn { background: #fff8c5; border-left: 4px solid #d4a72c; padding: .6rem .9rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
"""

STATUS_TEXT = {
    "removed": "removed",
    "would remove": "would remove <span class='dim'>(dry run)</span>",
    "kept": "kept <span class='dim'>(KEEP)</span>",
    "absent": "absent <span class='dim'>(nothing to do)</span>",
    "refused": "refused <span class='dim'>(not a plain directory in the bundle root)</span>",
    "failed": "FAILED",
}


def write_summary(rows, dry_run):
    """Write summaries/clean_up.html.

    ``rows`` is a list of ``(name, status, files, bytes, detail)`` tuples in
    TARGETS order, with the status as it actually resolved.
    """
    e = html.escape
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    gone = [r for r in rows if r[1] in ("removed", "would remove")]
    freed = sum(r[3] for r in gone)
    files = sum(r[2] for r in gone)

    H = ["<!doctype html><html><head><meta charset='utf-8'>",
         "<title>clean_up.py summary</title><style>%s</style></head><body>" % CSS]
    H.append("<h1>clean_up.py &mdash; step 1 intermediate directories</h1>")
    H.append("<p class='meta'>Generated %s. Final stage of the Step&nbsp;1 publications "
             "pipeline: the working directories that produced "
             "<code>gpu_bundle/experimental_ner/</code> are removed once the corpus "
             "is assembled.%s</p>"
             % (e(now), " <strong>Dry run &mdash; nothing was deleted.</strong>"
                if dry_run else ""))

    H.append("<div class='stat-grid'>")
    for k, v in (("directories %s" % ("that would go" if dry_run else "removed"), str(len(gone))),
                 ("files %s" % ("that would go" if dry_run else "deleted"), "{:,}".format(files)),
                 ("disk %s" % ("that would be freed" if dry_run else "freed"), human(freed)),
                 ("directories in the list", str(len(rows)))):
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>"
                 % (e(v), e(k)))
    H.append("</div>")

    H.append("<h2>Directories</h2>")
    H.append("<table><tr><th>directory</th><th>status</th><th class='num'>files</th>"
             "<th class='num'>size</th><th>note</th></tr>")
    for name, status, n_files, n_bytes, detail in rows:
        H.append("<tr><td><code>%s/</code></td><td>%s</td><td class='num'>%d</td>"
                 "<td class='num'>%s</td><td class='dim'>%s</td></tr>"
                 % (e(name), STATUS_TEXT.get(status, e(status)), n_files,
                    human(n_bytes), e(detail)))
    H.append("</table>")

    H.append("<div class='key'><strong>Not touched:</strong> "
             "<code>gpu_bundle/experimental_ner/</code> (the NER corpus handed to "
             "step&nbsp;2), <code>pmids/</code> (the query result and its resumable "
             "caches), <code>summaries/</code> (these reports) and every Python "
             "script. The removal list is a fixed literal of five names, each of "
             "which must be a real directory directly inside the bundle root.</div>")
    if any(r[1] == "failed" for r in rows):
        H.append("<div class='warn'>One or more directories could not be removed &mdash; "
                 "see the note column. A file held open by another process (an editor, "
                 "an indexer, a running stage) is the usual cause; re-run "
                 "<code>python clean_up.py</code> once it is released.</div>")
    H.append("<p class='dim'>Re-running earlier steps after a clean-up: the inputs of "
             "steps&nbsp;3&ndash;6b are gone, so resume from step&nbsp;2. Step&nbsp;6b is "
             "the exception &mdash; <code>from_archive.py</code> re-copies the hit set "
             "from the real archive (<code>ARCHIVE_DIR</code>), so "
             "<code>--only 6b</code> still works, it just pays the copy again. "
             "<code>KEEP=archive_xmls</code> preserves the staging copy.</p>")
    H.append("</body></html>")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(SUMMARY_HTML, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    print("[html] wrote %s" % os.path.relpath(SUMMARY_HTML, BASE_DIR))


def main(argv):
    args = argv[1:]
    dry_run = os.environ.get("DRY_RUN", "").strip() not in ("", "0")
    for a in args:
        if a in ("-h", "--help"):
            print(__doc__)
            return 0
        if a == "--dry-run":
            dry_run = True
        else:
            sys.exit("error: unknown argument %r (try --help)" % a)

    keep = {n.strip() for n in os.environ.get("KEEP", "").split(",") if n.strip()}
    unknown = keep - set(TARGETS)
    if unknown:
        sys.exit("error: KEEP names nothing this script removes: %s (valid: %s)"
                 % (", ".join(sorted(unknown)), ", ".join(TARGETS)))

    if dry_run:
        print("[clean_up] DRY RUN -- nothing will be deleted")

    rows = []
    failures = 0
    for name in TARGETS:
        status, path, n_files, n_bytes = classify(name, keep)
        detail = ""

        if status == "remove":
            if dry_run:
                status = "would remove"
            else:
                try:
                    shutil.rmtree(path)
                    status = "removed"
                except OSError as exc:
                    status = "failed"
                    detail = str(exc)
                    failures += 1
        elif status == "refused":
            detail = "not a plain directory inside %s" % BASE_DIR
        elif status == "kept":
            detail = "named in KEEP"

        rows.append((name, status, n_files, n_bytes, detail))
        print("%-18s %-13s %6d file(s)  %9s%s"
              % (name + "/", status, n_files, human(n_bytes),
                 "  -- " + detail if detail else ""))

    write_summary(rows, dry_run)

    gone = [r for r in rows if r[1] in ("removed", "would remove")]
    print("\n%s %d director%s, %d file(s), %s"
          % ("would remove" if dry_run else "removed", len(gone),
             "y" if len(gone) == 1 else "ies",
             sum(r[2] for r in gone), human(sum(r[3] for r in gone))))
    # A refusal is not a failure of the clean-up (nothing was supposed to be
    # deleted there) but it is anomalous -- one of the five names is present and
    # is not a plain directory in the bundle root -- so say so on stderr rather
    # than leaving it in the table only. The pipeline is not aborted for it.
    refused = [r[0] for r in rows if r[1] == "refused"]
    if refused:
        print("note: left alone, not a plain directory in the bundle root: %s"
              % ", ".join(refused), file=sys.stderr)
    if failures:
        print("%d director%s could not be removed (see above)"
              % (failures, "y" if failures == 1 else "ies"), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
