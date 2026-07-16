#!/usr/bin/env python3
"""subtract.py -- move directory_1 entries that also exist in directory_2 aside.

Compares the *names* of the entries directly inside two directories (read from
STDIN) and moves every entry of ``directory_1`` whose name also appears in
``directory_2`` into ``gpu_bundle/removed``. In set terms it takes
``directory_1 - directory_2`` in place: what stays in ``directory_1`` is the
part unique to it; the overlap is relocated (not deleted) so nothing is lost.

The comparison is by entry name only (basename, case-sensitive, extension
included) over the top level of each directory -- it does not recurse and does
not compare file contents. Both files and sub-directories count as entries.

Usage
-----
    python subtract.py
        -> prompts for the two directories, one per line:
             directory_1 (entries here may be moved out):
             directory_2 (the set to subtract):
    printf 'high_impact_xmls\ngrobid_xmls\n' | python subtract.py   # piped STDIN

Outputs
-------
* moves the overlapping entries into ``gpu_bundle/removed`` (created on demand),
* writes ``summaries/subtract_optional.html`` -- a summary of directory_1,
  directory_2 and exactly what was moved.

``gpu_bundle/removed`` and ``summaries/`` are resolved relative to this script's
directory (so both are stable no matter where you launch from). If a name
already sits in the destination from an earlier run, the incoming entry is given
a ``_1``/``_2``/... suffix rather than overwriting it.
"""

import os
import sys
import html
import shutil
import datetime

# This script's directory -- outputs are resolved against it so the pipeline can
# be launched from any working directory (matches the rest of the bundle, e.g.
# step_1_orchestrator.py).
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REMOVED_DIR = os.path.join(BASE_DIR, "gpu_bundle", "removed")
SUMMARY_DIR = os.path.join(BASE_DIR, "summaries")
SUMMARY_HTML = os.path.join(SUMMARY_DIR, "subtract_optional.html")


def prompt_dir(label):
    """Read one directory path from STDIN, prompting with ``label``.

    The prompt is flushed *before* blocking on STDIN so it shows even when STDIN
    is not a TTY (e.g. PyCharm's run console). A blank line or EOF exits with an
    error rather than silently proceeding.
    """
    sys.stdout.write(label)
    sys.stdout.flush()
    line = sys.stdin.readline()          # "" on EOF (closed/empty STDIN)
    path = line.strip().strip('"').strip("'")
    if not path:
        sys.exit("error: no directory provided on STDIN")
    if not os.path.isdir(path):
        sys.exit("error: not a directory: %s" % path)
    return path


def unique_destination(name):
    """Return a path in REMOVED_DIR for ``name`` that does not clobber an existing
    entry: ``name``, else ``name_1``, ``name_2``, ... (suffix before the extension).
    """
    dest = os.path.join(REMOVED_DIR, name)
    if not os.path.exists(dest):
        return dest
    stem, ext = os.path.splitext(name)
    i = 1
    while True:
        candidate = os.path.join(REMOVED_DIR, "%s_%d%s" % (stem, i, ext))
        if not os.path.exists(candidate):
            return candidate
        i += 1


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
"""


def write_summary(dir1, dir2, entries1, entries2, moved):
    """Write summaries/subtract_optional.html describing the comparison + the move.

    ``moved`` is the list of ``(name, src, dest)`` tuples actually relocated (with
    the real destination, so collision-renamed ``name_1`` entries show correctly).
    """
    e = html.escape
    kept = len(entries1) - len(moved)
    only2 = len(entries2) - len(moved)     # names in dir2 that had no match in dir1
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    H = ["<!doctype html><html><head><meta charset='utf-8'>",
         "<title>subtract.py summary</title><style>%s</style></head><body>" % CSS]
    H.append("<h1>subtract.py &mdash; directory_1 &minus; directory_2</h1>")
    H.append("<p class='meta'>Generated %s. Entries of <code>directory_1</code> whose "
             "<em>name</em> also appears in <code>directory_2</code> were moved to "
             "<code>%s</code>.</p>" % (e(now), e(os.path.relpath(REMOVED_DIR, BASE_DIR))))

    H.append("<div class='stat-grid'>")
    for k, v in (("directory_1 entries", len(entries1)),
                 ("directory_2 entries", len(entries2)),
                 ("moved (in both)", len(moved)),
                 ("kept in directory_1", kept)):
        H.append("<div class='stat'><div class='v'>%d</div><div class='k'>%s</div></div>"
                 % (v, e(k)))
    H.append("</div>")

    H.append("<h2>Inputs</h2>")
    H.append("<table><tr><th>role</th><th>path</th><th class='num'>entries</th></tr>")
    H.append("<tr><td>directory_1 <span class='dim'>(source)</span></td><td><code>%s</code></td>"
             "<td class='num'>%d</td></tr>" % (e(dir1), len(entries1)))
    H.append("<tr><td>directory_2 <span class='dim'>(subtracted)</span></td><td><code>%s</code></td>"
             "<td class='num'>%d</td></tr>" % (e(dir2), len(entries2)))
    H.append("</table>")
    H.append("<p class='dim'>%d name(s) unique to directory_1 stayed put; %d name(s) unique to "
             "directory_2 were irrelevant to the move.</p>" % (kept, only2))

    H.append("<h2>Moved to <code>%s</code> (%d)</h2>"
             % (e(os.path.relpath(REMOVED_DIR, BASE_DIR)), len(moved)))
    if moved:
        H.append("<table><tr><th>#</th><th>name</th><th>from (directory_1)</th>"
                 "<th>to (destination)</th></tr>")
        for i, (name, src, dest) in enumerate(moved, 1):
            H.append("<tr><td class='num'>%d</td><td><code>%s</code></td><td><code>%s</code></td>"
                     "<td><code>%s</code></td></tr>"
                     % (i, e(name), e(src), e(os.path.relpath(dest, BASE_DIR))))
        H.append("</table>")
    else:
        H.append("<div class='key'>Nothing was moved: no <code>directory_1</code> entry name "
                 "is also present in <code>directory_2</code>.</div>")

    H.append("</body></html>")

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    with open(SUMMARY_HTML, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    print("[html] wrote %s" % os.path.relpath(SUMMARY_HTML, BASE_DIR))


def main():
    dir1 = prompt_dir("directory_1 (entries here may be moved out): ")
    dir2 = prompt_dir("directory_2 (the set to subtract): ")

    # Names present in each directory; directory_2's set is what we subtract.
    entries1 = os.listdir(dir1)
    entries2 = os.listdir(dir2)
    names2 = set(entries2)
    overlap = sorted(name for name in entries1 if name in names2)

    moved = []
    if overlap:
        os.makedirs(REMOVED_DIR, exist_ok=True)
        for name in overlap:
            src = os.path.join(dir1, name)
            dest = unique_destination(name)
            shutil.move(src, dest)
            print("moved %s -> %s" % (src, os.path.relpath(dest, BASE_DIR)))
            moved.append((name, src, dest))

    write_summary(dir1, dir2, entries1, entries2, moved)

    if moved:
        print("\ndone: moved %d entr%s from %s into %s"
              % (len(moved), "y" if len(moved) == 1 else "ies", dir1,
                 os.path.relpath(REMOVED_DIR, BASE_DIR)))
    else:
        print("nothing to move: no directory_1 entry name is also in directory_2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
