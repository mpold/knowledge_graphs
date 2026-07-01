#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""zip_work.py -- bundle the whole run directory into a single kaggle_working.zip.

Final step of the gpu_bundle pipeline (gpu.py): once every earlier step has written
its outputs into the working dir, this packs that dir into one downloadable archive.

The working dir is /kaggle/working on Kaggle, or -- when launched by gpu.py, which
runs each step with cwd=work_root -- the staged work-root; standalone it falls back
to the current directory. The archive is built in a temp dir (so the growing zip
can't include itself) and then moved into the working dir. Inside a notebook it also
shows a clickable download link.
"""

import shutil
import tempfile
from pathlib import Path

# Dir to archive: /kaggle/working on Kaggle, else the current working dir
# (gpu.py runs this step with cwd = the staged work-root).
work_dir = Path("/kaggle/working") if Path("/kaggle/working").is_dir() else Path.cwd()

ARCHIVE = "kaggle_working"
target = work_dir / f"{ARCHIVE}.zip"

# Drop any archive from a previous run so it isn't bundled into the new one.
if target.exists():
    target.unlink()

# Build in a temp dir so the archive doesn't zip itself, then move it into the working dir.
built = shutil.make_archive(str(Path(tempfile.gettempdir()) / ARCHIVE), "zip", str(work_dir))
shutil.move(built, str(target))
print(f"bundled {work_dir} -> {target} ({target.stat().st_size:,} bytes)")

# Clickable download link when running inside a notebook (no-op as a plain script).
try:
    from IPython.display import display, FileLink
    display(FileLink(target.name if work_dir == Path.cwd() else str(target)))
except Exception:
    pass