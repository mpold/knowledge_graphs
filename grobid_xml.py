#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grobid_xml.py
=============
Convert the PDFs in ``ncbi_pdfs_grobid/`` (the full-text PDFs that ncbi_pdf.py
exported for the high-impact articles whose PMC XML had no <body>) into
structured TEI XML using GROBID, then write an HTML summary of the conversion.

Modeled on ``bacs/batch_grobid.py``: it drives the same GROBID
``processFulltextDocument`` service, here by POSTing each PDF to the GROBID REST
API with ``requests``. The TEI files GROBID produces (``<PMCID>.grobid.tei.xml``)
become the structured full text that was missing from the NCBI export, ready for
the same section analysis as the rest of the corpus.

Input  : ncbi_pdfs_grobid/PMC*.pdf
Output : grobid_xmls/PMC*.grobid.tei.xml
         summaries/grobid_xml_summary.html

NOT executed here -- run it yourself.

# =========================================================================== #
# PREREQUISITES
#   1. Docker installed and the engine running (the script can launch Docker
#      Desktop for you and starts the GROBID container itself).
#   2. Python package: requests (used for both the container health-check and the
#      PDF -> TEI REST calls -- no separate grobid client package is needed).
#   3. (optional) For the Chrome console step:  pip install selenium
#
# The script AUTOMATES the manual bring-up that bacs/batch_grobid.py left to the
# operator (set GROBID_AUTOLAUNCH=0 to skip it and use an already-running server):
#   1. Launches Docker Desktop and waits for the engine to be ready.
#   2. Starts the GROBID container DETACHED with Docker (pulling the image first
#      if missing):
#         docker run -d --rm --name grobid_xml -p 8070:8070 grobid/grobid:0.8.2
#      and waits for http://localhost:8070/api/isalive to return "true".
#   3. (optional) Opens Chrome at http://localhost:8070/ for visual confirmation.
#   4. Runs the batch PDF->TEI conversion by POSTing each PDF to the GROBID REST
#      API (/api/processFulltextDocument) with requests.
# =========================================================================== #

STRATEGY
--------
0. Bring-up automation (GROBID_AUTOLAUNCH, on by default). Launch Docker Desktop,
   start the GROBID container detached with Docker (pulling the image on first
   use), wait for the server to answer, then optionally open the Chrome console.
   Any step that cannot complete is logged and the script carries on.

1. Prerequisite: a running GROBID server (Docker, port 8070 -- see the comment
   block above). The script pings ``/api/isalive`` first and refuses to start the
   conversion if the server is unreachable (but still (re)writes the summary).

2. Convert per file over the GROBID REST API (after bacs/batch_grobid.py). Each
   PDF is POSTed to ``/api/processFulltextDocument`` with ``requests``, N at a time
   via a thread pool, writing one ``<PMCID>.grobid.tei.xml`` per PDF into
   ``grobid_xmls/``.

3. Resume-friendly & fault-tolerant. PDFs that already have a non-empty TEI are
   skipped, so a re-run only fills the gaps (set GROBID_FORCE=1 to reconvert all).
   Transient failures (read timeout / HTTP 503 busy / 5xx / dropped connections)
   are retried with linear backoff; permanent ones (HTTP 400, "failed to open PDF
   file", empty result) are logged and skipped so one bad PDF never aborts the run.

4. Summarise. After the batch, the script scans input PDFs vs output TEIs to report
   the conversion rate, then does a light TEI structure check (title, abstract,
   body sections, references) so the GROBID output can be compared with the NCBI
   XML analysis in summaries/ncbi_xml.html.
"""

import os
import re
import sys
import glob
import time
import random
import shutil
import threading
import subprocess
import webbrowser
import html as _html
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #
BASE         = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR    = os.path.join(BASE, "ncbi_pdfs_grobid")        # source PDFs
OUTPUT_DIR   = os.path.join(BASE, "grobid_xmls")             # TEI XML target
SUMMARY_DIR  = os.path.join(BASE, "summaries")
SUMMARY_HTML = os.path.join(SUMMARY_DIR, "grobid_xml.html")

# NB: use 127.0.0.1, not "localhost". On Windows + Docker Desktop, "localhost"
# resolves to IPv6 ::1 first, and Docker's IPv6 port-forward frequently does NOT
# service the published port -- so every health check (/api/isalive) and POST
# hangs until it times out even though GROBID is up and serving on IPv4. That
# made the script wrongly conclude "the server never came up". 127.0.0.1 forces
# IPv4 and connects instantly. Override with GROBID_SERVER if you really need a
# different host.
GROBID_SERVER = os.environ.get("GROBID_SERVER", "http://127.0.0.1:8070")
GROBID_URL    = GROBID_SERVER.rstrip("/") + "/"                # browser console
SERVICE       = "processFulltextDocument"
N_THREADS     = int(os.environ.get("GROBID_THREADS", "1"))
TIMEOUT       = int(os.environ.get("GROBID_TIMEOUT", "120"))
FORCE         = os.environ.get("GROBID_FORCE", "").strip() in ("1", "true", "yes")
TEI_SUFFIX    = ".grobid.tei.xml"

# --- Per-file fault tolerance ---------------------------------------------- #
# Transient HTTP statuses worth another attempt; everything else (e.g. 400 bad
# request / "failed to open PDF", 204 empty) is treated as permanent and skipped.
MAX_ATTEMPTS  = max(1, int(os.environ.get("GROBID_ATTEMPTS", "3")))   # tries per PDF
RETRY_SLEEP   = int(os.environ.get("GROBID_RETRY_SLEEP", "5"))        # s, linear backoff base
RETRY_STATUS  = (408, 429, 500, 502, 503, 504)

# --- Server saturation ("max connections exceeded") ------------------------ #
# GROBID keeps a fixed pool of parsing engines (its `concurrency` setting). When
# more requests are in flight than the pool holds -- including the warm-up window
# right after /api/isalive first answers, while models are still loading -- it
# rejects the surplus with HTTP 503 *before parsing*, body "the maximum number of
# concurrent connection has been reached". That request did no work, so it is
# ALWAYS safe to retry; we back off (with jitter, so the N threads desynchronise)
# and do NOT count it against MAX_ATTEMPTS, rather than logging it as a failure.
BUSY_STATUS    = (429, 503)
BUSY_TOKENS    = ("concurrent", "max connection", "saturat", "too many request",
                  "reached", "try again")
MAX_BUSY_WAITS = max(1, int(os.environ.get("GROBID_BUSY_RETRIES", "40")))   # per PDF
BUSY_WAIT_CAP  = int(os.environ.get("GROBID_BUSY_WAIT_CAP", "30"))          # s, backoff ceiling

# --- Server DOWN ("connection refused" / container crashed) ---------------- #
# A NewConnectionError / "max retries exceeded" / connection-reset means nothing
# is listening on port 8070: the GROBID container itself died -- typically an
# out-of-memory kill after a handful of memory-heavy fulltext parses -- rather
# than a transient HTTP error. The POST did no work, so it is ALWAYS safe to
# retry, but only AFTER the server is back. When a worker sees this it calls
# ensure_server_up() (which relaunches the container, under a lock so only one
# thread does it) and keeps retrying without spending the per-file attempt
# budget, so the batch PAUSES and RESUMES instead of mass-skipping every
# remaining PDF against a dead port.
MAX_CONN_WAITS = max(1, int(os.environ.get("GROBID_CONN_RETRIES", "60")))   # per PDF
CONN_WAIT_CAP  = int(os.environ.get("GROBID_CONN_WAIT_CAP", "30"))          # s, backoff ceiling
_SERVER_LOCK   = threading.Lock()    # serialises container recovery across workers

# Reuse keep-alive connections (pool sized to the thread count) instead of opening
# a fresh TCP connection per POST -- fewer connections held against the server.
_SESSION = requests.Session()
_adapter = requests.adapters.HTTPAdapter(pool_connections=N_THREADS,
                                         pool_maxsize=N_THREADS, max_retries=0)
_SESSION.mount("http://", _adapter)
_SESSION.mount("https://", _adapter)


def _is_busy(status, text):
    """True if GROBID rejected the request due to saturation (safe to retry)."""
    if status in BUSY_STATUS:
        return True
    t = (text or "").lower()
    return any(tok in t for tok in BUSY_TOKENS)

# --- Bring-up automation (set GROBID_AUTOLAUNCH=0 to disable) --------------- #
AUTOLAUNCH     = os.environ.get("GROBID_AUTOLAUNCH", "1").strip() not in ("0", "false", "no")
# Default to the CRF-only image (~540 MB): it starts Jetty in ~10s and uses a
# fraction of the RAM. The full deep-learning image "grobid/grobid:0.8.2"
# (~9.5 GB) stalls for minutes loading TensorFlow models on startup and is the
# memory hog that OOM-crashed the server after a handful of fulltext parses.
# Set GROBID_IMAGE=grobid/grobid:0.8.2 to use the full DL models instead.
DOCKER_IMAGE   = os.environ.get("GROBID_IMAGE", "grobid/grobid:0.8.2-crf")
CONTAINER_NAME = os.environ.get("GROBID_CONTAINER", "grobid_xml")
# The container is started DETACHED and managed by this script (no separate window).
DOCKER_RUN_CMD = "docker run -d --rm --name %s -p 8070:8070 %s" % (CONTAINER_NAME, DOCKER_IMAGE)
DOCKER_WAIT    = int(os.environ.get("GROBID_DOCKER_WAIT", "180"))   # s, engine ready
DOCKER_PULL_WAIT = int(os.environ.get("GROBID_PULL_WAIT", "1800"))  # s, first-time image pull
SERVER_WAIT    = int(os.environ.get("GROBID_SERVER_WAIT", "300"))   # s, /api/isalive

# Where Docker Desktop / Chrome usually live on Windows (first hit wins).
DOCKER_DESKTOP_PATHS = [
    os.environ.get("DOCKER_DESKTOP_EXE"),
    r"C:\Program Files\Docker\Docker\Docker Desktop.exe",
    os.path.expandvars(r"%LOCALAPPDATA%\Docker\Docker Desktop.exe"),
]
CHROME_PATHS = [
    os.environ.get("CHROME_EXE"),
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    os.path.expandvars(r"%LOCALAPPDATA%\Google\Chrome\Application\chrome.exe"),
]

PROMPT = (
    "design a strategy and write a script called 'grobid_xml.py' that converts pdfs in "
    "'ncbi_pdfs_grobid' into xmls; use '\\bacs\\batch_grobid.py' as an example; the prerequisites for "
    "this step are provided as comments (after #) in '\\bacs\\batch_grobid.py'; enable a summary file "
    "called 'summaries/grobid_xml_summary.html' as part of this script; do not execute 'grobid_xml.py'; "
    "include this prompt and the strategy in 'summaries/grobid_xml_summary.html'"
)


# --------------------------------------------------------------------------- #
# 1. GROBID server check + 2. batch conversion (after bacs/batch_grobid.py)
# --------------------------------------------------------------------------- #
def server_alive():
    """Return True if the GROBID server answers its health endpoint."""
    try:
        r = requests.get(GROBID_SERVER.rstrip("/") + "/api/isalive", timeout=5)
        return r.status_code == 200 and "true" in r.text.lower()
    except requests.RequestException:
        return False


# --------------------------------------------------------------------------- #
# 0. Bring-up automation: Docker Desktop -> GROBID container -> Chrome console
# --------------------------------------------------------------------------- #
def _find_exe(candidates):
    """First existing path among the candidates (skips None/empty)."""
    for p in candidates:
        if p and os.path.isfile(p):
            return p
    return None


def docker_engine_ready():
    """Return True if `docker info` succeeds (engine up and reachable)."""
    try:
        r = subprocess.run(["docker", "info"], stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, timeout=20)
        return r.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def launch_docker_desktop(wait=DOCKER_WAIT):
    """(1) Start Docker Desktop and block until the engine answers. Returns bool."""
    if docker_engine_ready():
        print("[docker] engine already running", file=sys.stderr)
        return True
    exe = _find_exe(DOCKER_DESKTOP_PATHS) or shutil.which("Docker Desktop")
    if exe:
        print("[docker] launching Docker Desktop: %s" % exe, file=sys.stderr)
        try:
            subprocess.Popen([exe], close_fds=True)
        except OSError as exc:
            print("[docker] could not launch Docker Desktop: %s" % exc, file=sys.stderr)
    else:
        print("[docker] Docker Desktop not found -- start it manually "
              "(or set DOCKER_DESKTOP_EXE).", file=sys.stderr)
    print("[docker] waiting up to %ds for the engine..." % wait, file=sys.stderr)
    deadline = time.time() + wait
    while time.time() < deadline:
        if docker_engine_ready():
            print("[docker] engine ready", file=sys.stderr)
            return True
        time.sleep(3)
    print("[docker] engine NOT ready after %ds" % wait, file=sys.stderr)
    return False


def _image_present(image):
    """True if the Docker image is already present locally."""
    try:
        r = subprocess.run(["docker", "image", "inspect", image],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
        return r.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def launch_grobid_server():
    """(2) Start the GROBID container with Docker (detached) and leave it running.

    Runs ``docker run -d`` directly (pulling the image first if missing) instead of
    popping a separate console window, so the container reliably comes up and this
    script can then talk to it on port 8070. Idempotent: if the server is already
    answering, or the container already exists, it is reused.
    """
    if server_alive():
        print("[grobid] server already up at %s" % GROBID_SERVER, file=sys.stderr)
        return

    # Pull the image on first use (it is ~1 GB, so this can take a few minutes).
    if not _image_present(DOCKER_IMAGE):
        print("[grobid] pulling image %s (first run -- may take several minutes)..."
              % DOCKER_IMAGE, file=sys.stderr)
        try:
            subprocess.run(["docker", "pull", DOCKER_IMAGE], timeout=DOCKER_PULL_WAIT)
        except (OSError, subprocess.SubprocessError) as exc:
            print("[grobid] docker pull failed: %s" % exc, file=sys.stderr)

    # Remove any stale container of the same name, then start fresh, detached.
    subprocess.run(["docker", "rm", "-f", CONTAINER_NAME],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    cmd = ["docker", "run", "-d", "--rm", "--name", CONTAINER_NAME,
           "-p", "8070:8070", DOCKER_IMAGE]
    print("[grobid] starting container: %s" % " ".join(cmd), file=sys.stderr)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            print("[grobid] docker run failed: %s"
                  % (r.stderr or r.stdout or "").strip(), file=sys.stderr)
        else:
            print("[grobid] container started: %s" % r.stdout.strip()[:12], file=sys.stderr)
    except (OSError, subprocess.SubprocessError) as exc:
        print("[grobid] could not start container: %s" % exc, file=sys.stderr)


def wait_for_server(wait=SERVER_WAIT):
    """Block until the GROBID server answers /api/isalive (or the timeout)."""
    print("[grobid] waiting up to %ds for %s/api/isalive..."
          % (wait, GROBID_SERVER.rstrip("/")), file=sys.stderr)
    deadline = time.time() + wait
    while time.time() < deadline:
        if server_alive():
            print("[grobid] server is up", file=sys.stderr)
            return True
        time.sleep(3)
    return False


def ensure_server_up(wait=SERVER_WAIT):
    """Make sure GROBID is answering again; (re)launch the container if it died.

    Called by a worker that hit a connection-refused error (the container
    crashed, usually OOM). Thread-safe: the first worker to grab the lock does
    the recovery while the others block, then everyone re-checks -- so the
    container is relaunched at most once per outage, not once per worker. With
    AUTOLAUNCH off we can only wait for the operator to bring it back. Returns
    True if the server is reachable on exit."""
    if server_alive():
        return True
    with _SERVER_LOCK:
        if server_alive():                 # another worker already revived it
            return True
        print("[grobid] server unreachable (container appears to have died) -- "
              "attempting recovery", file=sys.stderr)
        if AUTOLAUNCH:
            launch_docker_desktop()        # engine may also have stopped
            launch_grobid_server()         # docker run -d (no-op if already up)
            return wait_for_server(wait)
        print("[grobid] GROBID_AUTOLAUNCH disabled -- waiting for the server to "
              "come back (restart it manually)", file=sys.stderr)
        return wait_for_server(wait)


def _select_console_service():
    """Best-effort: drive Chrome with Selenium to click the 'TEI' tab and pick the
    'Process Fulltext Document' service. Returns True only if it fully succeeded.

    The console processes one PDF at a time, so this is for visual confirmation --
    the batch conversion itself runs through the GROBID REST API (same service)."""
    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import Select, WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
    except ImportError:
        return False
    try:
        opts = webdriver.ChromeOptions()
        opts.add_experimental_option("detach", True)   # leave the window open
        chrome = _find_exe(CHROME_PATHS)
        if chrome:
            opts.binary_location = chrome
        drv = webdriver.Chrome(options=opts)
        drv.get(GROBID_URL)
        wait = WebDriverWait(drv, 30)
        # Switch to the TEI tab, then choose the Process Fulltext Document service.
        wait.until(EC.element_to_be_clickable((By.LINK_TEXT, "TEI"))).click()
        sel = wait.until(EC.presence_of_element_located((By.ID, "selectedService")))
        Select(sel).select_by_value(SERVICE)            # processFulltextDocument
        print("[browser] Chrome console: 'TEI' tab + 'Process Fulltext Document' selected.",
              file=sys.stderr)
        return True
    except Exception as exc:                            # noqa: BLE001 - best effort
        print("[browser] Selenium console automation skipped: %s" % exc, file=sys.stderr)
        return False


def open_browser_ui():
    """(3) Open Chrome at the GROBID console; select TEI / Process Fulltext Document.

    Tries Selenium first (so the tab + service are actually selected); if Selenium
    is unavailable or fails, just opens Chrome at the console and tells the operator
    which tab/service to pick."""
    if _select_console_service():
        return
    chrome = _find_exe(CHROME_PATHS)
    try:
        if chrome:
            subprocess.Popen([chrome, GROBID_URL])
            print("[browser] opened Chrome at %s" % GROBID_URL, file=sys.stderr)
        else:
            webbrowser.open(GROBID_URL)
            print("[browser] opened default browser at %s (Chrome not found)"
                  % GROBID_URL, file=sys.stderr)
    except OSError as exc:
        print("[browser] could not open a browser: %s" % exc, file=sys.stderr)
    print("[browser] In the console, choose the 'TEI' tab and the "
          "'Process Fulltext Document' service -- the same service this script "
          "runs in batch via the GROBID REST API.", file=sys.stderr)


def _reason(status, text):
    """Human-readable explanation for a non-success GROBID HTTP status."""
    snippet = (text or "").strip().replace("\n", " ")[:160]
    names = {
        CONN_DOWN: "connection refused (GROBID server down / container crashed)",
        200: "ok",
        204: "empty TEI (no content extracted)",
        400: "bad request / failed to open PDF file",
        408: "read timed out",
        429: "too many requests",
        500: "server or connection error",
        502: "bad gateway",
        503: "server busy",
        504: "gateway timeout",
    }
    label = names.get(status, "HTTP %s" % status)
    return "%s%s" % (label, (" -- %s" % snippet) if snippet else "")


def process_fulltext(pdf):
    """POST one PDF to the GROBID REST API; return (status_code, tei_text).

    This is exactly the call grobid_client makes under the hood
    (POST <server>/api/processFulltextDocument with the PDF as multipart 'input'),
    done directly with ``requests`` so the script needs no extra client package.
    """
    url = GROBID_SERVER.rstrip("/") + "/api/" + SERVICE
    data = {
        "generateIDs": "0",
        "consolidateHeader": "0",
        "consolidateCitations": "0",
        "includeRawCitations": "0",
        "includeRawAffiliations": "0",
        "teiCoordinates": "0",
        "segmentSentences": "0",
    }
    with open(pdf, "rb") as fh:
        files = {"input": (os.path.basename(pdf), fh, "application/pdf")}
        resp = _SESSION.post(url, files=files, data=data, timeout=(10, TIMEOUT))
    return resp.status_code, resp.text


# Sentinel status used for "the server is not listening" (container crashed),
# as opposed to an HTTP error code returned by a live server.
CONN_DOWN = "conn-down"


def _convert_one(pdf, tei):
    """Convert one PDF to TEI, retrying transient failures. Returns (ok, status, reason).

    Writes ``tei`` only on a non-empty HTTP 200. Server saturation ("max
    connections exceeded": 503 / "concurrent connection" / busy) means the request
    was rejected before parsing, so it is retried indefinitely (up to
    MAX_BUSY_WAITS) with jittered, capped backoff and does NOT consume the attempt
    budget. Other transient statuses (read timeout, 5xx, dropped connections) are
    retried with linear backoff up to MAX_ATTEMPTS; permanent ones (400 / "failed
    to open PDF file" / empty result) return at once so one bad PDF never aborts
    the batch."""
    name = os.path.basename(pdf)
    status, text = None, ""
    attempt, busy_waits, conn_waits = 0, 0, 0
    while True:
        attempt += 1
        try:
            status, text = process_fulltext(pdf)    # processFulltextDocument via REST
        except requests.exceptions.ConnectionError as exc:
            # Nothing listening on the port: the container died (usually OOM).
            status, text = CONN_DOWN, str(exc)
        except Exception as exc:                    # noqa: BLE001 - keep the batch alive
            status, text = 500, str(exc)

        if status == 200 and text and text.strip():
            try:
                with open(tei, "w", encoding="utf-8") as fh:
                    fh.write(text)
            except OSError as exc:
                return (False, status, "write failed: %s" % exc)
            return (True, status, "ok")

        # Server DOWN (connection refused) -> the container crashed. Revive it
        # (once, under the lock) and keep retrying without spending an attempt,
        # so the whole batch pauses and resumes rather than skipping every
        # remaining PDF against a dead port.
        if status == CONN_DOWN:
            conn_waits += 1
            if conn_waits <= MAX_CONN_WAITS:
                back = min(RETRY_SLEEP * conn_waits, CONN_WAIT_CAP) + random.uniform(0, RETRY_SLEEP)
                print("[grobid]   %s server down %d/%d: %s -- recovering, retry in %.1fs"
                      % (name, conn_waits, MAX_CONN_WAITS,
                         _reason(status, text), back), file=sys.stderr)
                ensure_server_up()              # relaunch container if needed (thread-safe)
                time.sleep(back)
                attempt -= 1                    # connection refused did no work
                continue
            break

        # Server saturated -> request did no work, so always retry. Jittered,
        # capped backoff desynchronises the threads; this does not use up an attempt.
        if _is_busy(status, text):
            busy_waits += 1
            if busy_waits <= MAX_BUSY_WAITS:
                wait = min(RETRY_SLEEP * busy_waits, BUSY_WAIT_CAP) + random.uniform(0, RETRY_SLEEP)
                print("[grobid]   %s busy %d/%d: %s -- waiting %.1fs"
                      % (name, busy_waits, MAX_BUSY_WAITS, _reason(status, text), wait),
                      file=sys.stderr)
                time.sleep(wait)
                attempt -= 1
                continue
            break

        # Other transient errors: retry up to MAX_ATTEMPTS; permanent ones bail now.
        if status in RETRY_STATUS and attempt < MAX_ATTEMPTS:
            wait = RETRY_SLEEP * attempt
            print("[grobid]   %s attempt %d/%d: %s -- retrying in %ds"
                  % (name, attempt, MAX_ATTEMPTS, _reason(status, text), wait),
                  file=sys.stderr)
            time.sleep(wait)
            continue
        break

    return (False, status, _reason(status, text))


def run_grobid():
    """Convert PDFs in INPUT_DIR to TEI in OUTPUT_DIR, per file, with retries.

    (1) Skips any PDF that already has a non-empty TEI (resume -- never repeats
        completed work; set GROBID_FORCE=1 to reconvert everything).
    (2) Recovers from per-file errors: transient ones are retried (see
        _convert_one), permanent ones are logged and skipped. Returns True on run.
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # (1) Build the work list, skipping PDFs already converted.
    pdfs = sorted(glob.glob(os.path.join(INPUT_DIR, "*.pdf")))
    todo, skipped = [], 0
    for pdf in pdfs:
        stem = os.path.splitext(os.path.basename(pdf))[0]
        tei = os.path.join(OUTPUT_DIR, stem + TEI_SUFFIX)
        if not FORCE and os.path.exists(tei) and os.path.getsize(tei) > 0:
            skipped += 1
            continue
        todo.append((pdf, tei))

    print("[grobid] %d PDFs: %d already converted (skipped), %d to process "
          "(service=%s, threads=%d, attempts=%d, force=%s)"
          % (len(pdfs), skipped, len(todo), SERVICE, N_THREADS, MAX_ATTEMPTS, FORCE),
          file=sys.stderr)

    # (2) Convert the rest, N at a time; each task retries/recovers on its own.
    converted, failed, fail_list = 0, 0, []
    with ThreadPoolExecutor(max_workers=max(1, N_THREADS)) as pool:
        futures = {pool.submit(_convert_one, pdf, tei): pdf for pdf, tei in todo}
        for fut in as_completed(futures):
            pdf = futures[fut]
            name = os.path.basename(pdf)
            try:
                ok, _status, reason = fut.result()
            except Exception as exc:                # noqa: BLE001 - defensive
                ok, reason = False, "unexpected: %s" % exc
            if ok:
                converted += 1
                print("[grobid]   OK   %s" % name, file=sys.stderr)
            else:
                failed += 1
                fail_list.append((name, reason))
                print("[grobid]   SKIP %s (%s)" % (name, reason), file=sys.stderr)

    print("[grobid] done: %d converted, %d skipped (already done), %d failed/empty"
          % (converted, skipped, failed), file=sys.stderr)
    if fail_list:
        print("[grobid] %d not converted this pass (re-run to retry):" % len(fail_list),
              file=sys.stderr)
        for name, reason in fail_list:
            print("[grobid]     - %s: %s" % (name, reason), file=sys.stderr)
    return True


# --------------------------------------------------------------------------- #
# 4. TEI structure analysis (tolerant regex; GROBID TEI namespace)
# --------------------------------------------------------------------------- #
RE_TITLE   = re.compile(r'<title\b[^>]*type="main"[^>]*>(.*?)</title>', re.S)
RE_ABSTRACT= re.compile(r"<abstract\b[^>]*>(.*?)</abstract>", re.S)
RE_BODY    = re.compile(r"<body\b[^>]*>(.*?)</body>", re.S)
RE_DIV     = re.compile(r"<div\b")
RE_HEAD    = re.compile(r"<head\b[^>]*>(.*?)</head>", re.S)
RE_BIBL    = re.compile(r"<biblStruct\b")
RE_TAGS    = re.compile(r"<[^>]+>")


def _text(s):
    return RE_TAGS.sub(" ", s or "").strip()


def analyse_tei(text):
    """Return structure facts for one GROBID TEI document."""
    tm = RE_TITLE.search(text)
    has_title = bool(tm and _text(tm.group(1)))
    am = RE_ABSTRACT.search(text)
    has_abstract = bool(am and _text(am.group(1)))
    bm = RE_BODY.search(text)
    body = bm.group(1) if bm else ""
    n_divs = len(RE_DIV.findall(body))
    heads = [_text(h).lower() for h in RE_HEAD.findall(body) if _text(h)]
    n_refs = len(RE_BIBL.findall(text))
    return {
        "has_title": has_title,
        "has_abstract": has_abstract,
        "has_body": n_divs > 0,
        "n_divs": n_divs,
        "heads": heads,
        "has_refs": n_refs > 0,
        "n_refs": n_refs,
    }


# Map a GROBID <head> to a canonical IMRaD bucket (same spirit as xml_structure.py).
def classify_head(h):
    out = set()
    if "introduction" in h or h == "background" or h.startswith("background"):
        out.add("Introduction")
    if any(k in h for k in ("method", "materials", "experimental", "statistical")):
        out.add("Methods")
    if "result" in h or "finding" in h:
        out.add("Results")
    if "discussion" in h:
        out.add("Discussion")
    if "conclusion" in h or h == "summary":
        out.add("Conclusion")
    return out


# --------------------------------------------------------------------------- #
# Summary HTML
# --------------------------------------------------------------------------- #
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
 pre { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: .9rem 1rem; overflow-x: auto; font-size: .85rem; line-height: 1.4; white-space: pre-wrap; }
 .bar { display:inline-block; height:.72em; background:#3b7dd8; border-radius:2px; }
 .bar.g { background:#2da44e; } .bar.o { background:#bf8700; }
 .dim { color: #888; font-size: .85em; }
 .key  { background: #ddf4ff; border-left: 4px solid #0969da; padding: .6rem .9rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 .warn { background: #fff8c5; border-left: 4px solid #d4a72c; padding: .6rem .9rem; margin: 1rem 0; border-radius: 0 4px 4px 0; }
 ol.strategy > li { margin: .45rem 0; }
"""

STRATEGY_HTML = """
<ol class="strategy">
  <li><strong>Bring-up automation (<code>GROBID_AUTOLAUNCH</code>, on by default).</strong> The script
      (1) launches <strong>Docker Desktop</strong> and waits for the engine, (2) starts the GROBID
      container <strong>detached with Docker</strong>
      (<code>docker run -d --rm --name grobid_xml -p 8070:8070 grobid/grobid:0.8.2</code>, pulling the
      image first if missing) and waits for the server, and (3) optionally opens <strong>Chrome</strong>
      at <code>http://localhost:8070/</code> for visual confirmation. The batch itself converts each PDF
      by POSTing it to the GROBID REST API (<code>/api/processFulltextDocument</code>). Any step that
      cannot complete is logged and the script carries on.</li>
  <li><strong>Prerequisite &mdash; a running GROBID server.</strong> GROBID runs in Docker on port 8070
      (<code>docker run -t --rm -p 8070:8070 grobid/grobid:0.8.2</code>). The script pings
      <code>/api/isalive</code> and will not start the batch if the server is unreachable &mdash; but it
      still (re)writes this summary so the page always reflects the current state.</li>
  <li><strong>Convert per file over the GROBID REST API</strong> (after <code>bacs/batch_grobid.py</code>).
      Each PDF in <code>ncbi_pdfs_grobid/</code> is POSTed to
      <code>/api/processFulltextDocument</code> with <code>requests</code>, N at a time via a thread
      pool, writing one <code>&lt;PMCID&gt;.grobid.tei.xml</code> per PDF into
      <code>grobid_xmls/</code>.</li>
  <li><strong>Resume-friendly &amp; fault-tolerant.</strong> PDFs that already have a non-empty TEI are
      skipped, so a re-run only fills the gaps (set <code>GROBID_FORCE=1</code> to reconvert all).
      Transient failures (<em>read timed out</em> / HTTP&nbsp;503 busy / 5xx / dropped connections) are
      retried with linear backoff; permanent ones (<em>HTTP&nbsp;400</em>, <em>failed to open PDF file</em>,
      empty result) are logged and skipped so one bad PDF never aborts the run.</li>
  <li><strong>Summarise &amp; structure-check.</strong> Input PDFs are matched against output TEIs to
      report the conversion rate, then each TEI is lightly checked for title / abstract / body sections
      / references &mdash; the same structural lens used in
      <code>summaries/ncbi_xml.html</code>, so the recovered full text can be compared with the NCBI
      XML corpus.</li>
</ol>
"""


def _fmt(x):
    return format(x, ",")


def build_summary(server_up):
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    pdfs = sorted(glob.glob(os.path.join(INPUT_DIR, "*.pdf")))
    n_pdf = len(pdfs)

    converted = 0
    empty = 0
    tei_bytes = 0
    cov = {"Title": 0, "Abstract": 0, "Body sections": 0, "References": 0}
    sec_cov = {"Introduction": 0, "Methods": 0, "Results": 0, "Discussion": 0, "Conclusion": 0}
    full_imrad = 0
    div_total = 0

    for pdf in pdfs:
        stem = os.path.splitext(os.path.basename(pdf))[0]
        tei = os.path.join(OUTPUT_DIR, stem + TEI_SUFFIX)
        if not (os.path.exists(tei) and os.path.getsize(tei) > 0):
            continue
        converted += 1
        try:
            tei_bytes += os.path.getsize(tei)
            with open(tei, "r", encoding="utf-8", errors="ignore") as fh:
                text = fh.read()
        except OSError:
            continue
        r = analyse_tei(text)
        if not (r["has_title"] or r["has_abstract"] or r["has_body"]):
            empty += 1
        if r["has_title"]:
            cov["Title"] += 1
        if r["has_abstract"]:
            cov["Abstract"] += 1
        if r["has_body"]:
            cov["Body sections"] += 1
        if r["has_refs"]:
            cov["References"] += 1
        div_total += r["n_divs"]
        found = set()
        for h in r["heads"]:
            found |= classify_head(h)
        for s in found:
            sec_cov[s] += 1
        if {"Introduction", "Methods", "Results", "Discussion"} <= found:
            full_imrad += 1

    pending = n_pdf - converted
    rate = (converted / n_pdf * 100) if n_pdf else 0
    mb = tei_bytes / (1024 * 1024)
    mean_divs = (div_total / converted) if converted else 0
    pct = lambda c: (c / converted * 100) if converted else 0

    H = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<title>GROBID PDF&rarr;TEI conversion &mdash; summary</title>",
         "<style>" + CSS + "</style></head><body>"]
    H.append("<h1>GROBID PDF&rarr;TEI conversion &mdash; summary</h1>")
    H.append("<p class='meta'>Generated by <code>grobid_xml.py</code> (after "
             "<code>bacs/batch_grobid.py</code>) &middot; <code>ncbi_pdfs_grobid/*.pdf</code> &rarr; "
             "<code>grobid_xmls/*%s</code> &middot; service <code>%s</code></p>"
             % (TEI_SUFFIX, SERVICE))

    if not server_up:
        H.append("<div class='warn'><strong>GROBID server not reachable at <code>%s</code>.</strong> "
                 "Start it with <code>docker run -t --rm -p 8070:8070 grobid/grobid:0.8.2</code> and "
                 "re-run. Figures below reflect TEI files already on disk.</div>"
                 % _html.escape(GROBID_SERVER))

    cards = [
        (_fmt(n_pdf), "input PDFs"),
        (_fmt(converted), "TEI XMLs produced (%.1f%%)" % rate),
        (_fmt(pending), "not yet converted"),
        (_fmt(empty), "empty/near-empty TEI"),
        ("%.1f MB" % mb, "TEI on disk"),
        ("%.1f" % mean_divs, "mean body sections/doc"),
    ]
    H.append("<div class='stat-grid'>")
    for v, k in cards:
        H.append("<div class='stat'><div class='v'>%s</div><div class='k'>%s</div></div>" % (v, k))
    H.append("</div>")

    H.append("<div class='key'>GROBID converted <strong>%s</strong> of <strong>%s</strong> PDFs "
             "(<strong>%.1f%%</strong>) into TEI XML in <code>grobid_xmls/</code>. Of those, "
             "<strong>%s (%.1f%%)</strong> carry a structured body, <strong>%s (%.1f%%)</strong> an "
             "abstract and <strong>%s (%.1f%%)</strong> a reference list. Row data: the TEI files "
             "themselves.</div>"
             % (_fmt(converted), _fmt(n_pdf), rate, _fmt(cov["Body sections"]), pct(cov["Body sections"]),
                _fmt(cov["Abstract"]), pct(cov["Abstract"]), _fmt(cov["References"]), pct(cov["References"])))

    H.append("<h2>1. Strategy</h2>" + STRATEGY_HTML)

    H.append("<h2>2. Conversion outcome</h2>")
    H.append("<table><thead><tr><th>outcome</th><th class='num'>PDFs</th><th>dist.</th>"
             "<th class='num'>%</th></tr></thead><tbody>")
    omax = max(converted, pending) or 1
    for name, c, cls in (("converted to TEI", converted, "bar g"),
                         ("not yet converted", pending, "bar o")):
        w = max(2, round(c / omax * 260))
        p = (c / n_pdf * 100) if n_pdf else 0
        H.append("<tr><td>%s</td><td class='num'>%s</td>"
                 "<td><span class='%s' style='width:%dpx'></span></td>"
                 "<td class='num dim'>%.1f%%</td></tr>" % (name, _fmt(c), cls, w, p))
    H.append("</tbody></table>")
    H.append("<p class='dim'>'not yet converted' = no TEI on disk for that PDF (GROBID not run yet, "
             "still running, server was down, or the PDF failed to parse).</p>")

    H.append("<h2>3. TEI structure coverage (converted documents)</h2>")
    H.append("<table><thead><tr><th>component</th><th class='num'>documents</th><th>coverage</th>"
             "<th class='num'>%</th></tr></thead><tbody>")
    cmax = max(cov.values()) or 1
    for name in ("Title", "Abstract", "Body sections", "References"):
        c = cov[name]
        w = max(2, round(c / cmax * 260))
        H.append("<tr><td>%s</td><td class='num'>%s</td>"
                 "<td><span class='bar' style='width:%dpx'></span></td>"
                 "<td class='num dim'>%.1f%%</td></tr>" % (name, _fmt(c), w, pct(c)))
    H.append("</tbody></table>")

    H.append("<h2>4. IMRaD sections recovered (from TEI &lt;head&gt;s)</h2>")
    H.append("<table><thead><tr><th>section</th><th class='num'>documents</th><th>coverage</th>"
             "<th class='num'>%</th></tr></thead><tbody>")
    smax = max(sec_cov.values()) or 1
    for name in ("Introduction", "Methods", "Results", "Discussion", "Conclusion"):
        c = sec_cov[name]
        w = max(2, round(c / smax * 260))
        H.append("<tr><td>%s</td><td class='num'>%s</td>"
                 "<td><span class='bar g' style='width:%dpx'></span></td>"
                 "<td class='num dim'>%.1f%%</td></tr>" % (name, _fmt(c), w, pct(c)))
    H.append("</tbody></table>")
    H.append("<p class='dim'><strong>%s (%.1f%%)</strong> of converted documents expose all four IMRaD "
             "sections (Introduction, Methods, Results, Discussion) as TEI section heads.</p>"
             % (_fmt(full_imrad), pct(full_imrad)))

    H.append("""<h2>5. Caveats</h2><ul>
  <li><strong>Requires the GROBID server.</strong> No conversion happens unless the Docker server on
      port 8070 is running; this page then only reflects whatever TEI files already exist.</li>
  <li><strong>GROBID is a heuristic ML parser.</strong> Section heads and references are reconstructed
      from PDF layout, so coverage and labels are approximate &mdash; especially for scanned/OCR PDFs,
      which may yield sparse or empty TEI.</li>
  <li><strong>Counts reflect the files on disk at generation time.</strong> Re-run after the batch
      completes (or after starting the server) and this page refreshes; <code>force=False</code> means
      a re-run only fills the gaps.</li>
</ul>""")
    H.append("</body></html>")

    with open(SUMMARY_HTML, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    print("[summary] wrote %s (%d/%d converted)" % (SUMMARY_HTML, converted, n_pdf),
          file=sys.stderr)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    os.makedirs(INPUT_DIR, exist_ok=True)    # create input dir if missing
    if not os.path.isdir(INPUT_DIR):
        raise SystemExit("error: input dir not found: %s" % INPUT_DIR)
    n_pdf = len(glob.glob(os.path.join(INPUT_DIR, "*.pdf")))
    print("[grobid_xml] %d PDFs in %s" % (n_pdf, INPUT_DIR), file=sys.stderr)

    # 0. Bring-up: Docker Desktop -> GROBID container (detached) -> Chrome console.
    if AUTOLAUNCH:
        launch_docker_desktop()              # (1) launch Docker Desktop
        launch_grobid_server()               # (2) docker run -d (detached container)
        if wait_for_server():
            open_browser_ui()                # (3) Chrome at the TEI / Process Fulltext Document console
        else:
            print("[grobid_xml] GROBID server did not come up within %ds." % SERVER_WAIT,
                  file=sys.stderr)
    else:
        print("[grobid_xml] GROBID_AUTOLAUNCH disabled -- expecting an already-running server.",
              file=sys.stderr)

    up = server_alive()
    if not up:
        print("[grobid_xml] GROBID server NOT reachable at %s -- start the Docker server "
              "(see PREREQUISITES) and re-run. Writing summary from existing TEI only."
              % GROBID_SERVER, file=sys.stderr)
    else:
        try:
            run_grobid()
        except Exception as exc:                       # noqa: BLE001 - report & still summarise
            print("[grobid_xml] GROBID run error: %s" % exc, file=sys.stderr)

    # Always (re)write the summary so the page reflects the current state.
    build_summary(up)


if __name__ == "__main__":
    main()
