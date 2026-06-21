"""Shared helpers: downloading with a simple cache and progress, plus I/O."""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def get_json(url: str, *, timeout: int = 60, retries: int = 3, backoff: float = 2.0) -> dict:
    """GET a URL and parse JSON, retrying on transient errors."""
    req = urllib.request.Request(url, headers={"User-Agent": "pyCmmD/0.1",
                                               "Accept": "application/json"})
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            if attempt < retries - 1:
                time.sleep(min(backoff * (attempt + 1), 30.0))
    raise RuntimeError(f"GET failed after {retries} tries: {url}\n  {last}")


def download(url: str, dest: Path, *, force: bool = False, timeout: int = 120) -> Path:
    """Download ``url`` to ``dest``, skipping if the file already exists.

    Returns the path to the downloaded file. Set ``force=True`` to re-download.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        print(f"[cache] {dest.name} already present ({dest.stat().st_size:,} bytes)")
        return dest

    print(f"[download] {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "pyCmmD/0.1"})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as fh:
        total = int(resp.headers.get("Content-Length", 0))
        read = 0
        chunk = 1 << 20  # 1 MiB
        while True:
            buf = resp.read(chunk)
            if not buf:
                break
            fh.write(buf)
            read += len(buf)
            if total:
                pct = 100 * read / total
                print(f"\r  {read/1e6:7.1f} / {total/1e6:7.1f} MB ({pct:5.1f}%)",
                      end="", file=sys.stderr)
            else:
                print(f"\r  {read/1e6:7.1f} MB", end="", file=sys.stderr)
    print("", file=sys.stderr)
    tmp.replace(dest)
    print(f"[done] {dest.name} ({dest.stat().st_size:,} bytes)")
    return dest


def write_provenance(path: Path, meta: dict) -> None:
    """Write a JSON provenance sidecar next to a processed layer file."""
    meta = {"generated_utc": datetime.now(timezone.utc).isoformat(), **meta}
    path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"[meta] {path.name}")
