#!/usr/bin/env python3
import concurrent.futures
import os
import shutil
import ssl
import sys
import urllib.request
from pathlib import Path

url, target, total = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])
start = target.stat().st_size if target.exists() else 0
if start >= total:
    print(f"complete {target} {start}")
    raise SystemExit(0)

piece_dir = target.parent / (target.name + ".parts")
piece_dir.mkdir(exist_ok=True)
workers = min(12, max(1, (total - start + 4_999_999) // 5_000_000))
step = (total - start + workers - 1) // workers
ranges = [(i, start + i * step, min(total - 1, start + (i + 1) * step - 1)) for i in range(workers)]

ssl_context = ssl._create_unverified_context()

def fetch(item):
    i, a, b = item
    p = piece_dir / f"{i:03d}.part"
    expected = b - a + 1
    current = p.stat().st_size if p.exists() else 0
    if current > expected:
        raise RuntimeError(f"oversized piece {i}: {current} > {expected}")
    attempts = 0
    while current < expected and attempts < 8:
        req_start = a + current
        req = urllib.request.Request(url, headers={"Range": f"bytes={req_start}-{b}", "User-Agent": "Mozilla/5.0"})
        try:
            with urllib.request.urlopen(req, timeout=300, context=ssl_context) as r, p.open("ab") as f:
                if r.status != 206:
                    raise RuntimeError(f"range unsupported: HTTP {r.status}")
                cr = r.headers.get("Content-Range", "")
                if f"bytes {req_start}-{b}/" not in cr:
                    raise RuntimeError(f"unexpected Content-Range {cr}")
                shutil.copyfileobj(r, f, 1024 * 1024)
        except Exception:
            attempts += 1
        current = p.stat().st_size if p.exists() else 0
    if current != expected:
        raise RuntimeError(f"short piece {i}: {current} != {expected}")
    return p

with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
    pieces = list(ex.map(fetch, ranges))

assembled = target.with_suffix(target.suffix + ".assembled")
with assembled.open("wb") as out:
    if target.exists():
        with target.open("rb") as src:
            shutil.copyfileobj(src, out, 1024 * 1024)
    for p in sorted(pieces):
        with p.open("rb") as src:
            shutil.copyfileobj(src, out, 1024 * 1024)
if assembled.stat().st_size != total:
    raise RuntimeError(f"assembled size {assembled.stat().st_size} != {total}")
os.replace(assembled, target)
print(f"completed {target} {target.stat().st_size} bytes using {workers} ranges")
