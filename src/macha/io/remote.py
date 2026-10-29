"""Stream MS1 data out of remote mzML files without storing the raw file."""

from __future__ import annotations

import io
import sys
import time
import urllib.error
import urllib.request

from ..run import MSRun
from .mzml import read_mzml


class _CountingStream(io.RawIOBase):
    def __init__(self, raw):
        self.raw = raw
        self.n = 0

    def readinto(self, b):
        data = self.raw.read(len(b))
        self.n += len(data)
        b[: len(data)] = data
        return len(data)

    def readable(self):
        return True


def fetch_ms1(url: str, retries: int = 6, backoff: float = 30.0, verbose: bool = True, **read_kwargs) -> MSRun:
    """Download ``url`` as a stream and keep only MS1 centroids.

    Retries with linear backoff on HTTP 429/5xx (MassIVE rate limits).
    """
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "macha-fetch/0.2 (research)"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                total = resp.headers.get("Content-Length")
                counter = _CountingStream(resp)
                t0 = time.time()

                def progress(n_spectra):
                    if verbose:
                        mb = counter.n / 1e6
                        tot = f"/{int(total) / 1e6:.0f}" if total else ""
                        rate = mb / max(time.time() - t0, 1e-6)
                        print(f"\r  {n_spectra:>7} spectra  {mb:8.1f}{tot} MB  {rate:5.1f} MB/s", end="", file=sys.stderr)

                run = read_mzml(io.BufferedReader(counter, buffer_size=1 << 20), on_progress=progress, **read_kwargs)
                if verbose:
                    print(file=sys.stderr)
                run.meta.update(url=url, bytes_streamed=counter.n)
                return run
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (429, 500, 502, 503, 504):
                wait = backoff * (attempt + 1)
                if verbose:
                    print(f"  HTTP {e.code}; retrying in {wait:.0f}s", file=sys.stderr)
                time.sleep(wait)
                continue
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = e
            time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"failed to fetch {url}: {last}")
