"""Resolve alternate catalog names/aliases for a host star via SIMBAD's
identifier list, cached to disk by hostname so each host is only queried
once, ever.

Two separate problems this fixes (see GH issue #2, "Name crossmatching"):
  - A multi-target VizieR table whose per-row star-id column uses a
    catalog number SIMBAD already knows about (KOI, KIC, TIC, ...) -- not
    just the HD/HIP numbers the NEA export itself carries -- can't be
    filtered down to one host's rows without knowing that host's other
    designations too (see download.host_identifier_candidates).
  - DACE sometimes indexes a target under a different name than the
    host's own common name (e.g. a KIC/TIC id) -- trying every known
    alias before giving up can surface RVs a plain hostname lookup misses
    (see download.download_dace_table).
"""
import json
from pathlib import Path

from astroquery.simbad import Simbad


def load_cache(path):
    path = Path(path)
    return json.loads(path.read_text()) if path.exists() else {}


def save_cache(path, cache):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, indent=1))


def query_simbad_aliases(hostname):
    """Every alias string SIMBAD has on file for `hostname` (its own name
    plus catalog cross-ids: HD, HIP, KOI, KIC, TIC, Gaia DR*, 2MASS, ...),
    or [] if SIMBAD has no match for this name or the query fails."""
    try:
        result = Simbad.query_objectids(hostname)
    except Exception:
        return []
    if result is None:
        return []
    return [str(row[0]) for row in result]


def get_host_aliases(hostname, cache, cache_path="cache/simbad_aliases_cache.json"):
    """Cached wrapper around query_simbad_aliases. `cache` is the dict
    returned by load_cache; updated and persisted to `cache_path` on every
    new (uncached) lookup."""
    if hostname in cache:
        return cache[hostname]
    aliases = query_simbad_aliases(hostname)
    cache[hostname] = aliases
    save_cache(cache_path, cache)
    return aliases
