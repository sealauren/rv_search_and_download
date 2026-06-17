# rv_search_and_download

Given a NASA Exoplanet Archive (NEA) export, find and download the radial
velocity (RV) data behind each target star's planet masses.

For every host star in the input table, this pipeline:

1. **Searches VizieR** for an RV table associated with any paper the NEA
   cites for that system (matched by ADS bibcode).
2. **Falls back to ADS/arXiv** when a star has a literature Msini mass (proof
   an RV orbit was fit somewhere) but no VizieR table turned up -- reading
   the paper's abstract, and if needed its arXiv full text, to identify which
   spectrograph or archive (e.g. DACE, the California Legacy Survey) was
   actually used.
3. **Downloads the resolved data** -- a VizieR table, a redirect to the
   California Legacy Survey's table in Rosenthal et al. (2021), or a live
   query to the public DACE archive -- and saves it to
   `downloaded_rv_tables/<host>/`, with a commented header describing the
   source and every column's name/unit/meaning.

## Why

NEA mass/orbit columns cite *a* reference, but that reference is often a
paper that *reused* archival RVs rather than publishing its own VizieR
table (common for nearby, well-studied stars). This pipeline chases that
chain down automatically instead of requiring a manual literature search
per target.

## Setup

```bash
pip install -e .
```

This installs the package in editable mode and adds four commands to your
shell: `rv-search-and-download` (the full pipeline) and `rv-vizier-search`,
`rv-ads-search`, `rv-download-tables` for running any single stage on its
own.

Get a free ADS API token from your [ADS account settings](https://ui.adsabs.harvard.edu/user/settings/token)
and export it -- needed for stage 2 (the VizieR-only stage 1 and DACE
queries don't need one):

```bash
export ADS_API_TOKEN=your_token_here
```

## Quick start

A tiny 3-system sample (`sample_data/sample_data.csv`) is bundled so you
can try the pipeline immediately -- it deliberately covers all three
outcomes: a direct VizieR hit (K2-138), a DACE redirect (HD 20781), and an
unresolved case (tau Cet, whose RV paper isn't on VizieR or arXiv).

```bash
rv-search-and-download --catalog sample_data/sample_data.csv
```

Expected: 2 of the 3 hosts get a downloaded RV table; tau Cet is reported
as unresolved (with the implicated reference named, for manual follow-up).

### Running on your own data

Download a CSV from the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/)
(the "Planetary Systems Composite Parameters" table works well) and point
the pipeline at it. Both the raw archive export (with its leading `#`
metadata block) and a pre-stripped CSV are handled automatically.

```bash
rv-search-and-download --catalog PSCompPars_2026.01.14.csv --host "WASP-12,GJ 486"
```

Omit `--host` to process every host in the file -- for the full archive
(thousands of hosts) this means thousands of live VizieR/ADS requests, so
scoping to specific targets first is strongly recommended. Run
`rv-search-and-download --help` for all options.

### Using it as a library

```python
from rv_search_and_download import run_pipeline

result = run_pipeline(
    planet_catalog="sample_data/sample_data.csv",
    host="K2-138,HD 20781",
)
result["vizier_results"]      # DataFrame: one row per host/reference/matched-table
result["ads_instruments"]     # DataFrame: ADS/arXiv lookups (None if stage 2 wasn't needed)
result["downloaded_sources"]  # {hostname: [source dict, ...]} actually downloaded
result["unresolved_hosts"]    # hosts with a reported Msini but no downloadable source
```

## Project layout

```
rv_search_and_download/
├── pyproject.toml
├── src/rv_search_and_download/
│   ├── vizier_search.py    # stage 1
│   ├── ads_search.py       # stage 2
│   ├── download.py         # stage 3
│   └── pipeline.py         # run_pipeline() + the rv-search-and-download CLI
├── sample_data/
│   └── sample_data.csv     # tiny 3-system demo input
├── cache/                  # gitignored; bibcode -> lookup caches, created on first run
└── downloaded_rv_tables/   # gitignored; final output, created on first run
```

## How it works

| Stage | Module | Standalone command | What it does |
|---|---|---|---|
| 1 | `vizier_search.py` | `rv-vizier-search` | Parses every NEA reference column per host, queries VizieR by ADS bibcode (`-source=`), and flags tables whose title/catalog looks like a genuine per-target RV dataset (excluding VizieR's big homogeneous survey catalogs like Gaia, which aren't precision differential RVs). |
| 2 | `ads_search.py` | `rv-ads-search` | For hosts with an Msini-derived mass but no VizieR hit: checks the cited paper's ADS abstract, then (if needed) its arXiv full text via [ar5iv.org](https://ar5iv.labs.arxiv.org/), for a named RV instrument/survey, or a CLS/DACE mention. |
| 3 | `download.py` | `rv-download-tables` | Resolves each host's concrete source (VizieR table / Rosenthal+2021 CLS table / DACE target) and downloads it. |

Each stage caches its lookups by ADS bibcode under `cache/`
(`cache/vizier_bibcode_cache.json`, `cache/ads_abstract_cache.json`), so
reruns -- even against a different input catalog -- never re-query a paper
you've already resolved.

Each of the three modules above also works standalone via its own command
(see the table); `rv-search-and-download` just chains them together. See
the docstring at the top of each file (under `src/rv_search_and_download/`)
for its own `--flags`.

## Known limitations

- The "looks like RV" check is a keyword heuristic on VizieR table titles
  (`radial velocity`, `RV`, `spectroscopic orbit`), not a guarantee the
  table is the *exact* dataset behind a given mass measurement.
- ar5iv's LaTeXML full-text rendering occasionally fails for non-standard
  journal templates (Nature-family papers in particular render as
  near-empty stub pages) -- those papers fall back to "unresolved" even
  though they're real arXiv preprints.
- DACE and VizieR are queried in public/anonymous mode; private or
  not-yet-public datasets won't be found.
- Hosts where only a bare instrument name is found (no VizieR table, no
  CLS/DACE mention) aren't auto-downloaded -- there's no machine-readable
  location to fetch from without a manual archive lookup.

## License

MIT -- see [LICENSE](LICENSE).
