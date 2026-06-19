# rv_search_and_download

Given a NASA Exoplanet Archive (NEA) export, find and download the radial
velocity (RV) data behind each target star's planet masses.

For every host star in the input table, this pipeline:

1. **Searches VizieR** for an RV table associated with any paper the NEA
   cites for that system (matched by ADS bibcode). A candidate table must
   (a) belong to a journal-published catalog (`J/` class, not VizieR's own
   homogeneous survey catalogs like Gaia), (b) have a title or description
   matching RV keywords, and (c) actually contain an RV column -- verified
   by parsing the catalog's CDS ReadMe byte-by-byte column definitions,
   not just trusting the keyword match alone. Tables that textually
   identify as belonging to a *different* host in a shared multi-target
   catalog are also filtered out.

2. **Falls back to a three-tier ADS/arXiv search** for hosts that have a
   literature Msini mass (proof an RV orbit was fit somewhere) but whose
   provenance reference has no VizieR table:
   - **Tier 1**: the paper's ADS abstract.
   - **Tier 2**: if tier 1 finds nothing and the paper has an arXiv id,
     the full text via [ar5iv.org](https://ar5iv.labs.arxiv.org/) (LaTeXML
     HTML rendering, no PDF parsing needed).
   - **Tier 3**: if tiers 1--2 still find nothing -- including when there's
     no arXiv id, or ar5iv failed to render the paper -- ADS's own
     full-text search index, which is built from the publisher and isn't
     subject to arXiv/LaTeXML rendering failures.

   Each tier looks for named RV instruments/surveys (HARPS, HIRES, etc.),
   and two special-case redirects: a California Legacy Survey mention
   (the data lives in Rosenthal et al. 2021, not the citing paper) and a
   DACE mention or DACE-pipeline spectrograph name (ESPRESSO, HARPS,
   HARPS-N, CORALIE, CARMENES, SOPHIE).

3. **Downloads the resolved data** for every host with a reported Msini
   mass:
   - A direct VizieR table hit from stage 1.
   - Rosenthal et al. 2021's CLS table (J/ApJS/255/8/table6), filtered
     to this host's rows, when a CLS redirect was found in stage 2.
   - A live query to the public DACE archive. DACE is attempted for
     **every** Msini host -- not just those where stage 2 named DACE
     explicitly -- because its public API is a plain target-name lookup
     that needs no textual clue, and an aggregated dataset there is
     valuable even when another source already exists, even if there
     are some duplicate RVs.

   Each downloaded file is saved to `downloaded_rv_tables/<host>/` with
   a commented header describing the source and every column's
   name/unit/meaning.

4. **Aggregates each host's downloaded table(s)** into a single
   analysis-ready file: standardizes time/RV/RV-error into BJD, RV (m/s),
   and RV error (m/s); median-subtracts the RV column within each
   (table, instrument) group to remove pipeline-to-pipeline systemic
   offsets; and flags likely-duplicate RVs (same host and instrument,
   timestamps from different tables within 10 seconds of each other) for
   review, picking a default "preferred" row per a configurable policy.
   Written to `aggregated_rv_tables/<host>_rvs_aggregated.csv` -- this is
   the file meant to be fed straight into RV-fitting tools like `radvel`.

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

This installs the package in editable mode and adds five commands to your
shell: `rv-search-and-download` (the full pipeline) and `rv-vizier-search`,
`rv-ads-search`, `rv-download-tables`, `rv-aggregate-rvs` for running any
single stage on its own.

Get a free ADS API token from your [ADS account settings](https://ui.adsabs.harvard.edu/user/settings/token)
and export it -- needed for stage 2 (the VizieR-only stage 1 and DACE
queries don't need one):

```bash
export ADS_API_TOKEN=your_token_here
```

## Quick start

A tiny 3-system sample (`sample_data/sample_data.csv`) is bundled so you
can try the pipeline immediately -- it covers three different resolution
paths:

- **K2-138**: direct VizieR hit (Lopez et al. 2019, J/A+A/631/A90,
  194 rows) plus an independent DACE dataset (215 rows).
- **HD 20781**: no VizieR table; stage 2 finds CORALIE and HARPS named in
  the Udry et al. 2019 abstract (both are DACE-pipeline instruments), so
  stage 3 queries DACE directly (240 rows).
- **tau Cet**: no VizieR table; stage 2 finds HARPS in the Feng et al.
  2017 abstract, triggering a DACE query (20,202 rows -- tau Cet has
  extensive multi-decade HARPS coverage in the archive).

```bash
rv-search-and-download --catalog sample_data/sample_data.csv
```

Expected: all 3 hosts get at least one downloaded RV table (3/3), and a
matching `aggregated_rv_tables/<host>_rvs_aggregated.csv` for each.

You will see a harmless warning at startup about a missing `.dacerc` file
-- that file is only needed for authenticated/private DACE access; public
queries (which is all this pipeline does) work without it.

The `cache/` directory is populated on first run; subsequent runs against
the same or an overlapping catalog skip live network queries for any
bibcode already resolved.

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
result["aggregated_rvs"]      # {hostname: path} to each host's analysis-ready aggregated CSV
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
│   ├── aggregate.py        # stage 4
│   └── pipeline.py         # run_pipeline() + the rv-search-and-download CLI
├── sample_data/
│   └── sample_data.csv     # tiny 3-system demo input
├── cache/                  # gitignored; bibcode -> lookup caches, created on first run
├── downloaded_rv_tables/   # gitignored; stage 3 output, created on first run
├── aggregated_rv_tables/   # gitignored; stage 4 output, created on first run
└── rv_aggregate_config.json  # optional; local override of stage 4's duplicate-resolution policy
```

## How it works

| Stage | Module | Standalone command | What it does |
|---|---|---|---|
| 1 | `vizier_search.py` | `rv-vizier-search` | Parses every NEA reference column per host, queries VizieR by ADS bibcode (`-source=`), and flags tables that (a) are from a journal-published catalog (`J/` class), (b) match RV keywords in the title/description, and (c) have an actual RV column per the catalog's CDS ReadMe. Tables that positively name a different host (shared multi-target catalogs) are excluded. Writes `vizier_rv_results.csv`. |
| 2 | `ads_search.py` | `rv-ads-search` | For hosts with an Msini-derived mass but no VizieR hit: searches the cited paper's ADS abstract (tier 1), then its arXiv full text via ar5iv.org if available (tier 2), then ADS's own full-text index (tier 3), for a named RV instrument/survey or a CLS/DACE mention. Writes `ads_rv_instruments.csv`. Skipped entirely if stage 1 already resolved every Msini host. |
| 3 | `download.py` | `rv-download-tables` | Resolves each host's concrete source(s) and downloads: direct VizieR table hits, the Rosenthal+2021 CLS table (filtered per host), and DACE (queried for every Msini host regardless of whether stage 2 found an explicit DACE clue). |
| 4 | `aggregate.py` | `rv-aggregate-rvs` | Reads every downloaded table for a host (using the commented header `download.py` already wrote -- no re-querying VizieR/DACE), standardizes time/RV/RV-error into BJD/m/s/m/s, median-subtracts RV within each (table, instrument) group, flags likely-duplicate RVs across tables, and writes `aggregated_rv_tables/<host>_rvs_aggregated.csv`. |

Each stage caches its lookups by ADS bibcode under `cache/`
(`cache/vizier_bibcode_cache.json`, `cache/vizier_readme_cache.json`,
`cache/ads_abstract_cache.json`), so reruns -- even against a different
input catalog -- never re-query a paper you've already resolved.

Each of the four modules also works standalone via its own command
(see the table); `rv-search-and-download` just chains them together. See
the docstring at the top of each file (under `src/rv_search_and_download/`)
for its own `--flags`.

### Duplicate-RV resolution policy (stage 4)

Stage 4 flags RVs from different downloaded tables as likely duplicates
when they share a host and instrument and their timestamps are within 10
seconds of each other (see `duplicate_group` in the output -- every row is
kept for review, not silently dropped). Which row in a `duplicate_group`
is marked `is_preferred` (the one to use for analysis) is controlled by a
small JSON config:

```json
{
  "duplicate_default": "prefer_source",
  "preferred_source_order": ["dace", "vizier"]
}
```

- `"duplicate_default": "prefer_source"` (the default) always prefers
  whichever source appears earliest in `preferred_source_order`.
- `"duplicate_default": "most_recent"` instead prefers whichever row's
  publication is more recent (parsed from its ADS bibcode, or the
  retrieval date as a fallback).

Stage 4 looks for `rv_aggregate_config.json` in the current directory
first (a per-project/per-directory override), then falls back to
`~/.config/rv_search_and_download/config.json` (a global default), then
its built-in defaults shown above. Run `rv-aggregate-rvs --init-config
local` (or `--init-config global`) to write a starter file to edit.

## Known limitations

- The VizieR RV check is a two-step heuristic (keyword match on table
  title/description, then ReadMe column verification) and isn't a guarantee
  the table is the *exact* dataset behind a given mass measurement.
- The CDS ReadMe parser covers the common `J/` journal-catalog format; a
  ReadMe that can't be fetched or doesn't parse falls back to the title
  keyword match alone (noted in the output).
- The Rosenthal+2021 CLS filtering matches HD/GJ/HIP host names to the
  CPS identifier column; hosts with other designations (BD numbers, Giclas,
  lettered multiples) can't be isolated and are skipped rather than
  returning the whole-survey table.
- DACE and VizieR are queried in public/anonymous mode; private or
  not-yet-public datasets won't be found.
- Hosts where only a bare, non-DACE instrument name was found (no VizieR
  table, no CLS/DACE, and DACE returned no data) are reported as unresolved
  -- there's no machine-readable location to fetch from without a manual
  archive lookup.
- Stage 4's duplicate detection requires knowing each row's instrument.
  DACE rows carry a real per-row instrument name, but most VizieR RV
  tables don't have an instrument column at all; stage 4 falls back to a
  best-effort keyword match against known spectrograph names, and a
  per-table placeholder label when that fails -- which means a VizieR
  table and a DACE table from the *same* underlying instrument (e.g. both
  HARPS) won't always be recognized as overlapping, and won't be flagged
  as possible duplicates even though they might be redundant.

## License

MIT -- see [LICENSE](LICENSE).
