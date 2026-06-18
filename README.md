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
| 1 | `vizier_search.py` | `rv-vizier-search` | Parses every NEA reference column per host, queries VizieR by ADS bibcode (`-source=`), and flags tables that (a) are from a journal-published catalog (`J/` class), (b) match RV keywords in the title/description, and (c) have an actual RV column per the catalog's CDS ReadMe. Tables that positively name a different host (shared multi-target catalogs) are excluded. Writes `vizier_rv_results.csv`. |
| 2 | `ads_search.py` | `rv-ads-search` | For hosts with an Msini-derived mass but no VizieR hit: searches the cited paper's ADS abstract (tier 1), then its arXiv full text via ar5iv.org if available (tier 2), then ADS's own full-text index (tier 3), for a named RV instrument/survey or a CLS/DACE mention. Writes `ads_rv_instruments.csv`. Skipped entirely if stage 1 already resolved every Msini host. |
| 3 | `download.py` | `rv-download-tables` | Resolves each host's concrete source(s) and downloads: direct VizieR table hits, the Rosenthal+2021 CLS table (filtered per host), and DACE (queried for every Msini host regardless of whether stage 2 found an explicit DACE clue). |

Each stage caches its lookups by ADS bibcode under `cache/`
(`cache/vizier_bibcode_cache.json`, `cache/vizier_readme_cache.json`,
`cache/ads_abstract_cache.json`), so reruns -- even against a different
input catalog -- never re-query a paper you've already resolved.

Each of the three modules also works standalone via its own command
(see the table); `rv-search-and-download` just chains them together. See
the docstring at the top of each file (under `src/rv_search_and_download/`)
for its own `--flags`.

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

## License

MIT -- see [LICENSE](LICENSE).
