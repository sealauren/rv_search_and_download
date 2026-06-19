"""
rv_search_and_download.pipeline
================================

One-call wrapper around the four-stage RV-table pipeline:

  1. vizier_search -- for each target star ("host"), aggregate every
                       literature reference the NASA Exoplanet Archive
                       (NEA) cites for that system, and check VizieR
                       (by ADS bibcode) for a matching radial-velocity
                       (RV) data table.
  2. ads_search     -- for stars where step 1 found NO VizieR table,
                       but the NEA still reports an Msini mass (proof
                       an RV orbit was fit somewhere), read the paper's
                       ADS abstract -- and if needed, its arXiv full
                       text -- to figure out which RV instrument or
                       archive (e.g. "DACE", "California Legacy
                       Survey") was actually used.
  3. download       -- for every star where steps 1-2 found a concrete,
                       machine-readable RV source, download the actual
                       RV data and save it to disk under
                       downloaded_rv_tables/<host>/, with a commented
                       header describing where the data came from and
                       what each column means.
  4. aggregate      -- combine every downloaded table for a host into one
                       analysis-ready file: standardized BJD/RV(m/s)/RV
                       error(m/s) columns, per-table systemic RV offsets
                       removed, and likely-duplicate RVs (same host,
                       instrument, and timestamp across different tables)
                       flagged for review. Written to
                       aggregated_rv_tables/<host>_rvs_aggregated.csv.

Each stage writes (and re-reads) a small JSON cache keyed by ADS bibcode, so
re-running this on overlapping target lists -- or on a much bigger catalog
like the full PSCompPars export -- never re-asks VizieR/ADS for a paper it
has already resolved.

----------------------------------------------------------------------------
Use it as a library (API):

    from rv_search_and_download import run_pipeline

    result = run_pipeline(
        planet_catalog="PSCompPars_2026.01.14_08.42.07.csv",
        host="WASP-12,GJ 486",   # one host, a comma-separated list, or None for "every host in the file"
    )
    print(result["vizier_results"].head())   # a pandas DataFrame
    print(result["downloaded_sources"])      # {hostname: [source dicts]}

Use it from the command line (CLI), after `pip install -e .`:

    rv-search-and-download --catalog PSCompPars_2026.01.14_08.42.07.csv --host "WASP-12,GJ 486"
    rv-search-and-download --catalog sample_data/sample_data.csv --verbose
----------------------------------------------------------------------------
"""
import argparse
from pathlib import Path

# Each of these three modules also has its own console-script entry point
# (rv-vizier-search, rv-ads-search, rv-download-tables) and can be run
# standalone -- see the docstring at the top of each file. Importing them
# here just lets us call their logic directly instead of shelling out to
# separate processes.
from . import vizier_search, ads_search, download, aggregate


def run_pipeline(
    planet_catalog,
    host=None,
    output_dir=".",
    vizier_cache="cache/vizier_bibcode_cache.json",
    vizier_readme_cache="cache/vizier_readme_cache.json",
    ads_cache="cache/ads_abstract_cache.json",
    ads_token=None,
    downloaded_tables_dir="downloaded_rv_tables",
    aggregated_tables_dir="aggregated_rv_tables",
    vizier_sleep=1.0,
    ads_sleep=0.2,
    download_sleep=0.5,
    verbose=False,
):
    """Run all four pipeline stages back-to-back for one or more target stars.

    Parameters
    ----------
    planet_catalog : str
        Path to a NASA Exoplanet Archive CSV export (e.g. NEA_Multis.csv or
        a raw PSCompPars_*.csv download -- both formats are handled).
    host : str or None
        A single hostname, a comma-separated list of hostnames (e.g.
        "WASP-12,GJ 486"), or None to process every host in the catalog.
        Restricting to specific hosts is much faster, since stage 1 only
        queries VizieR for references actually cited by those stars.
    output_dir : str
        Where to write the two intermediate results CSVs
        (vizier_rv_results.csv, ads_rv_instruments.csv). The final RV
        tables always go under `downloaded_tables_dir`.
    vizier_cache, ads_cache : str
        Paths to the JSON bibcode caches stage 1 and stage 2 maintain.
        Shared/reused across runs and across different `planet_catalog`
        files, since a bibcode means the same paper everywhere.
    vizier_readme_cache : str
        Path to the JSON cache of catalog -> {table_stem: has_rv_column},
        parsed from each catalog's CDS ReadMe. Stage 1 uses this to confirm
        that a table whose title/description merely mentions "RV" actually
        has an RV column, instead of trusting the keyword match alone.
    ads_token : str or None
        ADS API token (free from https://ui.adsabs.harvard.edu/user/settings/token).
        If None, falls back to the $ADS_API_TOKEN environment variable
        (see ads_search.load_token).
    downloaded_tables_dir : str
        Where the final per-host RV data CSVs get written.
    aggregated_tables_dir : str
        Where stage 4's per-host analysis-ready RV CSVs get written.
    vizier_sleep, ads_sleep, download_sleep : float
        Seconds to pause between live requests to VizieR / ADS+ar5iv / VizieR+DACE
        respectively, to stay polite to those services.
    verbose : bool
        Print progress as each new reference is looked up.

    Returns
    -------
    dict with keys:
        "vizier_results"      -- DataFrame from stage 1 (one row per
                                  host/reference/matched-table combination).
        "ads_instruments"     -- DataFrame from stage 2 (one row per
                                  reference that needed an ADS/full-text
                                  lookup), or None if stage 1 already
                                  resolved every Msini host (stage 2 is
                                  skipped in that case).
        "downloaded_sources"  -- {hostname: [source dict, ...]} actually
                                  downloaded in stage 3.
        "aggregated_rvs"      -- {hostname: path} to stage 4's analysis-ready
                                  CSV, for every host that got at least one
                                  downloaded table.
        "unresolved_hosts"    -- hosts with a reported Msini planet but no
                                  downloadable RV source found by any stage.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    vizier_output_csv = output_dir / "vizier_rv_results.csv"
    ads_output_csv = output_dir / "ads_rv_instruments.csv"

    # --- Stage 1: search VizieR directly for each host's references ---
    print("=== Stage 1/4: searching VizieR for RV tables ===")
    vizier_results = vizier_search.search_vizier(
        input_csv=planet_catalog,
        output_csv=vizier_output_csv,
        cache_path=vizier_cache,
        host=host,
        sleep=vizier_sleep,
        verbose=verbose,
        readme_cache_path=vizier_readme_cache,
    )

    # --- Stage 2: only needed if some provenance reference still has no VizieR table ---
    # (search_ads_instruments reads vizier_output_csv itself and figures out
    # which hosts/references still need a lookup, so we just check here
    # whether that set is empty to decide whether to skip the stage. This is
    # checked per (hostname, bibcode), not per hostname: a host can cite many
    # papers, and a *different*, non-provenance reference having an RV table
    # on VizieR doesn't mean the paper actually credited with the reported
    # mass does too -- see vizier_search.unresolved_provenance_refs.)
    still_missing = set(vizier_search.unresolved_provenance_refs(vizier_results)["hostname"])

    ads_results = None
    if still_missing:
        print(f"\n=== Stage 2/4: {len(still_missing)} host(s) need an ADS/arXiv instrument lookup ===")
        ads_results = ads_search.search_ads_instruments(
            vizier_results_csv=vizier_output_csv,
            output_csv=ads_output_csv,
            cache_path=ads_cache,
            token=ads_token,
            sleep=ads_sleep,
            verbose=verbose,
        )
    else:
        print("\n=== Stage 2/4: skipped (every Msini host already has a VizieR RV table) ===")

    # --- Stage 3: download whatever concrete RV source(s) stages 1-2 found ---
    print("\n=== Stage 3/4: downloading resolved RV tables ===")
    if ads_results is None:
        # No stage-2 output exists yet (it was skipped), but download.download_tables
        # still needs *an* ads CSV to read -- write an empty placeholder with
        # the right columns so resolve_sources() finds no CLS/DACE redirects
        # and just falls back to whatever VizieR already resolved directly.
        import pandas as pd
        empty_ads = pd.DataFrame(columns=[
            "hostname", "bibcode", "display", "title", "arxiv_id", "tier",
            "rv_instruments", "mentions_archival_rv", "mentions_cls", "mentions_dace",
        ])
        empty_ads.to_csv(ads_output_csv, index=False)

    downloaded_sources, unresolved_hosts = download.download_tables(
        vizier_csv=vizier_output_csv,
        ads_csv=ads_output_csv,
        out_dir=downloaded_tables_dir,
        sleep=download_sleep,
    )

    # --- Stage 4: combine each host's downloaded table(s) into one analysis-ready file ---
    print("\n=== Stage 4/4: aggregating downloaded RV tables for analysis ===")
    aggregated_rvs = aggregate.aggregate_all(
        downloaded_dir=downloaded_tables_dir,
        out_dir=aggregated_tables_dir,
        host=host,
    )

    return {
        "vizier_results": vizier_results,
        "ads_instruments": ads_results,
        "downloaded_sources": downloaded_sources,
        "aggregated_rvs": aggregated_rvs,
        "unresolved_hosts": unresolved_hosts,
    }


def main():
    # CLI entry point, installed as the `rv-search-and-download` command.
    # Run `rv-search-and-download --help` to see this.
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalog", required=True,
                         help="Path to a NASA Exoplanet Archive CSV (e.g. NEA_Multis.csv or a PSCompPars_*.csv export)")
    parser.add_argument("--host", default=None,
                         help="One hostname, a comma-separated list (e.g. 'WASP-12,GJ 486'), "
                              "or omit to process every host in the catalog")
    parser.add_argument("--output-dir", default=".", help="Where to write the intermediate results CSVs")
    parser.add_argument("--downloaded-tables-dir", default="downloaded_rv_tables",
                         help="Where to write the final downloaded RV table CSVs")
    parser.add_argument("--aggregated-tables-dir", default="aggregated_rv_tables",
                         help="Where to write the per-host analysis-ready aggregated RV CSVs")
    parser.add_argument("--vizier-cache", default="cache/vizier_bibcode_cache.json")
    parser.add_argument("--vizier-readme-cache", default="cache/vizier_readme_cache.json")
    parser.add_argument("--ads-cache", default="cache/ads_abstract_cache.json")
    parser.add_argument("--ads-token", default=None,
                         help="ADS API token; defaults to the $ADS_API_TOKEN environment variable")
    parser.add_argument("--vizier-sleep", type=float, default=1.0)
    parser.add_argument("--ads-sleep", type=float, default=0.2)
    parser.add_argument("--download-sleep", type=float, default=0.5)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    result = run_pipeline(
        planet_catalog=args.catalog,
        host=args.host,
        output_dir=args.output_dir,
        vizier_cache=args.vizier_cache,
        vizier_readme_cache=args.vizier_readme_cache,
        ads_cache=args.ads_cache,
        ads_token=args.ads_token,
        downloaded_tables_dir=args.downloaded_tables_dir,
        aggregated_tables_dir=args.aggregated_tables_dir,
        vizier_sleep=args.vizier_sleep,
        ads_sleep=args.ads_sleep,
        download_sleep=args.download_sleep,
        verbose=args.verbose,
    )

    # Final wrap-up message for the CLI user.
    n_hosts = result["vizier_results"]["hostname"].nunique()
    n_downloaded = len(result["downloaded_sources"])
    n_aggregated = len(result["aggregated_rvs"])
    print(f"\nDone. {n_downloaded}/{n_hosts} host(s) processed got at least one downloaded RV table "
          f"({n_aggregated} aggregated for analysis).")
    if result["unresolved_hosts"]:
        print(f"{len(result['unresolved_hosts'])} host(s) still have no downloadable RV source: "
              f"{', '.join(result['unresolved_hosts'])}")


if __name__ == "__main__":
    main()
