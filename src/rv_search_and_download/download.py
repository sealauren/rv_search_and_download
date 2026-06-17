"""For every host with at least one planet whose mass is Msini-derived (i.e.
an RV source must exist in the literature), resolve a concrete, downloadable
RV data source and save it to downloaded_rv_tables/<host>/.

Combines the outputs of vizier_search.py and ads_search.py:
  - A direct VizieR RV table for that host, if one was found.
  - If a reference mentions the "California Legacy Survey" / "CLS", the
    actual RV table lives in Rosenthal et al. 2021's VizieR catalog
    (J/ApJS/255/8/table6), not the citing paper -- download that instead.
  - If a reference mentions "DACE", query the DACE archive directly for
    that host star's RV time series (public API, no auth needed for
    published data).

Hosts where only a bare instrument name was found (no VizieR table, no CLS/
DACE) are reported but not downloaded -- there's no machine-readable
location to fetch from without manual archive lookup.

Each saved CSV has a '#'-commented header describing the source, reference,
and column names/units/descriptions (pulled from VizieR's table metadata,
or a fixed key-column glossary for DACE).

Usage (after `pip install -e .`):
    rv-download-tables
    rv-download-tables --out-dir downloaded_rv_tables
"""
import argparse
import re
import time
from datetime import date
from pathlib import Path

import pandas as pd
from astroquery.vizier import Vizier
from dace_query.spectroscopy import Spectroscopy

ROSENTHAL_BIBCODE = "2021ApJS..255....8R"
ROSENTHAL_CATALOG = "J/ApJS/255/8"
ROSENTHAL_TABLE = "J/ApJS/255/8/table6"
ROSENTHAL_DISPLAY = "Rosenthal et al. 2021"

DACE_KEY_COLUMNS = {
    "rjd": "Reduced Julian Date (JD - 2400000)",
    "rv": "Radial velocity (m/s)",
    "rv_err": "RV uncertainty (m/s)",
    "instrument_name": "Spectrograph/instrument name",
    "ins_mode": "Instrument mode",
    "program_id": "Observing program ID",
    "pub_bibcode": "ADS bibcode of the publication this RV point came from",
    "drs_qc": "Data-reduction-software quality flag (True = passed QC)",
}


def sanitize(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")


def resolve_sources(vizier_csv, ads_csv):
    """Map hostname -> list of resolved RV source dicts, for hosts that
    have a reported Msini planet mass."""
    viz = pd.read_csv(vizier_csv)
    ads = pd.read_csv(ads_csv)
    msini_hosts = set(viz.loc[viz["rv_provenance"], "hostname"].unique())

    sources = {}

    direct = viz[viz["looks_like_rv"] & viz["hostname"].isin(msini_hosts)]
    for _, r in direct.iterrows():
        sources.setdefault(r["hostname"], []).append({
            "type": "vizier", "catalog": r["vizier_catalog"], "table_id": r["vizier_table"],
            "bibcode": r["bibcode"], "display": r["display"],
        })

    for _, r in ads[ads["mentions_cls"] & ads["hostname"].isin(msini_hosts)].iterrows():
        existing = sources.setdefault(r["hostname"], [])
        if not any(e["table_id"] == ROSENTHAL_TABLE for e in existing if e["type"] == "vizier"):
            existing.append({
                "type": "vizier", "catalog": ROSENTHAL_CATALOG, "table_id": ROSENTHAL_TABLE,
                "bibcode": ROSENTHAL_BIBCODE,
                "display": f"{ROSENTHAL_DISPLAY} (CLS table; redirected from {r['display']})",
            })

    for _, r in ads[ads["mentions_dace"] & ads["hostname"].isin(msini_hosts)].iterrows():
        existing = sources.setdefault(r["hostname"], [])
        if not any(e["type"] == "dace" for e in existing):
            existing.append({"type": "dace", "bibcode": r["bibcode"], "display": r["display"]})

    return sources, msini_hosts


def fetch_vizier_table(v, table_id):
    """Fetch one VizieR table by its full '<catalog>/<table>' id, retrying
    transient failures and falling back to querying the bare parent catalog
    id. VizieR's catalog service intermittently fails to resolve the full
    path for catalogs that only have one table (it works fine for
    multi-table catalogs like J/ApJS/255/8/table6), even though the data is
    there -- querying the parent id alone reliably finds it in that case.
    """
    for attempt in range(3):
        try:
            cats = v.get_catalogs(table_id)
        except Exception as exc:
            print(f"  ! attempt {attempt + 1} failed to fetch {table_id}: {exc}")
            cats = None
        if cats:
            return cats[0]
        time.sleep(3)

    parent_id = table_id.rsplit("/", 1)[0]
    if parent_id != table_id:
        try:
            cats = v.get_catalogs(parent_id)
        except Exception:
            cats = None
        for c in cats or []:
            if c.meta.get("name") == table_id:
                return c
        if cats and len(cats) == 1:
            return cats[0]
    return None


def download_vizier_table(hostname, source, out_dir):
    v = Vizier(columns=["**"], row_limit=-1, timeout=180)
    table = fetch_vizier_table(v, source["table_id"])
    if table is None:
        print(f"  ! no data returned for {source['table_id']}")
        return None

    lines = [
        f"# Host: {hostname}",
        f"# Source: VizieR catalog {source['catalog']} (table {source['table_id']})",
        f"# Reference: {source['display']} ({source['bibcode']})",
        f"# ADS: https://ui.adsabs.harvard.edu/abs/{source['bibcode']}/abstract",
        f"# Retrieved via astroquery.vizier on {date.today().isoformat()}",
        "# Columns:",
    ]
    for c in table.colnames:
        unit = table[c].unit
        desc = table[c].description or ""
        lines.append(f"#   {c} [{unit if unit else '-'}] - {desc}")

    out_path = out_dir / (sanitize(source["table_id"]) + ".csv")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    table.to_pandas().to_csv(out_path, mode="a", index=False)
    print(f"  wrote {out_path} ({len(table)} rows)")
    return out_path


def download_dace_table(hostname, source, out_dir):
    target = hostname.replace(" ", "")
    try:
        df = Spectroscopy.get_timeseries(target=target, output_format="pandas", sorted_by_instrument=False)
    except Exception as exc:
        print(f"  ! DACE query failed for {target}: {exc}")
        return None
    if df is None or len(df) == 0:
        print(f"  ! DACE returned no data for {target}")
        return None

    lines = [
        f"# Host: {hostname}",
        "# Source: DACE archive (public mode, dace_query.spectroscopy.Spectroscopy.get_timeseries)",
        f"# Target queried: {target}",
        f"# Triggered by: {source['display']} ({source['bibcode']}) mentioning DACE",
        f"# Retrieved via dace_query on {date.today().isoformat()}",
        f"# {len(df.columns)} total columns returned; key columns:",
    ]
    for c, desc in DACE_KEY_COLUMNS.items():
        if c in df.columns:
            lines.append(f"#   {c} - {desc}")
    lines.append("# See https://dace.unige.ch for the full column schema and additional ancillary columns.")

    out_path = out_dir / "dace.csv"
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    df.to_csv(out_path, mode="a", index=False)
    print(f"  wrote {out_path} ({len(df)} rows)")
    return out_path


def download_tables(vizier_csv, ads_csv, out_dir="downloaded_rv_tables", sleep=0.5):
    """Resolve and download RV tables for every Msini host, returning
    (sources, unresolved) where `sources` is the hostname -> resolved
    source list and `unresolved` is the list of hosts with no downloadable
    source. This is the importable entry point; main() is a CLI wrapper.

    Clears astroquery's on-disk HTTP cache first: a transient failure
    (timeout) on an earlier run can get cached as an empty result, which
    would otherwise silently poison every later retry of that same table.
    """
    Vizier.clear_cache()

    sources, msini_hosts = resolve_sources(vizier_csv, ads_csv)
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    for hostname in sorted(sources):
        host_dir = out_dir / sanitize(hostname)
        host_dir.mkdir(parents=True, exist_ok=True)
        print(hostname)
        for source in sources[hostname]:
            if source["type"] == "vizier":
                download_vizier_table(hostname, source, host_dir)
            elif source["type"] == "dace":
                download_dace_table(hostname, source, host_dir)
            time.sleep(sleep)

    unresolved = sorted(msini_hosts - set(sources))
    print(f"\nDownloaded RV source(s) for {len(sources)}/{len(msini_hosts)} hosts with a reported Msini.")
    if unresolved:
        print(f"{len(unresolved)} hosts still have no resolved, downloadable RV table source:")
        for h in unresolved:
            print(f"  {h}")

    return sources, unresolved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vizier-input", default="vizier_rv_results.csv")
    parser.add_argument("--ads-input", default="ads_rv_instruments.csv")
    parser.add_argument("--out-dir", default="downloaded_rv_tables")
    parser.add_argument("--sleep", type=float, default=0.5)
    args = parser.parse_args()

    download_tables(args.vizier_input, args.ads_input, args.out_dir, args.sleep)


if __name__ == "__main__":
    main()
