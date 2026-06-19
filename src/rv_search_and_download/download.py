"""For every host with at least one planet whose mass is Msini-derived (i.e.
an RV source must exist in the literature), resolve a concrete, downloadable
RV data source and save it to downloaded_rv_tables/<host>/.

Combines the outputs of vizier_search.py and ads_search.py:
  - A direct VizieR RV table for that host, if one was found.
  - If a reference mentions the "California Legacy Survey" / "CLS", the
    actual RV table lives in Rosenthal et al. 2021's VizieR catalog
    (J/ApJS/255/8/table6), not the citing paper -- download that instead.
  - If a reference mentions "DACE", or names a DACE-pipeline spectrograph
    (ESPRESSO, HARPS, HARPS-N, CORALIE, CARMENES, SOPHIE), query the DACE
    archive directly for that host star's RV time series (public API, no
    auth needed for published data).
  - DACE is also queried directly for *every* host with a reported Msini
    mass, regardless of whether a VizieR/CLS source was already resolved
    -- DACE's API needs no textual clue to try, and an independent dataset
    there is valuable even when another source already exists.

Hosts where only a bare, non-DACE instrument name was found (no VizieR
table, no CLS/DACE, and DACE itself returned no data) are reported but not
downloaded -- there's no machine-readable location to fetch from without
manual archive lookup.

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

from .ads_search import DACE_PIPELINE_INSTRUMENTS

ROSENTHAL_BIBCODE = "2021ApJS..255....8R"
ROSENTHAL_CATALOG = "J/ApJS/255/8"
ROSENTHAL_TABLE = "J/ApJS/255/8/table6"
ROSENTHAL_DISPLAY = "Rosenthal et al. 2021"

# table6 is a single ~104,720-row table spanning ~700 different California
# Legacy Survey targets, keyed by this "CPS" identifier column -- a bare HD
# number, or a lowercase "gl<n>"/"hip<n>" catalog id with no space. Without
# filtering by it, every host redirected here would get the *entire*
# survey's RVs, not just its own.
ROSENTHAL_ID_COLUMN = "CPS"

HD_RE = re.compile(r'^HD\s*(\d+)', re.IGNORECASE)
GJ_RE = re.compile(r'^G[JL]\s*(\d+)', re.IGNORECASE)
HIP_RE = re.compile(r'^HIP\s*(\d+)', re.IGNORECASE)

# Column names some VizieR RV tables use to identify which star each row
# belongs to, for tables that -- like Rosenthal+2021's CLS table6 above --
# combine many different targets' RVs into one table but (unlike table6)
# weren't specifically anticipated here. Checked in download_vizier_table
# whenever a fetched table isn't already known to need the CLS-specific
# CPS-column filtering above.
GENERIC_MULTI_TARGET_ID_COLUMNS = ["Name", "Star", "Target", "Object", "Source"]


def normalize_identifier(value):
    return re.sub(r'[^a-z0-9]', '', str(value).lower())


def host_identifier_candidates(hostname, hd_name=None, hip_name=None):
    """Normalized strings that might identify `hostname` in a VizieR
    table's per-row star-name column: its own name, plus its HD/HIP
    designations both with and without the catalog prefix (one survey's
    table might list "47 UMa", another "HD 95128", another a bare
    "95128") -- for filtering tables that combine many stars' RVs into
    one table down to just this host's rows.
    """
    candidates = {normalize_identifier(hostname)}
    for raw in (hd_name, hip_name):
        if pd.notna(raw):
            text = str(raw)
            candidates.add(normalize_identifier(text))
            m = re.search(r'(\d+)', text)
            if m:
                candidates.add(m.group(1).lstrip("0") or "0")
    return candidates


def detect_multi_target_column(colnames):
    for c in GENERIC_MULTI_TARGET_ID_COLUMNS:
        if c in colnames:
            return c
    return None


def host_to_cps_id(hostname):
    """Best-effort guess at this host's identifier in Rosenthal et al.
    2021's CLS table. Returns None if the hostname doesn't match one of
    the catalog-number patterns CPS uses for its bare/lowercase ids --
    CPS also has BD/Giclas/LHS designations and lettered multiple-star
    components we don't attempt to cover, so callers must treat None (or
    a guess that matches zero rows) as "could not isolate this host" and
    not fall back to handing over the unfiltered, all-stars table.
    """
    m = HD_RE.match(hostname)
    if m:
        return m.group(1)
    m = GJ_RE.match(hostname)
    if m:
        return f"gl{m.group(1)}"
    m = HIP_RE.match(hostname)
    if m:
        return f"hip{m.group(1)}"
    return None

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
            "hd_name": r.get("hd_name"), "hip_name": r.get("hip_name"),
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
            cited = {i.strip() for i in str(r.get("rv_instruments") or "").split(",") if i.strip()}
            dace_instruments = cited & DACE_PIPELINE_INSTRUMENTS
            clue = (f"names {', '.join(sorted(dace_instruments))} (DACE-pipeline instrument)"
                    if dace_instruments else "mentions DACE")
            existing.append({
                "type": "dace", "bibcode": r["bibcode"], "display": r["display"], "dace_clue": clue,
            })

    # Always try DACE directly for every Msini host, regardless of whether
    # a VizieR/CLS source was already resolved or an ads-search clue named
    # DACE explicitly. DACE's public API needs no textual clue to query --
    # it's just a target-name lookup -- and an independent RV dataset
    # there is valuable even when another source already exists (e.g. a
    # VizieR table from one paper, DACE RVs combining several others).
    for hostname in msini_hosts:
        existing = sources.setdefault(hostname, [])
        if not any(e["type"] == "dace" for e in existing):
            existing.append({
                "type": "dace", "bibcode": None, "display": None,
                "dace_clue": "queried directly (DACE is checked for every host with a reported Msini mass)",
            })

    # Multi-target survey tables (currently just Rosenthal et al. 2021's
    # CLS table6) need per-host row filtering no matter how a host ended
    # up pointed at them -- whether it directly cites Rosenthal+2021 itself
    # (the "direct" loop above) or got redirected here because some other
    # paper just said "California Legacy Survey" (the "mentions_cls" loop).
    # Annotated centrally so both paths get the same filtering treatment.
    for hostname, host_sources in sources.items():
        for s in host_sources:
            if s["type"] == "vizier" and s["table_id"] == ROSENTHAL_TABLE:
                s["multi_target"] = True
                s["target_filter"] = (ROSENTHAL_ID_COLUMN, host_to_cps_id(hostname))

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

    filter_note = None
    if source.get("multi_target"):
        n_total = len(table)
        print(f"  ! WARNING: {source['table_id']} is a survey-wide table covering "
              f"many different stars, not just {hostname} ({n_total} rows total) -- "
              "filtering to this host's rows only")
        col, target_id = source.get("target_filter") or (None, None)
        if col and target_id and col in table.colnames:
            table = table[[str(row_val).strip().lower() == target_id.lower() for row_val in table[col]]]
            print(f"    kept {len(table)}/{n_total} rows matching {col}={target_id}")
            filter_note = f"Filtered from the full {n_total}-row survey table to rows where {col}={target_id}"
        if not col or not target_id or len(table) == 0:
            print(f"    ! could not isolate {hostname}'s rows in this multi-target table "
                  "(no reliable id match) -- skipping rather than saving every star's RVs under this host")
            return None
    else:
        # Not a table we already know needs row filtering (e.g. Rosenthal's
        # CLS table6 above) -- but some other survey papers *also* publish
        # one VizieR table covering many different stars (e.g. a per-row
        # "Name" column listing several different HD numbers), and nothing
        # upstream can know that without actually looking at the fetched
        # data. Check here so those don't silently dump every star's RVs
        # under this one host's directory.
        id_col = detect_multi_target_column(table.colnames)
        if id_col:
            distinct = {str(v).strip() for v in table[id_col] if str(v).strip()}
            if len(distinct) > 1:
                n_total = len(table)
                print(f"  ! WARNING: {source['table_id']} combines RV data for {len(distinct)} different "
                      f"stars in one table, not just {hostname} ({n_total} rows total) -- "
                      "filtering to this host's rows only")
                candidates = host_identifier_candidates(hostname, source.get("hd_name"), source.get("hip_name"))
                table = table[[normalize_identifier(v) in candidates for v in table[id_col]]]
                if len(table) == 0:
                    print(f"    ! could not isolate {hostname}'s rows in this multi-target table "
                          f"(no value in column {id_col} matched this host's name/HD/HIP designation) -- "
                          "skipping rather than saving every star's RVs under this host")
                    return None
                print(f"    kept {len(table)}/{n_total} rows matching {id_col}")
                filter_note = f"Filtered from the full {n_total}-row table (column {id_col}) to this host's rows only"

    lines = [
        f"# Host: {hostname}",
        f"# Source: VizieR catalog {source['catalog']} (table {source['table_id']})",
        f"# Reference: {source['display']} ({source['bibcode']})",
        f"# ADS: https://ui.adsabs.harvard.edu/abs/{source['bibcode']}/abstract",
        f"# Retrieved via astroquery.vizier on {date.today().isoformat()}",
    ]
    if filter_note:
        lines.append(f"# NOTE: {filter_note}")
    lines.append("# Columns:")
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

    trigger_line = (f"# Triggered by: {source['display']} ({source['bibcode']}) -- {source['dace_clue']}"
                     if source.get("display") else f"# Triggered by: {source['dace_clue']}")
    lines = [
        f"# Host: {hostname}",
        "# Source: DACE archive (public mode, dace_query.spectroscopy.Spectroscopy.get_timeseries)",
        f"# Target queried: {target}",
        trigger_line,
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
