"""For each hostname in a NASA Exoplanet Archive multis CSV, aggregate the unique
literature references cited for that system and check VizieR (by ADS bibcode)
for an associated table that looks like a radial-velocity dataset.

Usage (after `pip install -e .`):
    rv-vizier-search --input NEA_Multis.csv --output vizier_rv_results.csv
    rv-vizier-search --host K2-138 --verbose
"""
import argparse
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

REF_ANCHOR_RE = re.compile(
    r'<a\s+refstr=(?P<refstr>\S+)\s+href=(?P<href>\S+)\s+target=\S+>(?P<text>[^<]*)</a>'
)
ADS_BIBCODE_RE = re.compile(r'adsabs\.harvard\.edu/abs/([^/]+)/abstract')
RV_KEYWORD_RE = re.compile(r'radial\s*veloc|spectroscopic\s*orbit|\bRVs?\b', re.IGNORECASE)

# NEA mass-provenance values that mean the mass was derived from an RV orbit
# (Msini), as opposed to e.g. a mass-radius relationship or TTVs.
MSINI_PROVENANCE = {"Msini", "Msin(i)/sin(i)"}

VIZIER_URL = "https://vizier.cds.unistra.fr/viz-bin/VizieR-3"
HEADERS = {"User-Agent": "vizier_rv_search.py (research script; contact lweiss4@nd.edu)"}


def read_nea_csv(path):
    """Read a NASA Exoplanet Archive CSV export, skipping any leading '#'
    metadata block. Doesn't use pandas' comment= kwarg because reference
    text can contain HTML numeric entities like '&#x131;', whose literal
    '#' would otherwise truncate the row mid-field.
    """
    with open(path) as f:
        skip = 0
        for line in f:
            if not line.startswith("#"):
                break
            skip += 1
    return pd.read_csv(path, skiprows=skip, low_memory=False)


def reference_columns(df):
    return [c for c in df.columns if c.endswith("_reflink")] + ["disc_refname"]


def parse_reference(raw):
    """Parse one NEA reference anchor string into refstr/bibcode/display text."""
    m = REF_ANCHOR_RE.match(raw.strip())
    if not m:
        return None
    href = m.group("href")
    bib_m = ADS_BIBCODE_RE.search(href)
    return {
        "refstr": m.group("refstr"),
        "bibcode": bib_m.group(1) if bib_m else None,
        "display": m.group("text").strip(),
        "href": href,
    }


def aggregate_references_by_host(df):
    """Map hostname -> list of unique parsed references (deduped by bibcode/refstr).

    Each reference also carries an `rv_provenance` flag: True if the NEA
    itself attributes an RV-derived quantity (an Msini mass, or a measured
    RV semi-amplitude) for some planet in this system to that reference.
    That's independent of whether VizieR happens to host a matching table --
    it's a clue that RV data for this host exists in the literature even if
    VizieR search comes up empty.
    """
    refcols = reference_columns(df)
    host_refs = {}
    for hostname, group in df.groupby("hostname"):
        seen = {}
        for col in refcols:
            for raw in group[col].dropna():
                ref = parse_reference(raw)
                if ref is None:
                    continue
                key = ref["bibcode"] or ref["refstr"]
                seen.setdefault(key, ref)
                seen[key].setdefault("rv_provenance", False)

        for _, row in group.iterrows():
            provenance_cols = []
            if row.get("pl_bmassprov") in MSINI_PROVENANCE:
                provenance_cols += ["pl_bmasse_reflink", "pl_bmassj_reflink"]
            if pd.notna(row.get("pl_rvamp")):
                provenance_cols.append("pl_rvamp_reflink")
            for col in provenance_cols:
                raw = row.get(col)
                if pd.isna(raw):
                    continue
                ref = parse_reference(raw)
                if ref is None:
                    continue
                key = ref["bibcode"] or ref["refstr"]
                if key in seen:
                    seen[key]["rv_provenance"] = True

        host_refs[hostname] = list(seen.values())
    return host_refs


def clean_table_title(raw_title, table_id):
    title = raw_title
    prefix = f"{table_id}: "
    if title.startswith(prefix):
        title = title[len(prefix):]
    title = re.split(r"\s*catContent", title)[0]
    return title.strip()


def query_vizier_bibcode(bibcode, session, sleep=1.0):
    """Look up VizieR catalogs/tables associated with an ADS bibcode."""
    params = f"-source={requests.utils.quote(bibcode)}&-meta.all&-out.max=50"
    resp = session.get(f"{VIZIER_URL}?{params}", headers=HEADERS, timeout=30)
    time.sleep(sleep)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "html.parser")
    select = soup.find("select", attrs={"name": "//tables"})
    if select is None:
        return []

    tables = []
    current_catalog, current_catalog_desc = None, None
    for opt in select.find_all("option"):
        value = opt.get("value")
        if value is None:
            current_catalog = opt.text.strip()
            current_catalog_desc = opt.get("title", "")
        else:
            tables.append({
                "catalog": current_catalog,
                "catalog_description": current_catalog_desc,
                "table_id": value,
                "table_title": clean_table_title(opt.get("title", ""), value),
            })
    return tables


def looks_like_rv(text):
    return bool(text) and bool(RV_KEYWORD_RE.search(text))


def is_journal_table_catalog(catalog_id):
    """True for VizieR catalogs published from a specific paper's tables
    (the 'J/' class). VizieR's Roman-numeral classes (I=astrometry,
    II=photometry, III=spectroscopy, ...) are its own large homogeneous
    surveys (Gaia, Hipparcos, 2MASS, ...); even though they can be linked
    to a paper's bibcode, they're not that paper's own per-target RV
    monitoring data and aren't useful for precision differential RVs.
    """
    return bool(catalog_id) and catalog_id.startswith("J/")


def load_cache(path):
    if path.exists():
        return json.loads(path.read_text())
    return {}


def save_cache(path, cache):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, indent=1))


def search_vizier(input_csv, output_csv=None, cache_path="cache/vizier_bibcode_cache.json",
                   host=None, sleep=1.0, verbose=False):
    """Run the full VizieR RV-table search and return the results as a DataFrame.

    This is the importable entry point (main() is just a thin CLI wrapper
    around this, installed as the `rv-vizier-search` command). Parameters
    mirror the CLI flags:

    input_csv    -- path to a NASA Exoplanet Archive CSV export.
    output_csv   -- if given, write the results table here.
    cache_path   -- JSON cache of bibcode -> VizieR tables, reused across runs.
    host         -- a hostname, or comma-separated list of hostnames, to
                    restrict the search to. None processes every host in the file.
    sleep        -- seconds to wait between live VizieR requests (politeness).
    verbose      -- print progress as each new bibcode is queried.
    """
    df = read_nea_csv(input_csv)
    if host:
        hosts = [h.strip() for h in host.split(",")]
        df = df[df["hostname"].isin(hosts)]

    host_refs = aggregate_references_by_host(df)
    cache_path = Path(cache_path)
    cache = load_cache(cache_path)

    session = requests.Session()
    rows = []
    for hostname, refs in host_refs.items():
        for ref in refs:
            bibcode = ref["bibcode"]
            rv_provenance = ref.get("rv_provenance", False)
            if bibcode is None:
                rows.append({
                    "hostname": hostname, "refstr": ref["refstr"], "bibcode": None,
                    "display": ref["display"], "href": ref["href"],
                    "vizier_catalog": None, "vizier_table": None,
                    "table_title": None, "looks_like_rv": False,
                    "rv_provenance": rv_provenance, "note": "no ADS bibcode",
                })
                continue

            if bibcode not in cache:
                if verbose:
                    print(f"querying VizieR for {bibcode} ({ref['display']})")
                try:
                    cache[bibcode] = query_vizier_bibcode(bibcode, session, sleep=sleep)
                except requests.RequestException as exc:
                    print(f"  failed: {exc}")
                    cache[bibcode] = []
                save_cache(cache_path, cache)

            tables = cache[bibcode]
            if not tables:
                rows.append({
                    "hostname": hostname, "refstr": ref["refstr"], "bibcode": bibcode,
                    "display": ref["display"], "href": ref["href"],
                    "vizier_catalog": None, "vizier_table": None,
                    "table_title": None, "looks_like_rv": False,
                    "rv_provenance": rv_provenance, "note": "no VizieR catalog found",
                })
            for t in tables:
                rv_match = is_journal_table_catalog(t["catalog"]) and (
                    looks_like_rv(t["table_title"]) or looks_like_rv(t["catalog_description"])
                )
                rows.append({
                    "hostname": hostname, "refstr": ref["refstr"], "bibcode": bibcode,
                    "display": ref["display"], "href": ref["href"],
                    "vizier_catalog": t["catalog"], "vizier_table": t["table_id"],
                    "table_title": t["table_title"],
                    "looks_like_rv": rv_match,
                    "rv_provenance": rv_provenance, "note": "",
                })

    out = pd.DataFrame(rows)
    if output_csv:
        out.to_csv(output_csv, index=False)
        print(f"wrote {len(out)} rows to {output_csv}")

    rv_hosts = set(out.loc[out["looks_like_rv"], "hostname"].unique())
    print(f"{len(rv_hosts)}/{out['hostname'].nunique()} hosts have at least one candidate RV table on VizieR")

    # Hosts where the NEA attributes an Msini mass or RV semi-amplitude to a
    # reference, but no VizieR table for that reference was flagged as RV --
    # the RV data probably exists in the literature, just isn't mirrored on
    # VizieR (or its table title didn't match our keyword heuristic).
    missing = out[out["rv_provenance"] & ~out["hostname"].isin(rv_hosts)]
    missing = missing.drop_duplicates(subset=["hostname", "bibcode"])
    if len(missing):
        print(f"\n{missing['hostname'].nunique()} hosts have an RV-derived mass/amplitude in the NEA "
              "but no matching VizieR table was found -- likely just not mirrored on VizieR:")
        for _, r in missing.sort_values("hostname").iterrows():
            print(f"  {r['hostname']}: {r['display']} ({r['bibcode']}) -- {r['note'] or 'catalog found, no RV-looking table'}")

    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="NEA_Multis.csv")
    parser.add_argument("--output", default="vizier_rv_results.csv")
    parser.add_argument("--cache", default="cache/vizier_bibcode_cache.json")
    parser.add_argument("--host", default=None,
                         help="Only process this hostname, or a comma-separated list of hostnames")
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds between VizieR requests")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    search_vizier(args.input, args.output, args.cache, args.host, args.sleep, args.verbose)


if __name__ == "__main__":
    main()
