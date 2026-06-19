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

# Always construct the results DataFrame with this fixed column set --
# otherwise an empty `rows` list (e.g. every --host value failed to match
# the catalog) produces a zero-column DataFrame, and the "looks_like_rv"
# lookup just below crashes with a KeyError instead of reporting 0 results.
RESULT_COLUMNS = [
    "hostname", "refstr", "bibcode", "display", "href",
    "vizier_catalog", "vizier_table", "table_title", "looks_like_rv",
    "rv_provenance", "note", "hd_name", "hip_name",
]

VIZIER_URL = "https://vizier.cds.unistra.fr/viz-bin/VizieR-3"
HEADERS = {"User-Agent": "vizier_rv_search.py (research script; contact lweiss4@nd.edu)"}

# CDS ReadMe files live at this fixed path for every "J/" journal-table
# catalog. We use them to check what a candidate table's columns actually
# are, because the keyword match against a VizieR table_title/catalog
# description alone is too easy to fool: a joint photometry+RV paper's
# catalog description (or a table title like "...joint photo-dynamical+RV
# fit") mentions "RV" even for tables that are pure photometry or fitted
# planet parameters, with no RV time series in them at all.
README_URL_TMPL = "https://cdsarc.cds.unistra.fr/ftp/{catalog}/ReadMe"

# Matches one column-definition line inside a ReadMe "Byte-by-byte
# Description of file: ..." block, e.g.:
#   "  18- 25  F7.2  m/s     RV        Radial velocity"
# Header/separator/continuation lines don't have the leading
# bytes-format-units-label structure and so don't match.
README_COLUMN_LINE_RE = re.compile(
    r'^\s*[\d\-\s]+\s+\S+\s+(?P<units>\S+)\s+(?P<label>\S+)\s+(?P<expl>.*)$'
)
README_RV_LABEL_RE = re.compile(r'^e?_?(RV|HRV|Vrad|VR)\d*$', re.IGNORECASE)
README_RV_UNIT_RE = re.compile(r'^k?c?m[./]s(-1)?$', re.IGNORECASE)
# Deliberately narrower than RV_KEYWORD_RE: a bare "RV" is exactly the
# substring that produces false positives in table titles (see above), so
# at the column-explanation level we only trust the unambiguous phrase.
README_RV_EXPLANATION_RE = re.compile(r'radial\s*veloc', re.IGNORECASE)


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


def fetch_readme(catalog_id, session, sleep=1.0):
    """Download a CDS ReadMe file for a "J/" journal-table catalog. Returns
    the raw text, or None if it couldn't be fetched (no ReadMe at that path,
    network error, etc.) -- callers should fall back to the title/
    description keyword heuristic in that case."""
    url = README_URL_TMPL.format(catalog=catalog_id)
    try:
        resp = session.get(url, headers=HEADERS, timeout=30)
        time.sleep(sleep)
        if resp.status_code != 200:
            return None
        return resp.text
    except requests.RequestException:
        return None


def parse_readme_rv_tables(text):
    """Parse a CDS ReadMe's "Byte-by-byte Description of file: ..." blocks
    and return {table_stem (lowercase, no extension): has_rv_column}.

    A file counts as having an RV column if some line in its block either
    (a) has an Explanations field matching "radial velocity" outright, or
    (b) has a Label like "RV"/"e_RV"/"Vrad" paired with a velocity Unit
    (km/s, m/s, ...). Both are deliberately narrower than the title/
    description keyword match -- they look at what the column actually is,
    not just whether "RV" appears somewhere nearby in prose.
    """
    blocks = re.split(r'\nByte-by-byte Description of file:\s*', text)
    results = {}
    for block in blocks[1:]:
        m = re.match(r'(.*?)\n-{10,}\n(?:.*?\n-{10,}\n)?(.*?)\n-{10,}', block, re.DOTALL)
        if not m:
            continue
        files_part, body = m.group(1), m.group(2)
        filenames = re.findall(r'(\S+)\.\w+', files_part)
        has_rv = False
        for line in body.splitlines():
            col = README_COLUMN_LINE_RE.match(line)
            if not col:
                continue
            if README_RV_EXPLANATION_RE.search(col.group("expl")):
                has_rv = True
                break
            if README_RV_LABEL_RE.match(col.group("label")) and README_RV_UNIT_RE.match(col.group("units")):
                has_rv = True
                break
        for stem in filenames:
            stem = stem.lower()
            results[stem] = results.get(stem, False) or has_rv
    return results


def get_readme_rv_tables(catalog_id, session, cache, sleep=1.0, verbose=False):
    """Cached wrapper around fetch_readme + parse_readme_rv_tables. `cache`
    is a dict (persisted alongside the bibcode cache) keyed by catalog_id,
    storing the parsed {stem: has_rv_column} map, or None if the ReadMe
    couldn't be fetched/parsed for that catalog."""
    if catalog_id in cache:
        return cache[catalog_id]
    if verbose:
        print(f"  fetching ReadMe for {catalog_id}")
    text = fetch_readme(catalog_id, session, sleep=sleep)
    parsed = parse_readme_rv_tables(text) if text else None
    cache[catalog_id] = parsed
    return parsed


def normalize_host_slug(name):
    return re.sub(r'[^a-z0-9]', '', name.lower())


HD_NUM_RE = re.compile(r'\bHD\s*0*(\d+)\b', re.IGNORECASE)
HIP_NUM_RE = re.compile(r'\bHIP\s*0*(\d+)\b', re.IGNORECASE)


def extract_catalog_numbers(text):
    """Pull out any HD/HIP catalog numbers mentioned in `text` (leading
    zeros stripped), e.g. {"HD": {"95128"}, "HIP": set()} for "HD 095128"."""
    text = text or ""
    return {"HD": set(HD_NUM_RE.findall(text)), "HIP": set(HIP_NUM_RE.findall(text))}


# Common stellar-designation prefixes that show up in VizieR table titles
# for *one specific star*, distinct from a survey-wide table that just
# describes its general content (e.g. "Stellar properties and RV data
# summary" -- no specific star named at all). Curated rather than a fully
# generic "letters+digits" regex, which would also match ordinary phrases
# like "Table 4" or "63 Kepler planet-hosting stars" (number-before-letters,
# which this deliberately excludes -- see DESIGNATION_TOKEN_RE).
DESIGNATION_PREFIXES = [
    "HD", "HIP", "GJ", "GL", "LHS", "LTT", "LP", "BD", "TOI", "KOI", "KIC",
    "TIC", "EPIC", "K2", "WASP", "HAT-P", "HATS", "KELT", "XO", "CoRoT",
    "NGTS", "TrES", "Kepler",
]
DESIGNATION_TOKEN_RE = re.compile(
    r'\b(' + '|'.join(re.escape(p) for p in DESIGNATION_PREFIXES) + r')[\s\-]*0*(\d+)\b',
    re.IGNORECASE,
)


def extract_designation_tokens(text):
    """Pull out every recognized stellar designation in `text` as a set of
    (prefix, number) tuples, e.g. {("gj", "1252")} for "GJ 1252". Used to
    tell whether a table's title names a *specific* star at all -- and if
    so, which one -- independent of whether that star happens to be a
    known host in this run's input catalog (see table_belongs_to_other_host).
    """
    text = text or ""
    return {(p.lower().replace("-", ""), n.lstrip("0") or "0") for p, n in DESIGNATION_TOKEN_RE.findall(text)}


def host_designation_tokens(hostname, hd_name=None, hip_name=None):
    """A host's own designation tokens (see extract_designation_tokens):
    parsed from its common NEA name (which for many hosts -- "TOI-1235",
    "Kepler-126", "GJ 1252" -- already *is* a recognized catalog
    designation) plus its hd_name/hip_name columns."""
    tokens = extract_designation_tokens(hostname)
    for raw in (hd_name, hip_name):
        if pd.notna(raw):
            tokens |= extract_designation_tokens(str(raw))
    return tokens


def build_host_designations(df):
    """hostname -> {"HD": {...}, "HIP": {...}} parsed from the NEA catalog's
    hd_name/hip_name columns. Many multi-target VizieR catalogs (e.g. an
    ELODIE/HARPS survey paper covering a dozen unrelated stars) title each
    per-target table with a bare HD/HIP number rather than the star's
    common name -- "HD 95128" rather than "47 UMa" -- so a pure hostname
    text match against the table title can't tell those tables apart.
    """
    designations = {}
    for hostname, group in df.groupby("hostname"):
        ids = {"HD": set(), "HIP": set()}
        for col, key in (("hd_name", "HD"), ("hip_name", "HIP")):
            if col in group.columns:
                for raw in group[col].dropna().unique():
                    ids[key] |= extract_catalog_numbers(str(raw))[key]
        designations[hostname] = ids
    return designations


def host_primary_names(df):
    """hostname -> (hd_name, hip_name) raw display strings (e.g. "HD 95128"),
    propagated through to the download stage so it can filter VizieR
    tables that combine many different stars' RVs into one table (see
    download.detect_multi_target_column) -- not just the Rosenthal+2021
    CLS table, which already had its own hardcoded filtering."""
    names = {}
    for hostname, group in df.groupby("hostname"):
        hd = next(iter(group["hd_name"].dropna().unique()), None) if "hd_name" in group else None
        hip = next(iter(group["hip_name"].dropna().unique()), None) if "hip_name" in group else None
        names[hostname] = (hd, hip)
    return names


def table_belongs_to_other_host(table, hostname, other_hostnames, host_designations=None,
                                 hd_name=None, hip_name=None):
    """True if a table's id/title textually identifies it as belonging to
    a *different* star than `hostname`. Some VizieR catalogs cover
    multiple distinct targets in one paper, with separate per-target
    tables (e.g. Lillo-Box et al. 2020's J/A+A/640/A48 has "k2-32rv" and
    "k2-233rv" for two unrelated multi-planet systems that both happen to
    cite that paper) -- without this check, every table under a shared
    bibcode gets attached to every host citing it.

    Two checks, in order:
      1. HD/HIP catalog numbers (via `host_designations`, see
         build_host_designations) and a plain hostname text match --
         needed because some surveys title each per-target table with a
         bare HD/HIP number (e.g. "HD 95128" for 47 UMa) that won't
         contain either host's common name as a literal substring. Only
         excludes when the title positively names a *different host
         that's also in this run's input catalog* (`other_hostnames`).
      2. A generic stellar-designation token (see extract_designation_tokens)
         that doesn't match this host's own tokens at all -- catches a
         table titled for some *other* star the title positively names
         (e.g. "GJ 1252 radial velocity curve"), even when that other star
         isn't itself a host in this run to check against in step 1 (e.g.
         a survey paper covering several stars, only one of which is in
         your NEA export).

    Returns False (not someone else's) whenever there's no hint either
    way, so single-target catalogs -- the common case -- are unaffected;
    a survey-wide table that doesn't name any specific star in its title
    (e.g. "Radial velocities from Keck-HIRES of 63 Kepler planet-hosting
    stars") also isn't excluded here -- isolating its rows to one host is
    download.py's job (see detect_multi_target_column), not this title-only
    check's.
    """
    raw_text = f"{table['table_id']} {table['table_title']}"
    text = normalize_host_slug(raw_text)
    own_slug = normalize_host_slug(hostname)
    host_designations = host_designations or {}

    table_nums = extract_catalog_numbers(raw_text)
    own_ids = host_designations.get(hostname, {"HD": set(), "HIP": set()})
    if (table_nums["HD"] & own_ids["HD"]) or (table_nums["HIP"] & own_ids["HIP"]):
        return False
    if own_slug and own_slug in text:
        return False

    for other in other_hostnames:
        if other == hostname:
            continue
        other_ids = host_designations.get(other, {"HD": set(), "HIP": set()})
        if (table_nums["HD"] & other_ids["HD"]) or (table_nums["HIP"] & other_ids["HIP"]):
            return True
        other_slug = normalize_host_slug(other)
        # Skip very short slugs (e.g. "K2-3" -> "k23") -- too easy to
        # match by coincidence inside unrelated text.
        if len(other_slug) >= 4 and other_slug in text:
            return True

    table_tokens = extract_designation_tokens(raw_text)
    if table_tokens and not (table_tokens & host_designation_tokens(hostname, hd_name, hip_name)):
        return True
    return False


def is_journal_table_catalog(catalog_id):
    """True for VizieR catalogs published from a specific paper's tables
    (the 'J/' class). VizieR's Roman-numeral classes (I=astrometry,
    II=photometry, III=spectroscopy, ...) are its own large homogeneous
    surveys (Gaia, Hipparcos, 2MASS, ...); even though they can be linked
    to a paper's bibcode, they're not that paper's own per-target RV
    monitoring data and aren't useful for precision differential RVs.
    """
    return bool(catalog_id) and catalog_id.startswith("J/")


def unresolved_provenance_refs(results_df):
    """Provenance references -- (hostname, bibcode) pairs the NEA credits
    with an Msini mass or RV semi-amplitude -- that don't themselves have a
    VizieR table flagged as RV, deduped to one row per (hostname, bibcode).

    Deliberately keyed on (hostname, bibcode), not hostname alone: a host
    can cite many papers, and a *different*, non-provenance reference for
    that same host having an RV-looking VizieR table (e.g. an earlier
    discovery paper using a different instrument) doesn't mean the
    provenance paper's own RV data was found.
    """
    resolved = set(
        results_df.loc[results_df["looks_like_rv"], ["hostname", "bibcode"]].itertuples(index=False, name=None)
    )
    is_unresolved = ~results_df.apply(lambda r: (r["hostname"], r["bibcode"]) in resolved, axis=1)
    missing = results_df[results_df["rv_provenance"] & is_unresolved]
    return missing.drop_duplicates(subset=["hostname", "bibcode"])


def load_cache(path):
    if path.exists():
        return json.loads(path.read_text())
    return {}


def save_cache(path, cache):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cache, indent=1))


def classify_bibcode_tables(tables, session, readme_cache, readme_cache_path, sleep, verbose):
    """Run the catalog-class + RV-keyword + ReadMe checks for every table
    returned for one bibcode, independent of which host is asking. Returns
    a list of {"table": t, "rv_match": bool, "note": str} dicts.

    Done once per bibcode (host-independent, cached by search_vizier in
    bibcode_classified) since the result doesn't depend on which host is
    currently being processed -- avoids re-running the ReadMe check for
    every host that happens to cite the same bibcode.
    """
    classified = []
    for t in tables:
        rv_match = is_journal_table_catalog(t["catalog"]) and (
            looks_like_rv(t["table_title"]) or looks_like_rv(t["catalog_description"])
        )
        note = ""
        # A title/description keyword match is easy to fool -- e.g. a
        # joint photometry+RV paper's table titled "...joint
        # photo-dynamical+RV fit" matches on the bare "RV" even
        # though that table holds fitted planet parameters, not RV
        # data. Cross-check against the catalog's actual ReadMe
        # column descriptions before trusting it.
        if rv_match:
            readme_tables = get_readme_rv_tables(
                t["catalog"], session, readme_cache, sleep=sleep, verbose=verbose
            )
            save_cache(readme_cache_path, readme_cache)
            stem = t["table_id"].rsplit("/", 1)[-1].lower()
            if readme_tables is None:
                note = "ReadMe unavailable; kept keyword match unverified"
            elif stem not in readme_tables:
                note = "table not found in ReadMe; kept keyword match unverified"
            elif not readme_tables[stem]:
                rv_match = False
                note = "ReadMe shows no RV column in this table; keyword match overridden"
        classified.append({"table": t, "rv_match": rv_match, "note": note})
    return classified


def search_vizier(input_csv, output_csv=None, cache_path="cache/vizier_bibcode_cache.json",
                   host=None, sleep=1.0, verbose=False,
                   readme_cache_path="cache/vizier_readme_cache.json"):
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
    readme_cache_path -- JSON cache of catalog -> {table_stem: has_rv_column},
                    parsed from each catalog's CDS ReadMe. Used to confirm a
                    title/description keyword match against the table's
                    actual columns (see get_readme_rv_tables).
    """
    df = read_nea_csv(input_csv)
    # Used for table_belongs_to_other_host below -- the full catalog's
    # hostnames (and their HD/HIP designations), not just the (possibly
    # --host-filtered) subset being processed, so a single-host run still
    # benefits from knowing about other targets that might share a
    # multi-target VizieR catalog.
    all_hostnames = set(df["hostname"].unique())
    all_designations = build_host_designations(df)
    all_names = host_primary_names(df)
    if host:
        hosts = [h.strip() for h in host.split(",")]
        df = df[df["hostname"].isin(hosts)]
        unmatched = [h for h in hosts if h not in all_hostnames]
        if unmatched:
            print(f"  ! WARNING: host(s) not found in {input_csv} (check spelling/spacing -- "
                  f"hostnames must match the catalog exactly): {', '.join(unmatched)}")

    host_refs = aggregate_references_by_host(df)
    cache_path = Path(cache_path)
    cache = load_cache(cache_path)
    readme_cache_path = Path(readme_cache_path)
    readme_cache = load_cache(readme_cache_path)

    session = requests.Session()
    rows = []
    failed_bibcodes = set()
    bibcode_classified = {}  # bibcode -> classify_bibcode_tables(...) result, reused across hosts sharing it
    for hostname, refs in host_refs.items():
        hd_name, hip_name = all_names.get(hostname, (None, None))
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
                    "hd_name": hd_name, "hip_name": hip_name,
                })
                continue

            if bibcode in cache:
                tables = cache[bibcode]
                query_failed = False
            elif bibcode in failed_bibcodes:
                # Already failed once this run -- don't hammer an
                # unresponsive endpoint again for a reference shared by
                # multiple hosts; just treat it as empty for this row.
                tables = []
                query_failed = True
            else:
                if verbose:
                    print(f"querying VizieR for {bibcode} ({ref['display']})")
                try:
                    tables = query_vizier_bibcode(bibcode, session, sleep=sleep)
                    query_failed = False
                except requests.RequestException as exc:
                    print(f"  failed: {exc}")
                    tables = []
                    query_failed = True
                    failed_bibcodes.add(bibcode)
                if not query_failed:
                    # Only persist genuine "VizieR has no catalog for this
                    # bibcode" results. Caching a transient network failure
                    # as [] would be indistinguishable from that on the next
                    # run and silently poison it forever.
                    cache[bibcode] = tables
                    save_cache(cache_path, cache)

            if not tables:
                rows.append({
                    "hostname": hostname, "refstr": ref["refstr"], "bibcode": bibcode,
                    "display": ref["display"], "href": ref["href"],
                    "vizier_catalog": None, "vizier_table": None,
                    "table_title": None, "looks_like_rv": False,
                    "rv_provenance": rv_provenance,
                    "note": "VizieR query failed; will retry next run" if query_failed else "no VizieR catalog found",
                    "hd_name": hd_name, "hip_name": hip_name,
                })
            if bibcode not in bibcode_classified:
                bibcode_classified[bibcode] = classify_bibcode_tables(
                    tables, session, readme_cache, readme_cache_path, sleep, verbose
                )
            for c in bibcode_classified[bibcode]:
                t, rv_match, note = c["table"], c["rv_match"], c["note"]
                if rv_match and table_belongs_to_other_host(t, hostname, all_hostnames, all_designations,
                                                              hd_name, hip_name):
                    rv_match = False
                    note = "table id/title names a different host -- shared multi-target catalog, not this host's data"
                rows.append({
                    "hostname": hostname, "refstr": ref["refstr"], "bibcode": bibcode,
                    "display": ref["display"], "href": ref["href"],
                    "vizier_catalog": t["catalog"], "vizier_table": t["table_id"],
                    "table_title": t["table_title"],
                    "looks_like_rv": rv_match,
                    "rv_provenance": rv_provenance, "note": note,
                    "hd_name": hd_name, "hip_name": hip_name,
                })

    out = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    if output_csv:
        out.to_csv(output_csv, index=False)
        print(f"wrote {len(out)} rows to {output_csv}")

    if out.empty:
        print("0/0 hosts have at least one candidate RV table on VizieR (no host matched the catalog -- "
              "see the WARNING above, if any)")
        return out

    rv_hosts = set(out.loc[out["looks_like_rv"], "hostname"].unique())
    print(f"{len(rv_hosts)}/{out['hostname'].nunique()} hosts have at least one candidate RV table on VizieR")

    # Hosts where the NEA attributes an Msini mass or RV semi-amplitude to a
    # reference, but that *specific* reference has no VizieR table flagged
    # as RV -- the RV data probably exists in the literature, just isn't
    # mirrored on VizieR (or its table title didn't match our keyword
    # heuristic). Checked per (hostname, bibcode), not just per hostname:
    # a host can cite many papers, and one of them (e.g. an earlier,
    # unrelated discovery paper) having an RV table on VizieR doesn't mean
    # the paper actually credited with the reported mass does too.
    missing = unresolved_provenance_refs(out)
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
    parser.add_argument("--readme-cache", default="cache/vizier_readme_cache.json")
    parser.add_argument("--host", default=None,
                         help="Only process this hostname, or a comma-separated list of hostnames")
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds between VizieR requests")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    search_vizier(args.input, args.output, args.cache, args.host, args.sleep, args.verbose,
                  readme_cache_path=args.readme_cache)


if __name__ == "__main__":
    main()
