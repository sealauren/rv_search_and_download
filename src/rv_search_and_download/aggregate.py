"""For every host with at least one downloaded RV table
(downloaded_rv_tables/<host>/*.csv, written by download.py), combine those
tables into a single analysis-ready file:

  - Standardizes each source's time/RV/RV-error columns into BJD, RV (m/s),
    and RV error (m/s), correctly accounting for any BJD0 offset (e.g.
    "BJD-2400000") and unit (e.g. km/s). Every other original column is
    carried over as-is -- no attempt is made to match these across
    different data sets (e.g. one table's "FWHM" and another's "ccf_fwhm"),
    per the project's column-naming inconsistency across sources.
  - Median-subtracts the RV column within each (source table, instrument)
    group, to remove pipeline-to-pipeline systemic RV offsets that would
    otherwise choke a multi-instrument orbital fit. Grouping by instrument
    *within* a table matters because a single DACE download can combine
    many physically distinct instruments (e.g. CORAVEL ~-16500 m/s vs.
    HAMILTON ~0 m/s for the same star) -- subtracting one whole-table
    median would leave most of those instruments' offsets untouched.
  - Flags likely-duplicate RVs: same host, same instrument, timestamps from
    different source tables within DUPLICATE_WINDOW_SECONDS of each other
    (transitively grouped, not strict pairwise matching). Each such group
    gets a shared duplicate_group id; exactly one row per group is marked
    is_preferred=True per a configurable default-source policy (see
    load_config), and the rest are kept (is_preferred=False) for review.
  - Writes aggregated_rv_tables/<host>_rvs_aggregated.csv.

Source type (VizieR vs. DACE) and column metadata are read directly back
out of the commented header download.py already writes at the top of each
file -- no need to re-query VizieR/DACE.

Usage (after `pip install -e .`):
    rv-aggregate-rvs
    rv-aggregate-rvs --host "tau Cet"
    rv-aggregate-rvs --init-config local   # write a starter rv_aggregate_config.json here
"""
import argparse
import json
import re
from datetime import date
from pathlib import Path

import pandas as pd

from .ads_search import RV_INSTRUMENTS
from .download import sanitize

DUPLICATE_WINDOW_SECONDS = 10.0

DEFAULT_CONFIG = {
    # "prefer_source": always pick the row whose source_type is earliest in
    # preferred_source_order. "most_recent": pick the row with the latest
    # `year` (parsed from its bibcode, falling back to retrieval date).
    "duplicate_default": "prefer_source",
    "preferred_source_order": ["dace", "vizier"],
}

LOCAL_CONFIG_NAME = "rv_aggregate_config.json"
GLOBAL_CONFIG_PATH = Path.home() / ".config" / "rv_search_and_download" / "config.json"

# Lines download.py writes look like "#   <name> [<unit>] - <description>"
# (vizier) or "#   <name> - <description>" (dace) -- always 3+ spaces after
# "#", which distinguishes them from single-space-indented narrative lines
# like "# Host: ..." or "# Retrieved via ... on ...".
COLUMN_LINE_RE = re.compile(r'^#\s{2,}(\S+)\s*(?:\[([^\]]*)\])?\s*-\s*(.*)$')
SOURCE_LINE_RE = re.compile(r'^# Source: (\S+)')
ADS_LINE_RE = re.compile(r'^# ADS: https://ui\.adsabs\.harvard\.edu/abs/(.+)/abstract$')
RETRIEVED_RE = re.compile(r'on (\d{4})-\d{2}-\d{2}$')
BIBCODE_IN_PARENS_RE = re.compile(r'\((\d{4}\S*)\)')

TIME_NAME_RE = re.compile(r'^(b?h?jd|time|rjd)$', re.IGNORECASE)
TIME_DESC_RE = re.compile(r'\b(b?h?jd|julian date)\b', re.IGNORECASE)
OFFSET_RE = re.compile(r'[bhm]?jd\s*-\s*(\d+(?:\.\d+)?)', re.IGNORECASE)

RV_NAME_RE = re.compile(r'^rv$|radial.?vel', re.IGNORECASE)
RV_ERR_NAME_RE = re.compile(r'^e_?rv$|rv_?err', re.IGNORECASE)
RV_DESC_RE = re.compile(r'radial veloc', re.IGNORECASE)
ERROR_DESC_RE = re.compile(r'uncertain|error|sigma', re.IGNORECASE)

# Some VizieR RV tables -- unlike the typical title-only case
# detect_vizier_instrument() handles -- do carry a per-row column naming
# which telescope/instrument took each RV, either as literal text (e.g. a
# "Tel" column with values like "HET") or as a numeric code explained right
# in the column's own description (e.g. "Instr [-] - [1,2] Instrument (1:
# ELODIE, 2: CORALIE)"). Recognizing these gives duplicate detection a real
# per-row instrument instead of one fallback label for the whole table.
INSTRUMENT_COLUMN_NAME_RE = re.compile(r'^(tel|telescope|instr|instrument|spectro|spectrograph|facility)$', re.IGNORECASE)
INSTRUMENT_LEGEND_RE = re.compile(r'\((?:\s*\d+\s*:\s*[^,()]+,?\s*)+\)')
LEGEND_CODE_PAIR_RE = re.compile(r'(\d+)\s*:\s*([^,()]+)')

# A handful of VizieR tables identify the *telescope* rather than the
# instrument/spectrograph that took each RV (e.g. Wittenmyer et al. 2009's
# "Tel" column lists "HET", the Hobby-Eberly Telescope -- whose only RV
# spectrograph at the time was HRS). Mapped here so those rows still get a
# real, cross-source-matchable instrument name instead of the bare
# telescope label.
TELESCOPE_TO_INSTRUMENT = {
    "HET": "HRS",
}


def parse_instrument_legend(desc):
    """Parse a column description like '[1,2] Instrument (1: ELODIE, 2:
    CORALIE)' into {'1': 'ELODIE', '2': 'CORALIE'}, or None if the
    description isn't an instrument/telescope code legend."""
    if not re.search(r'instrument|spectrograph|telescope', desc, re.IGNORECASE):
        return None
    m = INSTRUMENT_LEGEND_RE.search(desc)
    if not m:
        return None
    pairs = LEGEND_CODE_PAIR_RE.findall(m.group(0))
    return {code: label.strip() for code, label in pairs} if pairs else None


def find_instrument_row_source(columns):
    """Look for a per-row column identifying which telescope/instrument
    took each RV. Returns (column_name, legend_or_None), or (None, None)
    if no such column is named in the table's header metadata."""
    for name, info in columns.items():
        if INSTRUMENT_COLUMN_NAME_RE.match(name):
            return name, parse_instrument_legend(info["desc"])
    return None, None


def normalize_instrument_label(raw):
    raw = str(raw).strip()
    if raw in TELESCOPE_TO_INSTRUMENT:
        return TELESCOPE_TO_INSTRUMENT[raw]
    for inst in RV_INSTRUMENTS:
        if raw.lower() == inst.lower():
            return inst
    return raw


def load_config(start_dir="."):
    """Load duplicate-resolution settings: a local rv_aggregate_config.json
    in `start_dir` takes priority, then a global config under
    ~/.config/rv_search_and_download/, then the built-in defaults."""
    local_path = Path(start_dir) / LOCAL_CONFIG_NAME
    if local_path.exists():
        return {**DEFAULT_CONFIG, **json.loads(local_path.read_text())}, str(local_path)
    if GLOBAL_CONFIG_PATH.exists():
        return {**DEFAULT_CONFIG, **json.loads(GLOBAL_CONFIG_PATH.read_text())}, str(GLOBAL_CONFIG_PATH)
    return dict(DEFAULT_CONFIG), "built-in defaults"


def write_starter_config(scope):
    path = Path(LOCAL_CONFIG_NAME) if scope == "local" else GLOBAL_CONFIG_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(DEFAULT_CONFIG, indent=2) + "\n")
    return path


def parse_header(path):
    """Return the leading '#'-commented lines of a downloaded RV table."""
    comments = []
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            comments.append(line.rstrip("\n"))
    return comments


def parse_column_metadata(comments):
    columns = {}
    for line in comments:
        m = COLUMN_LINE_RE.match(line)
        if m:
            name, unit, desc = m.groups()
            columns[name] = {"unit": (unit or "").strip(), "desc": desc.strip()}
    return columns


def detect_source_type(comments):
    for line in comments:
        m = SOURCE_LINE_RE.match(line)
        if m:
            return "dace" if m.group(1).upper() == "DACE" else "vizier"
    return None


def extract_bibcode(comments, source_type):
    if source_type == "vizier":
        for line in comments:
            m = ADS_LINE_RE.match(line)
            if m:
                return m.group(1)
    elif source_type == "dace":
        for line in comments:
            if line.startswith("# Triggered by:"):
                m = BIBCODE_IN_PARENS_RE.search(line)
                if m:
                    return m.group(1)
    return None


def extract_retrieved_year(comments):
    for line in comments:
        m = RETRIEVED_RE.search(line)
        if m:
            return int(m.group(1))
    return None


def find_time_column(columns):
    for name, info in columns.items():
        if TIME_NAME_RE.match(name) or TIME_DESC_RE.search(info["desc"]):
            m = OFFSET_RE.search(info["desc"])
            return name, (float(m.group(1)) if m else 0.0)
    return None, None


def find_rv_columns(columns):
    rv_col, err_col = None, None
    for name, info in columns.items():
        is_rv_like = RV_NAME_RE.search(name) or RV_DESC_RE.search(info["desc"])
        if not is_rv_like:
            continue
        is_error = RV_ERR_NAME_RE.search(name) or ERROR_DESC_RE.search(info["desc"])
        if is_error and err_col is None:
            err_col = name
        elif not is_error and rv_col is None:
            rv_col = name
    return rv_col, err_col


def unit_to_mps_factor(unit):
    u = (unit or "").strip().lower().replace(" ", "")
    if u in ("km/s", "kms-1", "kms^-1"):
        return 1000.0
    return 1.0


def detect_vizier_instrument(columns, table_id):
    """Best-effort instrument label for a VizieR table, which usually has
    no per-row instrument column. Falls back to a stable per-table
    pseudo-label when no known RV instrument name shows up in the table's
    own column names/descriptions -- duplicate detection then simply won't
    match this table's rows against another source's, which is the correct
    (conservative) behavior when the instrument truly can't be determined.
    """
    haystack = table_id + " " + " ".join(f"{n} {i['desc']}" for n, i in columns.items())
    found = sorted({i for i in RV_INSTRUMENTS if re.search(rf'\b{re.escape(i)}\b', haystack, re.IGNORECASE)})
    return found[0] if len(found) == 1 else f"table:{table_id}"


def load_vizier_table(path, columns):
    time_col, offset = find_time_column(columns)
    rv_col, err_col = find_rv_columns(columns)
    df = pd.read_csv(path, comment="#")
    if not time_col or not rv_col or time_col not in df.columns or rv_col not in df.columns:
        return None, "could not identify time/RV columns from this table's header metadata"

    instr_col, legend = find_instrument_row_source(columns)
    if instr_col and instr_col in df.columns:
        raw = df[instr_col].astype(str).str.strip()
        if legend:
            raw = raw.map(lambda v: legend.get(v, v))
        instrument = raw.map(normalize_instrument_label)
    else:
        instrument = detect_vizier_instrument(columns, path.stem)

    out = pd.DataFrame({
        "time_bjd": df[time_col].astype(float) + offset,
        "rv_raw_mps": df[rv_col].astype(float) * unit_to_mps_factor(columns[rv_col]["unit"]),
        "rv_err_mps": (df[err_col].astype(float) * unit_to_mps_factor(columns[err_col]["unit"])
                       if err_col and err_col in df.columns else float("nan")),
        "instrument": instrument,
    })
    extras = df.drop(columns=[c for c in (time_col, rv_col, err_col) if c], errors="ignore")
    return pd.concat([out, extras.reset_index(drop=True)], axis=1), None


def load_dace_table(path):
    df = pd.read_csv(path, comment="#")
    missing = {"rjd", "rv", "rv_err"} - set(df.columns)
    if missing:
        return None, f"missing expected DACE column(s): {sorted(missing)}"

    out = pd.DataFrame({
        "time_bjd": df["rjd"].astype(float) + 2400000.0,
        "rv_raw_mps": df["rv"].astype(float),
        "rv_err_mps": df["rv_err"].astype(float),
        "instrument": df["instrument_name"] if "instrument_name" in df.columns else "unknown",
        "source_bibcode": df["pub_bibcode"] if "pub_bibcode" in df.columns else None,
    })
    drop = [c for c in ("rjd", "rv", "rv_err", "instrument_name", "pub_bibcode") if c in df.columns]
    extras = df.drop(columns=drop)
    return pd.concat([out, extras.reset_index(drop=True)], axis=1), None


def load_source_file(path):
    """Return (standardized DataFrame, warning) for one downloaded RV
    table; warning is None on success, or a string explaining why the file
    was skipped (df is then None)."""
    comments = parse_header(path)
    source_type = detect_source_type(comments)
    bibcode = extract_bibcode(comments, source_type)
    fallback_year = extract_retrieved_year(comments)

    if source_type == "dace":
        df, warning = load_dace_table(path)
    elif source_type == "vizier":
        df, warning = load_vizier_table(path, parse_column_metadata(comments))
    else:
        return None, "no recognized '# Source:' header (not produced by this pipeline's download stage?)"
    if df is None:
        return None, warning

    if "source_bibcode" not in df.columns:
        df["source_bibcode"] = bibcode
    else:
        df["source_bibcode"] = df["source_bibcode"].fillna(bibcode)
    df["source_table"] = path.name
    df["source_type"] = source_type
    df["_fallback_year"] = fallback_year
    return df, None


def find_duplicate_groups(df, window_seconds=DUPLICATE_WINDOW_SECONDS):
    """Union-find over rows sharing an instrument, sorted by time: any two
    rows from *different* source tables within `window_seconds` of each
    other are linked (transitively) into the same group. Returns a Series
    aligned to df.index: an integer group id for rows with at least one
    match, NA otherwise."""
    window_days = window_seconds / 86400.0
    n = len(df)
    parent = list(range(n))

    def find(i):
        root = i
        while parent[root] != root:
            root = parent[root]
        while parent[i] != root:
            parent[i], i = root, parent[i]
        return root

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj

    for _, idx in df.groupby("instrument").groups.items():
        order = sorted(idx, key=lambda i: df.at[i, "time_bjd"])
        for a in range(len(order)):
            i = order[a]
            b = a + 1
            while b < len(order) and df.at[order[b], "time_bjd"] - df.at[i, "time_bjd"] <= window_days:
                j = order[b]
                if df.at[i, "source_table"] != df.at[j, "source_table"]:
                    union(i, j)
                b += 1

    roots = pd.Series([find(i) for i in range(n)], index=df.index)
    counts = roots.value_counts()
    is_dup = roots.map(counts) > 1
    uniq = sorted(roots[is_dup].unique())
    relabel = {g: i + 1 for i, g in enumerate(uniq)}
    return roots.where(is_dup).map(relabel)


def choose_preferred(df, group_id, config):
    order = config.get("preferred_source_order", DEFAULT_CONFIG["preferred_source_order"])
    rank = {name: i for i, name in enumerate(order)}
    strategy = config.get("duplicate_default", DEFAULT_CONFIG["duplicate_default"])

    is_preferred = pd.Series(True, index=df.index)
    for _, idx in df.groupby(group_id, dropna=True).groups.items():
        idx = list(idx)
        if strategy == "most_recent":
            best = max(idx, key=lambda i: (
                df.at[i, "year"] if pd.notna(df.at[i, "year"]) else -1,
                -rank.get(df.at[i, "source_type"], len(order)),
            ))
        else:
            best = min(idx, key=lambda i: (rank.get(df.at[i, "source_type"], len(order)), df.at[i, "source_table"]))
        for i in idx:
            is_preferred.at[i] = (i == best)
    return is_preferred


FRONT_COLUMNS = [
    "host", "time_bjd", "rv_mps", "rv_err_mps", "rv_raw_mps", "rv_median_subtracted_mps",
    "instrument", "source_table", "source_type", "source_bibcode", "year",
    "duplicate_group", "is_preferred",
]


def aggregate_host(host_dir, hostname, config, window_seconds=DUPLICATE_WINDOW_SECONDS):
    """Combine every downloaded RV table for one host into a single
    analysis-ready DataFrame. Returns (DataFrame or None, warnings, source_files)."""
    files = sorted(p for p in host_dir.glob("*.csv") if not p.name.endswith("_rvs_aggregated.csv"))
    frames, warnings, used = [], [], []
    for path in files:
        df, warning = load_source_file(path)
        if warning:
            warnings.append(f"{path.name}: {warning}")
            continue
        df = df.dropna(subset=["time_bjd", "rv_raw_mps"])
        if len(df) == 0:
            warnings.append(f"{path.name}: no usable rows after parsing")
            continue
        frames.append(df)
        used.append(path.name)
    if not frames:
        return None, warnings, used

    combined = pd.concat(frames, ignore_index=True, sort=False)

    group_key = combined["source_table"].astype(str) + "::" + combined["instrument"].astype(str)
    medians = combined.groupby(group_key)["rv_raw_mps"].median()
    combined["rv_median_subtracted_mps"] = group_key.map(medians)
    combined["rv_mps"] = combined["rv_raw_mps"] - combined["rv_median_subtracted_mps"]

    combined["year"] = pd.to_numeric(combined["source_bibcode"].astype(str).str.slice(0, 4), errors="coerce")
    combined["year"] = combined["year"].fillna(combined["_fallback_year"])
    combined = combined.drop(columns=["_fallback_year"])

    combined["duplicate_group"] = find_duplicate_groups(combined, window_seconds=window_seconds)
    combined["is_preferred"] = choose_preferred(combined, combined["duplicate_group"], config)

    combined.insert(0, "host", hostname)
    combined = combined.sort_values("time_bjd").reset_index(drop=True)
    extra_cols = [c for c in combined.columns if c not in FRONT_COLUMNS]
    combined = combined[FRONT_COLUMNS + extra_cols]
    return combined, warnings, used


COLUMN_DOCS = [
    ("time_bjd", "d", "Barycentric (or similar) Julian Date, BJD0 offset restored"),
    ("rv_mps", "m/s", "Radial velocity, analysis-ready: median-subtracted per (source_table, instrument)"),
    ("rv_err_mps", "m/s", "Radial velocity uncertainty"),
    ("rv_raw_mps", "m/s", "Radial velocity before median subtraction"),
    ("rv_median_subtracted_mps", "m/s", "Per (source_table, instrument) median subtracted from rv_raw_mps"),
    ("instrument", "-", "Spectrograph/instrument name (best-effort 'table:<id>' fallback if undetermined)"),
    ("source_table", "-", "Downloaded file this row came from (under downloaded_rv_tables/<host>/)"),
    ("source_type", "-", "vizier or dace"),
    ("source_bibcode", "-", "ADS bibcode for this row's source/publication, if known"),
    ("year", "-", "Publication/retrieval year used for the 'most_recent' duplicate policy"),
    ("duplicate_group", "-", "Shared id for rows flagged as likely the same observation from different tables (blank = no duplicate found)"),
    ("is_preferred", "-", "True for the row to use for analysis within a duplicate_group; all non-duplicate rows are True"),
]


def write_aggregated(out_path, hostname, df, config, config_source, used_files, warnings):
    lines = [
        f"# Host: {hostname}",
        f"# Aggregated from {len(used_files)} source table(s): {', '.join(used_files)}",
        f"# Generated by rv_search_and_download.aggregate on {date.today().isoformat()}",
        f"# Duplicate-detection window: {DUPLICATE_WINDOW_SECONDS}s; "
        f"policy: {config.get('duplicate_default')} "
        f"(preferred_source_order={config.get('preferred_source_order')}) -- from {config_source}",
    ]
    for w in warnings:
        lines.append(f"# WARNING: skipped {w}")
    lines.append("# Columns:")
    for name, unit, desc in COLUMN_DOCS:
        lines.append(f"#   {name} [{unit}] - {desc}")
    lines.append("# All other columns are carried over verbatim from each source table and are NOT matched across tables.")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    df.to_csv(out_path, mode="a", index=False)


def aggregate_all(downloaded_dir="downloaded_rv_tables", out_dir="aggregated_rv_tables",
                   host=None, config_dir=".", window_seconds=DUPLICATE_WINDOW_SECONDS, verbose=False):
    """Aggregate every host (or a comma-separated subset via `host`) that
    has at least one downloaded RV table. Returns {hostname: out_path}."""
    downloaded_dir = Path(downloaded_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    config, config_source = load_config(config_dir)

    wanted = {h.strip() for h in host.split(",")} if host else None

    results = {}
    for host_dir in sorted(p for p in downloaded_dir.iterdir() if p.is_dir()):
        hostname = host_dir.name.replace("_", " ")
        if wanted and hostname not in wanted and host_dir.name not in wanted:
            continue
        print(hostname)
        df, warnings, used = aggregate_host(host_dir, hostname, config, window_seconds)
        for w in warnings:
            print(f"  ! {w}")
        if df is None:
            print("  ! no usable source tables -- skipping")
            continue
        n_dup = int(df["duplicate_group"].notna().sum())
        out_path = out_dir / f"{sanitize(hostname)}_rvs_aggregated.csv"
        write_aggregated(out_path, hostname, df, config, config_source, used, warnings)
        print(f"  wrote {out_path} ({len(df)} rows from {len(used)} table(s), {n_dup} flagged as possible duplicates)")
        results[hostname] = out_path

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--downloaded-dir", default="downloaded_rv_tables")
    parser.add_argument("--out-dir", default="aggregated_rv_tables")
    parser.add_argument("--host", default=None,
                         help="One hostname, a comma-separated list, or omit to process every downloaded host")
    parser.add_argument("--config-dir", default=".", help="Directory to look for a local rv_aggregate_config.json in")
    parser.add_argument("--window-seconds", type=float, default=DUPLICATE_WINDOW_SECONDS,
                         help="Max timestamp difference to treat two RVs as possible duplicates")
    parser.add_argument("--init-config", choices=["local", "global"], default=None,
                         help="Write a starter rv_aggregate_config.json (local: current dir; global: ~/.config/rv_search_and_download/) and exit")
    args = parser.parse_args()

    if args.init_config:
        path = write_starter_config(args.init_config)
        print(f"Wrote starter config to {path}")
        return

    aggregate_all(
        downloaded_dir=args.downloaded_dir,
        out_dir=args.out_dir,
        host=args.host,
        config_dir=args.config_dir,
        window_seconds=args.window_seconds,
    )


if __name__ == "__main__":
    main()
