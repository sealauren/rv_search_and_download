"""For hosts where vizier_search.py found an Msini/RV-amplitude provenance
clue but no matching VizieR table, identify which RV instrument/survey the
paper actually used -- a clue that the paper is reusing archival RVs (common
for nearby bright stars) rather than presenting its own VizieR-mirrored
dataset.

Three-tier search:
  1. ADS abstract (fast, works for any paper ADS indexes).
  2. If tier 1 finds nothing and the paper has an arXiv id, full text of
     the arXiv preprint via ar5iv.org (LaTeXML HTML rendering -- no PDF
     parsing needed). The arXiv id is resolved from ADS's `identifier`
     field, which works even when the paper's canonical bibcode is the
     published (non-arXiv) version. LaTeXML can fail to render a paper
     (common for Nature-formatted submissions) while still returning
     HTTP 200 -- that's detected and treated as "couldn't check", not as
     a confirmed empty result.
  3. If tiers 1-2 still find nothing -- including when there's no arXiv
     id at all, or ar5iv failed to render -- ADS's own full-text search
     index, which is built from the publisher and isn't subject to
     arXiv/LaTeXML rendering failures.

Two special-case redirects this also flags:
  - "California Legacy Survey" / "CLS" -> the actual RV table lives in
    Rosenthal et al. 2021 (2021ApJS..255...8R), not the citing paper.
  - "DACE", or a named DACE-pipeline spectrograph (ESPRESSO, HARPS,
    HARPS-N, CORALIE, CARMENES) -> the RVs likely live in the DACE
    archive and must be queried there rather than downloaded from a
    VizieR table.

Requires an ADS API token (free, from https://ui.adsabs.harvard.edu/user/settings/token)
set as the $ADS_API_TOKEN environment variable.

Usage (after `pip install -e .`):
    rv-ads-search --input vizier_rv_results.csv --output ads_rv_instruments.csv
"""
import argparse
import json
import os
import re
import time
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

from .vizier_search import unresolved_provenance_refs

ADS_URL = "https://api.adsabs.harvard.edu/v1/search/query"
AR5IV_URL = "https://ar5iv.labs.arxiv.org/html/{arxiv_id}"
HEADERS = {"User-Agent": "ads_rv_instrument_search.py (research script; contact lweiss4@nd.edu)"}

ROSENTHAL_2021_BIBCODE = "2021ApJS..255....8R"

# Specific RV spectrograph/survey names -- kept specific (not generic
# observatory names like "Keck" alone) to avoid false positives from papers
# that mention a facility for unrelated reasons (e.g. imaging).
RV_INSTRUMENTS = [
    "HARPS-N", "HARPS", "HIRES", "APF", "ESPRESSO", "CARMENES", "SOPHIE",
    "CORALIE", "PFS", "NEID", "MAROON-X", "EXPRES", "UVES", "FEROS", "TRES",
    "PARAS", "CHIRON", "ELODIE", "UCLES", "MINERVA", "SONG", "IRD", "GIANO",
    "iSHELL", "AAPS", "LCES", "SPIRou", "ESPaDOnS",
]
INSTRUMENT_RE = re.compile(r"\b(" + "|".join(re.escape(i) for i in RV_INSTRUMENTS) + r")\b")
ARCHIVAL_RE = re.compile(
    r"archival radial veloc|previously published (radial veloc|RV)|publicly available radial veloc",
    re.IGNORECASE,
)
CLS_RE = re.compile(r"California Legacy Survey|\bCLS\b")
DACE_RE = re.compile(r"\bDACE\b")

# These spectrographs' public RV time series are routinely archived in
# DACE even when a paper never says the word "DACE" -- a paper just naming
# one of them is itself a DACE clue (see find_rv_clues below). Originally
# just the Geneva Observatory pipeline instruments (ESPRESSO/HARPS/
# HARPS-N/CORALIE); CARMENES and SOPHIE added per direct confirmation that
# their public RVs are also mirrored there.
DACE_PIPELINE_INSTRUMENTS = {"ESPRESSO", "HARPS-N", "HARPS", "CORALIE", "CARMENES", "SOPHIE"}

ARXIV_IDENTIFIER_RE = re.compile(r"^arXiv:(.+)$")


def load_token():
    token = os.environ.get("ADS_API_TOKEN")
    if token:
        return token.strip()
    raise SystemExit(
        "No ADS API token found. Get a free token from "
        "https://ui.adsabs.harvard.edu/user/settings/token and `export ADS_API_TOKEN=...`"
    )


def fetch_ads_record(bibcode, token, session):
    """Return (title, abstract, arxiv_id) for a bibcode, or (None, None, None)
    if ADS has no record. Queries by `identifier:` rather than `bibcode:`,
    because preprint bibcodes (e.g. the NEA's reflink bibcode) are often
    superseded by the published bibcode once a paper comes out, and a direct
    `bibcode:` query doesn't follow that redirect -- `identifier:` matches
    against the full alternate-identifier list instead. arxiv_id is then
    resolved from that same `identifier` list.
    """
    resp = session.get(
        ADS_URL,
        params={"q": f"identifier:{bibcode}", "fl": "title,abstract,identifier"},
        headers={"Authorization": f"Bearer {token}"},
        timeout=30,
    )
    resp.raise_for_status()
    docs = resp.json()["response"]["docs"]
    if not docs:
        return None, None, None
    doc = docs[0]
    title = doc.get("title", [None])[0]
    abstract = doc.get("abstract")
    arxiv_id = None
    for ident in doc.get("identifier", []):
        m = ARXIV_IDENTIFIER_RE.match(ident)
        if m:
            arxiv_id = m.group(1)
            break
    return title, abstract, arxiv_id


# LaTeXML (the engine behind ar5iv) can fail to render a paper's LaTeX --
# common for Nature-formatted submissions -- and still return an HTTP 200
# page saying so, rather than an error status. Treating that page's
# (near-empty) text as "the full text, which mentions no RV instrument"
# would be wrong: we never actually read the paper.
AR5IV_FAILURE_MARKER = "Conversion to HTML had a Fatal error"


def fetch_ar5iv_fulltext(arxiv_id, session):
    resp = session.get(AR5IV_URL.format(arxiv_id=arxiv_id), headers=HEADERS, timeout=60)
    if resp.status_code != 200:
        return None
    soup = BeautifulSoup(resp.text, "html.parser")
    text = soup.get_text(" ")
    if AR5IV_FAILURE_MARKER in text:
        return None
    return text


def find_rv_clues(text):
    if not text:
        return {"instruments": [], "archival": False, "cls": False, "dace": False}
    instruments = sorted(set(INSTRUMENT_RE.findall(text)))
    return {
        "instruments": instruments,
        "archival": bool(ARCHIVAL_RE.search(text)),
        "cls": bool(CLS_RE.search(text)),
        "dace": bool(DACE_RE.search(text)) or bool(DACE_PIPELINE_INSTRUMENTS & set(instruments)),
    }


def fetch_ads_fulltext_clues(bibcode, token, session):
    """Query ADS's own full-text search index (built from the publisher,
    not arXiv/LaTeXML) for RV instrument/CLS/DACE keywords -- a fallback
    for papers ar5iv can't render, or that have no arXiv id at all. ADS
    indexes full text independent of arXiv's HTML conversion pipeline, so
    it isn't subject to the same rendering failures.

    Returns a clues dict in the same shape as find_rv_clues, or None if the
    ADS query itself fails (network error) -- as opposed to a confirmed
    "zero matches" result, which is a legitimate, cacheable answer and
    returned normally.

    Does one combined existence query first (cheap), and only follows up
    with per-instrument queries -- to find out *which* instrument(s) hit --
    for the rare paper where that combined query found something.
    """
    def fulltext_hit(query):
        resp = session.get(
            ADS_URL,
            params={"q": f"bibcode:{bibcode} AND full:({query})", "fl": "bibcode"},
            headers={"Authorization": f"Bearer {token}"},
            timeout=30,
        )
        resp.raise_for_status()
        return resp.json()["response"]["numFound"] > 0

    all_terms = RV_INSTRUMENTS + ["DACE", "California Legacy Survey"]
    try:
        if not fulltext_hit(" OR ".join(f'"{t}"' for t in all_terms)):
            return {"instruments": [], "archival": False, "cls": False, "dace": False}
        instruments = [i for i in RV_INSTRUMENTS if fulltext_hit(f'"{i}"')]
        cls_hit = fulltext_hit('"California Legacy Survey" OR "CLS"')
        dace_hit = fulltext_hit('"DACE"')
    except requests.RequestException:
        return None

    return {
        "instruments": sorted(instruments),
        "archival": False,
        "cls": cls_hit,
        "dace": dace_hit or bool(DACE_PIPELINE_INSTRUMENTS & set(instruments)),
    }


def has_any_clue(clues):
    return bool(clues["instruments"]) or clues["archival"] or clues["cls"] or clues["dace"]


def missing_rv_references(results_path):
    """Reproduce the same 'RV provenance but no VizieR table' set that
    vizier_search.py reports at the end of its run."""
    df = pd.read_csv(results_path)
    return unresolved_provenance_refs(df)


def search_ads_instruments(vizier_results_csv, output_csv=None, cache_path="cache/ads_abstract_cache.json",
                            token=None, sleep=0.2, verbose=False):
    """Run the three-tier ADS/ar5iv/ADS-fulltext RV-instrument search and return the
    results as a DataFrame. This is the importable entry point
    (main() is just a thin CLI wrapper around this, installed as the
    `rv-ads-search` command).

    vizier_results_csv -- the output of vizier_search.search_vizier().
    output_csv          -- if given, write the results table here.
    cache_path           -- JSON cache of bibcode -> {abstract, arxiv_id,
                             fulltext_clues}, reused across runs.
    token                -- ADS API token; defaults to load_token() ($ADS_API_TOKEN).
    sleep                 -- seconds to wait between live ADS/ar5iv requests.
    verbose               -- print progress as each new bibcode is fetched.
    """
    token = token or load_token()
    missing = missing_rv_references(vizier_results_csv)

    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}

    session = requests.Session()
    rows = []
    for _, r in missing.sort_values("hostname").iterrows():
        bibcode = r["bibcode"]
        if bibcode not in cache or "arxiv_id" not in cache[bibcode]:
            if verbose:
                print(f"fetching ADS record for {bibcode} ({r['display']})")
            title, abstract, arxiv_id = fetch_ads_record(bibcode, token, session)
            cache[bibcode] = {"title": title, "abstract": abstract, "arxiv_id": arxiv_id}
            cache_path.write_text(json.dumps(cache, indent=1))
            time.sleep(sleep)

        entry = cache[bibcode]
        title, abstract, arxiv_id = entry["title"], entry["abstract"], entry["arxiv_id"]

        # Tier 1: does the abstract alone name an instrument/survey/CLS/DACE?
        clues = find_rv_clues(abstract)
        tier = "abstract"

        # Tier 2: fall back to the arXiv full text only if tier 1 came up empty.
        if not has_any_clue(clues) and arxiv_id:
            if "fulltext_clues" not in entry:
                if verbose:
                    print(f"  no abstract match -- fetching ar5iv full text for arXiv:{arxiv_id}")
                fulltext = fetch_ar5iv_fulltext(arxiv_id, session)
                # Only cache a real "read the text, found nothing" result.
                # fulltext is None both on a fetch error and on a LaTeXML
                # conversion failure (see fetch_ar5iv_fulltext) -- caching
                # that as a confirmed-empty result would permanently hide
                # a paper we never actually got to read.
                if fulltext is not None:
                    entry["fulltext_clues"] = find_rv_clues(fulltext)
                    cache_path.write_text(json.dumps(cache, indent=1))
                time.sleep(sleep)
            if "fulltext_clues" in entry and has_any_clue(entry["fulltext_clues"]):
                clues = entry["fulltext_clues"]
                tier = "fulltext"

        # Tier 3: ADS's own full-text index. Tried whenever tiers 1-2 found
        # nothing -- including when there's no arXiv id to even attempt
        # tier 2, or when ar5iv failed to render the paper (common for
        # Nature-formatted LaTeX) -- since ADS indexes full text
        # independent of arXiv's HTML conversion pipeline.
        if not has_any_clue(clues):
            if "ads_fulltext_clues" not in entry:
                if verbose:
                    print(f"  no clue from abstract/ar5iv -- trying ADS's own full-text index")
                ads_clues = fetch_ads_fulltext_clues(bibcode, token, session)
                if ads_clues is not None:
                    entry["ads_fulltext_clues"] = ads_clues
                    cache_path.write_text(json.dumps(cache, indent=1))
                time.sleep(sleep)
            if "ads_fulltext_clues" in entry and has_any_clue(entry["ads_fulltext_clues"]):
                clues = entry["ads_fulltext_clues"]
                tier = "ADS fulltext"

        if not has_any_clue(clues):
            tier = "none (checked abstract/ar5iv/ADS fulltext)"

        rows.append({
            "hostname": r["hostname"], "bibcode": bibcode, "display": r["display"],
            "title": title, "arxiv_id": arxiv_id, "tier": tier,
            "rv_instruments": ", ".join(clues["instruments"]),
            "mentions_archival_rv": clues["archival"],
            "mentions_cls": clues["cls"],
            "mentions_dace": clues["dace"],
        })

    out = pd.DataFrame(rows)
    if output_csv:
        out.to_csv(output_csv, index=False)
        print(f"wrote {len(out)} rows to {output_csv}")

    with_clue = out[
        out["rv_instruments"].astype(bool) | out["mentions_archival_rv"]
        | out["mentions_cls"] | out["mentions_dace"]
    ]
    print(f"{len(with_clue)}/{len(out)} references give an RV-source clue:")
    for _, r in with_clue.iterrows():
        tags = []
        if r["rv_instruments"]:
            tags.append(r["rv_instruments"])
        if r["mentions_cls"]:
            tags.append(f"CLS -> see Rosenthal et al. 2021 ({ROSENTHAL_2021_BIBCODE})")
        if r["mentions_dace"]:
            tags.append("DACE -> query DACE archive directly")
        if r["mentions_archival_rv"] and not tags:
            tags.append("archival RV language, no named instrument")
        print(f"  [{r['tier']}] {r['hostname']}: {r['display']} -- {'; '.join(tags)}")

    no_clue = out[~out.index.isin(with_clue.index)]
    if len(no_clue):
        print(f"\n{len(no_clue)} references still have no RV-source clue:")
        for _, r in no_clue.iterrows():
            print(f"  [{r['tier']}] {r['hostname']}: {r['display']}")

    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="vizier_rv_results.csv")
    parser.add_argument("--output", default="ads_rv_instruments.csv")
    parser.add_argument("--cache", default="cache/ads_abstract_cache.json")
    parser.add_argument("--sleep", type=float, default=0.2)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    search_ads_instruments(args.input, args.output, args.cache, sleep=args.sleep, verbose=args.verbose)


if __name__ == "__main__":
    main()
