#!/usr/bin/env python3
"""
OnsiteObs_scrape.py
====================
Download sargassum beaching observations from Sargassum Monitoring's
uMap-hosted maps, filter by island bounding box,
and export to CSVs.

How it works
------------
Each yearly map at sargassummonitoring.com is just an embedded uMap
(OpenStreetMap France). uMap stores every map's content as one or more
"datalayers", each of which is served as a public GeoJSON file.
Two endpoints are all we need:

  1. Map metadata (lists the datalayers):
       GET https://umap.openstreetmap.fr/en/map/<map_id>/geojson/
  2. One datalayer's full GeoJSON (the actual points):
       GET https://umap.openstreetmap.fr/en/datalayer/<map_id>/<uuid>/

The script walks every year in MAP_IDS, fetches every datalayer,
parses each Point feature, pulls a date out of `name`+`description`
(the website warns: dates are MM/DD/YYYY), and filters by the per-island
bounding boxes in ISLAND_BBOXES.

Data ownership
--------------
This script is a technical aid for retrieving publicly served JSON.
Before using the data in a publication, I have contacted Sargassum Monitoring
(https://sargassummonitoring.com/en/contact/) for permission, and I will credit
them in my final writings.

Output layout
-------------
  data/
    raw/<year>/<datalayer_uuid>.geojson      # untouched raw GeoJSON per layer
    raw/<year>/_map_meta.json                # map metadata (list of datalayers)
    by_year/<island>_<year>.csv              # one island, one year
    <island>_all.csv                         # one island, all years combined
    all_observations.csv                     # every island, every year

Run
---
    python3 OnsiteObs_scrape.py
    # subsequent runs skip raw files that already exist on disk

Dependencies: standard library only (no pandas / requests required).
Tested with Python 3.9+.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
import sys
import time
from datetime import datetime
from pathlib import Path
from urllib.parse import urljoin
import urllib.request
import urllib.error

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

# Map IDs harvested from the iframe URLs at sargassummonitoring.com
# Verify before running by visiting e.g. sargassummonitoring.com/en/official-map-2025/
# and looking at the iframe URL: umap.openstreetmap.fr/.../map/<slug>_<map_id>
MAP_IDS = {
    2018: 951751,
    2019: 952044,
    2020: 952063,
    2021: 952070,
    2022: 952071,
    2023: 953077,
    2024: 1003076,
    2025: 1154887,
    2026: 1336079,
}

UMAP_BASE = "https://umap.openstreetmap.fr"

# Bounding boxes per island: (lat_min, lat_max, lon_min, lon_max)
# Generous around the coast to catch shore-adjacent observations.
ISLAND_BBOXES = {
    "bonaire":  (11.97, 12.36, -68.48, -68.13),
    "aruba":    (12.35, 12.68, -70.12, -69.81),
    "barbados": (12.99, 13.39, -59.70, -59.36),
}

# Be polite. uMap is a free community service.
USER_AGENT = "Master-thesis-research-scraper/1.0 (contact: 755056xl@eur.nl)"
REQUEST_DELAY_SEC = 1.0
REQUEST_TIMEOUT_SEC = 60
MAX_RETRIES = 3

OUTPUT_DIR = Path("/Users/arthliams/Documents/2026 Thesis/Data&Code/OBSdata")
RAW_DIR = OUTPUT_DIR / "raw"
BY_YEAR_DIR = OUTPUT_DIR / "by_year"

CSV_COLUMNS = [
    "year", "island", "date", "datetime_raw", "longitude", "latitude",
    "name", "description", "photo_url", "map_id", "datalayer_uuid",
    "feature_id",
]


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------

def http_get_json(url: str) -> dict:
    """GET a URL and parse JSON. Retries on transient errors. Handles gzip."""
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "application/json, application/geo+json, */*",
            "Accept-Encoding": "gzip",
        },
    )
    last_err = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT_SEC) as resp:
                raw = resp.read()
                if resp.headers.get("Content-Encoding") == "gzip":
                    raw = gzip.decompress(raw)
                return json.loads(raw.decode("utf-8"))
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_err = e
            wait = 2 ** attempt
            print(f"  ! request failed (attempt {attempt}/{MAX_RETRIES}): {e} -- sleeping {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"giving up on {url}: {last_err}")


# ---------------------------------------------------------------------------
# uMap-specific helpers
# ---------------------------------------------------------------------------

def map_geojson_url(map_id: int) -> str:
    return f"{UMAP_BASE}/en/map/{map_id}/geojson/"


def datalayer_url(map_id: int, layer_uuid: str) -> str:
    return f"{UMAP_BASE}/en/datalayer/{map_id}/{layer_uuid}/"


def extract_datalayer_ids(map_meta: dict) -> list[str]:
    """Return a list of datalayer UUIDs from a map's metadata response.

    uMap stores the list under properties.datalayers. Each entry has an `id`
    that is a UUID string (modern uMap) or an integer (legacy maps). The
    function returns them as-is so they can be plugged into the URL.
    """
    props = map_meta.get("properties", {}) or {}
    layers = props.get("datalayers", []) or []
    out = []
    for layer in layers:
        # Older versions had children under "datalayers"; newer ones nest groups.
        if isinstance(layer, dict):
            if "id" in layer:
                out.append(str(layer["id"]))
            # Some maps group layers; recurse on nested 'datalayers'
            for sub in layer.get("datalayers", []) or []:
                if isinstance(sub, dict) and "id" in sub:
                    out.append(str(sub["id"]))
    return out


# ---------------------------------------------------------------------------
# Feature parsing
# ---------------------------------------------------------------------------

HTML_TAG_RE = re.compile(r"<[^>]+>")
MD_LINK_RE = re.compile(r"\[([^\]]*)\]\(([^)]+)\)")
IMG_URL_RE = re.compile(r"https?://\S+\.(?:jpe?g|png|webp|gif|mp4)", re.IGNORECASE)

# Date patterns. The site says MM/DD/YYYY but we also accept DD/MM/YYYY
# (some volunteer entries get it backwards) and ISO YYYY-MM-DD.
DATE_PATTERNS = [
    # MM/DD/YYYY or M/D/YYYY  (preferred per site warning)
    (re.compile(r"\b(\d{1,2})/(\d{1,2})/(\d{4})\b"), "mdy"),
    # YYYY-MM-DD
    (re.compile(r"\b(\d{4})-(\d{1,2})-(\d{1,2})\b"), "ymd"),
    # YYYY/MM/DD
    (re.compile(r"\b(\d{4})/(\d{1,2})/(\d{1,2})\b"), "ymd"),
]


def clean_text(s: str | None) -> str:
    if not s:
        return ""
    # markdown links: keep the label
    s = MD_LINK_RE.sub(r"\1", s)
    # strip HTML
    s = HTML_TAG_RE.sub(" ", s)
    # collapse whitespace
    return re.sub(r"\s+", " ", s).strip()


def find_photo_url(text: str) -> str:
    if not text:
        return ""
    m = IMG_URL_RE.search(text)
    return m.group(0) if m else ""


def parse_date(text: str, fallback_year: int) -> tuple[str, str]:
    """Return (ISO date string or '', original matched string).

    Strategy: scan with all patterns, pick the first that yields a plausible
    date in or near fallback_year (the year the map represents). Plausible =
    fallback_year-1 to fallback_year+1 (some observations posted on the rim
    of a calendar year end up on the next year's map).
    """
    if not text:
        return "", ""
    candidates: list[tuple[datetime, str]] = []
    for pattern, kind in DATE_PATTERNS:
        for m in pattern.finditer(text):
            try:
                if kind == "mdy":
                    a, b, y = int(m.group(1)), int(m.group(2)), int(m.group(3))
                    # try MM/DD first
                    for month, day in ((a, b), (b, a)):
                        try:
                            dt = datetime(y, month, day)
                        except ValueError:
                            continue
                        if fallback_year - 1 <= dt.year <= fallback_year + 1:
                            candidates.append((dt, m.group(0)))
                            break
                else:  # ymd
                    y, month, day = int(m.group(1)), int(m.group(2)), int(m.group(3))
                    dt = datetime(y, month, day)
                    if fallback_year - 1 <= dt.year <= fallback_year + 1:
                        candidates.append((dt, m.group(0)))
            except (ValueError, IndexError):
                continue
    if not candidates:
        return "", ""
    # Earliest valid candidate tends to be the observation date
    candidates.sort(key=lambda x: x[0])
    dt, raw = candidates[0]
    return dt.date().isoformat(), raw


def feature_centroid(geom: dict) -> tuple[float | None, float | None]:
    """Return (lon, lat) of a feature's representative point, or (None, None)."""
    if not geom:
        return None, None
    gtype = geom.get("type")
    coords = geom.get("coordinates")
    try:
        if gtype == "Point":
            return float(coords[0]), float(coords[1])
        if gtype == "MultiPoint" or gtype == "LineString":
            pts = coords
        elif gtype == "MultiLineString" or gtype == "Polygon":
            pts = [p for ring in coords for p in ring]
        elif gtype == "MultiPolygon":
            pts = [p for poly in coords for ring in poly for p in ring]
        else:
            return None, None
        if not pts:
            return None, None
        xs = [float(p[0]) for p in pts]
        ys = [float(p[1]) for p in pts]
        return sum(xs) / len(xs), sum(ys) / len(ys)
    except (TypeError, ValueError, IndexError):
        return None, None


def in_bbox(lon: float, lat: float, bbox: tuple[float, float, float, float]) -> bool:
    lat_min, lat_max, lon_min, lon_max = bbox
    return (lat_min <= lat <= lat_max) and (lon_min <= lon <= lon_max)


def which_island(lon: float, lat: float) -> str | None:
    for name, bbox in ISLAND_BBOXES.items():
        if in_bbox(lon, lat, bbox):
            return name
    return None


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def fetch_year(year: int, map_id: int) -> list[dict]:
    """Fetch every feature for one yearly map. Returns list of dicts."""
    year_raw_dir = RAW_DIR / str(year)
    year_raw_dir.mkdir(parents=True, exist_ok=True)

    # 1. Map metadata (datalayer list)
    meta_path = year_raw_dir / "_map_meta.json"
    if meta_path.exists():
        print(f"[{year}] using cached map metadata")
        map_meta = json.loads(meta_path.read_text(encoding="utf-8"))
    else:
        url = map_geojson_url(map_id)
        print(f"[{year}] GET {url}")
        map_meta = http_get_json(url)
        meta_path.write_text(json.dumps(map_meta, ensure_ascii=False, indent=2), encoding="utf-8")
        time.sleep(REQUEST_DELAY_SEC)

    layer_ids = extract_datalayer_ids(map_meta)
    print(f"[{year}] {len(layer_ids)} datalayers")

    rows: list[dict] = []
    for i, uuid_str in enumerate(layer_ids, 1):
        layer_path = year_raw_dir / f"{uuid_str}.geojson"
        if layer_path.exists():
            print(f"  ({i}/{len(layer_ids)}) cached {uuid_str}")
            fc = json.loads(layer_path.read_text(encoding="utf-8"))
        else:
            url = datalayer_url(map_id, uuid_str)
            print(f"  ({i}/{len(layer_ids)}) GET {url}")
            fc = http_get_json(url)
            layer_path.write_text(json.dumps(fc, ensure_ascii=False), encoding="utf-8")
            time.sleep(REQUEST_DELAY_SEC)

        for feat in fc.get("features", []) or []:
            lon, lat = feature_centroid(feat.get("geometry") or {})
            if lon is None or lat is None:
                continue
            island = which_island(lon, lat)
            if not island:
                continue
            props = feat.get("properties", {}) or {}
            name = clean_text(props.get("name"))
            desc_raw = props.get("description") or ""
            desc = clean_text(desc_raw)
            # photo URL: scan both raw description (in case markdown was stripped)
            # and any properties matching image-like keys
            photo = find_photo_url(desc_raw) or find_photo_url(json.dumps(props))
            iso_date, raw_date = parse_date(f"{name} {desc}", fallback_year=year)
            rows.append({
                "year": year,
                "island": island,
                "date": iso_date,
                "datetime_raw": raw_date,
                "longitude": lon,
                "latitude": lat,
                "name": name,
                "description": desc,
                "photo_url": photo,
                "map_id": map_id,
                "datalayer_uuid": uuid_str,
                "feature_id": str(feat.get("id", "")),
            })
    print(f"[{year}] kept {len(rows)} features inside island bboxes")
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in CSV_COLUMNS})
    print(f"  wrote {len(rows)} rows -> {path}")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    BY_YEAR_DIR.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []
    for year in sorted(MAP_IDS):
        map_id = MAP_IDS[year]
        try:
            year_rows = fetch_year(year, map_id)
        except Exception as e:
            print(f"[{year}] FAILED: {e}", file=sys.stderr)
            continue
        all_rows.extend(year_rows)

        # Per-island per-year CSV
        for island in ISLAND_BBOXES:
            island_year_rows = [r for r in year_rows if r["island"] == island]
            if island_year_rows:
                write_csv(BY_YEAR_DIR / f"{island}_{year}.csv", island_year_rows)

    # Per-island combined CSV
    for island in ISLAND_BBOXES:
        island_rows = [r for r in all_rows if r["island"] == island]
        island_rows.sort(key=lambda r: (r["date"] or "", r["year"]))
        write_csv(OUTPUT_DIR / f"{island}_all.csv", island_rows)

    # Master CSV
    all_rows.sort(key=lambda r: (r["year"], r["island"], r["date"] or ""))
    write_csv(OUTPUT_DIR / "all_observations.csv", all_rows)

    # Quick summary
    print("\nSummary:")
    for island in ISLAND_BBOXES:
        n = sum(1 for r in all_rows if r["island"] == island)
        print(f"  {island:>9s}: {n} observations")
    print(f"  {'TOTAL':>9s}: {len(all_rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#Summary:
#    bonaire: 161 observations
#      aruba: 25 observations
#   barbados: 4725 observations
#      TOTAL: 4911