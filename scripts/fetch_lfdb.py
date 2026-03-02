#!/usr/bin/env python3
"""
fetch_lfdb.py — Download and process LfDB Ka values from GlyCosmos.

LfDB (Lectin Frontier Database) provides quantitative Ka values from 
frontal affinity chromatography (FAC) at 25°C for lectin-glycan pairs.

Data source: https://download.glycosmos.org/glycosmos_lectins_lfdb_list.csv
Units: Ka in 10⁴ M⁻¹ (confirmed from Hirabayashi et al. 2002 Biochim Biophys Acta)

Run from glycan-scoring root:
    python scripts/fetch_lfdb.py
"""

import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

try:
    import urllib.request
except ImportError:
    print("ERROR: urllib required")
    sys.exit(1)

# Project lectins with LfDB IDs
PROJECT_LECTINS = {
    "LfDB0163": {"name": "WGA", "full": "Wheat Germ Agglutinin", "spec": "GlcNAc", "uniprot": "P10969"},
    "LfDB0172": {"name": "PNA", "full": "Peanut Agglutinin", "spec": "Gal", "uniprot": "P02872"},
    "LfDB0056": {"name": "Galectin-3", "full": "Galectin-3 (Human)", "spec": "Gal", "uniprot": "P17931"},
    "LfDB0065": {"name": "Galectin-3 CRD", "full": "Galectin-3 CRD (Human)", "spec": "Gal", "uniprot": "P17931"},
    "LfDB0166": {"name": "SBA", "full": "Soybean Agglutinin", "spec": "GalNAc", "uniprot": "P05046"},
    "LfDB0159": {"name": "Conarva", "full": "Conarva (ConA-like)", "spec": "Man/Glc", "uniprot": "Q9FVA0"},
}

LFDB_URL = "https://download.glycosmos.org/glycosmos_lectins_lfdb_list.csv"
LECTINS_URL = "https://download.glycosmos.org/glycosmos_lectins_list.csv"


def download_csv(url: str, dest: str) -> int:
    """Download CSV from URL. Returns row count."""
    print(f"  Downloading {url}...")
    urllib.request.urlretrieve(url, dest)
    with open(dest) as f:
        count = sum(1 for _ in f) - 1  # minus header
    print(f"  → {count} rows saved to {dest}")
    return count


def process_lfdb(raw_path: str, outdir: str) -> dict:
    """Process raw LfDB CSV into per-lectin files. Returns counts."""
    os.makedirs(outdir, exist_ok=True)
    
    lfdb = defaultdict(list)
    with open(raw_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            lfdb_id = row.get("LfDB ID", "")
            lfdb[lfdb_id].append(row)
    
    counts = {}
    for lfdb_id, info in PROJECT_LECTINS.items():
        entries = lfdb.get(lfdb_id, [])
        if not entries:
            print(f"  WARNING: No data for {info['name']} ({lfdb_id})")
            counts[info["name"]] = 0
            continue
        
        fname = f"{info['name'].replace(' ', '_')}_Ka.csv"
        fpath = os.path.join(outdir, fname)
        with open(fpath, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["lfdb_id", "lectin", "specificity", "glycotoucan_id",
                         "Ka_x10e4_per_M", "source"])
            for e in sorted(entries, key=lambda x: float(x.get("Value", 0) or 0), reverse=True):
                w.writerow([lfdb_id, info["name"], info["spec"],
                            e.get("GlyTouCan ID", ""),
                            e.get("Value", ""),
                            "LfDB/FAC/GlyCosmos"])
        print(f"  {fname}: {len(entries)} Ka values")
        counts[info["name"]] = len(entries)
    
    return counts


def build_cross_lectin_matrix(raw_path: str, outpath: str) -> int:
    """Build cross-lectin selectivity matrix for glycans tested in ≥3 lectins."""
    lfdb = defaultdict(list)
    with open(raw_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            lfdb[row.get("LfDB ID", "")].append(row)
    
    glycan_to_lectins = defaultdict(dict)
    for lfdb_id, info in PROJECT_LECTINS.items():
        for entry in lfdb.get(lfdb_id, []):
            gid = entry.get("GlyTouCan ID", "")
            ka = float(entry.get("Value", 0) or 0)
            if ka > 0 and gid:
                glycan_to_lectins[gid][info["name"]] = ka
    
    shared = {gid: lectins for gid, lectins in glycan_to_lectins.items()
              if len(lectins) >= 3}
    
    lectin_names = sorted(set(n for lectins in shared.values() for n in lectins))
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["glycotoucan_id"] + [f"Ka_{name}_x10e4_per_M" for name in lectin_names])
        for gid in sorted(shared.keys()):
            row = [gid] + [shared[gid].get(name, "") for name in lectin_names]
            w.writerow(row)
    
    print(f"  cross_lectin_selectivity.csv: {len(shared)} glycans × {len(lectin_names)} lectins")
    return len(shared)


def main():
    root = Path(__file__).parent.parent
    data = root / "data"
    lfdb_dir = data / "lfdb"
    
    print("=" * 60)
    print("LfDB Data Fetch — GlyCosmos Lectin Frontier Database")
    print("=" * 60)
    
    # 1. Download raw data
    print("\n[1/4] Downloading LfDB Ka values...")
    raw_path = str(lfdb_dir / "lfdb_raw_all.csv")
    os.makedirs(str(lfdb_dir), exist_ok=True)
    n_raw = download_csv(LFDB_URL, raw_path)
    
    # 2. Download lectins cross-reference
    print("\n[2/4] Downloading GlyCosmos lectins list...")
    lectins_path = str(data / "glycosmos_lectins_list.csv")
    download_csv(LECTINS_URL, lectins_path)
    
    # 3. Process per-lectin CSVs
    print("\n[3/4] Processing per-lectin Ka files...")
    counts = process_lfdb(raw_path, str(lfdb_dir))
    
    # 4. Build cross-lectin matrix
    print("\n[4/4] Building cross-lectin selectivity matrix...")
    n_shared = build_cross_lectin_matrix(raw_path, 
        str(data / "answer_keys" / "cross_lectin_selectivity.csv"))
    
    # Summary
    total = sum(counts.values())
    print("\n" + "=" * 60)
    print(f"DONE: {total} Ka values for {sum(1 for c in counts.values() if c > 0)} lectins")
    print(f"Raw LfDB: {n_raw} entries across all lectins in database")
    print(f"Cross-selectivity: {n_shared} glycans with ≥3 lectin measurements")
    print()
    print("Ka units: 10⁴ M⁻¹ (FAC at 25°C)")
    print("Reference: Hirabayashi J et al. (2002) Biochim Biophys Acta 1572:232-254")
    print()
    print("STILL NEEDED (no LfDB data):")
    print("  ConA (LfDB0170): no FAC data in LfDB")
    print("    → Use Dam & Brewer 2002 Chem Rev 102:387 Tables 1-3")
    print("    → Or Chervenak & Toone 1995 Biochemistry 34:5685 (ITC)")
    print("  Lysozyme: not in LfDB (not a lectin)")
    print("    → Use Kumagai 1993 Eur J Biochem 212:151")
    print("=" * 60)


if __name__ == "__main__":
    main()
