"""
ingest_and_fit.py — Data ingestion + parameter determination pipeline.

Students run this after retrieving paper tables. It:
  1. Reads filled CSV templates from data/
  2. Validates data format and completeness
  3. Runs Phase G1-G3 parameter determination
  4. Saves locked parameters to glycan_params.yaml
  5. Re-runs predictions and reports statistics

Usage:
  python ingest_and_fit.py              # full pipeline
  python ingest_and_fit.py --validate   # validate data only, don't fit
  python ingest_and_fit.py --diagnose   # show per-OH diagnostics
"""

import sys
import os
import csv
import json
import numpy as np
from typing import Dict, List, Optional, Tuple

from glycan_scorer import GlycanParams
from parameter_fitting import (
    fit_desolvation_from_schwarz,
    fit_NAc_desolvation,
    fit_conformational_entropy,
    fit_ch_pi_from_synthetic_hosts,
    assemble_params,
)


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


# ═══════════════════════════════════════════════════════════════════════════
# Data readers
# ═══════════════════════════════════════════════════════════════════════════

def read_schwarz_csv(filepath: str) -> Dict[str, float]:
    """Read dissolution calorimetry data from schwarz_dissolution.csv.

    Returns dict: compound_name → ΔH_sol (kJ/mol) at 25°C.
    Skips rows with empty dH_sol_kJ_mol.
    """
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            if parts[0] == 'compound':  # header
                continue
            compound = parts[0].strip()
            dH_str = parts[1].strip()
            if dH_str:
                try:
                    data[compound] = float(dH_str)
                except ValueError:
                    print(f"  WARNING: Cannot parse dH_sol for '{compound}': '{dH_str}'")
    return data


def read_jasra_csv(filepath: str) -> Dict[str, float]:
    """Read polyol dissolution data from jasra_polyols.csv.

    Same format as Schwarz. Returns compound → ΔH_sol (kJ/mol).
    """
    return read_schwarz_csv(filepath)  # same CSV format


def read_laughrey_csv(filepath: str) -> List[Dict]:
    """Read CH-π host-guest data from laughrey_ch_pi.csv.

    Returns list of dicts for fit_ch_pi_from_synthetic_hosts().
    """
    entries = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith('#')),
        )
        for row in reader:
            try:
                entry = {
                    "host": row.get("host", "").strip(),
                    "guest": row.get("guest", "").strip(),
                    "dG_measured": float(row["dG_kJ_mol"]),
                    "dG_other_terms": float(row.get("dG_other_kJ_mol", 0)),
                    "n_CH_pi_contacts": int(row.get("n_ch_pi", 0)),
                }
                if entry["n_CH_pi_contacts"] > 0:
                    entries.append(entry)
            except (ValueError, KeyError):
                continue
    return entries


# ═══════════════════════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════════════════════

def validate_data():
    """Check which data files have actual values (not empty templates)."""
    print("=" * 60)
    print("DATA VALIDATION")
    print("=" * 60)

    files_status = {}

    # Schwarz dissolution
    schwarz_path = os.path.join(DATA_DIR, "schwarz_dissolution.csv")
    if os.path.exists(schwarz_path):
        data = read_schwarz_csv(schwarz_path)
        n = len(data)
        files_status["schwarz"] = n
        status = "✓ READY" if n >= 4 else f"⚠ Only {n} values (need ≥4)"
        if n == 0:
            status = "✗ EMPTY TEMPLATE"
        print(f"  schwarz_dissolution.csv: {n} compounds with ΔH_sol — {status}")
    else:
        files_status["schwarz"] = 0
        print(f"  schwarz_dissolution.csv: MISSING")

    # Jasra polyols
    jasra_path = os.path.join(DATA_DIR, "jasra_polyols.csv")
    if os.path.exists(jasra_path):
        data = read_jasra_csv(jasra_path)
        n = len(data)
        files_status["jasra"] = n
        status = "✓ READY" if n >= 2 else f"⚠ Only {n} values"
        if n == 0:
            status = "✗ EMPTY TEMPLATE"
        print(f"  jasra_polyols.csv:       {n} compounds with ΔH_sol — {status}")
    else:
        files_status["jasra"] = 0
        print(f"  jasra_polyols.csv:       MISSING")

    # Laughrey CH-π
    laughrey_path = os.path.join(DATA_DIR, "laughrey_ch_pi.csv")
    if os.path.exists(laughrey_path):
        data = read_laughrey_csv(laughrey_path)
        n = len(data)
        files_status["laughrey"] = n
        status = "✓ READY" if n >= 2 else f"⚠ Only {n} entries"
        if n == 0:
            status = "✗ EMPTY TEMPLATE"
        print(f"  laughrey_ch_pi.csv:      {n} host-guest entries — {status}")
    else:
        files_status["laughrey"] = 0
        print(f"  laughrey_ch_pi.csv:      MISSING")

    # Answer keys
    ak_dir = os.path.join(DATA_DIR, "answer_keys")
    if os.path.isdir(ak_dir):
        for fname in sorted(os.listdir(ak_dir)):
            if fname.endswith('.csv'):
                fpath = os.path.join(ak_dir, fname)
                n_total, n_filled = 0, 0
                with open(fpath, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#') or 'VERIFIED' in line:
                            continue
                        parts = line.split(',')
                        n_total += 1
                        # Check if Ka or dG column has data
                        if len(parts) > 1 and parts[1].strip():
                            n_filled += 1
                status = f"{n_filled}/{n_total} entries"
                print(f"  answer_keys/{fname}: {status}")

    # Summary
    print("\n" + "-" * 60)
    can_fit_G1 = files_status.get("schwarz", 0) >= 4
    can_fit_G3 = files_status.get("laughrey", 0) >= 2
    print(f"  Phase G1 (desolvation): {'READY' if can_fit_G1 else 'BLOCKED — need Schwarz table values'}")
    print(f"  Phase G2 (conf entropy): BLOCKED — need GLYCAM06 QM profiles")
    print(f"  Phase G3 (CH-π):         {'READY' if can_fit_G3 else 'BLOCKED — need Laughrey table values'}")

    return files_status


# ═══════════════════════════════════════════════════════════════════════════
# Fitting pipeline
# ═══════════════════════════════════════════════════════════════════════════

def run_fitting(diagnose: bool = False) -> Optional[GlycanParams]:
    """Run Phase G1-G3 parameter determination from available data."""
    print("\n" + "=" * 60)
    print("PARAMETER FITTING")
    print("=" * 60)

    # --- Phase G1: Desolvation ---
    print("\n--- Phase G1: Polyol desolvation ---")
    schwarz_path = os.path.join(DATA_DIR, "schwarz_dissolution.csv")
    schwarz_data = read_schwarz_csv(schwarz_path) if os.path.exists(schwarz_path) else {}

    if len(schwarz_data) >= 4:
        desolv_result = fit_desolvation_from_schwarz(schwarz_data)
        print(f"  k_desolv_eq = {desolv_result['k_desolv_eq']:.2f} kJ/mol "
              f"(from {desolv_result['n_eq_points']} equatorial points)")
        print(f"  k_desolv_ax = {desolv_result['k_desolv_ax']:.2f} kJ/mol "
              f"(from {desolv_result['n_ax_points']} axial points)")

        if diagnose and desolv_result.get("per_position"):
            print("\n  Per-position diagnostics:")
            for pos, info in desolv_result["per_position"].items():
                print(f"    {pos:10s}: ΔΔH = {info['ddH_kJ']:+.2f} kJ/mol ({info['stereo']})")

        # NAc from Jasra
        jasra_path = os.path.join(DATA_DIR, "jasra_polyols.csv")
        jasra_data = read_jasra_csv(jasra_path) if os.path.exists(jasra_path) else {}
        k_NAc = fit_NAc_desolvation(jasra_data)
        desolv_result["k_desolv_NAc"] = k_NAc
        print(f"  k_desolv_NAc = {k_NAc:.2f} kJ/mol")
    else:
        print(f"  SKIPPED — only {len(schwarz_data)} Schwarz values (need ≥4)")
        print(f"  Using placeholder: k_desolv_eq=8.0, k_desolv_ax=5.0, k_desolv_NAc=6.0")
        desolv_result = {"k_desolv_eq": 8.0, "k_desolv_ax": 5.0, "k_desolv_NAc": 6.0}

    # --- Phase G2: Conformational entropy ---
    print("\n--- Phase G2: Conformational entropy ---")
    # QM torsion profiles not available via web — use literature consensus
    print("  Using literature consensus: eps_conf ≈ 3.5 kJ/mol per φ/ψ pair")
    print("  (Imberty 1999, Kirschner 2008 GLYCAM06)")
    conf_result = {"eps_conf": 3.5}

    # --- Phase G3: CH-π stacking ---
    print("\n--- Phase G3: CH-π stacking ---")
    laughrey_path = os.path.join(DATA_DIR, "laughrey_ch_pi.csv")
    laughrey_data = read_laughrey_csv(laughrey_path) if os.path.exists(laughrey_path) else []

    if len(laughrey_data) >= 2:
        ch_pi_result = fit_ch_pi_from_synthetic_hosts(laughrey_data)
        print(f"  eps_CH_pi = {ch_pi_result['eps_CH_pi']:.2f} kJ/mol per contact "
              f"(from {ch_pi_result['n_entries']} entries)")
        if ch_pi_result.get("range"):
            print(f"  Range: {ch_pi_result['range'][0]:.2f} to {ch_pi_result['range'][1]:.2f}")
    else:
        print(f"  SKIPPED — only {len(laughrey_data)} Laughrey entries (need ≥2)")
        print(f"  Using Nishio consensus: eps_CH_pi ≈ -2.5 kJ/mol")
        ch_pi_result = {"eps_CH_pi": -2.5}

    # --- Assemble ---
    print("\n--- Assembling parameter set ---")
    params = assemble_params(
        desolv_result=desolv_result,
        conf_result=conf_result,
        ch_pi_result=ch_pi_result,
        eps_HB=-6.0,
        dG_0=12.0,  # placeholder; will be calibrated below
    )
    params.eps_water_bridge = -4.0  # kJ/mol per conserved water (literature estimate)

    # Fit dG_0 from ConA MeαMan if answer key available
    from contact_maps_predictions import PREDICTION_2_MAPS
    from glycan_scorer import predict_dG
    meas_man = -22.38  # kJ/mol (MeαMan, Chervenak/Dam consensus)
    pred_man = predict_dG(params, PREDICTION_2_MAPS["MeαMan"]).dG_predicted
    offset_correction = meas_man - (pred_man - params.dG_0)
    params.dG_0 = round(offset_correction, 2)
    print(f"  dG_0 calibrated to ConA:MeαMan = {params.dG_0:.2f} kJ/mol")

    print(f"\n  FINAL PARAMETERS:")
    print(f"    k_desolv_eq  = {params.k_desolv_eq:.2f} kJ/mol")
    print(f"    k_desolv_ax  = {params.k_desolv_ax:.2f} kJ/mol")
    print(f"    k_desolv_NAc = {params.k_desolv_NAc:.2f} kJ/mol")
    print(f"    eps_CH_pi    = {params.eps_CH_pi:.2f} kJ/mol")
    print(f"    eps_HB       = {params.eps_HB:.2f} kJ/mol")
    print(f"    eps_conf     = {params.eps_conf:.2f} kJ/mol")
    print(f"    eps_water    = {params.eps_water_bridge:.2f} kJ/mol")
    print(f"    dG_0         = {params.dG_0:.2f} kJ/mol")

    return params


def save_params(params: GlycanParams, filepath: str = "glycan_params.yaml"):
    """Save locked parameters to YAML."""
    content = f"""# Glycan Scoring Parameters — LOCKED
# Generated by ingest_and_fit.py
# DO NOT EDIT manually after Phase G1-G3 lock.
#
k_desolv_eq: {params.k_desolv_eq}
k_desolv_ax: {params.k_desolv_ax}
k_desolv_NAc: {params.k_desolv_NAc}
eps_CH_pi: {params.eps_CH_pi}
eps_HB: {params.eps_HB}
eps_conf: {params.eps_conf}
eps_water_bridge: {params.eps_water_bridge}
dG_0: {params.dG_0}
"""
    with open(filepath, 'w') as f:
        f.write(content)
    print(f"\n  Parameters saved to {filepath}")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    validate_only = "--validate" in sys.argv
    diagnose = "--diagnose" in sys.argv

    status = validate_data()

    if validate_only:
        return

    params = run_fitting(diagnose=diagnose)
    if params is None:
        print("\nFitting failed. Check data files.")
        return

    save_params(params)

    # Re-run predictions with fitted params
    print("\n" + "=" * 60)
    print("RE-RUNNING PREDICTIONS WITH FITTED PARAMETERS")
    print("=" * 60)
    from run_predictions import (
        run_prediction_1, run_prediction_2, run_prediction_3,
        run_prediction_4, run_prediction_6, generate_figures,
    )
    run_prediction_1(params)
    run_prediction_2(params)
    run_prediction_3(params)
    run_prediction_4(params)
    run_prediction_6(params)

    print("\n  Generating figures...")
    generate_figures(params, output_dir="figures")

    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()