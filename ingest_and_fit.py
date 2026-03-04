#!/usr/bin/env python3
"""
MABE Glycan Module — ingest_and_fit.py
Phase 1+2: Score all systems with locked v2.2 parameters.
No fitting to biological data. All parameters from pure chemistry sources.
Outputs: predictions table, residuals, R², sanity checks.
"""

import math
import csv
import sys
from parameters_v22 import score, EPS_HB, EPS_CH_PI, K_DESOLV_EQ, K_DESOLV_AX, K_DESOLV_C6, K_DESOLV_NAC, BETA_CONTEXT
from contact_maps_v1 import ALL_SYSTEMS

RT = 2.479  # kJ/mol at 298 K (8.314e-3 × 298)

def Ka_from_dG(dG_kJ):
    """Convert ΔG (kJ/mol) to Ka (M⁻¹)."""
    return math.exp(-dG_kJ / RT)

def dG_from_Ka(Ka):
    """Convert Ka (M⁻¹) to ΔG (kJ/mol)."""
    return -RT * math.log(Ka)

# ============================================================
# PHASE 1: SANITY CHECKS ON PARAMETERS
# ============================================================
print("=" * 60)
print("PHASE 1: PARAMETER SANITY CHECKS")
print("=" * 60)

gates = []

# Gate G1: k_desolv_eq in expected range
assert 1.5 <= K_DESOLV_EQ <= 5.0, f"k_desolv_eq={K_DESOLV_EQ} out of range [1.5, 5.0]"
print(f"[PASS] k_desolv_eq = {K_DESOLV_EQ:.1f} kJ/mol  (range 1.5–5.0)")

# Gate G2: k_desolv_ax > k_desolv_eq (axial more costly)
assert K_DESOLV_AX > K_DESOLV_EQ, "k_desolv_ax must exceed k_desolv_eq"
print(f"[PASS] k_desolv_ax = {K_DESOLV_AX:.1f} > k_desolv_eq = {K_DESOLV_EQ:.1f}")

# Gate G3: eps_CH_pi in expected range
assert -4.0 <= EPS_CH_PI <= -1.0, f"eps_CH_pi={EPS_CH_PI} out of range [-4.0, -1.0]"
print(f"[PASS] eps_CH_pi   = {EPS_CH_PI:.1f} kJ/mol  (range -4.0 to -1.0)")

# Gate G4: eps_HB in expected range
assert -8.0 <= EPS_HB <= -3.0, f"eps_HB={EPS_HB} out of range"
print(f"[PASS] eps_HB      = {EPS_HB:.1f} kJ/mol  (range -8.0 to -3.0)")

# Gate G5: effective HB after beta_context
eff_HB = EPS_HB * BETA_CONTEXT
print(f"[INFO] Effective HB energy (×β_context): {eff_HB:.2f} kJ/mol")
assert -4.0 <= eff_HB <= -1.0, f"Effective HB {eff_HB} out of expected range"
print(f"[PASS] Effective HB in range [-4.0, -1.0]")

print()

# ============================================================
# PHASE 2: PREDICTIONS
# ============================================================
print("=" * 60)
print("PHASE 2: PREDICTIONS vs OBSERVED")
print("=" * 60)
print(f"{'System':<38} {'ΔG_pred':>8} {'ΔG_obs':>8} {'Resid':>7} {'Ka_pred':>10} {'Ka_obs':>10} {'Conf'}")
print("-" * 95)

results = []
for sys in ALL_SYSTEMS:
    dG_pred = score(sys["n_HB"], sys["n_CHP"], sys["buried_ohs"])
    dG_obs  = sys["dG_obs"]
    Ka_pred = Ka_from_dG(dG_pred)
    Ka_obs  = sys["Ka_obs"]
    resid   = dG_pred - dG_obs

    results.append({
        "name":     sys["name"],
        "dG_pred":  dG_pred,
        "dG_obs":   dG_obs,
        "resid":    resid,
        "Ka_pred":  Ka_pred,
        "Ka_obs":   Ka_obs,
        "conf":     sys["confidence"],
        "n_HB":     sys["n_HB"],
        "n_CHP":    sys["n_CHP"],
        "source":   sys["source"],
    })
    print(f"{sys['name']:<38} {dG_pred:>8.1f} {dG_obs:>8.1f} {resid:>+7.1f} {Ka_pred:>10.0f} {Ka_obs:>10} {sys['confidence']}")

print()

# ============================================================
# STATISTICS
# ============================================================
print("=" * 60)
print("STATISTICS")
print("=" * 60)

# Separate by confidence
high_med = [r for r in results if r["conf"] in ("HIGH", "MEDIUM")]
high_only = [r for r in results if r["conf"] == "HIGH"]

def stats(rset, label):
    n = len(rset)
    if n < 2:
        print(f"{label}: n={n}, insufficient for R²")
        return
    obs   = [r["dG_obs"]  for r in rset]
    pred  = [r["dG_pred"] for r in rset]
    resids = [r["resid"]  for r in rset]

    mean_obs  = sum(obs) / n
    ss_tot    = sum((o - mean_obs)**2 for o in obs)
    ss_res    = sum(r**2 for r in resids)
    r2        = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    mae       = sum(abs(r) for r in resids) / n
    rmse      = math.sqrt(ss_res / n)

    print(f"{label} (n={n}):")
    print(f"  R²   = {r2:.3f}")
    print(f"  MAE  = {mae:.2f} kJ/mol")
    print(f"  RMSE = {rmse:.2f} kJ/mol")
    return r2, mae, rmse

stats(high_only, "HIGH confidence only")
print()
stats(high_med,  "HIGH + MEDIUM confidence")

print()

# ============================================================
# SELECTIVITY RATIOS (key qualitative tests)
# ============================================================
print("=" * 60)
print("SELECTIVITY RATIOS")
print("=" * 60)

def get(name):
    for r in results:
        if r["name"] == name:
            return r
    return None

conA_man = get("ConA + αMeMan")
conA_glu = get("ConA + αMeGlu")
davis_glu = get("Davis hexaurea + β-D-Glu")
davis_gal = get("Davis hexaurea + β-D-Gal")
davis_man = get("Davis hexaurea + β-D-Man")
davis_2dg = get("Davis hexaurea + 2-deoxy-D-Glu")

print("ConA selectivity (should be Man > Glu):")
if conA_man and conA_glu:
    pred_ratio  = conA_man["Ka_pred"] / conA_glu["Ka_pred"]
    obs_ratio   = conA_man["Ka_obs"]  / conA_glu["Ka_obs"]
    print(f"  Predicted Man/Glu ratio: {pred_ratio:.1f}×  |  Observed: {obs_ratio:.1f}×")
    print(f"  {'PASS ✓' if pred_ratio > 1.0 else 'FAIL ✗'} — Man predicted {'>' if pred_ratio > 1 else '<'} Glu")

print()
print("Davis selectivity (should be Glu >> Man ≈ Gal — OPPOSITE from ConA):")
if davis_glu and davis_gal and davis_man:
    pred_glu_gal = davis_glu["Ka_pred"] / davis_gal["Ka_pred"]
    pred_glu_man = davis_glu["Ka_pred"] / davis_man["Ka_pred"]
    obs_glu_gal  = davis_glu["Ka_obs"]  / davis_gal["Ka_obs"]
    obs_glu_man  = davis_glu["Ka_obs"]  / davis_man["Ka_obs"]
    print(f"  Predicted Glu/Gal: {pred_glu_gal:.0f}×  |  Observed: {obs_glu_gal:.0f}×")
    print(f"  Predicted Glu/Man: {pred_glu_man:.0f}×  |  Observed: {obs_glu_man:.0f}×")
    scaffold_pass = pred_glu_gal > 1.0 and pred_glu_man > 1.0
    print(f"  {'PASS ✓' if scaffold_pass else 'FAIL ✗'} — scaffold-independence test")

print()
print("Davis 2dGlu vs Glu (C2-OH contribution):")
if davis_glu and davis_2dg:
    pred_ratio = davis_glu["Ka_pred"] / davis_2dg["Ka_pred"]
    obs_ratio  = davis_glu["Ka_obs"]  / davis_2dg["Ka_obs"]
    print(f"  Predicted Glu/2dGlu: {pred_ratio:.1f}×  |  Observed: {obs_ratio:.1f}×")

print()

# ============================================================
# CSV OUTPUT
# ============================================================
outfile = "/home/claude/glycan_scorer/phase2_predictions.csv"
with open(outfile, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["name","n_HB","n_CHP","dG_pred","dG_obs","resid","Ka_pred","Ka_obs","conf","source"])
    writer.writeheader()
    for r in results:
        writer.writerow({k: (f"{r[k]:.2f}" if isinstance(r[k], float) else r[k]) for k in writer.fieldnames})

print(f"Results written to: {outfile}")
print()
print("=" * 60)
print("PHASE 1+2 COMPLETE")
print("=" * 60)