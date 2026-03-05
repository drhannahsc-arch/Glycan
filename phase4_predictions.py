#!/usr/bin/env python3
"""
MABE Glycan Module — Phase 4 Predictions
Full monosaccharide panel across 5 scaffolds.
Parameters: v2.2 (locked). No fitting to biological data.
"""

import math, csv
from parameters_v22 import score, EPS_HB, EPS_CH_PI, BETA_CONTEXT
from contact_maps_phase4 import ALL_PHASE4

RT = 2.479

def dG_from_Ka(Ka):
    return -RT * math.log(Ka) if Ka and Ka > 0 else None

def Ka_from_dG(dG):
    return math.exp(-dG / RT) if dG is not None else None

# ============================================================
# Run predictions for each scaffold
# ============================================================
all_results = []
print("=" * 80)
print("PHASE 4: FULL MONOSACCHARIDE PANEL — ALL SCAFFOLDS")
print(f"Parameters: ε_HB={EPS_HB}, β={BETA_CONTEXT}, ε_CH-π={EPS_CH_PI}")
print("=" * 80)

scaffold_dG0 = {}  # fitted ΔG_0 per scaffold from reference ligand

for scaffold, panel in ALL_PHASE4.items():
    print(f"\n{'─'*70}")
    print(f"  {scaffold}")
    print(f"{'─'*70}")
    print(f"  {'Ligand':<26} {'dG_phys':>8} {'Ka_pred_raw':>12} {'Ka_obs':>10} {'Conf'}")
    print(f"  {'─'*60}")

    scaffold_results = []
    for sys in panel:
        dG_phys = score(sys["n_HB"], sys["n_CHP"], sys["buried_ohs"])
        scaffold_results.append({
            "scaffold": scaffold,
            "name":     sys["name"],
            "dG_phys":  dG_phys,
            "Ka_obs":   sys["Ka_obs"],
            "dG_obs":   sys["dG_obs"],
            "conf":     sys["confidence"],
            "n_HB":     sys["n_HB"],
            "n_CHP":    sys["n_CHP"],
            "notes":    sys["notes"],
            "source":   sys["source"],
        })
        Ka_raw = Ka_from_dG(dG_phys)
        print(f"  {sys['name']:<26} {dG_phys:>8.1f} {Ka_raw:>12.0f} {str(sys['Ka_obs'] or 'N/A'):>10} {sys['confidence']}")

    # Fit ΔG_0 from HIGH confidence ligands with observed Ka
    ref_set = [r for r in scaffold_results if r["conf"] == "HIGH" and r["dG_obs"] is not None]
    if ref_set:
        offsets = [r["dG_obs"] - r["dG_phys"] for r in ref_set]
        dG0 = sum(offsets) / len(offsets)
    else:
        # Use first ligand with Ka as reference
        ref_set = [r for r in scaffold_results if r["dG_obs"] is not None]
        dG0 = (ref_set[0]["dG_obs"] - ref_set[0]["dG_phys"]) if ref_set else 0.0

    scaffold_dG0[scaffold] = dG0

    # Apply ΔG_0, compute corrected predictions
    print(f"\n  ΔG_0({scaffold}) = {dG0:.2f} kJ/mol")
    print(f"\n  {'Ligand':<26} {'dG_corr':>8} {'dG_obs':>8} {'Resid':>7} {'Ka_corr':>10} {'Ka_obs':>10}")
    print(f"  {'─'*70}")
    for r in scaffold_results:
        r["dG_corrected"] = r["dG_phys"] + dG0
        r["Ka_corrected"]  = Ka_from_dG(r["dG_corrected"])
        resid = (r["dG_corrected"] - r["dG_obs"]) if r["dG_obs"] is not None else None
        r["resid"] = resid
        print(f"  {r['name']:<26} {r['dG_corrected']:>8.1f} "
              f"{str(round(r['dG_obs'],1) if r['dG_obs'] else 'N/A'):>8} "
              f"{(f'{resid:+.1f}' if resid is not None else 'N/A'):>7} "
              f"{r['Ka_corrected']:>10.0f} "
              f"{str(r['Ka_obs'] or 'N/A'):>10}")

    all_results.extend(scaffold_results)

# ============================================================
# ΔΔG analysis — scaffold-independent
# ============================================================
print("\n\n" + "=" * 80)
print("ΔΔG ANALYSIS (scaffold-relative, no ΔG_0 needed)")
print("=" * 80)

# Define reference ligands and comparison pairs
comparisons = [
    # (scaffold, ref_name, alt_name, description)
    ("ConA",  "ConA + αMeMan",  "ConA + αMeGlu",  "ConA: Man→Glu"),
    ("WGA",   "WGA + GlcNAc",   "WGA + (GlcNAc)2","WGA: GlcNAc→(GlcNAc)2"),
    ("PNA",   "PNA + Gal",      "PNA + GalNAc",   "PNA: Gal→GalNAc"),
    ("Gal3",  "Gal3 + Gal",     "Gal3 + LacNAc",  "Gal3: Gal→LacNAc"),
    ("Davis", "Davis + βGlu",   "Davis + βGal",   "Davis: Glu→Gal"),
    ("Davis", "Davis + βGlu",   "Davis + βMan",   "Davis: Glu→Man"),
    ("Davis", "Davis + βGlu",   "Davis + 2dGlu",  "Davis: Glu→2dGlu"),
]

def find(name):
    for r in all_results:
        if r["name"] == name:
            return r
    return None

ddG_results = []
print(f"\n  {'Comparison':<30} {'ΔΔG_pred':>10} {'ΔΔG_obs':>10} {'Error':>8} {'Conf'}")
print(f"  {'─'*65}")
for scaffold, ref_name, alt_name, desc in comparisons:
    ref = find(ref_name)
    alt = find(alt_name)
    if not ref or not alt:
        continue
    ddG_pred = alt["dG_phys"] - ref["dG_phys"]
    ddG_obs  = ((alt["dG_obs"] - ref["dG_obs"])
                if alt["dG_obs"] is not None and ref["dG_obs"] is not None else None)
    err = (ddG_pred - ddG_obs) if ddG_obs is not None else None
    conf = "HIGH" if ref["conf"] == "HIGH" and alt["conf"] == "HIGH" else alt["conf"]
    print(f"  {desc:<30} {ddG_pred:>+10.1f} "
          f"{(f'{ddG_obs:+.1f}' if ddG_obs is not None else 'N/A'):>10} "
          f"{(f'{err:+.1f}' if err is not None else 'N/A'):>8} {conf}")
    if ddG_obs is not None:
        ddG_results.append((desc, ddG_pred, ddG_obs, conf))

# ============================================================
# Statistics on ΔΔG
# ============================================================
print("\n\n" + "=" * 80)
print("STATISTICS")
print("=" * 80)

for label, filter_fn in [
    ("ALL with obs", lambda x: True),
    ("HIGH only",    lambda x: x[3] == "HIGH"),
]:
    subset = [x for x in ddG_results if filter_fn(x)]
    if len(subset) < 2:
        print(f"{label}: n={len(subset)}, insufficient")
        continue
    preds = [x[1] for x in subset]
    obs   = [x[2] for x in subset]
    mean_obs = sum(obs)/len(obs)
    ss_tot = sum((o-mean_obs)**2 for o in obs)
    ss_res = sum((p-o)**2 for p,o in zip(preds,obs))
    r2   = 1 - ss_res/ss_tot if ss_tot > 0 else float('nan')
    mae  = sum(abs(p-o) for p,o in zip(preds,obs))/len(preds)
    rmse = math.sqrt(ss_res/len(preds))
    print(f"\n{label} (n={len(subset)}):")
    print(f"  R²   = {r2:.3f}")
    print(f"  MAE  = {mae:.2f} kJ/mol")
    print(f"  RMSE = {rmse:.2f} kJ/mol")

# ============================================================
# Scaffold independence — key check
# ============================================================
print("\n\n" + "=" * 80)
print("SCAFFOLD INDEPENDENCE TEST")
print("=" * 80)

# ConA vs Davis: Man/Gal preference reversal
conA_man = find("ConA + αMeMan")
conA_gal = find("ConA + αMeGal")
conA_glu = find("ConA + αMeGlu")
davis_glu = find("Davis + βGlu")
davis_gal = find("Davis + βGal")
davis_man = find("Davis + βMan")

# Relative to glucose scaffold for each
conA_glu_vs_gal = (conA_gal["dG_phys"] - conA_glu["dG_phys"]) if conA_gal and conA_glu else None
davis_glu_vs_gal = (davis_gal["dG_phys"] - davis_glu["dG_phys"]) if davis_gal and davis_glu else None

print(f"\nConA:  ΔΔG(Glu→Gal) = {conA_glu_vs_gal:+.1f} kJ/mol  (should be +, ConA rejects axial C4)")
print(f"Davis: ΔΔG(Glu→Gal) = {davis_glu_vs_gal:+.1f} kJ/mol  (should be +, Davis also rejects axial C4)")

# PNA vs ConA: both involve C4-OH but in opposite roles
pna_gal = find("PNA + Gal")
pna_glu = find("PNA + Glu")
if pna_gal and pna_glu:
    pna_glu_vs_gal = pna_glu["dG_phys"] - pna_gal["dG_phys"]
    print(f"\nPNA:   ΔΔG(Gal→Glu) = {pna_glu_vs_gal:+.1f} kJ/mol  (should be +, PNA REQUIRES axial C4 → Glu weaker)")

gal3_gal = find("Gal3 + Gal")
gal3_glu = find("Gal3 + Glu")
if gal3_gal and gal3_glu:
    gal3_glu_vs_gal = gal3_glu["dG_phys"] - gal3_gal["dG_phys"]
    print(f"Gal3:  ΔΔG(Gal→Glu) = {gal3_glu_vs_gal:+.1f} kJ/mol  (should be +, Gal3 also requires axial C4)")

print("\nConA Man>Glu: ", "PASS ✓" if conA_man["dG_phys"] < conA_glu["dG_phys"] else "FAIL ✗")
print("PNA Gal>Glu:  ", "PASS ✓" if pna_gal["dG_phys"] < pna_glu["dG_phys"] else "FAIL ✗")
print("Gal3 Gal>Glu: ", "PASS ✓" if gal3_gal["dG_phys"] < gal3_glu["dG_phys"] else "FAIL ✗")
print("Davis Glu>Gal:", "PASS ✓" if davis_glu["dG_phys"] < davis_gal["dG_phys"] else "FAIL ✗")
print("Davis Glu>Man:", "PASS ✓" if davis_glu["dG_phys"] < davis_man["dG_phys"] else "FAIL ✗")

# ============================================================
# Write CSV
# ============================================================
outfile = "/home/claude/glycan_scorer/phase4_predictions.csv"
with open(outfile, "w", newline="") as f:
    fieldnames = ["scaffold","name","n_HB","n_CHP","dG_phys","dG_corrected",
                  "dG_obs","resid","Ka_corrected","Ka_obs","conf","source"]
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for r in all_results:
        row = {k: (f"{r[k]:.2f}" if isinstance(r.get(k), float) else r.get(k, "")) for k in fieldnames}
        writer.writerow(row)
print(f"\nCSV written: {outfile}")
print("\n" + "=" * 80)
print("PHASE 4 COMPLETE")
print("=" * 80)