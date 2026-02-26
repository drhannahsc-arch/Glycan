"""
Tests for PDB lectin reference data and parameter integration.
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.pdb_lectin_reference import (
    LECTIN_REFERENCE_SET, LECTIN_MULTIVALENCY, WATER_BRIDGE_STATISTICS,
    BOUND_STATE_TORSIONS, get_water_bridge_norm, get_inter_site_distances,
)
from glycan_params_integration import (
    compute_params_from_tier1, GlycanParams,
    score_monosaccharide, validate_against_reference,
)

PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS: {name}")
    else:
        FAIL += 1
        print(f"  FAIL: {name} — {detail}")


print("=" * 60)
print("TEST SUITE: PDB Reference + Integration (Tier 1)")
print("=" * 60)


# ─── PDB REFERENCE DATA ────────────────────────────────────

print("\n--- PDB Lectin Reference Set ---")

check("Reference set has entries",
      len(LECTIN_REFERENCE_SET) >= 7,
      f"got {len(LECTIN_REFERENCE_SET)}")

for ls in LECTIN_REFERENCE_SET:
    check(f"{ls.pdb_id} has valid PDB ID",
          len(ls.pdb_id) >= 3 and len(ls.pdb_id) <= 4,
          f"'{ls.pdb_id}' is {len(ls.pdb_id)} chars")

for ls in LECTIN_REFERENCE_SET:
    check(f"{ls.pdb_id} resolution <= 2.0 A",
          ls.resolution_A <= 2.0,
          f"got {ls.resolution_A}")

pdb_ids = [ls.pdb_id for ls in LECTIN_REFERENCE_SET]
for required in ["5CNA", "3ZSJ", "2UVO"]:
    check(f"Required PDB {required} present",
          required in pdb_ids)

cona = [ls for ls in LECTIN_REFERENCE_SET if ls.pdb_id == "5CNA"][0]
check("ConA has ITC dG",
      cona.sites[0].dG_exp_kJ is not None and cona.sites[0].dG_exp_kJ < 0)
check("ConA has metal ions",
      len(cona.sites[0].metal_ions) >= 2)


print("\n--- Water Bridge Statistics ---")

check("Water bridge data has entries",
      len(WATER_BRIDGE_STATISTICS) >= 5)

for pdb, data in WATER_BRIDGE_STATISTICS.items():
    check(f"{pdb} water density in range",
          0.003 < data["density"] < 0.03,
          f"got {data['density']}")

p9 = get_water_bridge_norm()
check("P9 mean > 0 and < 0.02",
      0 < p9["P9_mean_waters_per_A2"] < 0.02,
      f"got {p9['P9_mean_waters_per_A2']}")


print("\n--- Multivalency Geometry ---")

check("Multivalency data has entries",
      len(LECTIN_MULTIVALENCY) >= 5)


print("\n--- Parameter Integration ---")

params = compute_params_from_tier1()

check("P4 computed (not 0)",
      params.P4_eps_glycosidic_kJ != 0.0,
      f"got {params.P4_eps_glycosidic_kJ}")
check("P4 positive (entropy cost)",
      params.P4_eps_glycosidic_kJ > 0,
      f"got {params.P4_eps_glycosidic_kJ}")
check("P4 in 2-40 kJ/mol (avg across diverse linkages)",
      2 < params.P4_eps_glycosidic_kJ < 40,
      f"got {params.P4_eps_glycosidic_kJ}")
check("P6 negative",
      params.P6_eps_CH_pi_pyranose_kJ < 0)
check("P9 from PDB",
      params.P9_n_water_norm > 0.005)
check("P12 computed",
      params.P12_d_optimal_multi_A > 10)
check("All 12 params have provenance",
      len(params.provenance()) == 12)


print("\n--- Monosaccharide Scoring ---")

r = score_monosaccharide("alpha-Man", "5CNA")
check("ConA scoring succeeds", "error" not in r)
check("ConA dG negative", r["dG_predicted_kJ"] < 0)
check("ConA has components", "components" in r)
check("ConA exp dG present", r["dG_exp_kJ"] is not None and r["dG_exp_kJ"] < 0)

r2 = score_monosaccharide("beta-Gal", "3ZSJ")
check("Galectin-3 succeeds", "error" not in r2)
check("Galectin-3 CH-pi < 0", r2["components"]["CH_pi_kJ"] < 0)
check("Galectin-3 no metal", r2["components"]["metal_coord_kJ"] == 0.0)

r3 = score_monosaccharide("alpha-Man", "1SL5")
check("DC-SIGN metal < 0", r3["components"]["metal_coord_kJ"] < 0)
check("DC-SIGN no CH-pi", r3["components"]["CH_pi_kJ"] == 0.0)


print("\n--- Validation ---")

val = validate_against_reference()
check("Validation scored >= 3 pairs", val["n_pairs"] >= 3)
check("MAE finite", val["MAE_kJ"] < 100)


print()
print("=" * 60)
print(f"RESULTS: {PASS} passed, {FAIL} failed out of {PASS + FAIL}")
print("=" * 60)

if FAIL > 0:
    sys.exit(1)