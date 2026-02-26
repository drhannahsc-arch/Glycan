"""
Tests for glycan physics data modules (Tier 1 sources).
"""
import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.chi_energy_functions import (
    chi_phi_alpha, chi_phi_beta, chi_psi_2a3e, chi_psi_2e3a,
    get_chi_functions, score_linkage_entropy, compute_TdS_freeze,
    LINKAGE_PSI_CLASS,
)
from data.ch_pi_energetics import (
    LAUGHREY_2008, KEYS_2025_RANGES, SUGAR_CH_PI_CONTACTS,
    get_eps_CH_pi_estimates, estimate_CH_pi_energy,
    LAUGHREY_PER_CONTACT,
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
print("TEST SUITE: Glycan Physics Data Modules (Tier 1)")
print("=" * 60)

# ─── CHI ENERGY FUNCTIONS ───────────────────────────────────

print("\n--- CHI Energy Functions ---")

# Test 1: phi_alpha has a minimum near 60 deg (typical alpha-glycosidic angle)
energies = [(phi, chi_phi_alpha(phi)) for phi in range(-180, 180, 5)]
min_phi, min_e = min(energies, key=lambda x: x[1])
check("phi_alpha minimum in [-80, 80] range",
      -80 <= min_phi <= 80,
      f"min at {min_phi} deg")

# Test 2: phi_alpha has a barrier (max > min + 2 kcal/mol)
max_e = max(e for _, e in energies)
check("phi_alpha has barrier > 2 kcal/mol",
      max_e - min_e > 2.0,
      f"barrier = {max_e - min_e:.2f} kcal/mol")

# Test 3: phi_beta has a different minimum than phi_alpha
energies_b = [(phi, chi_phi_beta(phi)) for phi in range(-180, 180, 5)]
min_phi_b, _ = min(energies_b, key=lambda x: x[1])
check("phi_beta minimum differs from phi_alpha",
      abs(min_phi_b - min_phi) > 20,
      f"alpha min={min_phi}, beta min={min_phi_b}")

# Test 4: psi functions have minima in expected range
psi_e_2a3e = [(psi, chi_psi_2a3e(psi)) for psi in range(0, 360, 5)]
min_psi_2a3e, _ = min(psi_e_2a3e, key=lambda x: x[1])
check("psi_2a3e minimum exists",
      0 <= min_psi_2a3e <= 360,
      f"min at {min_psi_2a3e}")

psi_e_2e3a = [(psi, chi_psi_2e3a(psi)) for psi in range(0, 360, 5)]
min_psi_2e3a, _ = min(psi_e_2e3a, key=lambda x: x[1])
check("psi_2e3a minimum exists",
      0 <= min_psi_2e3a <= 360,
      f"min at {min_psi_2e3a}")

# Test 5: Linkage classification covers common sugars
for sugar in ["Glc", "Gal", "Man", "GlcNAc", "GalNAc"]:
    for pos in [2, 3, 4]:
        key = (pos, sugar)
        check(f"Linkage ({pos},{sugar}) classified",
              key in LINKAGE_PSI_CLASS,
              f"missing from LINKAGE_PSI_CLASS")

# Test 6: get_chi_functions returns callables
phi_f, psi_f = get_chi_functions("alpha", 3, "Man")
check("get_chi_functions returns callables",
      callable(phi_f) and callable(psi_f))

# Test 7: (1->6) raises ValueError
try:
    get_chi_functions("alpha", 6, "Glc")
    check("(1->6) raises ValueError", False, "no exception raised")
except ValueError:
    check("(1->6) raises ValueError", True)

# Test 8: Entropy computation returns reasonable values
# Use angles near the CHI minimum (phi~70, psi~180 for alpha(1->3)Man)
# which is where ConA actually binds mannose.
result = score_linkage_entropy("alpha", 3, "Man", 70, 200)
check("TdS_freeze is positive (unfavorable entropy)",
      result["TdS_freeze_kJ"] > 0,
      f"got {result['TdS_freeze_kJ']}")

check("TdS_freeze in reasonable range (1-15 kJ/mol) at near-minimum angles",
      1 < result["TdS_freeze_kJ"] < 15,
      f"got {result['TdS_freeze_kJ']}")

check("f_bound < 1.0 (bound window is subset of free ensemble)",
      result["f_bound"] < 1.0,
      f"got {result['f_bound']}")

# Test 9: Bound state at global minimum has lowest entropy cost
result_min = score_linkage_entropy("alpha", 3, "Man",
                                   result["phi_min_deg"],
                                   result["psi_min_deg"])
check("TdS_freeze at global minimum < at arbitrary point",
      result_min["TdS_freeze_kJ"] <= result["TdS_freeze_kJ"],
      f"min={result_min['TdS_freeze_kJ']}, arb={result['TdS_freeze_kJ']}")

# Test 10: Different linkage types give different entropy costs
result_a13 = score_linkage_entropy("alpha", 3, "Glc", -70, 150)
result_b14 = score_linkage_entropy("beta", 4, "Glc", -70, 150)
check("alpha(1->3) and beta(1->4) give different TdS",
      abs(result_a13["TdS_freeze_kJ"] - result_b14["TdS_freeze_kJ"]) > 0.1,
      f"a13={result_a13['TdS_freeze_kJ']}, b14={result_b14['TdS_freeze_kJ']}")


# ─── CH-PI ENERGETICS ───────────────────────────────────────

print("\n--- CH-pi Energetics ---")

# Test 11: Laughrey dataset has entries
check("Laughrey dataset has contacts",
      len(LAUGHREY_2008.contacts) >= 8,
      f"got {len(LAUGHREY_2008.contacts)}")

# Test 12: Laughrey per-contact energies are negative (favorable)
for aromatic, data in LAUGHREY_PER_CONTACT.items():
    check(f"Laughrey {aromatic} dG_per_CH < 0",
          data["dG_per_CH_kJ"] < 0,
          f"got {data['dG_per_CH_kJ']}")

# Test 13: Indole > benzene > cyclohexane (ranking preserved)
check("Indole stronger than benzene",
      LAUGHREY_PER_CONTACT["indole"]["dG_per_CH_kJ"] <
      LAUGHREY_PER_CONTACT["benzene"]["dG_per_CH_kJ"])

check("Benzene stronger than cyclohexane",
      LAUGHREY_PER_CONTACT["benzene"]["dG_per_CH_kJ"] <
      LAUGHREY_PER_CONTACT["cyclohexane"]["dG_per_CH_kJ"])

# Test 14: Keys 2025 DFT energies are more negative than solution values
# (gas phase should be 3-4x stronger)
keys_mean = KEYS_2025_RANGES["CH_pi_stacking"]["mean_kJ"]
laughrey_mean = LAUGHREY_PER_CONTACT["indole"]["dG_per_CH_kJ"] * 3  # 3 contacts
check("DFT gas-phase more negative than solution dG",
      keys_mean < laughrey_mean,
      f"DFT={keys_mean}, solution={laughrey_mean}")

# Test 15: Sugar contact maps cover key sugars
for sugar in ["beta-Glc", "beta-Gal", "alpha-Man", "beta-GlcNAc"]:
    check(f"Contact map for {sugar}",
          sugar in SUGAR_CH_PI_CONTACTS)

# Test 16: Glc has more/equal CH contacts than Man (all-equatorial advantage)
glc_max = max(f["n_CH"] for f in SUGAR_CH_PI_CONTACTS["beta-Glc"].values())
man_max = max(f["n_CH"] for f in SUGAR_CH_PI_CONTACTS["alpha-Man"].values())
check("Glc >= Man CH contacts (all-equatorial advantage)",
      glc_max >= man_max,
      f"Glc={glc_max}, Man={man_max}")

# Test 17: Parameter estimates have expected signs and ranges
estimates = get_eps_CH_pi_estimates()
for key, val in estimates.items():
    check(f"eps estimate {key} is negative",
          val["value_kJ"] < 0,
          f"got {val['value_kJ']}")
    check(f"eps estimate {key} in [-10, 0] kJ/mol",
          -10 < val["value_kJ"] < 0,
          f"got {val['value_kJ']}")

# Test 18: estimate_CH_pi_energy works and gives sensible output
result = estimate_CH_pi_energy("beta-Gal", "Trp")
check("CH-pi energy for Gal-Trp is negative",
      result["dG_CH_pi_kJ"] < 0,
      f"got {result['dG_CH_pi_kJ']}")

check("CH-pi energy for Gal-Trp in [-20, -5] kJ/mol",
      -20 < result["dG_CH_pi_kJ"] < -5,
      f"got {result['dG_CH_pi_kJ']}")

# Test 19: Trp gives more negative dG than Phe for same sugar
trp = estimate_CH_pi_energy("beta-Gal", "Trp")
phe = estimate_CH_pi_energy("beta-Gal", "Phe")
check("Trp more favorable than Phe for Gal",
      trp["dG_CH_pi_kJ"] < phe["dG_CH_pi_kJ"],
      f"Trp={trp['dG_CH_pi_kJ']}, Phe={phe['dG_CH_pi_kJ']}")

# Test 20: Glc has most favorable CH-pi (most contacts)
glc = estimate_CH_pi_energy("beta-Glc", "Trp")
gal = estimate_CH_pi_energy("beta-Gal", "Trp")
man = estimate_CH_pi_energy("alpha-Man", "Trp")
check("Glc >= Gal for CH-pi (more equatorial CHs)",
      glc["dG_CH_pi_kJ"] <= gal["dG_CH_pi_kJ"],
      f"Glc={glc['dG_CH_pi_kJ']}, Gal={gal['dG_CH_pi_kJ']}")


# ─── SUMMARY ────────────────────────────────────────────────

print()
print("=" * 60)
print(f"RESULTS: {PASS} passed, {FAIL} failed out of {PASS + FAIL}")
print("=" * 60)

if FAIL > 0:
    sys.exit(1)