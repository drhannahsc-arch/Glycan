"""
CH-pi interaction energetics for carbohydrate-aromatic contacts.

Data compiled from multiple sources for MABE glycan parameter P6 (eps_CH_pi_pyranose)
and P7 (eps_CH_pi_furanose) calibration.

SOURCE HIERARCHY (by parameter isolation quality):
  1. Laughrey 2008 (JACS) — absolute dG per contact in designed system
  2. Jimenez-Moreno 2015 (Chem. Sci.) — 117-complex relative ddG in water
  3. Keys et al. 2025 (Chem. Sci.) — 774 DFT contacts from PDB structures
  4. Hudson et al. 2015 (JACS) — PDB statistics + NMR Hammett analysis
  5. Fernandez-Alonso 2005/2009 — calorimetric dH in organic solvent

IMPORTANT: Gas-phase/DFT interaction energies (3-8 kcal/mol) are 3-4x larger
than solution dG values (0.5-2 kcal/mol) due to desolvation penalty and entropy.
MABE needs SOLUTION dG, not gas-phase interaction energy.

Units: All energies stored in kJ/mol unless noted.
"""

from dataclasses import dataclass, field
from typing import List, Optional

KCAL_TO_KJ = 4.184


@dataclass
class CHPiContact:
    """A single CH-pi interaction measurement."""
    sugar: str                # e.g., "beta-Gal", "Glc", "MeAlphaMan"
    aromatic: str             # e.g., "Trp", "Tyr", "Phe", "indole", "benzene"
    value_kJ: float           # Interaction energy or dG in kJ/mol
    value_type: str           # "dG_solution", "dH_solution", "E_DFT_gas", "ddG_relative"
    n_CH_contacts: Optional[int] = None  # Number of CH groups engaged
    face: Optional[str] = None           # "alpha", "beta", or None
    source: str = ""
    doi: str = ""
    notes: str = ""


@dataclass
class CHPiDataset:
    """Collection of CH-pi measurements from one source."""
    name: str
    description: str
    contacts: List[CHPiContact] = field(default_factory=list)
    data_type: str = ""       # "solution_dG", "DFT_gas", "calorimetric_dH"
    solvent: str = "water"
    temperature_K: float = 298.15


# =====================================================================
# DATASET 1: Laughrey et al. 2008 (JACS 130, 14625-14633)
# Designed beta-hairpin peptide system. Sugar at diagonal position
# from aromatic. 1:1 geometry. dG from thermal denaturation.
# =====================================================================

LAUGHREY_2008 = CHPiDataset(
    name="Laughrey 2008",
    description="Carbohydrate-pi interactions in designed beta-hairpin peptides. "
                "Absolute dG per CH-pi contact from thermal denaturation.",
    data_type="solution_dG",
    solvent="water (phosphate buffer pH 7)",
    contacts=[
        # Ac4Glc (peracetylated glucose, enhanced CH polarization)
        CHPiContact("Ac4Glc", "indole",      -3.35, "dG_solution", n_CH_contacts=3,
                    source="Laughrey 2008", doi="10.1021/ja804166f",
                    notes="Most favorable pair. -0.8 kcal/mol"),
        CHPiContact("Ac4Glc", "naphthalene",  -2.93, "dG_solution", n_CH_contacts=3,
                    source="Laughrey 2008", doi="10.1021/ja804166f",
                    notes="-0.7 kcal/mol"),
        CHPiContact("Ac4Glc", "benzene",      -2.51, "dG_solution", n_CH_contacts=3,
                    source="Laughrey 2008", doi="10.1021/ja804166f",
                    notes="-0.6 kcal/mol"),
        CHPiContact("Ac4Glc", "cyclohexane",  -0.42, "dG_solution", n_CH_contacts=3,
                    source="Laughrey 2008", doi="10.1021/ja804166f",
                    notes="Control: -0.1 kcal/mol. Hydrophobic only, no pi."),
        # Glc (unmodified glucose)
        CHPiContact("Glc", "indole",           -2.93, "dG_solution", n_CH_contacts=3,
                    source="Laughrey 2008", doi="10.1021/ja804166f",
                    notes="-0.7 kcal/mol"),
        CHPiContact("Glc", "naphthalene",      -2.51, "dG_solution", n_CH_contacts=3,
                    source="Laughrey 2008", doi="10.1021/ja804166f",
                    notes="-0.6 kcal/mol"),
        CHPiContact("Glc", "benzene",          -2.09, "dG_solution", n_CH_contacts=3,
                    source="Laughrey 2008", doi="10.1021/ja804166f",
                    notes="-0.5 kcal/mol"),
        # Me4Glc (permethylated glucose)
        CHPiContact("Me4Glc", "indole",        -2.51, "dG_solution", n_CH_contacts=3,
                    source="Laughrey 2008", doi="10.1021/ja804166f",
                    notes="-0.6 kcal/mol. Reduced vs Glc — OH->OMe reduces polarization"),
    ],
)

# Key derived values from Laughrey:
# CH-pi per contact (indole): ~(-2.93 - (-0.42)) / 3 = -0.84 kcal/mol = -3.5 kJ/mol
# CH-pi per contact (benzene): ~(-2.09 - (-0.42)) / 3 = -0.56 kcal/mol = -2.3 kJ/mol
# Pure hydrophobic (cyclohexane control): -0.14 kcal/mol per contact = -0.6 kJ/mol

LAUGHREY_PER_CONTACT = {
    "indole":      {"dG_per_CH_kJ": -3.5, "notes": "(Glc-indole - Glc-cyclohexane) / 3"},
    "naphthalene": {"dG_per_CH_kJ": -2.8, "notes": "estimated from naphthalene series"},
    "benzene":     {"dG_per_CH_kJ": -2.3, "notes": "(Glc-benzene - Glc-cyclohexane) / 3"},
    "cyclohexane": {"dG_per_CH_kJ": -0.6, "notes": "hydrophobic only, no pi component"},
}


# =====================================================================
# DATASET 2: Keys et al. 2025 (Chem. Sci. 16, 1746-1761)
# 774 galactose-aromatic DFT contacts from PDB structures.
# B3LYP-D3 interaction energies. GAS PHASE — not solution dG.
# =====================================================================

KEYS_2025_SUMMARY = CHPiDataset(
    name="Keys 2025",
    description="DFT (B3LYP-D3) interaction energies for 774 galactose-aromatic "
                "close contacts from PDB protein-carbohydrate structures. "
                "Gas-phase electronic energies, NOT solution free energies.",
    data_type="DFT_gas",
    solvent="gas phase (implicit protein environment)",
    contacts=[
        # Summary statistics by amino acid type (Table S4 in paper)
        CHPiContact("beta-Gal", "Trp_CH_pi",  -25.5, "E_DFT_gas",
                    source="Keys 2025", doi="10.1039/d4sc06246a",
                    notes="Average CH-pi stacking: -6.1 kcal/mol. n=many"),
        CHPiContact("beta-Gal", "Tyr_CH_pi",  -25.1, "E_DFT_gas",
                    source="Keys 2025", doi="10.1039/d4sc06246a",
                    notes="Tyr CH-pi slightly weaker than Trp"),
        CHPiContact("beta-Gal", "Phe_CH_pi",  -24.7, "E_DFT_gas",
                    source="Keys 2025", doi="10.1039/d4sc06246a",
                    notes="Phe CH-pi weakest of three aromatics"),
        CHPiContact("beta-Gal", "any_Hbond",  -18.4, "E_DFT_gas",
                    source="Keys 2025", doi="10.1039/d4sc06246a",
                    notes="Average H-bond: -4.4 kcal/mol"),
        CHPiContact("beta-Gal", "any_other",  -13.4, "E_DFT_gas",
                    source="Keys 2025", doi="10.1039/d4sc06246a",
                    notes="Average other contact: -3.2 kcal/mol"),
    ],
)

# Key findings from Keys 2025:
KEYS_2025_RANGES = {
    "CH_pi_stacking": {
        "mean_kcal": -6.1,
        "mean_kJ": -25.5,
        "range_kcal": (-10.1, -0.6),
        "range_kJ": (-42.3, -2.5),
        "n_contacts": 774,  # total close contacts analyzed
    },
    "H_bonding": {
        "mean_kcal": -4.4,
        "mean_kJ": -18.4,
    },
    "other_contacts": {
        "mean_kcal": -3.2,
        "mean_kJ": -13.4,
    },
    # Amino acid ranking:
    "Trp_vs_Tyr_vs_Phe": {
        "Trp": {"mean_kcal": -6.1, "notes": "bicyclic, largest pi system"},
        "Tyr": {"mean_kcal": -5.8, "notes": "unicyclic + OH, intermediate"},
        "Phe": {"mean_kcal": -5.5, "notes": "unicyclic, smallest pi system"},
    },
    # CRITICAL SCALING FACTOR: gas-phase E_int vs solution dG
    # NMR solution dG: 1-2 kcal/mol (Laughrey, Hudson)
    # DFT gas-phase:   3-8 kcal/mol (Keys, Tsuzuki)
    # Ratio: ~3-4x attenuation from gas to solution
    "gas_to_solution_factor": 0.25,  # approximate: dG_soln ~ 0.25 * E_DFT_gas
}


# =====================================================================
# DATASET 3: Jimenez-Moreno et al. 2015 (Chem. Sci. 6, 6076-6085)
# 117 carbohydrate/aromatic complexes. Relative ddG from dynamic
# combinatorial chemistry in water. Best dataset for RANKING sugars.
# =====================================================================

# Key findings (paraphrased, not full 117-entry dataset):
JIMENEZ_MORENO_2015_KEY_FINDINGS = {
    "description": "117 carbohydrate/aromatic complexes in water. "
                   "Dynamic combinatorial imine/hemiaminal approach. "
                   "Relative stabilities from NMR competition experiments.",
    "doi": "10.1039/c5sc02108a",
    "n_complexes": 117,
    "key_results": {
        "aromatic_ranking": "indole > naphthalene > benzene (consistent with Laughrey)",
        "sugar_face_preference": {
            "Glc_beta_face": "best (all-equatorial CHs for stacking)",
            "Gal_C3456_patch": "good (axial C4-OH creates electropositive CH patch)",
            "Man_alpha_face": "moderate",
        },
        "axial_OH_disruption": {
            "equatorial_OMe_on_stacking_face": "0.38-0.64 kcal/mol penalty",
            "notes": "Axial polar groups on the interacting face disrupt CH-pi",
        },
        "electrostatics_dominate": (
            "Electrostatic and polarization contributions correlate better with "
            "experimental dG than dispersion. Despite dispersion being dominant in "
            "gas phase, solution selectivity is driven by electrostatics."
        ),
    },
    # These relative ddG values need the full supplementary table from the paper.
    # Students should extract Table S1 for the complete 117-entry dataset.
    "data_status": "SUMMARY_ONLY — full 117-entry table requires paper access",
}


# =====================================================================
# DATASET 4: Hudson et al. 2015 (JACS 137, 15152-15160)
# PDB-wide survey + NMR Hammett analysis.
# =====================================================================

HUDSON_2015_SUMMARY = {
    "description": "PDB survey of carbohydrate-aromatic contacts across all "
                   "high-resolution structures. NMR chemical shift perturbation "
                   "of sugar protons by indole. Hammett analysis with substituted indoles.",
    "doi": "10.1021/jacs.5b08424",
    "key_results": {
        "Trp_prevalence": "9x higher than expected at sugar binding sites",
        "Gal_stacking_face": "C4-C5-C6 (electropositive patch from axial C4-OH)",
        "Glc_stacking_face": "All-equatorial CHs, both faces accessible",
        "Hammett": {
            "electron_rich_indoles": "stronger CH-pi (more negative delta-delta)",
            "electron_poor_indoles": "weaker CH-pi",
            "5_nitro_indole": "abolishes interaction entirely",
            "correlation": "linear free energy relationship with sigma_para",
        },
    },
    "data_status": "QUALITATIVE — NMR gives relative shifts, not absolute dG",
}


# =====================================================================
# DATASET 5: Fernandez-Alonso / Asensio calorimetry
# =====================================================================

FERNANDEZ_ALONSO_CALORIMETRY = CHPiDataset(
    name="Fernandez-Alonso 2005/2009",
    description="Calvet microcalorimetry: enthalpies of solvation of permethylated "
                "sugar derivatives in benzene. Difference between sugars = CH-pi enthalpy.",
    data_type="calorimetric_dH",
    solvent="benzene (as both solvent and pi-partner)",
    contacts=[
        CHPiContact("alpha-Me5Man", "benzene", -78.8 * KCAL_TO_KJ / 78.8,  # normalized
                    "dH_solution",
                    source="Fernandez-Alonso 2005", doi="10.1021/ja057453r",
                    notes="dH_solv = -78.8 kJ/mol in benzene"),
        CHPiContact("alpha-Me5Gal", "benzene", -88.7 * KCAL_TO_KJ / 88.7,
                    "dH_solution",
                    source="Fernandez-Alonso 2005",
                    notes="dH_solv = -88.7 kJ/mol. 9.9 kJ/mol more favorable than Man"),
        CHPiContact("beta-Me5Gal", "benzene", -88.7 * KCAL_TO_KJ / 88.7,
                    "dH_solution",
                    source="Fernandez-Alonso 2005",
                    notes="Same as alpha-Me5Gal within error"),
        CHPiContact("beta-Ac5Gal", "benzene", -132.5 * KCAL_TO_KJ / 132.5,
                    "dH_solution",
                    source="Fernandez-Alonso 2005",
                    notes="dH_solv = -132.5 kJ/mol. Acetyl groups enhance interaction"),
    ],
)

# The CH-pi face-face enthalpy isolated by Fernandez-Alonso:
# dH(beta-Me5Gal) - dH(alpha-Me5Man) = -88.7 - (-78.8) = -9.9 kJ/mol
# This is the ENTHALPY difference between Gal and Man stacking faces.
# Caveat: measured in organic solvent, not aqueous.

FERNANDEZ_ALONSO_FACE_DH = {
    "Gal_minus_Man_dH_kJ": -9.9,
    "notes": "In benzene solvent. Gal stacks 9.9 kJ/mol more favorably than Man. "
             "Consistent with Gal having more electropositive CH patch.",
}


# =====================================================================
# CONSENSUS VALUES FOR MABE PARAMETER ESTIMATION
# =====================================================================

def get_eps_CH_pi_estimates():
    """Return consensus estimates for eps_CH_pi parameters.

    These are starting values derived from the literature survey.
    The final values come from back-solve fitting against the
    Jimenez-Moreno dataset (when available) or Laughrey data.

    Returns:
        Dict with parameter estimates and confidence.
    """
    return {
        "eps_CH_pi_Trp": {
            "value_kJ": -3.5,
            "range_kJ": (-4.5, -2.5),
            "confidence": "HIGH",
            "basis": "Laughrey 2008 per-contact dG, Glc-indole minus cyclohexane control",
        },
        "eps_CH_pi_Tyr": {
            "value_kJ": -2.8,
            "range_kJ": (-3.5, -2.0),
            "confidence": "MEDIUM",
            "basis": "Interpolated from Laughrey benzene/naphthalene + Keys Tyr/Trp ratio",
        },
        "eps_CH_pi_Phe": {
            "value_kJ": -2.3,
            "range_kJ": (-3.0, -1.5),
            "confidence": "MEDIUM",
            "basis": "Laughrey benzene contact + Keys Phe/Trp ratio",
        },
        "eps_CH_pi_pyranose_avg": {
            "value_kJ": -3.0,
            "range_kJ": (-4.0, -2.0),
            "confidence": "HIGH",
            "basis": "Weighted average across Trp/Tyr/Phe prevalence in lectin sites",
            "notes": "Use when amino acid identity unknown. Trp-weighted since Trp is 9x enriched.",
        },
        "eps_CH_pi_furanose": {
            "value_kJ": -2.0,
            "range_kJ": (-3.5, -1.0),
            "confidence": "LOW",
            "basis": "Extrapolated from pyranose. Fewer CH contacts on furanose ring. "
                     "Only gas-phase QM (Tsuzuki fucose-benzene) available.",
        },
    }


# =====================================================================
# MAPPING: Sugar stereochemistry -> CH-pi contact quality
# Which face, which CHs, how many effective contacts
# =====================================================================

SUGAR_CH_PI_CONTACTS = {
    # sugar: {face: {n_CH_contacts, quality, interacting_carbons}}
    "beta-Glc": {
        "alpha_face": {"n_CH": 5, "quality": "excellent", "carbons": "C1,C2,C3,C4,C5",
                       "notes": "All-equatorial. Best stacking geometry."},
        "beta_face":  {"n_CH": 5, "quality": "excellent", "carbons": "C1,C2,C3,C4,C5"},
    },
    "alpha-Glc": {
        "alpha_face": {"n_CH": 4, "quality": "good", "carbons": "C2,C3,C4,C5",
                       "notes": "Anomeric OH blocks C1 on alpha face."},
        "beta_face":  {"n_CH": 5, "quality": "excellent", "carbons": "C1,C2,C3,C4,C5"},
    },
    "beta-Gal": {
        "alpha_face": {"n_CH": 3, "quality": "good", "carbons": "C3,C5,C6",
                       "notes": "Axial C4-OH disrupts this face."},
        "beta_face":  {"n_CH": 4, "quality": "excellent", "carbons": "C3,C4,C5,C6",
                       "notes": "C4-OH axial creates electropositive CH patch on beta face. "
                                "Preferred stacking face per Hudson 2015."},
    },
    "alpha-Man": {
        "alpha_face": {"n_CH": 3, "quality": "moderate", "carbons": "C1,C3,C5",
                       "notes": "Axial C2-OH disrupts. Fewer aligned CHs."},
        "beta_face":  {"n_CH": 4, "quality": "good", "carbons": "C2,C3,C4,C5"},
    },
    "beta-GlcNAc": {
        "alpha_face": {"n_CH": 4, "quality": "good", "carbons": "C1,C3,C4,C5",
                       "notes": "NAc at C2 reduces stacking at that position."},
        "beta_face":  {"n_CH": 4, "quality": "good", "carbons": "C1,C3,C4,C5"},
    },
    "alpha-Fuc": {
        "alpha_face": {"n_CH": 3, "quality": "good", "carbons": "C1,C3,C5",
                       "notes": "6-deoxy sugar. Methyl at C6 adds hydrophobic contact."},
        "beta_face":  {"n_CH": 4, "quality": "good", "carbons": "C2,C3,C4,C5"},
    },
}


def estimate_CH_pi_energy(sugar: str, aromatic: str = "Trp",
                          face: Optional[str] = None) -> dict:
    """Estimate CH-pi stacking energy for a sugar-aromatic pair.

    Uses per-contact energy from Laughrey data scaled by number of
    effective contacts from the sugar's stacking geometry.

    Args:
        sugar: Key from SUGAR_CH_PI_CONTACTS (e.g., "beta-Gal")
        aromatic: "Trp", "Tyr", or "Phe"
        face: "alpha_face" or "beta_face". If None, uses best face.

    Returns:
        Dict with estimated dG and components.
    """
    estimates = get_eps_CH_pi_estimates()

    aromatic_map = {"Trp": "eps_CH_pi_Trp", "Tyr": "eps_CH_pi_Tyr", "Phe": "eps_CH_pi_Phe"}
    eps_key = aromatic_map.get(aromatic, "eps_CH_pi_pyranose_avg")
    eps = estimates[eps_key]["value_kJ"]

    sugar_info = SUGAR_CH_PI_CONTACTS.get(sugar)
    if sugar_info is None:
        return {"dG_CH_pi_kJ": eps * 3, "n_CH": 3, "notes": f"Unknown sugar {sugar}, default 3 contacts"}

    if face is None:
        # Pick best face
        best_face = max(sugar_info.keys(), key=lambda f: sugar_info[f]["n_CH"])
        face = best_face

    face_info = sugar_info.get(face, list(sugar_info.values())[0])
    n_ch = face_info["n_CH"]

    # Not all CHs contribute equally — effective contacts < geometric contacts
    # Typical: 2-3 effective contacts even for 5-CH face
    n_effective = min(n_ch, 3.5)  # cap effective contacts

    dG = eps * n_effective

    return {
        "dG_CH_pi_kJ": round(dG, 2),
        "eps_per_contact_kJ": eps,
        "n_geometric_CH": n_ch,
        "n_effective_CH": n_effective,
        "face": face,
        "aromatic": aromatic,
        "sugar": sugar,
    }


if __name__ == "__main__":
    print("CH-pi Energetics - Parameter Estimates")
    print("=" * 50)

    estimates = get_eps_CH_pi_estimates()
    for key, val in estimates.items():
        print(f"  {key}: {val['value_kJ']:.1f} kJ/mol  [{val['confidence']}]")

    print()
    print("Sugar-specific estimates (Trp):")
    for sugar in ["beta-Glc", "beta-Gal", "alpha-Man", "beta-GlcNAc"]:
        result = estimate_CH_pi_energy(sugar, "Trp")
        print(f"  {sugar:15s}: dG = {result['dG_CH_pi_kJ']:6.1f} kJ/mol "
              f"({result['n_effective_CH']:.1f} eff. contacts on {result['face']})")