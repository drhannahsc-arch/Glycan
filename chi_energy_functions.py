"""
CHI Energy Functions for glycosidic torsion scoring.

Source: Nivedha AK, Makeneni S, Foley BL, Tessier MB, Woods RJ.
        J. Comput. Chem. 2014, 35, 526-539. DOI: 10.1002/jcc.23517
        https://glycam.org/docs/forcefield/chi-energy-functions/

These functions score oligosaccharide conformations based on glycosidic
linkage torsion angles (phi, psi). Derived from QM (B3LYP) calculations
on tetrahydropyran-based disaccharide models.

General form: E(theta) = sum_i [ a_i * exp(-(theta - b_i)^2 / (2 * c_i^2)) ] + d

Units: kcal/mol (original publication). Converted to kJ/mol where noted.

MABE integration:
  - Constrains P4 (eps_glycosidic_freeze) without fitting.
  - TdS_freeze = -RT * ln(p_bound / p_free) computed from Boltzmann
    integration over these potentials.
"""

import math
import numpy as np
from typing import Tuple, Dict, Optional

KCAL_TO_KJ = 4.184
R_KCAL = 1.9872e-3  # kcal/(mol*K)
R_KJ = 8.314e-3     # kJ/(mol*K)


# =====================================================================
# CHI EQUATION PARAMETERS (from glycam.org summary tables)
# Each entry: (a_i, b_i, c_i^2) where c_i^2 = 2*c_i^2 in the Gaussian
# NOTE: The website gives the denominator as the FULL 2*c^2 value directly.
# =====================================================================

# E(phi|alpha): phi for axial C1-O1 bond (alpha-linked sugars)
_PHI_ALPHA_PARAMS = [
    ( 2.977,     -199.49,    677.81),
    ( 102.25,     170.6,    1696.8),
    ( 10.745,    -105.31,   4724.6),
    ( 3.6735,      6.2012,  1347.7),
    ( 2.061,      91.655,   1500.0),
    ( 6.1939,    -22.979,   2122.3),
    (-2.1115,     83.602,   1254.1),
    (-98.001,    170.01,    1598.7),
]
_PHI_ALPHA_D = 0.0

# E(phi|beta): phi for equatorial C1-O1 bond (beta-linked sugars)
_PHI_BETA_PARAMS = [
    ( 450.54,    -330.77,   4449.8),
    ( 23.712,     304.63,   8375.2),
    ( 5.9353,    -152.08,   6049.8),
    ( 22.467,     -23.516,   606.90),
    ( 10.036,     120.96,   4038.0),
    (-18.141,     -24.268,   543.05),
    ( 5.8823,      19.632,   897.93),
]
_PHI_BETA_D = -2.1283

# E(psi|2a3e): psi where C2'/C4' OH is axial or C3' OH is equatorial
_PSI_2A3E_PARAMS = [
    ( 4.6237,      5.0456,  5005.8),
    ( 4.6139,    362.49,    2090.6),
    ( 4.9419,    121.2,     2093.8),
    ( 0.4029,    241.43,     456.38),
    ( 0.79888,    68.425,    678.81),
    ( 0.22299,   192.93,     347.25),
]
_PSI_2A3E_D = -0.12565

# E(psi|2e3a): psi where C2'/C4' OH is equatorial or C3' OH is axial
_PSI_2E3A_PARAMS = [
    ( 4.4681,      1e-30,   1279.6),
    ( 4.382,     357.77,    6050.1),
    ( 284.95,    146.64,    1551.8),
    ( 4.7613,    220.68,    5892.9),
    (-169.2,     147.37,    1742.5),
    (-118.44,    146.06,    1359.8),
]
_PSI_2E3A_D = 1.0220


def _eval_chi(angle_deg: float, params: list, d: float) -> float:
    """Evaluate a CHI energy function at a given angle.

    Args:
        angle_deg: Torsion angle in degrees.
        params: List of (a, b, c_sq) tuples.
        d: Constant offset.

    Returns:
        Energy in kcal/mol.
    """
    E = d
    for a, b, c_sq in params:
        E += a * math.exp(-((angle_deg - b) ** 2) / c_sq)
    return E


def chi_phi_alpha(phi_deg: float) -> float:
    """CHI energy for phi torsion of alpha-linked disaccharides (axial C1-O1).

    Args:
        phi_deg: Phi angle in degrees (-180 to 180).

    Returns:
        Energy in kcal/mol.
    """
    return _eval_chi(phi_deg, _PHI_ALPHA_PARAMS, _PHI_ALPHA_D)


def chi_phi_beta(phi_deg: float) -> float:
    """CHI energy for phi torsion of beta-linked disaccharides (equatorial C1-O1).

    Args:
        phi_deg: Phi angle in degrees (-180 to 180).

    Returns:
        Energy in kcal/mol.
    """
    return _eval_chi(phi_deg, _PHI_BETA_PARAMS, _PHI_BETA_D)


def chi_psi_2a3e(psi_deg: float) -> float:
    """CHI energy for psi torsion (C2'/C4' axial or C3' equatorial).

    Args:
        psi_deg: Psi angle in degrees (0 to 360).

    Returns:
        Energy in kcal/mol.
    """
    return _eval_chi(psi_deg, _PSI_2A3E_PARAMS, _PSI_2A3E_D)


def chi_psi_2e3a(psi_deg: float) -> float:
    """CHI energy for psi torsion (C2'/C4' equatorial or C3' axial).

    Args:
        psi_deg: Psi angle in degrees (0 to 360).

    Returns:
        Energy in kcal/mol.
    """
    return _eval_chi(psi_deg, _PSI_2E3A_PARAMS, _PSI_2E3A_D)


# =====================================================================
# LINKAGE TYPE CLASSIFICATION
# Maps common glycosidic linkages to the appropriate CHI equation pair.
# =====================================================================

# For phi: depends on anomeric configuration (alpha = axial, beta = equatorial)
# For psi: depends on the orientation of substituents at the reducing-end carbon.
#
# Common linkage types and their psi class:
# (1->2): C2' OH equatorial in Glc/Man/Gal -> 2e3a
# (1->3): C3' OH equatorial in Glc, axial in nothing common -> 2a3e for most
# (1->4): C4' OH equatorial in Glc/Man, axial in Gal -> depends on sugar
# (1->6): Not covered by two-bond CHI (three-bond linkage with omega)

LINKAGE_PSI_CLASS = {
    # (linkage_position, reducing_sugar): psi_class
    # Glucose-based (all equatorial)
    (2, "Glc"): "2e3a",
    (3, "Glc"): "2a3e",   # C3-OH is equatorial in Glc -> class 2a3e
    (4, "Glc"): "2a3e",   # C4-OH is equatorial in Glc
    # Galactose-based (C4-OH axial)
    (2, "Gal"): "2e3a",
    (3, "Gal"): "2a3e",
    (4, "Gal"): "2e3a",   # C4-OH is axial in Gal -> class 2e3a
    # Mannose-based (C2-OH axial)
    (2, "Man"): "2a3e",   # C2-OH is axial in Man
    (3, "Man"): "2a3e",
    (4, "Man"): "2a3e",
    # GlcNAc (same ring as Glc)
    (2, "GlcNAc"): "2e3a",
    (3, "GlcNAc"): "2a3e",
    (4, "GlcNAc"): "2a3e",
    # GalNAc (same ring as Gal)
    (2, "GalNAc"): "2e3a",
    (3, "GalNAc"): "2a3e",
    (4, "GalNAc"): "2e3a",
}


def get_chi_functions(anomer: str, linkage_pos: int,
                      reducing_sugar: str = "Glc"):
    """Return the appropriate (phi_func, psi_func) for a given linkage.

    Args:
        anomer: "alpha" or "beta"
        linkage_pos: Position on reducing sugar (2, 3, or 4).
            Position 6 linkages have an extra omega torsion not covered by CHI.
        reducing_sugar: Three-letter code of the reducing-end sugar.

    Returns:
        Tuple of (phi_function, psi_function) or raises ValueError.
    """
    if linkage_pos == 6:
        raise ValueError(
            "(1->6) linkages have three torsion angles (phi, psi, omega). "
            "CHI covers only two-bond linkages. Use GLYCAM06 for omega."
        )

    phi_func = chi_phi_alpha if anomer == "alpha" else chi_phi_beta

    key = (linkage_pos, reducing_sugar)
    psi_class = LINKAGE_PSI_CLASS.get(key)
    if psi_class is None:
        raise ValueError(
            f"No CHI psi classification for linkage ({linkage_pos}, {reducing_sugar}). "
            f"Known: {sorted(LINKAGE_PSI_CLASS.keys())}"
        )

    psi_func = chi_psi_2a3e if psi_class == "2a3e" else chi_psi_2e3a

    return phi_func, psi_func


# =====================================================================
# CONFORMATIONAL ENTROPY FROM CHI POTENTIALS
# =====================================================================

def compute_TdS_freeze(phi_func, psi_func,
                       phi_bound_deg: float, psi_bound_deg: float,
                       T_K: float = 298.15,
                       resolution_deg: float = 1.0,
                       libration_window_deg: float = 30.0) -> dict:
    """Compute the entropy cost of freezing a glycosidic linkage upon binding.

    Uses Boltzmann integration over the CHI potential. The "free" state
    samples all conformations weighted by energy. The "bound" state is
    modeled as a librational window (typically +/-30 deg) around the
    crystal structure torsion angles, reflecting thermal fluctuations
    in the binding pocket.

    TdS_freeze = -RT * ln(Z_bound / Z_free)

    where Z_bound integrates over the librational window and Z_free
    integrates over all angles.

    Args:
        phi_func: CHI phi energy function.
        psi_func: CHI psi energy function.
        phi_bound_deg: Bound-state phi from crystal structure.
        psi_bound_deg: Bound-state psi from crystal structure.
        T_K: Temperature in Kelvin.
        resolution_deg: Integration grid spacing in degrees.
        libration_window_deg: Half-width of bound-state librational window.

    Returns:
        Dict with:
            TdS_freeze_kcal: Entropy cost in kcal/mol (positive = unfavorable)
            TdS_freeze_kJ: Same in kJ/mol
            E_bound_kcal: CHI energy at the bound conformation
            phi_min_deg: Phi at the global energy minimum
            psi_min_deg: Psi at the global energy minimum
    """
    beta = 1.0 / (R_KCAL * T_K)
    hw = libration_window_deg

    # Determine phi range
    if phi_func is chi_phi_alpha or phi_func is chi_phi_beta:
        phi_range = np.arange(-180, 180, resolution_deg)
    else:
        phi_range = np.arange(0, 360, resolution_deg)

    psi_range = np.arange(0, 360, resolution_deg)

    # Compute 1D partition functions (phi and psi are independent in CHI)
    phi_energies = np.array([phi_func(p) for p in phi_range])
    phi_min_e = phi_energies.min()
    phi_energies_shifted = phi_energies - phi_min_e
    phi_weights = np.exp(-beta * phi_energies_shifted)
    Z_phi_free = phi_weights.sum()

    psi_energies = np.array([psi_func(p) for p in psi_range])
    psi_min_e = psi_energies.min()
    psi_energies_shifted = psi_energies - psi_min_e
    psi_weights = np.exp(-beta * psi_energies_shifted)
    Z_psi_free = psi_weights.sum()

    # Bound-state partition function: integrate over librational window
    def _in_window(angles, center, half_width):
        """Check if angles are within half_width of center (circular)."""
        diff = np.abs(angles - center)
        diff = np.minimum(diff, 360 - diff)
        return diff <= half_width

    phi_bound_mask = _in_window(phi_range, phi_bound_deg, hw)
    Z_phi_bound = phi_weights[phi_bound_mask].sum()

    psi_bound_mask = _in_window(psi_range, psi_bound_deg, hw)
    Z_psi_bound = psi_weights[psi_bound_mask].sum()

    # Fraction of partition function in bound window
    f_phi = Z_phi_bound / Z_phi_free if Z_phi_free > 0 else 1e-10
    f_psi = Z_psi_bound / Z_psi_free if Z_psi_free > 0 else 1e-10
    f_total = f_phi * f_psi

    # Entropy cost
    if f_total > 0 and f_total < 1:
        TdS_freeze_kcal = -R_KCAL * T_K * math.log(f_total)
    elif f_total >= 1:
        TdS_freeze_kcal = 0.0
    else:
        TdS_freeze_kcal = 15.0  # cap

    # Energy at bound point
    E_phi_bound = phi_func(phi_bound_deg)
    E_psi_bound = psi_func(psi_bound_deg)

    # Find global minimum
    phi_min_idx = np.argmin(phi_energies)
    psi_min_idx = np.argmin(psi_energies)

    return {
        "TdS_freeze_kcal": round(TdS_freeze_kcal, 3),
        "TdS_freeze_kJ": round(TdS_freeze_kcal * KCAL_TO_KJ, 3),
        "E_bound_kcal": round(E_phi_bound + E_psi_bound, 3),
        "E_bound_kJ": round((E_phi_bound + E_psi_bound) * KCAL_TO_KJ, 3),
        "phi_min_deg": round(float(phi_range[phi_min_idx]), 1),
        "psi_min_deg": round(float(psi_range[psi_min_idx]), 1),
        "f_bound": round(f_total, 6),
    }


def score_linkage_entropy(anomer: str, linkage_pos: int,
                          reducing_sugar: str,
                          phi_bound: float, psi_bound: float,
                          T_K: float = 298.15) -> dict:
    """Convenience: score conformational entropy for a specific linkage.

    Args:
        anomer: "alpha" or "beta"
        linkage_pos: 2, 3, or 4
        reducing_sugar: "Glc", "Gal", "Man", "GlcNAc", "GalNAc"
        phi_bound: Bound-state phi from crystal structure (degrees)
        psi_bound: Bound-state psi from crystal structure (degrees)
        T_K: Temperature

    Returns:
        Dict with TdS_freeze and related values.
    """
    phi_func, psi_func = get_chi_functions(anomer, linkage_pos, reducing_sugar)
    return compute_TdS_freeze(phi_func, psi_func, phi_bound, psi_bound, T_K)


# =====================================================================
# QUICK SANITY CHECK
# =====================================================================

if __name__ == "__main__":
    print("CHI Energy Functions - Sanity Check")
    print("=" * 50)

    # Check phi_alpha at known minimum (~60 deg for alpha linkages)
    for phi in [-60, 0, 60, 120, 180]:
        e = chi_phi_alpha(phi)
        print(f"  phi_alpha({phi:4d}) = {e:8.3f} kcal/mol")

    print()
    for psi in [0, 60, 120, 180, 240, 300]:
        e = chi_psi_2a3e(psi)
        print(f"  psi_2a3e({psi:4d}) = {e:8.3f} kcal/mol")

    print()
    # ConA binds MeAlphaMan: typical phi~-90, psi~150 for alpha(1->3) linkage
    result = score_linkage_entropy(
        anomer="alpha", linkage_pos=3, reducing_sugar="Man",
        phi_bound=-90, psi_bound=150
    )
    print(f"  ConA alpha(1->3)Man: TdS_freeze = {result['TdS_freeze_kJ']:.2f} kJ/mol")
    print(f"    f_bound = {result['f_bound']:.6f}")
    print(f"    E_bound = {result['E_bound_kJ']:.2f} kJ/mol")