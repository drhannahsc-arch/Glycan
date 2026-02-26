"""
parameter_fitting.py — Determine glycan scoring parameters from pure chemistry.

Phase G1: Polyol desolvation from dissolution calorimetry (Schwarz 1996, Jasra 1982)
Phase G2: Conformational entropy from QM torsion potentials (GLYCAM06)
Phase G3: CH-π stacking from synthetic host-guest data (Laughrey 2008)

No protein-sugar binding data used anywhere in this module.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from glycan_scorer import GlycanParams


# ═══════════════════════════════════════════════════════════════════════════
# Phase G1: Polyol desolvation from dissolution calorimetry
# ═══════════════════════════════════════════════════════════════════════════
#
# Logic: ΔH_sol(sugar) − ΔH_sol(deoxy-sugar) = hydration cost of that OH.
# Schwarz (1996) measured ΔH_sol for glucose, mannose, and their deoxy
# derivatives. The per-hydroxyl ΔΔH, averaged by stereochemistry
# (equatorial vs axial), gives k_desolv_eq and k_desolv_ax.
#
# This is the SAME approach as the metal module: isolate a single physics
# term by comparing systematically varied compounds.

def fit_desolvation_from_schwarz(schwarz_data: Dict[str, float]) -> Dict[str, float]:
    """Extract per-hydroxyl desolvation costs from dissolution calorimetry.

    Args:
        schwarz_data: Dict mapping compound name → ΔH_sol (kJ/mol).
            Expected keys include parent sugars and their deoxy derivatives.
            e.g. {"D-glucose": -11.0, "2-deoxy-D-glucose": -8.5, ...}

    Returns:
        Dict with keys: k_desolv_eq, k_desolv_ax, k_desolv_NAc,
        plus per-position values for diagnostics.
    """
    # Placeholder structure — populate with actual Schwarz table values
    # in Phase G1 execution.
    #
    # For each deoxy derivative:
    #   ΔΔH = ΔH_sol(parent) − ΔH_sol(deoxy) = hydration cost of removed OH
    #
    # Then classify by stereochemistry:
    #   C2-OH in glucose = equatorial → contributes to k_desolv_eq
    #   C2-OH in mannose = axial → contributes to k_desolv_ax

    # Position → (parent, deoxy_derivative, OH_stereochemistry)
    DEOXY_PAIRS = {
        # Glucose series
        "Glc_C2": ("D-glucose", "2-deoxy-D-glucose", "eq"),
        "Glc_C3": ("D-glucose", "3-deoxy-D-glucose", "eq"),
        "Glc_C4": ("D-glucose", "4-deoxy-D-glucose", "eq"),
        "Glc_C6": ("D-glucose", "6-deoxy-D-glucose", "eq"),  # 6-deoxy = quinovose
        # Mannose series
        "Man_C2": ("D-mannose", "2-deoxy-D-mannose", "ax"),
        "Man_C3": ("D-mannose", "3-deoxy-D-mannose", "eq"),
        "Man_C4": ("D-mannose", "4-deoxy-D-mannose", "eq"),
        "Man_C6": ("D-mannose", "6-deoxy-D-mannose", "eq"),  # 6-deoxy = rhamnose
        # Galactose series (if available)
        "Gal_C4": ("D-galactose", "4-deoxy-D-galactose", "ax"),
    }

    eq_values = []
    ax_values = []
    per_position = {}

    for label, (parent, deoxy, stereo) in DEOXY_PAIRS.items():
        if parent in schwarz_data and deoxy in schwarz_data:
            ddH = schwarz_data[parent] - schwarz_data[deoxy]
            per_position[label] = {"ddH_kJ": ddH, "stereo": stereo}
            if stereo == "eq":
                eq_values.append(ddH)
            elif stereo == "ax":
                ax_values.append(ddH)

    k_eq = np.mean(eq_values) if eq_values else 0.0
    k_ax = np.mean(ax_values) if ax_values else 0.0

    # Sanity check: equatorial OH should be MORE hydrated (higher desolvation cost)
    # because equatorial OH makes more water contacts than axial
    if k_eq > 0 and k_ax > 0 and k_eq < k_ax:
        print("WARNING: k_desolv_eq < k_desolv_ax — expected eq > ax. "
              "Check data or stereochemistry assignments.")

    return {
        "k_desolv_eq": round(k_eq, 2),
        "k_desolv_ax": round(k_ax, 2),
        "k_desolv_NAc": 0.0,  # From Jasra: GlcNAc vs Glc comparison (Phase G1)
        "per_position": per_position,
        "n_eq_points": len(eq_values),
        "n_ax_points": len(ax_values),
    }


def fit_NAc_desolvation(jasra_data: Dict[str, float]) -> float:
    """Extract NAc desolvation cost from Jasra (1982) polyol data.

    Compare ΔH_sol(GlcNAc) vs ΔH_sol(Glc) to isolate
    the N-acetyl hydration contribution.
    """
    if "GlcNAc" in jasra_data and "D-glucose" in jasra_data:
        return round(jasra_data["D-glucose"] - jasra_data["GlcNAc"], 2)
    return 0.0


# ═══════════════════════════════════════════════════════════════════════════
# Phase G2: Conformational entropy from QM torsion potentials
# ═══════════════════════════════════════════════════════════════════════════
#
# Logic: When a glycosidic bond freezes upon binding, the conformational
# entropy cost = RT × Σ pᵢ ln(pᵢ) over the Boltzmann-populated torsion states.
#
# Source: GLYCAM06 QM torsion potential energy surfaces (Kirschner 2008).
# Or: Woods group published φ/ψ population maps.

def compute_torsion_entropy(
    torsion_energies_kJ: np.ndarray,
    angles_deg: np.ndarray,
    T_K: float = 298.15
) -> float:
    """Compute conformational entropy of a torsion from its energy profile.

    Args:
        torsion_energies_kJ: Energy at each torsion angle (kJ/mol)
        angles_deg: Torsion angles (degrees), same length as energies
        T_K: Temperature (K)

    Returns:
        TΔS in kJ/mol (positive = entropy cost of freezing)
    """
    R = 8.314e-3  # kJ/(mol·K)
    RT = R * T_K

    # Boltzmann populations
    E_rel = torsion_energies_kJ - np.min(torsion_energies_kJ)
    weights = np.exp(-E_rel / RT)
    Z = np.sum(weights)
    populations = weights / Z

    # Shannon entropy: S = -R Σ pᵢ ln(pᵢ)
    # Filter zero populations to avoid log(0)
    mask = populations > 1e-15
    S_conf = -R * np.sum(populations[mask] * np.log(populations[mask]))

    # TΔS_freeze = T × S_conf (this is the cost of freezing the torsion)
    return round(T_K * S_conf, 2)


def fit_conformational_entropy(
    linkage_torsion_profiles: Dict[str, Tuple[np.ndarray, np.ndarray]]
) -> Dict[str, float]:
    """Compute ε_conf from QM torsion profiles for common linkage types.

    Args:
        linkage_torsion_profiles: Dict mapping linkage type (e.g. "1->4", "1->6")
            to (angles_deg, energies_kJ) arrays.

    Returns:
        Dict with per-linkage TΔS values and the mean ε_conf.
    """
    per_linkage = {}
    all_TdS = []

    for linkage, (angles, energies) in linkage_torsion_profiles.items():
        TdS = compute_torsion_entropy(energies, angles)
        per_linkage[linkage] = TdS
        all_TdS.append(TdS)

    eps_conf = np.mean(all_TdS) if all_TdS else 0.0

    return {
        "eps_conf": round(eps_conf, 2),
        "per_linkage": per_linkage,
        "n_linkages": len(all_TdS),
    }


# ═══════════════════════════════════════════════════════════════════════════
# Phase G3: CH-π stacking from synthetic host-guest data
# ═══════════════════════════════════════════════════════════════════════════
#
# Logic: Synthetic aromatic hosts (Laughrey 2008, Diederich cyclophanes,
# Asensio model systems) bind guests via CH-π contacts. Measure ΔG_bind,
# subtract known contributions (HG framework handles SASA, H-bond, shape),
# divide residual by n_CH-π contacts.

def fit_ch_pi_from_synthetic_hosts(
    host_guest_data: List[Dict]
) -> Dict[str, float]:
    """Extract ε_CH-π per contact from synthetic host-guest binding data.

    Args:
        host_guest_data: List of dicts, each with:
            - "host": str (host name)
            - "guest": str (guest name)
            - "dG_measured": float (kJ/mol, measured binding ΔG)
            - "dG_other_terms": float (kJ/mol, sum of non-CH-π terms
              from locked HG framework)
            - "n_CH_pi_contacts": int (counted from structure/model)

    Returns:
        Dict with eps_CH_pi and diagnostics.
    """
    per_contact_values = []

    for entry in host_guest_data:
        n = entry.get("n_CH_pi_contacts", 0)
        if n == 0:
            continue
        residual = entry["dG_measured"] - entry["dG_other_terms"]
        per_contact = residual / n
        per_contact_values.append(per_contact)

    if not per_contact_values:
        return {"eps_CH_pi": 0.0, "n_entries": 0}

    eps_CH_pi = np.mean(per_contact_values)

    # Sanity check: Nishio consensus for CH-π in water is -1.5 to -3.0 kJ/mol
    if eps_CH_pi > 0:
        print("WARNING: ε_CH-π is positive (repulsive). Expected negative.")
    if eps_CH_pi < -4.0:
        print("WARNING: ε_CH-π < -4.0 kJ/mol — stronger than Nishio consensus. "
              "Check for double-counting with hydrophobic term.")

    return {
        "eps_CH_pi": round(eps_CH_pi, 2),
        "std_dev": round(np.std(per_contact_values), 2),
        "n_entries": len(per_contact_values),
        "range": (round(min(per_contact_values), 2),
                  round(max(per_contact_values), 2)),
    }


# ═══════════════════════════════════════════════════════════════════════════
# Assemble full parameter set
# ═══════════════════════════════════════════════════════════════════════════

def assemble_params(
    desolv_result: Dict,
    conf_result: Dict,
    ch_pi_result: Dict,
    eps_HB: float = -6.0,
    dG_0: float = 0.0,
) -> GlycanParams:
    """Assemble a GlycanParams object from Phase G1-G3 fitting results."""
    return GlycanParams(
        k_desolv_eq=desolv_result["k_desolv_eq"],
        k_desolv_ax=desolv_result["k_desolv_ax"],
        k_desolv_NAc=desolv_result["k_desolv_NAc"],
        eps_CH_pi=ch_pi_result["eps_CH_pi"],
        eps_HB=eps_HB,
        eps_conf=conf_result["eps_conf"],
        dG_0=dG_0,
    )