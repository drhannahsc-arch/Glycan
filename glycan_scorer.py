"""
glycan_scorer.py — Predict lectin-sugar binding from pure chemistry parameters.

All parameters derived from non-biological sources:
  - Polyol desolvation:    NIST dissolution calorimetry (Schwarz 1996, Jasra 1982)
  - CH-π stacking:         Synthetic host-guest data (Laughrey 2008, Asensio 2013)
  - Conformational entropy: QM torsion potentials (GLYCAM06, Kirschner 2008)
  - H-bond energy:         Literature consensus (no fitting)

No protein-sugar binding data used in parameter determination.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional


# ═══════════════════════════════════════════════════════════════════════════
# Data structures
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class GlycanParams:
    """All parameters from pure chemistry sources.

    Phase G1: Desolvation (from Schwarz/Jasra dissolution calorimetry)
    Phase G2: Conformational entropy (from QM torsion potentials)
    Phase G3: CH-π stacking (from synthetic host-guest data)
    Phase G5: Structural water (from lectin mutant series — stretch goal)
    """
    # Phase G1: Polyol desolvation
    k_desolv_eq: float = 0.0    # kJ/mol per equatorial OH buried
    k_desolv_ax: float = 0.0    # kJ/mol per axial OH buried
    k_desolv_NAc: float = 0.0   # kJ/mol per NAc group buried

    # Phase G3: CH-π stacking
    eps_CH_pi: float = 0.0      # kJ/mol per CH-π contact

    # H-bond (literature consensus, not fitted)
    eps_HB: float = -6.0        # kJ/mol per H-bond in water (range: -4 to -8)

    # Phase G2: Conformational entropy
    eps_conf: float = 0.0       # kJ/mol per frozen torsion pair (φ/ψ)

    # Phase G5 (stretch): Structural water
    eps_water_bridge: float = 0.0  # kJ/mol per conserved structural water

    # Offset (translational entropy loss, etc.)
    dG_0: float = 0.0


@dataclass
class ContactMap:
    """Contacts extracted from a crystal structure (PDB + PLIP).

    One ContactMap per lectin-sugar pair. Built manually from
    PyMOL visualization + PLIP automated contact analysis.
    """
    lectin: str                      # e.g. "ConA", "WGA", "PNA", "Davis"
    sugar: str                       # e.g. "mannose", "glucose", "2-deoxy-mannose"
    pdb_id: str                      # e.g. "5CNA"

    # Contact counts
    n_HB: int = 0                    # H-bonds (sugar OH ↔ protein, d < 3.2 Å)
    n_CH_pi: int = 0                 # CH-π contacts (sugar CH ↔ aromatic ring)

    # Buried hydroxyls by type
    n_OH_eq_buried: int = 0          # equatorial OH groups buried on binding
    n_OH_ax_buried: int = 0          # axial OH groups buried on binding
    n_NAc_buried: int = 0            # N-acetyl groups buried on binding

    # Conformational cost
    n_frozen_torsions: int = 0       # glycosidic torsion pairs frozen

    # Structural water (Phase G5 stretch)
    n_water_bridges: int = 0         # conserved bridging waters in interface

    # Metadata
    resolution_A: float = 0.0        # crystal structure resolution (Å)
    notes: str = ""


@dataclass
class SugarPropertyCard:
    """Intrinsic properties of a monosaccharide — no receptor information.

    Built in Phase G0 from pure chemistry data (ring geometry,
    hydroxyl stereochemistry, SASA, Cremer-Pople parameters).
    """
    name: str                        # e.g. "D-mannose", "D-glucose"
    abbreviation: str                # e.g. "Man", "Glc", "Gal"

    # Hydroxyl inventory
    n_OH_equatorial: int = 0
    n_OH_axial: int = 0
    has_NAc: bool = False

    # Ring geometry (Cremer-Pople)
    puckering_Q: float = 0.0         # total puckering amplitude
    puckering_theta: float = 0.0     # polar angle (°)
    puckering_phi: float = 0.0       # azimuthal angle (°)

    # Solvent-accessible surface
    sasa_total: float = 0.0          # Å²
    sasa_apolar: float = 0.0         # Å² (CH face area)
    sasa_polar: float = 0.0          # Å² (OH-rich face area)

    # CH-π donor potential
    n_axial_CH: int = 0              # axial C-H bonds (α-face CH-π donors)

    # Dissolution enthalpy (Phase G1 source data)
    dH_sol_kJ: Optional[float] = None        # kJ/mol, from Schwarz/Jasra
    dH_sol_deoxy_series: Dict[str, float] = field(default_factory=dict)
    # e.g. {"2-deoxy": -12.3, "3-deoxy": -11.1, ...}


@dataclass
class PredictionResult:
    """Output of the scoring function for one lectin-sugar pair."""
    lectin: str
    sugar: str
    dG_predicted: float              # kJ/mol
    dG_measured: Optional[float]     # kJ/mol (from ITC, if available)

    # Decomposition
    dG_HB: float = 0.0
    dG_CH_pi: float = 0.0
    dG_desolv: float = 0.0
    dG_conf: float = 0.0
    dG_water: float = 0.0
    dG_offset: float = 0.0


# ═══════════════════════════════════════════════════════════════════════════
# Scoring function
# ═══════════════════════════════════════════════════════════════════════════

def predict_dG(params: GlycanParams, contacts: ContactMap) -> PredictionResult:
    """Score a lectin-sugar interaction from locked pure-chemistry parameters.

    Sign convention (all terms summed directly):
      dG_HB     < 0  (favorable: H-bonds stabilize)
      dG_CH_pi  < 0  (favorable: CH-π stacking stabilizes)
      dG_desolv > 0  (unfavorable: stripping water from OH costs energy)
      dG_conf   > 0  (unfavorable: freezing torsions costs entropy)
      dG_water  < 0  (favorable: structural water bridges stabilize)
      dG_0      = offset (translational entropy loss, etc.)
    """
    # Term 1: H-bond contribution (favorable)
    dG_HB = params.eps_HB * contacts.n_HB

    # Term 2: CH-π stacking (favorable)
    dG_CH_pi = params.eps_CH_pi * contacts.n_CH_pi

    # Term 3: Desolvation penalty (unfavorable — stripping water costs energy)
    # Positive k_desolv × n_buried → positive dG_desolv → weakens binding
    dG_desolv = (
        params.k_desolv_eq * contacts.n_OH_eq_buried
        + params.k_desolv_ax * contacts.n_OH_ax_buried
        + params.k_desolv_NAc * contacts.n_NAc_buried
    )

    # Term 4: Conformational entropy cost (unfavorable — freezing torsions costs entropy)
    # Positive eps_conf × n_frozen → positive dG_conf → weakens binding
    dG_conf = params.eps_conf * contacts.n_frozen_torsions

    # Term 5: Structural water bridges (favorable, stretch goal)
    dG_water = params.eps_water_bridge * contacts.n_water_bridges

    # Term 6: Offset
    dG_offset = params.dG_0

    dG_predicted = dG_HB + dG_CH_pi + dG_desolv + dG_conf + dG_water + dG_offset

    return PredictionResult(
        lectin=contacts.lectin,
        sugar=contacts.sugar,
        dG_predicted=dG_predicted,
        dG_measured=None,
        dG_HB=dG_HB,
        dG_CH_pi=dG_CH_pi,
        dG_desolv=dG_desolv,
        dG_conf=dG_conf,
        dG_water=dG_water,
        dG_offset=dG_offset,
    )


def predict_ddG(params: GlycanParams,
                contacts_wt: ContactMap,
                contacts_variant: ContactMap) -> float:
    """Predict ΔΔG between wildtype sugar and a variant (deoxy, epimer, mutant).

    Used for Prediction 1 (ConA deoxy series) and Prediction 3 (lysozyme mutant).
    """
    result_wt = predict_dG(params, contacts_wt)
    result_var = predict_dG(params, contacts_variant)
    return result_var.dG_predicted - result_wt.dG_predicted


# ═══════════════════════════════════════════════════════════════════════════
# Batch scoring and statistics
# ═══════════════════════════════════════════════════════════════════════════

def score_panel(params: GlycanParams,
                contact_maps: List[ContactMap],
                measured_dG: Optional[Dict[str, float]] = None
                ) -> List[PredictionResult]:
    """Score a panel of sugar-receptor pairs. Attach measured values if provided."""
    results = []
    for cm in contact_maps:
        result = predict_dG(params, cm)
        if measured_dG is not None:
            key = f"{cm.lectin}:{cm.sugar}"
            result.dG_measured = measured_dG.get(key)
        results.append(result)
    return results


def compute_statistics(results: List[PredictionResult]) -> Dict[str, float]:
    """Compute R², MAE, RMSE, Spearman ρ for a set of predictions vs measurements."""
    from scipy import stats as sp_stats

    paired = [(r.dG_predicted, r.dG_measured) for r in results
              if r.dG_measured is not None]
    if len(paired) < 3:
        return {"n": len(paired), "R2": float("nan"), "MAE": float("nan"),
                "RMSE": float("nan"), "spearman_rho": float("nan")}

    pred = np.array([p[0] for p in paired])
    meas = np.array([p[1] for p in paired])

    ss_res = np.sum((meas - pred) ** 2)
    ss_tot = np.sum((meas - np.mean(meas)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    mae = np.mean(np.abs(meas - pred))
    rmse = np.sqrt(np.mean((meas - pred) ** 2))
    rho, p_val = sp_stats.spearmanr(pred, meas)

    return {
        "n": len(paired),
        "R2": round(r2, 3),
        "MAE": round(mae, 2),
        "RMSE": round(rmse, 2),
        "spearman_rho": round(rho, 3),
        "spearman_p": round(p_val, 4),
    }


# ═══════════════════════════════════════════════════════════════════════════
# Quick demo / sanity check
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # Placeholder parameters — will be replaced in Phase G1-G3
    params = GlycanParams(
        k_desolv_eq=8.0,    # kJ/mol (placeholder)
        k_desolv_ax=5.0,    # kJ/mol (placeholder)
        k_desolv_NAc=6.0,   # kJ/mol (placeholder)
        eps_CH_pi=-2.0,     # kJ/mol per contact (Nishio consensus range)
        eps_HB=-6.0,        # kJ/mol per H-bond (literature consensus)
        eps_conf=3.0,        # kJ/mol per frozen torsion pair (placeholder)
        dG_0=10.0,          # kJ/mol offset (placeholder)
    )

    # Example: ConA binding mannose (placeholder contacts)
    conA_mannose = ContactMap(
        lectin="ConA",
        sugar="mannose",
        pdb_id="5CNA",
        n_HB=4,
        n_CH_pi=1,
        n_OH_eq_buried=3,
        n_OH_ax_buried=1,
        n_frozen_torsions=0,  # monosaccharide, no glycosidic bond
    )

    result = predict_dG(params, conA_mannose)
    print(f"ConA + mannose (placeholder params):")
    print(f"  ΔG_predicted = {result.dG_predicted:.1f} kJ/mol")
    print(f"    H-bond:      {result.dG_HB:.1f}")
    print(f"    CH-π:        {result.dG_CH_pi:.1f}")
    print(f"    Desolvation: {result.dG_desolv:.1f}")
    print(f"    Conf. ent.:  {result.dG_conf:.1f}")
    print(f"    Offset:      {result.dG_offset:.1f}")