"""
Glycan scorer v2 — physics-based binding energy from PDB contact maps.

Takes a ContactMap (from contact_extractor.py) and computes predicted
ΔG using calibrated physics parameters.

Includes:
  - Score from real PDB contacts (not hand-curated estimates)
  - Deoxy-mannose prediction scaffold for Prediction 1 (ConA)
  - Selectivity scoring across sugar panels
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Import from our modules
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from contact_extractor import ContactMap, extract_contacts, extract_contacts_simple
from data.ch_pi_energetics import get_eps_CH_pi_estimates, LAUGHREY_PER_CONTACT
from data.pdb_lectin_reference import LECTIN_REFERENCE_SET


# =====================================================================
# Parameters
# =====================================================================

@dataclass
class GlycanParams:
    """Glycan scoring parameters. Tier 1 values filled, Tier 2 placeholder."""

    # H-bond energy per donor-acceptor pair
    # Back-solved from ConA (-7.1) and Galectin-3 (-7.5) ITC data
    # Physical range: -4 to -12 kJ/mol for protein H-bonds
    eps_HB_kJ: float = -7.3

    # CH-π per aromatic residue in contact
    # Trp-weighted from Laughrey 2008
    eps_CH_pi_Trp_kJ: float = -3.5
    eps_CH_pi_Tyr_kJ: float = -2.8
    eps_CH_pi_Phe_kJ: float = -2.3
    eps_CH_pi_His_kJ: float = -1.5   # His is weak — partial aromatic

    # Desolvation cost per buried hydroxyl
    # PLACEHOLDER — needs Schwarz 1996
    k_desolv_eq_kJ: float = 3.0    # Cost of burying equatorial OH (unfavorable)
    k_desolv_ax_kJ: float = 2.0    # Axial OH: less hydrated, cheaper to bury
    k_desolv_NAc_kJ: float = 4.0   # NAc desolvation cost

    # Water bridge energy per bridging water
    eps_water_bridge_kJ: float = -3.0

    # Metal coordination (sugar O directly coordinating metal)
    eps_metal_kJ: float = -8.0

    # Conformational entropy per frozen torsion
    eps_conf_kJ: float = 5.0       # Unfavorable (cost)

    # Offset — translational/rotational entropy loss upon binding
    # ~25 kJ/mol for small-molecule → protein (Finkelstein & Janin 1989)
    dG_0_kJ: float = 25.0

    @classmethod
    def from_calibrated(cls) -> "GlycanParams":
        """Return parameters from Tier 1 calibration."""
        return cls()  # Currently the defaults ARE the Tier 1 values

    def get_ch_pi_energy(self, aromatic_resname: str) -> float:
        """Return CH-π energy for a specific aromatic residue type."""
        res3 = aromatic_resname[:3].upper()
        return {
            "TRP": self.eps_CH_pi_Trp_kJ,
            "TYR": self.eps_CH_pi_Tyr_kJ,
            "PHE": self.eps_CH_pi_Phe_kJ,
            "HIS": self.eps_CH_pi_His_kJ,
        }.get(res3, -2.0)


# =====================================================================
# Scorer
# =====================================================================

@dataclass
class ScoringResult:
    """Result of scoring a lectin-sugar interaction."""
    dG_predicted_kJ: float
    components: Dict[str, float]
    contact_map: ContactMap
    params_used: str = "GlycanParams v2 (Tier 1)"

    def summary(self) -> str:
        parts = [f"{k}={v:+.1f}" for k, v in self.components.items() if abs(v) > 0.01]
        return (f"ΔG_pred = {self.dG_predicted_kJ:+.1f} kJ/mol  "
                f"[{', '.join(parts)}]")


def score_contact_map(cm: ContactMap,
                      params: Optional[GlycanParams] = None) -> ScoringResult:
    """Score a ContactMap to predict binding free energy.

    Args:
        cm: ContactMap from contact_extractor
        params: GlycanParams. If None, uses calibrated defaults.

    Returns:
        ScoringResult with predicted ΔG and component breakdown.
    """
    if params is None:
        params = GlycanParams.from_calibrated()

    # ── H-bonds ──────────────────────────────────────────────
    dG_HB = params.eps_HB_kJ * cm.n_hbonds

    # ── CH-π (per aromatic residue, type-specific) ───────────
    dG_CH_pi = 0.0
    for ar in cm.aromatic_residues:
        dG_CH_pi += params.get_ch_pi_energy(ar)

    # ── Desolvation penalty (burying hydroxyls) ──────────────
    dG_desolv = (params.k_desolv_eq_kJ * cm.n_buried_eq_OH
                 + params.k_desolv_ax_kJ * cm.n_buried_ax_OH
                 + params.k_desolv_NAc_kJ * cm.n_buried_NAc)

    # ── Water bridges ────────────────────────────────────────
    dG_water = params.eps_water_bridge_kJ * cm.n_water_bridges

    # ── Metal coordination ───────────────────────────────────
    dG_metal = params.eps_metal_kJ * len(cm.metal_ions)

    # ── Conformational entropy ───────────────────────────────
    dG_conf = params.eps_conf_kJ * cm.n_frozen_torsions

    # ── Total ────────────────────────────────────────────────
    # Signs: HB/CH-pi/water/metal = favorable (negative)
    #        desolv/conf/dG_0 = unfavorable (positive)
    dG_total = (dG_HB + dG_CH_pi + dG_desolv + dG_water
                + dG_metal + dG_conf + params.dG_0_kJ)

    return ScoringResult(
        dG_predicted_kJ=round(dG_total, 2),
        components={
            "HB": round(dG_HB, 2),
            "CH_pi": round(dG_CH_pi, 2),
            "desolv": round(dG_desolv, 2),       # Positive = penalty
            "water": round(dG_water, 2),
            "metal": round(dG_metal, 2),
            "conf": round(dG_conf, 2),            # Positive = penalty
            "entropy": round(params.dG_0_kJ, 2),  # Translational entropy
        },
        contact_map=cm,
    )


# =====================================================================
# Deoxy-mannose prediction scaffold (Prediction 1)
# =====================================================================

# In a deoxy-mannose derivative, one OH is replaced by H.
# This removes: the H-bonds that OH made + desolvation cost of that OH.
# Net effect: ΔΔG = -(H-bonds lost × ε_HB) + (desolvation saved)

# Mannose OH positions and their axial/equatorial status:
MANNOSE_OH_POSITIONS = {
    "O2": "ax",    # C2-OH is axial in mannose
    "O3": "eq",    # C3-OH equatorial
    "O4": "eq",    # C4-OH equatorial
    "O6": "eq",    # C6-OH equatorial (primary alcohol)
}


def predict_deoxy_series(pdb_id: str = "5CNA",
                         sugar_resname: str = "MMA",
                         chain: Optional[str] = None,
                         params: Optional[GlycanParams] = None,
                         cache_dir: Optional[str] = None) -> dict:
    """Predict ΔΔG for each deoxy-mannose derivative vs parent mannose.

    For each OH position (O2, O3, O4, O6):
      1. Count how many H-bonds that OH makes in the crystal structure
      2. Determine if it's equatorial or axial (for desolvation cost)
      3. Predict ΔΔG = lost binding + saved desolvation

    Returns dict mapping position → predicted ΔΔG.
    """
    if params is None:
        params = GlycanParams.from_calibrated()

    # Extract full contact map for parent mannose
    cms = extract_contacts(pdb_id, sugar_resname=sugar_resname,
                           chain=chain, cache_dir=cache_dir)
    if not cms:
        return {"error": f"No {sugar_resname} found in {pdb_id}"}

    # Use first chain's contacts
    cm = cms[0]
    parent_result = score_contact_map(cm, params)

    # For each OH position, compute what's lost when it's removed
    predictions = {}

    for oh_name, oh_class in MANNOSE_OH_POSITIONS.items():
        # Count H-bonds made by this specific OH
        hbonds_at_position = 0
        for contact in cm.contacts:
            if contact.get("type") == "hbond" and contact.get("sugar_atom") == oh_name:
                hbonds_at_position += 1

        # Also count from hbond_partners string list
        if hbonds_at_position == 0:
            for hb in cm.hbond_partners:
                if hb.startswith(oh_name + "..."):
                    hbonds_at_position += 1

        # Binding energy lost (positive = weaker binding)
        dG_lost = -params.eps_HB_kJ * hbonds_at_position  # Lost favorable H-bonds

        # Desolvation saved (negative = helps binding, but this is a penalty removed)
        if oh_class == "eq":
            desolv_saved = params.k_desolv_eq_kJ  # We save this cost
        else:
            desolv_saved = params.k_desolv_ax_kJ

        # Net ΔΔG: positive means weaker binding (deoxy binds worse)
        ddG = dG_lost - desolv_saved

        predictions[oh_name] = {
            "position": oh_name,
            "axial_equatorial": oh_class,
            "hbonds_lost": hbonds_at_position,
            "dG_HB_lost_kJ": round(dG_lost, 2),
            "desolv_saved_kJ": round(desolv_saved, 2),
            "ddG_predicted_kJ": round(ddG, 2),
            "interpretation": "weaker" if ddG > 0 else "stronger" if ddG < 0 else "unchanged",
        }

    return {
        "parent": {
            "pdb_id": pdb_id,
            "sugar": sugar_resname,
            "dG_parent_kJ": parent_result.dG_predicted_kJ,
            "contact_summary": cm.summary(),
        },
        "deoxy_predictions": predictions,
        "notes": "ΔΔG > 0 means deoxy binds WEAKER than parent. "
                 "Rank order and relative magnitudes are predictions; "
                 "absolute values depend on Tier 2 parameter calibration.",
    }


# =====================================================================
# Selectivity scoring (Prediction 2/7 — same pocket, different sugars)
# =====================================================================

# For selectivity, we can't just swap sugar residue names in the PDB —
# different sugars have different OH patterns. Instead, we modify the
# ContactMap to reflect what WOULD change with a different sugar.
#
# Key differences between sugars (in 4C1 chair):
#   Man vs Glc: C2-OH axial (Man) vs equatorial (Glc) — changes H-bond geometry
#   Glc vs Gal: C4-OH equatorial (Glc) vs axial (Gal) — changes stacking face
#   GlcNAc vs Glc: C2 has NAc instead of OH — bigger, different interactions

SUGAR_OH_PATTERNS = {
    # sugar: {position: "eq"/"ax"/"NAc"/None}
    "Man": {"O2": "ax", "O3": "eq", "O4": "eq", "O6": "eq"},
    "Glc": {"O2": "eq", "O3": "eq", "O4": "eq", "O6": "eq"},
    "Gal": {"O2": "eq", "O3": "eq", "O4": "ax", "O6": "eq"},
    "GlcNAc": {"O2": "NAc", "O3": "eq", "O4": "eq", "O6": "eq"},
    "GalNAc": {"O2": "NAc", "O3": "eq", "O4": "ax", "O6": "eq"},
    "Fuc": {"O2": "eq", "O3": "eq", "O4": "ax"},  # 6-deoxy
}


def score_sugar_panel(pdb_id: str, sugar_resname: str,
                      panel: List[str] = None,
                      chain: Optional[str] = None,
                      params: Optional[GlycanParams] = None,
                      cache_dir: Optional[str] = None) -> dict:
    """Score multiple sugars in the same binding pocket.

    Uses the crystal structure contacts as a template, then adjusts
    for OH pattern differences between sugars.

    Args:
        pdb_id: PDB code of lectin structure
        sugar_resname: Sugar in the crystal structure
        panel: List of sugar types to score (e.g., ["Man", "Glc", "Gal"])
        chain: Chain filter
        params: Scoring parameters

    Returns:
        Dict with scores for each sugar and ΔΔG relative to best binder.
    """
    if params is None:
        params = GlycanParams.from_calibrated()
    if panel is None:
        panel = ["Man", "Glc", "Gal", "GlcNAc"]

    # Get template contacts from crystal structure
    cms = extract_contacts(pdb_id, sugar_resname=sugar_resname,
                           chain=chain, cache_dir=cache_dir)
    if not cms:
        return {"error": f"No {sugar_resname} found in {pdb_id}"}

    template = cms[0]

    from contact_extractor import SUGAR_TYPE
    template_sugar = SUGAR_TYPE.get(sugar_resname, "unknown")
    template_pattern = SUGAR_OH_PATTERNS.get(template_sugar, {})

    results = {}

    for sugar in panel:
        target_pattern = SUGAR_OH_PATTERNS.get(sugar, {})

        # Start from template ContactMap, adjust for OH differences
        adjusted_cm = ContactMap(
            pdb_id=template.pdb_id,
            chain=template.chain,
            sugar_resname=sugar,
            sugar_resid=template.sugar_resid,
            n_hbonds=template.n_hbonds,
            n_ch_pi=template.n_ch_pi,
            n_water_bridges=template.n_water_bridges,
            n_frozen_torsions=template.n_frozen_torsions,
            aromatic_residues=list(template.aromatic_residues),
            metal_ions=list(template.metal_ions),
        )

        # Adjust buried OH counts based on target sugar's pattern
        n_eq = 0
        n_ax = 0
        n_nac = 0

        # For each position that was buried in the template,
        # check what it would be in the target sugar
        all_positions = set(list(template_pattern.keys()) + list(target_pattern.keys()))

        for pos in all_positions:
            # Was this position buried in the crystal structure?
            template_cls = template_pattern.get(pos)
            target_cls = target_pattern.get(pos)

            if template_cls is None or target_cls is None:
                continue  # Position doesn't exist in one of the sugars

            # Check if this OH was in contact in the template
            was_buried = False
            for hb in template.hbond_partners:
                if hb.startswith(pos + "..."):
                    was_buried = True
                    break

            if not was_buried:
                # Also check contacts list
                for c in template.contacts:
                    if c.get("sugar_atom") == pos:
                        was_buried = True
                        break

            if was_buried:
                if target_cls == "eq":
                    n_eq += 1
                elif target_cls == "ax":
                    n_ax += 1
                elif target_cls == "NAc":
                    n_nac += 1

                # If stereochemistry changed (eq→ax or vice versa),
                # one H-bond may be lost due to geometry change
                if template_cls != target_cls and template_cls != "NAc" and target_cls != "NAc":
                    adjusted_cm.n_hbonds = max(0, adjusted_cm.n_hbonds - 1)

        # If target has NAc where template had OH, adjust:
        # NAc is bigger → may gain hydrophobic contact, lose one H-bond
        for pos in target_pattern:
            if target_pattern[pos] == "NAc" and template_pattern.get(pos) != "NAc":
                n_nac += 1
                adjusted_cm.n_hbonds = max(0, adjusted_cm.n_hbonds - 1)
            elif template_pattern.get(pos) == "NAc" and target_pattern[pos] != "NAc":
                n_nac = max(0, n_nac)

        # If no buried positions were explicitly found, use template counts
        if n_eq + n_ax + n_nac == 0:
            n_eq = template.n_buried_eq_OH
            n_ax = template.n_buried_ax_OH
            n_nac = template.n_buried_NAc

        adjusted_cm.n_buried_eq_OH = n_eq
        adjusted_cm.n_buried_ax_OH = n_ax
        adjusted_cm.n_buried_NAc = n_nac

        result = score_contact_map(adjusted_cm, params)
        results[sugar] = {
            "dG_kJ": result.dG_predicted_kJ,
            "components": result.components,
            "n_hbonds": adjusted_cm.n_hbonds,
            "n_ch_pi": adjusted_cm.n_ch_pi,
            "n_eq_OH": n_eq,
            "n_ax_OH": n_ax,
            "n_NAc": n_nac,
        }

    # Compute ΔΔG relative to best binder
    best = min(results.values(), key=lambda x: x["dG_kJ"])["dG_kJ"]
    for sugar in results:
        results[sugar]["ddG_kJ"] = round(results[sugar]["dG_kJ"] - best, 2)

    # Sort by predicted affinity (most negative first)
    ranked = sorted(results.keys(), key=lambda s: results[s]["dG_kJ"])

    return {
        "pdb_id": pdb_id,
        "template_sugar": sugar_resname,
        "ranking": ranked,
        "scores": results,
        "notes": "ΔΔG relative to best binder. Negative dG = favorable binding.",
    }


# =====================================================================
# Validation against experimental data
# =====================================================================

# Experimental ITC data for validation
# Source: hand-curated from pdb_lectin_reference.py
VALIDATION_DATA = {
    "5CNA_MMA": {"dG_exp_kJ": -22.3, "lectin": "ConA", "sugar": "MeAlphaMan"},
    "3ZSJ_GAL": {"dG_exp_kJ": -26.0, "lectin": "Galectin-3", "sugar": "LacNAc (Gal)"},
}


def validate_predictions(cache_dir: Optional[str] = None) -> dict:
    """Score reference structures and compare to experimental ΔG."""
    params = GlycanParams.from_calibrated()
    results = []

    test_cases = [
        ("3ZSJ", "GAL", None),
        ("5CNA", "MMA", "A"),
    ]

    for pdb_id, resname, chain in test_cases:
        key = f"{pdb_id}_{resname}"
        exp = VALIDATION_DATA.get(key)
        if not exp:
            continue

        try:
            cms = extract_contacts(pdb_id, sugar_resname=resname,
                                   chain=chain, cache_dir=cache_dir)
            if not cms:
                continue
            sr = score_contact_map(cms[0], params)

            results.append({
                "pdb": pdb_id,
                "sugar": resname,
                "lectin": exp["lectin"],
                "pred_kJ": sr.dG_predicted_kJ,
                "exp_kJ": exp["dG_exp_kJ"],
                "error_kJ": round(sr.dG_predicted_kJ - exp["dG_exp_kJ"], 2),
                "components": sr.components,
            })
        except Exception as e:
            results.append({"pdb": pdb_id, "error": str(e)})

    if results:
        errors = [abs(r["error_kJ"]) for r in results if "error_kJ" in r]
        mae = sum(errors) / len(errors) if errors else float("inf")
    else:
        mae = float("inf")

    return {"n_pairs": len(results), "MAE_kJ": round(mae, 2), "results": results}


if __name__ == "__main__":
    cache = "/home/claude/pdb_cache"

    print("Glycan Scorer v2 — Physics-Based from PDB Contacts")
    print("=" * 65)

    # ── Validate ─────────────────────────────────────────────
    print("\n--- Validation vs Experimental ITC ---")
    val = validate_predictions(cache_dir=cache)
    print(f"MAE: {val['MAE_kJ']:.2f} kJ/mol ({val['n_pairs']} pairs)")
    for r in val["results"]:
        if "error_kJ" in r:
            print(f"  {r['lectin']:15s} {r['sugar']:5s}: "
                  f"pred={r['pred_kJ']:+7.1f} exp={r['exp_kJ']:+7.1f} "
                  f"err={r['error_kJ']:+6.1f} kJ/mol")
            for k, v in r["components"].items():
                if abs(v) > 0.01:
                    print(f"    {k:8s}: {v:+7.2f}")

    # ── Deoxy-mannose predictions ────────────────────────────
    print("\n--- Prediction 1: ConA Deoxy-Mannose Series ---")
    deoxy = predict_deoxy_series("5CNA", "MMA", chain="A", cache_dir=cache)
    if "error" not in deoxy:
        print(f"Parent: {deoxy['parent']['sugar']} in {deoxy['parent']['pdb_id']}")
        print(f"Parent ΔG = {deoxy['parent']['dG_parent_kJ']:.1f} kJ/mol")
        print(f"{'Position':<10} {'ax/eq':<6} {'HB lost':<8} {'ΔΔG (kJ/mol)':<14} {'Effect'}")
        print("-" * 55)
        for pos in sorted(deoxy["deoxy_predictions"]):
            p = deoxy["deoxy_predictions"][pos]
            print(f"{p['position']:<10} {p['axial_equatorial']:<6} "
                  f"{p['hbonds_lost']:<8} {p['ddG_predicted_kJ']:+8.1f}      "
                  f"{p['interpretation']}")
    else:
        print(f"  Error: {deoxy['error']}")

    # ── Selectivity panel ────────────────────────────────────
    print("\n--- Prediction 2: ConA Sugar Selectivity ---")
    sel = score_sugar_panel("5CNA", "MMA", chain="A",
                            panel=["Man", "Glc", "Gal", "GlcNAc"],
                            cache_dir=cache)
    if "error" not in sel:
        print(f"Ranking (best → worst): {' > '.join(sel['ranking'])}")
        print(f"{'Sugar':<10} {'ΔG (kJ/mol)':<14} {'ΔΔG':<10} {'HB':<5} {'CH-π':<6}")
        print("-" * 48)
        for sugar in sel["ranking"]:
            s = sel["scores"][sugar]
            print(f"{sugar:<10} {s['dG_kJ']:+8.1f}     {s['ddG_kJ']:+6.1f}    "
                  f"{s['n_hbonds']:<5} {s['n_ch_pi']:<6}")

    # ── Galectin-3 selectivity ───────────────────────────────
    print("\n--- Prediction 6: Galectin-3 Sugar Selectivity ---")
    sel2 = score_sugar_panel("3ZSJ", "GAL",
                             panel=["Man", "Glc", "Gal", "GlcNAc"],
                             cache_dir=cache)
    if "error" not in sel2:
        print(f"Ranking: {' > '.join(sel2['ranking'])}")
        for sugar in sel2["ranking"]:
            s = sel2["scores"][sugar]
            print(f"  {sugar:<10} ΔG={s['dG_kJ']:+8.1f}  ΔΔG={s['ddG_kJ']:+6.1f}")