"""
Glycan parameter integration module.

Computes MABE glycan scoring parameters from Tier 1 data sources:
  - CHI Energy Functions (glycam.org)    → P4 (eps_glycosidic_freeze)
  - CH-pi energetics (literature survey) → P6 (eps_CH_pi_pyranose), P7 (furanose)
  - PDB lectin structures                → P9 (water_bridge_norm), P12 (multivalency)

PARAMETER MAP:
  P1  k_desolv_eq       — equatorial OH desolvation    [PAYWALLED: Schwarz 1996]
  P2  k_desolv_ax       — axial OH desolvation          [PAYWALLED: Schwarz 1996]
  P3  k_desolv_NAc      — NAc group desolvation         [HIGH RISK: single source]
  P4  eps_glycosidic     — conformational entropy cost   [THIS MODULE: from CHI]
  P5  TdS_ring_pucker    — ring pucker entropy           [PAYWALLED: Peric-Hassler]
  P6  eps_CH_pi_pyranose — pyranose CH-pi per contact    [THIS MODULE: from literature]
  P7  eps_CH_pi_furanose — furanose CH-pi per contact    [LOW CONFIDENCE]
  P8  eps_water_bridge   — per-water bridge energy       [PAYWALLED: mutant ITC]
  P9  n_water_bridge_norm— waters per interface area     [THIS MODULE: from PDB]
  P10 eps_Ca_coord       — Ca2+ coordination energy      [residual from metal scorer]
  P11 alpha_multivalent  — multivalent enhancement exp   [PAYWALLED: Dam & Brewer]
  P12 d_optimal_multi    — optimal inter-site distance   [THIS MODULE: from PDB]
"""

from typing import Dict, Optional, Tuple
from dataclasses import dataclass

from data.chi_energy_functions import (
    score_linkage_entropy, get_chi_functions, LINKAGE_PSI_CLASS,
)
from data.ch_pi_energetics import (
    get_eps_CH_pi_estimates, estimate_CH_pi_energy, SUGAR_CH_PI_CONTACTS,
    LAUGHREY_PER_CONTACT,
)
from data.pdb_lectin_reference import (
    get_water_bridge_norm, get_inter_site_distances, LECTIN_REFERENCE_SET,
    BOUND_STATE_TORSIONS, WATER_BRIDGE_STATISTICS,
)

KCAL_TO_KJ = 4.184


@dataclass
class GlycanParams:
    """All glycan scoring parameters with provenance tracking."""
    # Desolvation (Tier 2 — paywalled, placeholder values)
    P1_k_desolv_eq_kJ: float = -5.0       # kJ/mol per equatorial OH buried
    P2_k_desolv_ax_kJ: float = -3.5       # kJ/mol per axial OH buried
    P3_k_desolv_NAc_kJ: float = -4.0      # kJ/mol per NAc buried

    # Conformational entropy (Tier 1 — from CHI)
    P4_eps_glycosidic_kJ: float = 0.0      # Computed per-linkage
    P5_TdS_ring_pucker_kJ: float = 2.5     # Default: ~2.5 kJ/mol per ring

    # CH-pi (Tier 1 — from literature)
    P6_eps_CH_pi_pyranose_kJ: float = -3.0 # kJ/mol per effective contact
    P7_eps_CH_pi_furanose_kJ: float = -2.0 # Low confidence

    # Water bridges (Tier 1/2 hybrid)
    P8_eps_water_bridge_kJ: float = -5.0   # kJ/mol per bridging water
    P9_n_water_norm: float = 0.01          # waters per Å² of interface

    # Metal coordination
    P10_eps_Ca_coord_kJ: float = -8.0      # kJ/mol per Ca-OH coordination

    # Multivalency
    P11_alpha_multivalent: float = 1.5     # Enhancement exponent
    P12_d_optimal_multi_A: float = 40.0    # Optimal inter-site distance

    # Offset
    dG_0_kJ: float = 0.0

    def provenance(self) -> Dict[str, str]:
        """Return data source for each parameter."""
        return {
            "P1": "PLACEHOLDER — needs Schwarz 1996 (paywalled)",
            "P2": "PLACEHOLDER — needs Schwarz 1996 (paywalled)",
            "P3": "PLACEHOLDER — needs verification (HIGH RISK)",
            "P4": "CHI Energy Functions (glycam.org), Nivedha et al. 2014",
            "P5": "DEFAULT — needs Peric-Hassler 2010 (paywalled)",
            "P6": "Laughrey 2008 + Keys 2025 + Hudson 2015",
            "P7": "EXTRAPOLATED from P6 (LOW CONFIDENCE)",
            "P8": "PLACEHOLDER — needs lectin mutant ITC data",
            "P9": "PDB lectin structures (UniLectin3D reference set)",
            "P10": "PLACEHOLDER — residual from MABE metal scorer",
            "P11": "PLACEHOLDER — needs Dam & Brewer 2002 (paywalled)",
            "P12": "PDB lectin oligomer geometry",
        }


def compute_params_from_tier1() -> GlycanParams:
    """Compute all parameters derivable from Tier 1 (free) data sources.

    Returns GlycanParams with Tier 1 values filled in and placeholders
    for Tier 2+ parameters.
    """
    params = GlycanParams()

    # ── P4: Conformational entropy from CHI ──────────────────
    # Compute average TdS_freeze across common linkage/bound-state pairs
    tds_values = []
    for linkage_key, entries in BOUND_STATE_TORSIONS.items():
        linkage_str, sugar = linkage_key
        # Parse linkage: "alpha-1,3" -> anomer="alpha", pos=3
        parts = linkage_str.split("-")
        anomer = parts[0]
        pos = int(parts[1].split(",")[1])

        # Map sugar to reducing-end sugar for CHI lookup
        sugar_map = {
            "Man": "Man", "Glc": "Glc", "Gal": "Gal",
            "GlcNAc": "GlcNAc", "GalNAc": "GalNAc",
        }
        reducing = sugar_map.get(sugar, "Glc")

        for phi, psi, pdb_id, lectin in entries:
            try:
                result = score_linkage_entropy(anomer, pos, reducing, phi, psi)
                tds_values.append(result["TdS_freeze_kJ"])
            except (ValueError, KeyError):
                pass

    if tds_values:
        params.P4_eps_glycosidic_kJ = round(sum(tds_values) / len(tds_values), 2)

    # ── P6: CH-pi from literature consensus ──────────────────
    estimates = get_eps_CH_pi_estimates()
    # Use Trp-weighted average since Trp is 9x enriched at lectin sites
    params.P6_eps_CH_pi_pyranose_kJ = estimates["eps_CH_pi_pyranose_avg"]["value_kJ"]
    params.P7_eps_CH_pi_furanose_kJ = estimates["eps_CH_pi_furanose"]["value_kJ"]

    # ── P9: Water bridge normalization from PDB ──────────────
    p9_data = get_water_bridge_norm()
    params.P9_n_water_norm = p9_data["P9_mean_waters_per_A2"]

    # ── P12: Multivalency distances from PDB ─────────────────
    distances = get_inter_site_distances()
    all_d = []
    for lectin, data in distances.items():
        all_d.extend(data["distances_A"])
    if all_d:
        params.P12_d_optimal_multi_A = round(sum(all_d) / len(all_d), 1)

    return params


def score_monosaccharide(sugar: str, lectin_pdb: str,
                         params: Optional[GlycanParams] = None) -> dict:
    """Score a monosaccharide binding to a lectin binding site.

    This is a simplified scorer for monosaccharide-lectin pairs.
    For oligosaccharides, use the full glycan_scorer pipeline.

    Args:
        sugar: Sugar identity (e.g., "alpha-Man", "beta-Gal")
        lectin_pdb: PDB ID of lectin structure
        params: GlycanParams to use. If None, computes from Tier 1.

    Returns:
        Dict with predicted dG and component breakdown.
    """
    if params is None:
        params = compute_params_from_tier1()

    # Look up lectin structure info
    lectin_info = None
    for ls in LECTIN_REFERENCE_SET:
        if ls.pdb_id == lectin_pdb:
            lectin_info = ls
            break

    if lectin_info is None:
        return {"error": f"PDB {lectin_pdb} not in reference set"}

    site = lectin_info.sites[0] if lectin_info.sites else None
    if site is None:
        return {"error": "No binding site data"}

    # ── Component 1: CH-pi stacking ──────────────────────────
    dG_CH_pi = 0.0
    if site.aromatic_residues:
        # Pick dominant aromatic
        dominant = "Trp" if any("Trp" in r for r in site.aromatic_residues) else \
                   "Tyr" if any("Tyr" in r for r in site.aromatic_residues) else "Phe"

        ch_pi = estimate_CH_pi_energy(sugar, dominant)
        dG_CH_pi = ch_pi["dG_CH_pi_kJ"]

    # ── Component 2: Water bridges ───────────────────────────
    wb_data = WATER_BRIDGE_STATISTICS.get(lectin_pdb, {})
    n_waters = wb_data.get("n_waters", 2)  # default 2
    dG_water = params.P8_eps_water_bridge_kJ * n_waters

    # ── Component 3: Metal coordination ──────────────────────
    dG_metal = 0.0
    if site.metal_ions:
        n_metal_bridges = len(site.metal_ions)
        dG_metal = params.P10_eps_Ca_coord_kJ * n_metal_bridges

    # ── Component 4: Desolvation penalty ─────────────────────
    # Estimate number of buried OH groups (simplified)
    sugar_contacts = SUGAR_CH_PI_CONTACTS.get(sugar, {})
    best_face = max(sugar_contacts.values(), key=lambda x: x["n_CH"]) if sugar_contacts else {"n_CH": 3}
    n_buried_OH = max(1, 5 - best_face["n_CH"])  # rough: non-CH positions have OH
    dG_desolv = params.P1_k_desolv_eq_kJ * n_buried_OH  # simplified: all equatorial

    # ── Total ────────────────────────────────────────────────
    dG_total = dG_CH_pi + dG_water + dG_metal + dG_desolv + params.dG_0_kJ

    return {
        "dG_predicted_kJ": round(dG_total, 2),
        "components": {
            "CH_pi_kJ": round(dG_CH_pi, 2),
            "water_bridge_kJ": round(dG_water, 2),
            "metal_coord_kJ": round(dG_metal, 2),
            "desolvation_kJ": round(dG_desolv, 2),
        },
        "dG_exp_kJ": site.dG_exp_kJ,
        "sugar": sugar,
        "lectin": lectin_info.name,
        "pdb_id": lectin_pdb,
        "params_provenance": "Tier 1 (CHI + CH-pi literature + PDB)",
    }


def validate_against_reference() -> dict:
    """Score all reference lectin-sugar pairs and compare to experimental dG.

    Returns prediction accuracy metrics.
    """
    params = compute_params_from_tier1()
    results = []

    test_pairs = [
        ("alpha-Man", "5CNA"),   # ConA + mannose
        ("alpha-Man", "1CVN"),   # ConA + MeAlphaMan
        ("beta-Gal", "3ZSJ"),   # Galectin-3 + LacNAc (Gal component)
        ("alpha-Man", "1SL5"),   # DC-SIGN + Man
    ]

    for sugar, pdb in test_pairs:
        r = score_monosaccharide(sugar, pdb, params)
        if "error" not in r and r["dG_exp_kJ"] is not None:
            pred = r["dG_predicted_kJ"]
            exp = r["dG_exp_kJ"]
            error = pred - exp
            results.append({
                "sugar": sugar,
                "pdb": pdb,
                "pred_kJ": pred,
                "exp_kJ": exp,
                "error_kJ": round(error, 2),
            })

    if results:
        errors = [abs(r["error_kJ"]) for r in results]
        mae = sum(errors) / len(errors)
    else:
        mae = float("inf")

    return {
        "n_pairs": len(results),
        "MAE_kJ": round(mae, 2),
        "results": results,
        "params_used": {
            "P4": params.P4_eps_glycosidic_kJ,
            "P6": params.P6_eps_CH_pi_pyranose_kJ,
            "P9": params.P9_n_water_norm,
            "P12": params.P12_d_optimal_multi_A,
        },
        "notes": "Pre-fitting validation. P1/P2/P3/P8 are placeholders. "
                 "MAE will improve once paywalled data fills those parameters."
    }


if __name__ == "__main__":
    print("Glycan Parameter Integration — Tier 1 Data")
    print("=" * 60)

    params = compute_params_from_tier1()
    print("\nComputed Parameters:")
    print(f"  P4  eps_glycosidic:      {params.P4_eps_glycosidic_kJ:7.2f} kJ/mol")
    print(f"  P6  eps_CH_pi_pyranose:  {params.P6_eps_CH_pi_pyranose_kJ:7.2f} kJ/mol")
    print(f"  P7  eps_CH_pi_furanose:  {params.P7_eps_CH_pi_furanose_kJ:7.2f} kJ/mol")
    print(f"  P9  n_water_norm:        {params.P9_n_water_norm:7.4f} waters/Å²")
    print(f"  P12 d_optimal_multi:     {params.P12_d_optimal_multi_A:7.1f} Å")

    print("\nPlaceholder Parameters (need Tier 2 data):")
    print(f"  P1  k_desolv_eq:         {params.P1_k_desolv_eq_kJ:7.2f} kJ/mol")
    print(f"  P2  k_desolv_ax:         {params.P2_k_desolv_ax_kJ:7.2f} kJ/mol")
    print(f"  P3  k_desolv_NAc:        {params.P3_k_desolv_NAc_kJ:7.2f} kJ/mol")
    print(f"  P8  eps_water_bridge:    {params.P8_eps_water_bridge_kJ:7.2f} kJ/mol")
    print(f"  P10 eps_Ca_coord:        {params.P10_eps_Ca_coord_kJ:7.2f} kJ/mol")

    print("\n\nProvenance:")
    for k, v in params.provenance().items():
        print(f"  {k}: {v}")

    print("\n\nValidation against reference structures:")
    val = validate_against_reference()
    print(f"  Pairs scored: {val['n_pairs']}")
    print(f"  MAE: {val['MAE_kJ']:.2f} kJ/mol")
    for r in val["results"]:
        print(f"    {r['sugar']:15s} in {r['pdb']}: pred={r['pred_kJ']:7.2f} "
              f"exp={r['exp_kJ']:7.2f} error={r['error_kJ']:+7.2f} kJ/mol")