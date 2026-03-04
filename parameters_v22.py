"""
MABE Glycan Module — Parameter Set v2.2
Sources documented inline. All values in kJ/mol.
Anti-fabrication: every value has explicit provenance.
"""

# ============================================================
# ENERGY PARAMETERS (kJ/mol)
# ============================================================

EPS_HB = -5.0
# Source: Literature consensus (Fersht 1985, Pace 2014).
# Status: LOCKED v2.1

EPS_CH_PI = -2.5
# Source: Asensio 2013 (Acc. Chem. Res. 46:946) glycopeptide model = -3.3 kJ/mol;
#         Diehl 2024 (JACS Au) Gal3 W181 ITC differential ~-1.5 kJ/mol per contact;
#         MMGBSA estimate ~-3.6 kJ/mol per contact (known to overestimate).
#         -2.5 sits between ITC lower bound and MMGBSA upper bound.
# Status: CONFIRMED v2.2

BETA_CONTEXT = 0.45
# Source: GLYCAM06 QM torsion analysis (Kirschner 2008, J. Comput. Chem. 29:622)
# Accounts for partial solvation of OH in binding pocket vs bulk solvent.
# Status: LOCKED v2.1

# ============================================================
# DESOLVATION PARAMETERS (kJ/mol, cost of burying one OH)
# ============================================================

K_DESOLV_EQ = 2.4
# Source: Schwarz 1996, J. Solution Chem. 25:471, Table I, DMG buffer pH 6.95.
#         Mean of C1 (+2.83) and C2 (+2.0) equatorial positions.
#         C3 anomalous (negative, excluded); C4 absent (only 4F-Glu present).
# Status: LOCKED v2.2; uncertainty ±0.6 kJ/mol

K_DESOLV_AX = 6.3
# Source: Jasra & Ahluwalia 1982, J. Solution Chem. 11:325, Table I.
#         Derived from Gal - 2-deoxy-Gal: 17.20 - 10.91 = +6.29 kJ/mol.
#         Axial C2-OH of galactose removed. Conservative (raw Jasra, no offset correction).
#         Uncertainty ±1.5 kJ/mol (method offset sign-inconsistent across positions).
# Axial/equatorial ratio: 6.3/2.4 = 2.6× (physically reasonable).
# Status: PROVISIONAL v2.2; use Schwarz-corrected 7.0 as sensitivity check.

K_DESOLV_C6 = 11.2
# Source: Schwarz 1996, Table I. Primary –CH₂OH group.
#         ΔΔH = αGlu (12.2) - 6HGlu (0.98) = +11.22 kJ/mol.
#         C6 is conformationally flexible — use as separate pool, not pooled with ring OHs.
# Status: LOCKED v2.2

K_DESOLV_NAC = 8.5
# Source: NOT DIRECTLY MEASURED. Estimate: K_DESOLV_EQ + acetamide desolvation penalty.
#         ΔH_sol(acetamide) ≈ -6 to -8 kJ/mol (NIST WebBook); net desolv cost higher.
#         Bracketed estimate: 8-10 kJ/mol. Central value 8.5 used.
# Status: ESTIMATE — flag all GlcNAc predictions as LOW CONFIDENCE until locked.

# Per-position map: which k_desolv to use
# Equatorial ring OH: K_DESOLV_EQ
# Axial ring OH: K_DESOLV_AX
# Primary C6-OH (–CH₂OH): K_DESOLV_C6
# N-acetyl group (C2-NHAc): K_DESOLV_NAC

# ============================================================
# SCORING FORMULA
# ============================================================
# ΔG_pred = Σ_HB [EPS_HB × BETA_CONTEXT] 
#          + Σ_CH-π [EPS_CH_PI]
#          + Σ_buried_OH [k_desolv_i]   ← desolvation PENALTY (positive)
#
# Note: HB term is negative (favorable); desolvation is positive (unfavorable).
# Net: ΔG_pred = favorable contacts − desolvation cost
#
# Final: ΔG_pred = n_HB × EPS_HB × BETA_CONTEXT + n_CHP × EPS_CH_PI + Σ k_desolv_i

def score(n_HB, n_CHP, buried_ohs):
    """
    n_HB: number of H-bonds formed
    n_CHP: number of CH-pi contacts
    buried_ohs: list of k_desolv values for each buried OH
                e.g. [K_DESOLV_EQ, K_DESOLV_EQ, K_DESOLV_C6]
    Returns ΔG_pred in kJ/mol (negative = favorable)
    """
    hb_term  = n_HB  * EPS_HB  * BETA_CONTEXT
    chp_term = n_CHP * EPS_CH_PI
    desolv   = sum(buried_ohs)
    return hb_term + chp_term + desolv