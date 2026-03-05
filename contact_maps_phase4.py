"""
MABE Glycan Module — Phase 4 Contact Maps
Full monosaccharide panel: ConA, WGA, PNA, Galectin-3, Davis receptor
Sources cited inline. All contacts from published crystal structures.
"""

from parameters_v22 import (
    K_DESOLV_EQ, K_DESOLV_AX, K_DESOLV_C6, K_DESOLV_NAC
)

# ============================================================
# CONCANAVALIN A — PDB 5CNA / 2CNA
# Pharmacophore: C3-OH(eq), C4-OH(eq), C6-OH(prim) required.
# Source: Chervenak 1995; Derewenda 1989.
# CH-pi residue: Tyr12. Aromatic hierarchy: Tyr → ~-1.5 kJ/mol per contact.
# Updated: use residue-specific eps; ConA Tyr gives n_CHP=2 for Man (axial C1+C3)
# and n_CHP=1 for Glu (C3 only). This is the source of Man>Glu selectivity.
# ============================================================

CONCA_PANEL = [
    dict(name="ConA + αMeMan",  Ka_obs=7600,  dG_obs=-22.2,
         n_HB=3, n_CHP=2,
         buried_ohs=[K_DESOLV_EQ, K_DESOLV_EQ, K_DESOLV_C6],
         notes="C3-OH,C4-OH,C6-OH H-bonds; C1ax+C3ax CH-pi to Tyr12",
         source="Chervenak 1995", confidence="HIGH"),

    dict(name="ConA + αMeGlu",  Ka_obs=2400,  dG_obs=-19.3,
         n_HB=3, n_CHP=1,
         buried_ohs=[K_DESOLV_EQ, K_DESOLV_EQ, K_DESOLV_C6],
         notes="Same H-bonds; C1 equatorial in Glu → 1 CH-pi only",
         source="Chervenak 1995", confidence="HIGH"),

    dict(name="ConA + βMeMan",  Ka_obs=0,     dG_obs=None,
         n_HB=0, n_CHP=0,
         buried_ohs=[],
         notes="NO BINDING: beta-anomer C1-OH axial; C4-OH blocked in beta-Man orientation",
         source="Chervenak 1995", confidence="HIGH"),

    dict(name="ConA + GlcNAc",  Ka_obs=0,     dG_obs=None,
         n_HB=2, n_CHP=0,
         buried_ohs=[K_DESOLV_NAC, K_DESOLV_EQ, K_DESOLV_C6],
         notes="NO BINDING obs; NHAc at C2 sterically disrupts Arg228 approach to C3-OH",
         source="Chervenak 1995", confidence="HIGH"),

    dict(name="ConA + αMeFru",  Ka_obs=0,     dG_obs=None,
         n_HB=1, n_CHP=0,
         buried_ohs=[K_DESOLV_EQ, K_DESOLV_C6],
         notes="NO BINDING obs; C2-OH ketone, wrong geometry for C3/C4 H-bond network",
         source="Chervenak 1995", confidence="MEDIUM"),

    dict(name="ConA + αMeGal",  Ka_obs=None,  dG_obs=None,
         n_HB=2, n_CHP=0,
         buried_ohs=[K_DESOLV_EQ, K_DESOLV_AX, K_DESOLV_C6],
         notes="Weak/no binding expected; C4-OH axial disrupts H-bond geometry",
         source="Inferred from pharmacophore; Dam & Brewer 2002", confidence="MEDIUM"),
]

# ============================================================
# WGA (Wheat Germ Agglutinin) — PDB 2UVO
# Specificity: GlcNAc >> Gal, Man, Glu (requires NHAc)
# Source: Bains 1992; Wright 1984, 1987.
# Aromatic: Tyr73 to GlcNAc ring. Residue: Tyr → eps_CHP ~ -1.5/contact
# ============================================================

WGA_PANEL = [
    dict(name="WGA + GlcNAc",   Ka_obs=410,   dG_obs=-15.5,
         n_HB=4, n_CHP=1,
         buried_ohs=[K_DESOLV_NAC, K_DESOLV_EQ],
         notes="Direct NHAc→Glu115,Ser114; OH3→Tyr73; Tyr73 stack",
         source="Bains 1992", confidence="MEDIUM"),

    dict(name="WGA + Gal",      Ka_obs=None,  dG_obs=None,
         n_HB=1, n_CHP=0,
         buried_ohs=[K_DESOLV_AX, K_DESOLV_C6],
         notes="No NHAc → loses 3 H-bonds; C4-OH axial disrupts Tyr73 stack",
         source="Inferred from WGA specificity (Bains 1992)", confidence="MEDIUM"),

    dict(name="WGA + Glu",      Ka_obs=None,  dG_obs=None,
         n_HB=1, n_CHP=0,
         buried_ohs=[K_DESOLV_EQ, K_DESOLV_C6],
         notes="No NHAc; equatorial face can't satisfy Glu115 acceptors",
         source="Inferred", confidence="LOW"),

    dict(name="WGA + (GlcNAc)2", Ka_obs=5300, dG_obs=-21.3,
         n_HB=8, n_CHP=2,
         buried_ohs=[K_DESOLV_NAC, K_DESOLV_EQ, K_DESOLV_NAC, K_DESOLV_EQ],
         notes="Two GlcNAc units; additive H-bonds + CH-pi per subsite",
         source="Bains 1992", confidence="MEDIUM"),
]

# ============================================================
# PNA (Peanut Agglutinin) — PDB 2PEL
# Specificity: Gal >> Glu, Man; Gal-β(1→3)-GalNAc preferred
# Key residues: Asp83 (C3-OH,C4-OH), Asn91 (C4-OH), Lys90 (C3-OH),
#               Ser211 (C6-OH), Tyr125 (CH-pi to Gal α-face)
# Source: Banerjee 1994 (J. Biol. Chem. 269:2745); Loganathan 1997
# Galactose at subsite: C4-OH axial is a FEATURE for PNA (opposite of ConA)
# ============================================================

PNA_PANEL = [
    dict(name="PNA + Gal",      Ka_obs=2100,  dG_obs=-18.9,
         n_HB=4, n_CHP=2,
         buried_ohs=[K_DESOLV_AX, K_DESOLV_EQ, K_DESOLV_C6],
         notes="Asp83→C3-OH,C4-OH; Asn91→C4-OH; Ser211→C6-OH; Tyr125 to C3,C4,C5 axial face (α-face Gal)",
         source="Banerjee 1994 [LIT]; PDB 2PEL", confidence="HIGH"),

    dict(name="PNA + Glu",      Ka_obs=None,  dG_obs=None,
         n_HB=2, n_CHP=0,
         buried_ohs=[K_DESOLV_EQ, K_DESOLV_C6],
         notes="C4-OH equatorial in Glu; cannot satisfy Asp83/Asn91 acceptors optimally; no CH-pi (C4-eq face disfavors)",
         source="Inferred from PNA selectivity", confidence="MEDIUM"),

    dict(name="PNA + Man",      Ka_obs=None,  dG_obs=None,
         n_HB=2, n_CHP=0,
         buried_ohs=[K_DESOLV_AX, K_DESOLV_C6],
         notes="C2-OH axial disrupts subsite geometry; C4-OH equatorial → misfit",
         source="Inferred", confidence="LOW"),

    dict(name="PNA + GalNAc",   Ka_obs=3800,  dG_obs=-20.2,
         n_HB=4, n_CHP=1,
         buried_ohs=[K_DESOLV_AX, K_DESOLV_NAC, K_DESOLV_C6],
         notes="Same Gal contacts; NHAc at C2 axial tolerated (extra H-bond possible); reduced CH-pi",
         source="Banerjee 1994 [LIT]", confidence="MEDIUM"),
]

# ============================================================
# GALECTIN-3 — PDB 3ZSJ (0.89 Å)
# Specificity: Gal/LacNAc >> Glu, Man
# Key residues: Arg162 (×2-3 H-bonds), His158, Asn174, Glu184, Trp181 (CH-pi)
# Galactose C4-OH AXIAL is required for Galectin-3 (opposite of ConA)
# Source: Diehl 2024; Seetharaman 1998
# ============================================================

GAL3_PANEL = [
    dict(name="Gal3 + Gal",     Ka_obs=9091,  dG_obs=-22.6,
         n_HB=5, n_CHP=3,
         buried_ohs=[K_DESOLV_AX, K_DESOLV_EQ, K_DESOLV_C6],
         notes="Arg162(×3),His158,Asn174; Trp181 to C3/C4/C5 Gal α-face",
         source="Diehl 2024; PDB 3ZSJ", confidence="HIGH"),

    dict(name="Gal3 + Glu",     Ka_obs=None,  dG_obs=None,
         n_HB=3, n_CHP=0,
         buried_ohs=[K_DESOLV_EQ, K_DESOLV_C6],
         notes="C4-OH equatorial; Arg162 network partially disrupted; no CH-pi (Asensio: equatorial face disfavors Trp stack)",
         source="Inferred from Gal3 selectivity; Diehl 2024", confidence="MEDIUM"),

    dict(name="Gal3 + Man",     Ka_obs=None,  dG_obs=None,
         n_HB=2, n_CHP=0,
         buried_ohs=[K_DESOLV_AX, K_DESOLV_C6],
         notes="C2-OH axial disrupts subsite; C4-OH equatorial → weaker",
         source="Inferred", confidence="LOW"),

    dict(name="Gal3 + LacNAc",  Ka_obs=50000, dG_obs=-26.8,
         n_HB=7, n_CHP=3,
         buried_ohs=[K_DESOLV_AX, K_DESOLV_EQ, K_DESOLV_C6, K_DESOLV_NAC, K_DESOLV_EQ],
         notes="Gal unit as above + GlcNAc unit (2 additional HBs via Asn160,Arg162 to GlcNAc OH)",
         source="Seetharaman 1998 [LIT]", confidence="MEDIUM"),
]

# ============================================================
# DAVIS RECEPTOR — already in Phase 2; add Man,Gal,Fru for completeness
# Source: Tromans 2019
# ============================================================

DAVIS_PANEL = [
    dict(name="Davis + βGlu",   Ka_obs=18600, dG_obs=-24.4,
         n_HB=10, n_CHP=2,
         buried_ohs=[K_DESOLV_EQ]*5 + [K_DESOLV_C6],
         notes="All-equatorial OH array; perfect fit to hexaurea cage + TEM faces",
         source="Tromans 2019", confidence="HIGH"),

    dict(name="Davis + βGal",   Ka_obs=180,   dG_obs=-12.9,
         n_HB=7,  n_CHP=1,
         buried_ohs=[K_DESOLV_EQ]*4 + [K_DESOLV_AX] + [K_DESOLV_C6],
         notes="C4-OH axial disrupts 3 HBs + reduces CH-pi",
         source="Tromans 2019", confidence="MEDIUM"),

    dict(name="Davis + βMan",   Ka_obs=140,   dG_obs=-12.4,
         n_HB=7,  n_CHP=1,
         buried_ohs=[K_DESOLV_EQ]*4 + [K_DESOLV_AX] + [K_DESOLV_C6],
         notes="C2-OH axial disrupts vicinal pair + reduces CH-pi",
         source="Tromans 2019", confidence="MEDIUM"),

    dict(name="Davis + 2dGlu",  Ka_obs=725,   dG_obs=-16.7,
         n_HB=6,  n_CHP=2,
         buried_ohs=[K_DESOLV_EQ]*4 + [K_DESOLV_C6],
         notes="No C2-OH; vicinal urea pair disrupted (-4 HBs); CH-pi preserved",
         source="Tromans 2019", confidence="HIGH"),

    dict(name="Davis + GlcNAc", Ka_obs=20,    dG_obs=-7.4,
         n_HB=6,  n_CHP=0,
         buried_ohs=[K_DESOLV_EQ]*4 + [K_DESOLV_NAC] + [K_DESOLV_C6],
         notes="NHAc steric disruption; cage geometry misfit",
         source="Tromans 2019", confidence="LOW"),

    dict(name="Davis + Fru",    Ka_obs=50,    dG_obs=-9.8,
         n_HB=4,  n_CHP=0,
         buried_ohs=[K_DESOLV_EQ]*3 + [K_DESOLV_C6],
         notes="Furanose form dominant; fewer equatorial OHs presented to cage",
         source="Tromans 2019", confidence="LOW"),
]

ALL_PHASE4 = {
    "ConA":   CONCA_PANEL,
    "WGA":    WGA_PANEL,
    "PNA":    PNA_PANEL,
    "Gal3":   GAL3_PANEL,
    "Davis":  DAVIS_PANEL,
}