"""
contact_maps_predictions.py — Prediction-ready ContactMap objects built from PDB analysis.

Each ContactMap represents one lectin-sugar pair with consensus contacts
averaged across crystallographic chains. These are the INPUTS to the scorer.

Sources:
  - contact_maps_auto.csv (BioPython analysis with geometry-verified CH-π)
  - Manual hydroxyl classification from sugar_properties.py

No parameter values here. No fitting. Just structural contact counts.
"""

from glycan_scorer import ContactMap


# ═══════════════════════════════════════════════════════════════════════════
# PREDICTION 1: ConA deoxy-mannose series
# PDB: 5CNA (ConA + MeαMan, 2.0 Å)
# Contact consensus across 4 chains: 10.5 HB, 1 CH-π (TYR12), 0.75 waters
# ═══════════════════════════════════════════════════════════════════════════

# Parent: MeαMan in ConA monosaccharide site
# H-bond assignment per OH from 5CNA crystal contacts:
#   C2-OH(ax): Asp208 Oδ (1 HB)
#   C3-OH(eq): Asn14 Nδ, Leu99 backbone O (2 HB)
#   C4-OH(eq): Asp208 Oδ, Tyr100 OH (2 HB)
#   C6-OH:     Arg228 NH, Asp208 Oδ (2 HB)
#   Non-OH HB: ring O5, OMe, backbone contacts (3 HB)
# Total: 1+2+2+2+3 = 10 HB ✓
#
# CRITICAL RULE: when removing a deoxy OH, reduce BOTH n_OH_*_buried
# AND n_HB by the per-OH HB count. Otherwise desolvation saving
# without H-bond loss yields unphysical ΔΔG < 0.

CONA_MEAMAN = ContactMap(
    lectin="ConA", sugar="MeαMan", pdb_id="5CNA",
    n_HB=10,           # consensus across 4 chains (10,10,12,10)
    n_CH_pi=1,         # TYR12, all 4 chains
    n_OH_eq_buried=3,  # C3-OH(eq), C4-OH(eq), C6-OH(eq)
    n_OH_ax_buried=1,  # C2-OH(ax) — H-bond to Asp208
    n_NAc_buried=0,
    n_frozen_torsions=0,
    n_water_bridges=1,
    resolution_A=2.0,
    notes="Consensus 4 chains. C1-OMe not counted as OH.",
)

# Deoxy derivatives: remove OH → lose BOTH its H-bonds AND its desolvation
CONA_2DEOXY_MAN = ContactMap(
    lectin="ConA", sugar="2-deoxy-Man", pdb_id="5CNA",
    n_HB=9,            # lose 1 HB (C2-OH → Asp208)
    n_CH_pi=1,
    n_OH_eq_buried=3,  # unchanged
    n_OH_ax_buried=0,  # C2-OH(ax) removed
    n_NAc_buried=0,
    n_frozen_torsions=0,
    n_water_bridges=1,
    resolution_A=2.0,
    notes="C2-OH removed. Lose 1 HB + 1 ax desolv. Net ΔΔG small (C2 not essential).",
)

CONA_3DEOXY_MAN = ContactMap(
    lectin="ConA", sugar="3-deoxy-Man", pdb_id="5CNA",
    n_HB=8,            # lose 2 HB (C3-OH → Asn14, Leu99)
    n_CH_pi=1,
    n_OH_eq_buried=2,  # C3-OH(eq) removed
    n_OH_ax_buried=1,
    n_NAc_buried=0,
    n_frozen_torsions=0,
    n_water_bridges=1,
    resolution_A=2.0,
    notes="C3-OH removed. Lose 2 HB + 1 eq desolv. Essential contact → large ΔΔG.",
)

CONA_4DEOXY_MAN = ContactMap(
    lectin="ConA", sugar="4-deoxy-Man", pdb_id="5CNA",
    n_HB=8,            # lose 2 HB (C4-OH → Asp208, Tyr100)
    n_CH_pi=1,
    n_OH_eq_buried=2,  # C4-OH(eq) removed
    n_OH_ax_buried=1,
    n_NAc_buried=0,
    n_frozen_torsions=0,
    n_water_bridges=0,  # lose 1 water bridge (C4-OH mediated)
    resolution_A=2.0,
    notes="C4-OH removed. Lose 2 HB + 1 eq desolv + 1 water. Largest ΔΔG expected.",
)

CONA_6DEOXY_MAN = ContactMap(
    lectin="ConA", sugar="6-deoxy-Man", pdb_id="5CNA",
    n_HB=8,            # lose 2 HB (C6-OH → Arg228, Asp208)
    n_CH_pi=1,
    n_OH_eq_buried=2,  # C6-OH removed
    n_OH_ax_buried=1,
    n_NAc_buried=0,
    n_frozen_torsions=0,
    n_water_bridges=1,
    resolution_A=2.0,
    notes="C6-OH removed (= rhamnose). Lose 2 HB + 1 eq desolv. Essential contact.",
)

PREDICTION_1_MAPS = {
    "MeαMan": CONA_MEAMAN,
    "2-deoxy-Man": CONA_2DEOXY_MAN,
    "3-deoxy-Man": CONA_3DEOXY_MAN,
    "4-deoxy-Man": CONA_4DEOXY_MAN,
    "6-deoxy-Man": CONA_6DEOXY_MAN,
}


# ═══════════════════════════════════════════════════════════════════════════
# PREDICTION 2: ConA epimer selectivity (Man vs Glc vs Gal)
# ═══════════════════════════════════════════════════════════════════════════

CONA_MEAGLC = ContactMap(
    lectin="ConA", sugar="MeαGlc", pdb_id="5CNA",
    n_HB=10,          # same H-bond network; C2-OH eq vs ax doesn't kill binding
    n_CH_pi=1,        # TYR12 still stacks
    n_OH_eq_buried=4, # C2(eq), C3(eq), C4(eq), C6(eq) — all equatorial in Glc
    n_OH_ax_buried=0, # no axial OH buried (C2 flipped to eq vs Man)
    n_NAc_buried=0,
    n_frozen_torsions=0,
    n_water_bridges=1,
    resolution_A=2.0,
    notes="Glucose: C2-OH eq (Man has ax). ~4× weaker Ka.",
)

CONA_MEAGAL = ContactMap(
    lectin="ConA", sugar="MeαGal", pdb_id="5CNA",
    n_HB=6,           # C4-OH axial causes steric clash, loses contacts
    n_CH_pi=0,        # steric disruption may displace sugar from stacking
    n_OH_eq_buried=2, # C3(eq), C6(eq) only
    n_OH_ax_buried=2, # C2(eq in Gal actually), C4(ax) — C4-ax causes clash
    n_NAc_buried=0,
    n_frozen_torsions=0,
    n_water_bridges=0,
    resolution_A=2.0,
    notes="Galactose: C4-OH axial disrupts pocket. No detectable binding.",
)

PREDICTION_2_MAPS = {
    "MeαMan": CONA_MEAMAN,
    "MeαGlc": CONA_MEAGLC,
    "MeαGal": CONA_MEAGAL,
}


# ═══════════════════════════════════════════════════════════════════════════
# PREDICTION 3: Lysozyme Trp62 mutant ΔΔG
# PDB: 1LZB (HEWL + chitobiose, 1.5 Å)
# ═══════════════════════════════════════════════════════════════════════════

LYSO_WT = ContactMap(
    lectin="Lysozyme_WT", sugar="(GlcNAc)2", pdb_id="1LZB",
    n_HB=6,           # NAG1: 5 HB, NAG2: 1 HB
    n_CH_pi=2,        # TRP108 (NAG1, 3.88Å/34°) + TRP62 (NAG2, 3.77Å/5.3°)
    n_OH_eq_buried=4, # 2 GlcNAc × 2 eq OH each (C3,C4,C6 eq; C2=NAc)
    n_OH_ax_buried=0,
    n_NAc_buried=2,   # 2 NAc groups
    n_frozen_torsions=1,  # one β1-4 glycosidic bond
    n_water_bridges=5,    # NAG1: 4, NAG2: 1
    resolution_A=1.5,
    notes="WT HEWL. Trp62 makes geometry-perfect CH-π with NAG2 C5.",
)

LYSO_W62Y = ContactMap(
    lectin="Lysozyme_W62Y", sugar="(GlcNAc)2", pdb_id="1LZG",
    n_HB=6,
    n_CH_pi=1,        # lose Trp62 CH-π; TRP108 remains; Tyr62 partial stacking
    n_OH_eq_buried=4,
    n_OH_ax_buried=0,
    n_NAc_buried=2,
    n_frozen_torsions=1,
    n_water_bridges=4,  # lose 1 water bridge at site B
    resolution_A=1.8,
    notes="W62Y: Tyr62 smaller ring, weaker stacking. Binds in 2 modes (A-B-C + B-C-D).",
)

LYSO_W62F = ContactMap(
    lectin="Lysozyme_W62F", sugar="(GlcNAc)2", pdb_id="1LZC",
    n_HB=5,
    n_CH_pi=1,        # Phe62 can stack but no H-bond capacity; TRP108 remains
    n_OH_eq_buried=4,
    n_OH_ax_buried=0,
    n_NAc_buried=2,
    n_frozen_torsions=1,
    n_water_bridges=3,
    resolution_A=1.8,
    notes="W62F: Phe cannot H-bond. Even weaker B-site.",
)

LYSO_W62H = ContactMap(
    lectin="Lysozyme_W62H", sugar="(GlcNAc)2", pdb_id="1LZB",
    n_HB=4,
    n_CH_pi=1,        # His62 imidazole: very different geometry; TRP108 remains
    n_OH_eq_buried=4,
    n_OH_ax_buried=0,
    n_NAc_buried=2,
    n_frozen_torsions=1,
    n_water_bridges=2,
    resolution_A=1.5,
    notes="W62H: drastic. A-E/B-F ratio shifts from 4:1 to 0.8:1.",
)

PREDICTION_3_MAPS = {
    "WT": LYSO_WT,
    "W62Y": LYSO_W62Y,
    "W62F": LYSO_W62F,
    "W62H": LYSO_W62H,
}


# ═══════════════════════════════════════════════════════════════════════════
# PREDICTION 4: ConA mono→di→tri→penta chain extension
# PDB: 1CVN (ConA + trimannoside, 2.3 Å)
# ═══════════════════════════════════════════════════════════════════════════

# Trimannoside consensus from 4 chains in 1CVN
CONA_TRIMANNOSE = ContactMap(
    lectin="ConA", sugar="trimannoside", pdb_id="1CVN",
    n_HB=19,          # sum across 3 mannose residues (consensus ~6+6+7)
    n_CH_pi=1,        # TYR12 contacts terminal mannose (MAN3)
    n_OH_eq_buried=8, # ~3 eq per mannose × 3 minus those not buried
    n_OH_ax_buried=2, # C2-OH(ax) on 2 mannose residues in binding site
    n_NAc_buried=0,
    n_frozen_torsions=2,  # 2 glycosidic bonds (α1-3 and α1-6)
    n_water_bridges=2,
    resolution_A=2.3,
    notes="Consensus 4 chains. α(1,6)Man in monosaccharide site, α(1,3)Man in extended site.",
)

PREDICTION_4_MAPS = {
    "MeαMan (mono)": CONA_MEAMAN,
    "trimannoside": CONA_TRIMANNOSE,
    # disaccharides and pentasaccharide: need additional structures or modeling
}


# ═══════════════════════════════════════════════════════════════════════════
# PREDICTION 6: Cross-lectin selectivity
# ═══════════════════════════════════════════════════════════════════════════

# WGA + GlcNAc: consensus from 2UVO
WGA_GLCNAC = ContactMap(
    lectin="WGA", sugar="GlcNAc", pdb_id="2UVO",
    n_HB=5,           # consensus across primary binding sites
    n_CH_pi=2,        # TYR66+TYR73 tandem pair (or TYR23+TYR30 in alt site)
    n_OH_eq_buried=3, # C3(eq), C4(eq), C6(eq)
    n_OH_ax_buried=0,
    n_NAc_buried=1,
    n_frozen_torsions=0,
    n_water_bridges=1,
    resolution_A=0.95,
    notes="Ultra-high resolution. Primary site with Tyr tandem CH-π.",
)

# PNA + Gal: consensus from 2PEL
PNA_GAL = ContactMap(
    lectin="PNA", sugar="galactose", pdb_id="2PEL",
    n_HB=9,           # consensus 4 chains (9,9,8,9)
    n_CH_pi=1,        # TYR125
    n_OH_eq_buried=3, # C2(eq), C3(eq), C6(eq)
    n_OH_ax_buried=1, # C4(ax) — Gal-specific recognition
    n_NAc_buried=0,
    n_frozen_torsions=0,
    n_water_bridges=2,
    resolution_A=2.0,
    notes="PNA strongly prefers Gal. Does NOT bind GalNAc (E129 steric block).",
)

# Galectin-3 + LacNAc: from 1A3K
GAL3_LACNAC = ContactMap(
    lectin="Galectin-3", sugar="LacNAc", pdb_id="1A3K",
    n_HB=12,          # NAG1: 5 HB, GAL2: 7 HB
    n_CH_pi=2,        # TRP181 (3.68Å, 5.7°) + HIS158 (3.98Å, 46.4°)
    n_OH_eq_buried=5,
    n_OH_ax_buried=1, # Gal C4-OH(ax)
    n_NAc_buried=1,
    n_frozen_torsions=1,  # one glycosidic bond
    n_water_bridges=4,    # 2+2 per residue
    resolution_A=2.1,
    notes="Canonical β-galactoside binding. TRP181 perfect stacking geometry.",
)

PREDICTION_6_MAPS = {
    "ConA:MeαMan": CONA_MEAMAN,
    "ConA:MeαGlc": CONA_MEAGLC,
    "WGA:GlcNAc": WGA_GLCNAC,
    "PNA:Gal": PNA_GAL,
    "Galectin-3:LacNAc": GAL3_LACNAC,
}


# ═══════════════════════════════════════════════════════════════════════════
# Utility: get all maps for a given prediction
# ═══════════════════════════════════════════════════════════════════════════

ALL_PREDICTIONS = {
    1: PREDICTION_1_MAPS,
    2: PREDICTION_2_MAPS,
    3: PREDICTION_3_MAPS,
    4: PREDICTION_4_MAPS,
    6: PREDICTION_6_MAPS,
}


def get_prediction_maps(prediction_number: int) -> dict:
    """Return the ContactMap dict for a given prediction number."""
    return ALL_PREDICTIONS.get(prediction_number, {})


if __name__ == "__main__":
    print("Prediction-ready ContactMaps")
    print("=" * 60)
    for pred_num, maps in sorted(ALL_PREDICTIONS.items()):
        print(f"\nPrediction {pred_num}: {len(maps)} lectin-sugar pairs")
        for label, cm in maps.items():
            print(f"  {label:25s}  HB={cm.n_HB:2d}  CH-π={cm.n_CH_pi}  "
                  f"OH_eq={cm.n_OH_eq_buried}  OH_ax={cm.n_OH_ax_buried}  "
                  f"NAc={cm.n_NAc_buried}  tors={cm.n_frozen_torsions}  "
                  f"water={cm.n_water_bridges}")