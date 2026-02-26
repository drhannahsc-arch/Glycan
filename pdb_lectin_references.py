"""
Curated lectin-glycan PDB structures for MABE glycan parameter extraction.

Used for:
  P9  — water bridge normalization (n_water / interface_area)
  P12 — multivalency inter-site distances
  P4/P5 validation — bound-state glycosidic torsion angles

Sources:
  - UniLectin3D (2709 structures, 1670 with glycan, as of 2025)
  - Manual curation from literature for highest-quality ITC-validated entries
  - GFDB (glycanstructure.org) for bound-state torsion statistics

SELECTION CRITERIA for this reference set:
  1. Resolution < 2.0 Å (for reliable water positions)
  2. Has ITC binding data in literature (for dG validation)
  3. Clear electron density for sugar and crystallographic waters
  4. Representative of major lectin families
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict


@dataclass
class LectinBindingSite:
    """One binding site in a lectin structure."""
    site_id: str                    # e.g., "A" for chain A site
    sugar: str                      # Bound sugar identity
    dG_exp_kJ: Optional[float] = None   # Experimental binding dG (ITC) in kJ/mol
    dH_exp_kJ: Optional[float] = None   # Experimental binding dH in kJ/mol
    Ka_M: Optional[float] = None        # Association constant
    n_crystal_waters: Optional[int] = None # Crystallographic waters in binding site
    aromatic_residues: List[str] = field(default_factory=list)  # CH-pi partners
    metal_ions: List[str] = field(default_factory=list)
    phi_bound: Optional[float] = None   # Glycosidic phi if disaccharide
    psi_bound: Optional[float] = None   # Glycosidic psi if disaccharide
    notes: str = ""


@dataclass
class LectinStructure:
    """A PDB lectin structure with curated binding information."""
    pdb_id: str
    name: str
    organism: str
    resolution_A: float
    lectin_class: str          # e.g., "legume", "C-type", "galectin", "siglec"
    fold: str                  # e.g., "beta-sandwich", "beta-trefoil"
    oligomeric_state: str      # "monomer", "dimer", "tetramer", "hexamer"
    n_binding_sites: int       # Per biological assembly
    inter_site_distance_A: Optional[float]  # For multivalent lectins
    sites: List[LectinBindingSite] = field(default_factory=list)
    doi_itc: str = ""          # Paper with ITC data
    doi_crystal: str = ""      # Paper with crystal structure
    notes: str = ""


# =====================================================================
# TIER 1 REFERENCE SET: High-resolution structures with ITC data
# These are the gold-standard entries for parameter extraction.
# =====================================================================

LECTIN_REFERENCE_SET = [

    # ─── CONCANAVALIN A (ConA) — the anchor system ──────────────
    LectinStructure(
        pdb_id="5CNA",
        name="Concanavalin A + trimannoside",
        organism="Canavalia ensiformis (jack bean)",
        resolution_A=1.8,
        lectin_class="legume",
        fold="beta-sandwich (jelly roll)",
        oligomeric_state="tetramer",
        n_binding_sites=4,
        inter_site_distance_A=65.0,  # ~65 Å between adjacent sites in tetramer
        sites=[
            LectinBindingSite(
                site_id="A",
                sugar="alpha-D-Man",
                dG_exp_kJ=-24.7,     # Ka ~ 8200 M-1 for MeAlphaMan
                dH_exp_kJ=-36.4,     # Chervenak & Toone 1995
                Ka_M=8200,
                n_crystal_waters=4,   # 3-4 conserved waters in binding site
                aromatic_residues=["Tyr12", "Tyr100"],
                metal_ions=["Ca2+", "Mn2+"],
                notes="Reference system for deoxy series. Chervenak & Toone 1995."
            ),
        ],
        doi_itc="10.1021/bi00009a001",   # Chervenak & Toone 1995
        doi_crystal="10.1016/0022-2836(76)90244-X",  # Hardman & Ainsworth 1976
        notes="ConA is THE validation anchor. Deoxy-mannose series gives P1/P2. "
              "Ca2+ and Mn2+ are structural, not directly involved in sugar binding."
    ),

    LectinStructure(
        pdb_id="1CVN",
        name="Concanavalin A + methyl alpha-D-mannopyranoside",
        organism="Canavalia ensiformis",
        resolution_A=1.6,
        lectin_class="legume",
        fold="beta-sandwich",
        oligomeric_state="tetramer",
        n_binding_sites=4,
        inter_site_distance_A=65.0,
        sites=[
            LectinBindingSite(
                site_id="A",
                sugar="MeAlphaMan",
                dG_exp_kJ=-22.3,
                Ka_M=8200,
                n_crystal_waters=3,
                aromatic_residues=["Tyr12", "Tyr100"],
                metal_ions=["Ca2+", "Mn2+"],
            ),
        ],
        doi_crystal="",
        notes="Higher resolution ConA-mannose for water position extraction."
    ),

    # ─── WHEAT GERM AGGLUTININ (WGA) — GlcNAc binder ──────────
    LectinStructure(
        pdb_id="2UVO",
        name="Wheat Germ Agglutinin + GlcNAc",
        organism="Triticum aestivum (wheat)",
        resolution_A=1.5,
        lectin_class="hevein domain",
        fold="hevein",
        oligomeric_state="dimer",
        n_binding_sites=8,           # 4 sites per monomer
        inter_site_distance_A=14.0,  # Adjacent sites ~14 Å apart
        sites=[
            LectinBindingSite(
                site_id="primary",
                sugar="GlcNAc",
                dG_exp_kJ=-18.0,     # ~Ka 1500 M-1
                Ka_M=1500,
                n_crystal_waters=2,
                aromatic_residues=["Trp", "Tyr"],
                notes="Multiple Trp stack on GlcNAc. NAc group contributes."
            ),
        ],
        doi_crystal="10.1107/S0907444905023838",
        notes="WGA tests P3 (NAc desolvation) and P6 (Trp CH-pi). "
              "Dense binding sites test multivalency (P11/P12)."
    ),

    # ─── PEANUT AGGLUTININ (PNA) — Gal binder ─────────────────
    LectinStructure(
        pdb_id="2PEL",
        name="Peanut Agglutinin + T-antigen (Gal-beta1,3-GalNAc)",
        organism="Arachis hypogaea (peanut)",
        resolution_A=1.95,
        lectin_class="legume",
        fold="beta-sandwich",
        oligomeric_state="tetramer",
        n_binding_sites=4,
        inter_site_distance_A=75.0,
        sites=[
            LectinBindingSite(
                site_id="A",
                sugar="Gal-beta1,3-GalNAc",
                dG_exp_kJ=-21.8,
                Ka_M=6700,
                n_crystal_waters=3,
                aromatic_residues=["Tyr125"],
                metal_ions=[],
                notes="Gal binding with Tyr CH-pi. No metal requirement."
            ),
        ],
        doi_crystal="10.1016/S0022-2836(96)90537-X",
        notes="Tests Gal selectivity over Glc. No metal — isolates sugar terms."
    ),

    # ─── GALECTIN-3 — beta-Gal binder, no metal ───────────────
    LectinStructure(
        pdb_id="3ZSJ",
        name="Galectin-3 CRD + LacNAc",
        organism="Homo sapiens",
        resolution_A=1.4,
        lectin_class="galectin",
        fold="beta-sandwich (S-type)",
        oligomeric_state="monomer (CRD)",
        n_binding_sites=1,
        inter_site_distance_A=None,
        sites=[
            LectinBindingSite(
                site_id="A",
                sugar="LacNAc (Gal-beta1,4-GlcNAc)",
                dG_exp_kJ=-26.0,     # ~Ka ~3e4 M-1
                Ka_M=30000,
                n_crystal_waters=2,
                aromatic_residues=["Trp181"],
                metal_ions=[],
                notes="Classic Trp-Gal CH-pi stacking. No metal. Clean system."
            ),
        ],
        doi_crystal="10.1074/jbc.M112.423459",
        notes="Galectin-3 is key for CH-pi validation. Trp181 stacks on Gal beta face. "
              "No metal ions — pure sugar recognition. High resolution for waters."
    ),

    # ─── CHOLERA TOXIN B — GM1 ganglioside binder ─────────────
    LectinStructure(
        pdb_id="3CHB",
        name="Cholera Toxin B pentamer + GM1 pentasaccharide",
        organism="Vibrio cholerae",
        resolution_A=1.25,
        lectin_class="AB5 toxin",
        fold="beta-pentamer (OB fold)",
        oligomeric_state="pentamer",
        n_binding_sites=5,
        inter_site_distance_A=35.0,  # Sites around pentamer ring
        sites=[
            LectinBindingSite(
                site_id="A",
                sugar="GM1 (Gal-GalNAc-[NeuAc]-Gal-Glc)",
                dG_exp_kJ=-40.6,     # Very tight for a lectin
                Ka_M=1.2e6,
                n_crystal_waters=5,
                aromatic_residues=["Trp88"],
                metal_ions=[],
                notes="High-affinity ganglioside binding. Multiple sugars contribute."
            ),
        ],
        doi_crystal="10.1126/science.1707152",
        notes="Tests multivalent display (5 sites). Extremely high resolution. "
              "Complex oligosaccharide ligand."
    ),

    # ─── DC-SIGN — C-type lectin, Ca2+-dependent ──────────────
    LectinStructure(
        pdb_id="1SL5",
        name="DC-SIGN CRD + Man4",
        organism="Homo sapiens",
        resolution_A=1.55,
        lectin_class="C-type",
        fold="C-type lectin fold",
        oligomeric_state="tetramer (neck + CRD)",
        n_binding_sites=4,
        inter_site_distance_A=40.0,
        sites=[
            LectinBindingSite(
                site_id="A",
                sugar="alpha-D-Man",
                dG_exp_kJ=-22.0,
                Ka_M=7000,
                n_crystal_waters=3,
                aromatic_residues=[],       # No aromatic in binding site!
                metal_ions=["Ca2+"],        # Ca2+ directly coordinates sugar OH
                notes="Ca2+ bridges 3-OH and 4-OH of mannose. Tests P10."
            ),
        ],
        doi_crystal="10.1074/jbc.M305636200",
        notes="C-type lectin: Ca2+ is essential for binding. No CH-pi. "
              "Tests metal coordination (P10) isolated from aromatic terms."
    ),

    # ─── RICIN B — beta-trefoil, Gal binder ──────────────────
    LectinStructure(
        pdb_id="2AAI",
        name="Ricin B chain + galactose",
        organism="Ricinus communis",
        resolution_A=1.8,
        lectin_class="R-type",
        fold="beta-trefoil",
        oligomeric_state="monomer",
        n_binding_sites=2,
        inter_site_distance_A=25.0,
        sites=[
            LectinBindingSite(
                site_id="1alpha",
                sugar="Gal",
                dG_exp_kJ=-19.0,
                Ka_M=2000,
                n_crystal_waters=2,
                aromatic_residues=["Trp37"],
                metal_ions=[],
            ),
            LectinBindingSite(
                site_id="2gamma",
                sugar="Gal",
                dG_exp_kJ=-17.0,
                Ka_M=1000,
                n_crystal_waters=2,
                aromatic_residues=[],
                metal_ions=[],
                notes="Weaker site, fewer contacts."
            ),
        ],
        doi_crystal="",
        notes="Two non-equivalent Gal sites: tests site-specific scoring."
    ),

    # ─── CYANOVIRIN-N — high-mannose binder ───────────────────
    LectinStructure(
        pdb_id="3GXZ",
        name="Cyanovirin-N + Man-alpha-1,2-Man",
        organism="Nostoc ellipsosporum",
        resolution_A=1.5,
        lectin_class="cyanovirin",
        fold="domain-swapped dimer",
        oligomeric_state="monomer",
        n_binding_sites=2,
        inter_site_distance_A=35.0,
        sites=[
            LectinBindingSite(
                site_id="A",
                sugar="Man-alpha-1,2-Man",
                dG_exp_kJ=-29.0,
                Ka_M=1e5,
                n_crystal_waters=4,
                aromatic_residues=[],
                metal_ions=[],
                phi_bound=-40.0,    # Alpha(1->2) linkage
                psi_bound=170.0,
                notes="Tests conformational entropy (P4/P5) for alpha(1->2) linkage."
            ),
        ],
        doi_crystal="10.1021/ja805696u",
        notes="High affinity for dimannose. Good for P4 validation via bound torsions."
    ),
]


# =====================================================================
# MULTIVALENCY GEOMETRY (P12): Inter-site distances for common lectins
# =====================================================================

LECTIN_MULTIVALENCY = {
    # lectin: {oligomeric_state, n_sites, inter_site_distances_A}
    "ConA": {
        "state": "tetramer",
        "n_sites": 4,
        "distances_A": [65.0, 65.0, 40.0],  # adjacent, adjacent, diagonal
        "geometry": "rectangular tetramer",
        "notes": "Sites on corners of roughly rectangular arrangement.",
    },
    "WGA": {
        "state": "dimer",
        "n_sites": 8,   # 4 per monomer
        "distances_A": [14.0, 14.0, 14.0, 52.0],  # adjacent (within monomer), cross-dimer
        "geometry": "linear array (4 per chain, paired)",
        "notes": "Very dense sites. Ideal for chitooligosaccharide bridging.",
    },
    "Cholera_Toxin_B": {
        "state": "pentamer",
        "n_sites": 5,
        "distances_A": [35.0],  # between adjacent sites around ring
        "geometry": "pentagonal ring",
        "notes": "All sites equivalent. 5-fold symmetric.",
    },
    "DC-SIGN": {
        "state": "tetramer",
        "n_sites": 4,
        "distances_A": [40.0],
        "geometry": "extended neck, flexible CRDs",
        "notes": "CRDs at end of coiled-coil necks. Flexible inter-site distance.",
    },
    "Galectin-1": {
        "state": "dimer",
        "n_sites": 2,
        "distances_A": [50.0],
        "geometry": "back-to-back dimer",
        "notes": "Sites point in opposite directions. Cross-links cell surfaces.",
    },
    "Galectin-3": {
        "state": "monomer (pentamerizes at high conc)",
        "n_sites": 1,
        "distances_A": [],
        "geometry": "monomer with N-terminal oligomerization",
        "notes": "CRD is monomeric. Multivalency via N-terminal collagen-like domain.",
    },
    "PNA": {
        "state": "tetramer",
        "n_sites": 4,
        "distances_A": [75.0],
        "geometry": "tetrahedral-like",
        "notes": "Wide spacing between sites.",
    },
    "Ricin_B": {
        "state": "monomer",
        "n_sites": 2,
        "distances_A": [25.0],
        "geometry": "beta-trefoil lobes",
        "notes": "Two non-equivalent sites within single domain.",
    },
}


# =====================================================================
# WATER BRIDGE STATISTICS from PDB lectin structures
# Used for P9 normalization: n_water_bridge / interface_area
# =====================================================================

# Compiled from Clarke 2001, Kadirvelraj 2008, and PDB survey.
# These are literature-validated counts of conserved crystallographic
# waters that bridge protein and sugar via hydrogen bonds.

WATER_BRIDGE_STATISTICS = {
    # pdb_id: {n_waters_bridging, interface_area_A2, n_waters_per_A2}
    "5CNA": {
        "n_waters": 4,
        "interface_area_A2": 350,
        "density": 0.0114,
        "sugar": "trimannoside",
        "notes": "3-4 highly conserved waters. One bridges Ca2+ to sugar.",
    },
    "1CVN": {
        "n_waters": 3,
        "interface_area_A2": 280,
        "density": 0.0107,
        "sugar": "MeAlphaMan",
    },
    "3ZSJ": {
        "n_waters": 2,
        "interface_area_A2": 320,
        "density": 0.0063,
        "sugar": "LacNAc",
        "notes": "Fewer waters than ConA. Galectin uses extensive direct H-bonds.",
    },
    "3CHB": {
        "n_waters": 5,
        "interface_area_A2": 520,
        "density": 0.0096,
        "sugar": "GM1 pentasaccharide",
        "notes": "Larger interface, more waters. High resolution (1.25 Å).",
    },
    "1SL5": {
        "n_waters": 3,
        "interface_area_A2": 290,
        "density": 0.0103,
        "sugar": "Man4",
        "notes": "Ca2+-mediated binding. Some waters coordinate Ca2+.",
    },
    "2UVO": {
        "n_waters": 2,
        "interface_area_A2": 240,
        "density": 0.0083,
        "sugar": "GlcNAc",
        "notes": "Small binding site, compact. Dominated by Trp stacking.",
    },
}


def get_water_bridge_norm() -> dict:
    """Compute P9 normalization value from PDB statistics.

    P9 = mean(n_waters / interface_area) across reference structures.

    Returns:
        Dict with mean density, range, and per-structure values.
    """
    densities = [v["density"] for v in WATER_BRIDGE_STATISTICS.values()]
    mean_d = sum(densities) / len(densities)

    return {
        "P9_mean_waters_per_A2": round(mean_d, 4),
        "P9_range": (round(min(densities), 4), round(max(densities), 4)),
        "n_structures": len(densities),
        "notes": "From curated PDB lectin structures with resolution < 2.0 Å. "
                 "Use to normalize water bridge count by interface size.",
    }


def get_inter_site_distances() -> Dict[str, dict]:
    """Return multivalency geometry data for P12 calibration."""
    return LECTIN_MULTIVALENCY


# =====================================================================
# BOUND-STATE TORSION ANGLES for P4/P5 validation
# From crystal structures and GFDB statistics
# =====================================================================

# These are the most common bound-state glycosidic angles observed
# in lectin complexes. Compare against CHI free-state minima to
# compute TdS_freeze.

BOUND_STATE_TORSIONS = {
    # (linkage_type, sugar): [(phi, psi, pdb_id, lectin), ...]
    ("alpha-1,3", "Man"): [
        (-85, 155, "5CNA", "ConA"),
        (-75, 160, "1CVN", "ConA"),
    ],
    ("alpha-1,2", "Man"): [
        (-40, 170, "3GXZ", "Cyanovirin-N"),
    ],
    ("beta-1,4", "GlcNAc"): [
        (-90, 140, "2UVO", "WGA"),
    ],
    ("beta-1,3", "Gal"): [
        (-70, 130, "2PEL", "PNA"),
    ],
    ("beta-1,4", "Gal"): [
        (-65, 125, "3ZSJ", "Galectin-3"),
    ],
}


if __name__ == "__main__":
    print("PDB Lectin Reference Dataset")
    print("=" * 60)

    print(f"\n{len(LECTIN_REFERENCE_SET)} reference structures:")
    for ls in LECTIN_REFERENCE_SET:
        sites_str = ", ".join(s.sugar for s in ls.sites)
        print(f"  {ls.pdb_id} | {ls.name[:40]:40s} | {ls.resolution_A:.1f} Å | {sites_str}")

    print(f"\nWater bridge normalization (P9):")
    p9 = get_water_bridge_norm()
    print(f"  Mean: {p9['P9_mean_waters_per_A2']:.4f} waters/Å²")
    print(f"  Range: {p9['P9_range']}")
    print(f"  From {p9['n_structures']} structures")

    print(f"\nMultivalency geometry (P12):")
    for name, data in LECTIN_MULTIVALENCY.items():
        d_str = ", ".join(f"{d:.0f}" for d in data["distances_A"]) if data["distances_A"] else "N/A"
        print(f"  {name:20s}: {data['state']:12s} {data['n_sites']} sites, d={d_str} Å")