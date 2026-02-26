"""
sugar_properties.py — Phase G0: Sugar property extraction from pure chemistry.

Builds a SugarPropertyCard for every monosaccharide and common disaccharide.
No binding data. No receptor information. Pure molecular descriptors.

Sources:
  - Ring geometry: Cremer-Pople puckering (Cremer & Pople 1975, JACS 97:1354)
  - Hydroxyl inventory: standard nomenclature (equatorial/axial assignment)
  - SASA: RDKit or FreeSASA calculation
  - Dissolution enthalpy: Schwarz 1996, Jasra 1982
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
from glycan_scorer import SugarPropertyCard


# ═══════════════════════════════════════════════════════════════════════════
# Monosaccharide definitions — hydroxyl stereochemistry
# ═══════════════════════════════════════════════════════════════════════════
#
# Convention: 4C1 chair conformation (standard for D-sugars).
# Positions C1–C6. OH orientation: eq = equatorial, ax = axial.
#
# Reference: Any carbohydrate chemistry textbook, e.g. Varki et al.
# "Essentials of Glycobiology" (free at NCBI Bookshelf).

MONOSACCHARIDE_HYDROXYL_MAP = {
    # name: {position: "eq"|"ax"} for each OH in 4C1
    # α-anomers shown; β shifts C1-OH from ax→eq
    "alpha-D-mannose": {
        "C1": "ax", "C2": "ax", "C3": "eq", "C4": "eq", "C6": "eq",
    },
    "beta-D-mannose": {
        "C1": "eq", "C2": "ax", "C3": "eq", "C4": "eq", "C6": "eq",
    },
    "alpha-D-glucose": {
        "C1": "ax", "C2": "eq", "C3": "eq", "C4": "eq", "C6": "eq",
    },
    "beta-D-glucose": {
        "C1": "eq", "C2": "eq", "C3": "eq", "C4": "eq", "C6": "eq",
    },
    "alpha-D-galactose": {
        "C1": "ax", "C2": "eq", "C3": "eq", "C4": "ax", "C6": "eq",
    },
    "beta-D-galactose": {
        "C1": "eq", "C2": "eq", "C3": "eq", "C4": "ax", "C6": "eq",
    },
    "alpha-D-GlcNAc": {
        "C1": "ax", "C2": "NAc", "C3": "eq", "C4": "eq", "C6": "eq",
    },
    "beta-D-GlcNAc": {
        "C1": "eq", "C2": "NAc", "C3": "eq", "C4": "eq", "C6": "eq",
    },
    "alpha-D-GalNAc": {
        "C1": "ax", "C2": "NAc", "C3": "eq", "C4": "ax", "C6": "eq",
    },
    "beta-D-GalNAc": {
        "C1": "eq", "C2": "NAc", "C3": "eq", "C4": "ax", "C6": "eq",
    },
    "alpha-L-fucose": {
        "C1": "ax", "C2": "eq", "C3": "eq", "C4": "eq",
        # C6 is CH3 (6-deoxy), no OH
    },
    "alpha-D-xylose": {
        "C1": "ax", "C2": "eq", "C3": "eq", "C4": "eq",
        # No C6 (pentose)
    },
    # Deoxy series for Prediction 1 (ConA deoxy-mannose)
    "2-deoxy-D-mannose": {
        "C1": "ax", "C3": "eq", "C4": "eq", "C6": "eq",
        # C2-OH removed
    },
    "3-deoxy-D-mannose": {
        "C1": "ax", "C2": "ax", "C4": "eq", "C6": "eq",
    },
    "4-deoxy-D-mannose": {
        "C1": "ax", "C2": "ax", "C3": "eq", "C6": "eq",
    },
    "6-deoxy-D-mannose": {  # = rhamnose
        "C1": "ax", "C2": "ax", "C3": "eq", "C4": "eq",
    },
}

# CH-π donor inventory: count of axial C-H bonds on the α-face.
# These are the primary donors for aromatic stacking (Asensio 2013).
AXIAL_CH_COUNT = {
    "alpha-D-mannose":    3,  # C1-H(ax), C2-H(ax), C5-H — strong α-face
    "beta-D-mannose":     2,  # C2-H(ax), C5-H
    "alpha-D-glucose":    2,  # C1-H(ax), C5-H
    "beta-D-glucose":     1,  # C5-H only — weakest CH-π donor
    "alpha-D-galactose":  3,  # C1-H(ax), C4-H(ax), C5-H
    "beta-D-galactose":   2,  # C4-H(ax), C5-H
    "alpha-D-GlcNAc":    2,
    "beta-D-GlcNAc":     1,
    "alpha-D-GalNAc":    3,
    "beta-D-GalNAc":     2,
    "alpha-L-fucose":     3,  # strong CH-π donor (selectin recognition)
    "alpha-D-xylose":     2,
}


def build_property_card(name: str) -> SugarPropertyCard:
    """Build a SugarPropertyCard from the hydroxyl map.

    SASA and Cremer-Pople values require RDKit/FreeSASA and are
    left as 0.0 until Phase G0 implementation.
    """
    oh_map = MONOSACCHARIDE_HYDROXYL_MAP.get(name, {})

    n_eq = sum(1 for v in oh_map.values() if v == "eq")
    n_ax = sum(1 for v in oh_map.values() if v == "ax")
    has_nac = any(v == "NAc" for v in oh_map.values())

    # Abbreviation extraction (rough)
    abbrev_map = {
        "mannose": "Man", "glucose": "Glc", "galactose": "Gal",
        "GlcNAc": "GlcNAc", "GalNAc": "GalNAc",
        "fucose": "Fuc", "xylose": "Xyl",
    }
    abbrev = "?"
    for key, val in abbrev_map.items():
        if key in name:
            abbrev = val
            break

    return SugarPropertyCard(
        name=name,
        abbreviation=abbrev,
        n_OH_equatorial=n_eq,
        n_OH_axial=n_ax,
        has_NAc=has_nac,
        n_axial_CH=AXIAL_CH_COUNT.get(name, 0),
    )


def build_all_cards() -> Dict[str, SugarPropertyCard]:
    """Build property cards for all defined monosaccharides."""
    return {name: build_property_card(name) for name in MONOSACCHARIDE_HYDROXYL_MAP}


# ═══════════════════════════════════════════════════════════════════════════
# Display
# ═══════════════════════════════════════════════════════════════════════════

def print_card_table(cards: Dict[str, SugarPropertyCard]) -> None:
    """Print a summary table of all sugar property cards."""
    print(f"{'Sugar':<28} {'Abbr':<8} {'OH_eq':>5} {'OH_ax':>5} {'NAc':>4} {'ax_CH':>5}")
    print("-" * 60)
    for name, card in sorted(cards.items()):
        print(f"{card.name:<28} {card.abbreviation:<8} "
              f"{card.n_OH_equatorial:>5} {card.n_OH_axial:>5} "
              f"{'Y' if card.has_NAc else 'N':>4} {card.n_axial_CH:>5}")


if __name__ == "__main__":
    cards = build_all_cards()
    print(f"Sugar Property Cards — Phase G0")
    print(f"Total: {len(cards)} monosaccharides\n")
    print_card_table(cards)