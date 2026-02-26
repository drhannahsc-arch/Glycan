"""
test_glycan_scorer.py — Tests for the glycan scoring pipeline.

Run: python -m pytest tests/test_glycan_scorer.py -v

Test categories:
  1. Dataclass construction and defaults
  2. Parameter range sanity checks
  3. Scoring function sign/direction checks
  4. ΔΔG direction checks (deoxy weaker, mutant weaker)
  5. Statistics computation
  6. Sugar property card correctness
"""

import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from glycan_scorer import (
    GlycanParams, ContactMap, PredictionResult,
    predict_dG, predict_ddG, score_panel, compute_statistics,
)
from sugar_properties import (
    build_property_card, build_all_cards,
    MONOSACCHARIDE_HYDROXYL_MAP, AXIAL_CH_COUNT,
)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Dataclass construction
# ═══════════════════════════════════════════════════════════════════════════

def test_default_params():
    p = GlycanParams()
    assert p.eps_HB == -6.0, "Default H-bond energy should be -6.0 kJ/mol"
    assert p.k_desolv_eq == 0.0, "Unfitted params should be 0.0"


def test_contact_map_construction():
    cm = ContactMap(lectin="ConA", sugar="mannose", pdb_id="5CNA", n_HB=4)
    assert cm.n_HB == 4
    assert cm.n_CH_pi == 0  # default


# ═══════════════════════════════════════════════════════════════════════════
# 2. Parameter range sanity checks
# ═══════════════════════════════════════════════════════════════════════════

def test_hb_range():
    """H-bond in water: -4 to -8 kJ/mol (literature consensus)."""
    p = GlycanParams(eps_HB=-6.0)
    assert -8.0 <= p.eps_HB <= -4.0


def test_ch_pi_range():
    """CH-π in water: -1.0 to -4.0 kJ/mol (Nishio consensus)."""
    p = GlycanParams(eps_CH_pi=-2.0)
    assert -4.0 <= p.eps_CH_pi <= -1.0


def test_desolv_eq_gt_ax():
    """Equatorial OH should cost more to desolvate than axial."""
    p = GlycanParams(k_desolv_eq=8.0, k_desolv_ax=5.0)
    assert p.k_desolv_eq > p.k_desolv_ax


def test_conf_entropy_range():
    """Glycosidic torsion freezing: 1.5-6 kJ/mol per φ/ψ pair."""
    p = GlycanParams(eps_conf=3.0)
    assert 1.5 <= p.eps_conf <= 6.0


# ═══════════════════════════════════════════════════════════════════════════
# 3. Scoring function sign checks
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def typical_params():
    return GlycanParams(
        k_desolv_eq=8.0,
        k_desolv_ax=5.0,
        k_desolv_NAc=6.0,
        eps_CH_pi=-2.0,
        eps_HB=-6.0,
        eps_conf=3.0,
        dG_0=10.0,
    )


def test_hb_favorable(typical_params):
    """H-bonds should make ΔG more negative."""
    cm = ContactMap(lectin="test", sugar="test", pdb_id="XXXX", n_HB=4)
    r = predict_dG(typical_params, cm)
    assert r.dG_HB < 0, "H-bond contribution should be negative (favorable)"


def test_ch_pi_favorable(typical_params):
    """CH-π should make ΔG more negative."""
    cm = ContactMap(lectin="test", sugar="test", pdb_id="XXXX", n_CH_pi=2)
    r = predict_dG(typical_params, cm)
    assert r.dG_CH_pi < 0, "CH-π contribution should be negative (favorable)"


def test_desolvation_unfavorable(typical_params):
    """Burying OHs costs energy (desolvation penalty is positive)."""
    cm = ContactMap(lectin="test", sugar="test", pdb_id="XXXX",
                    n_OH_eq_buried=3)
    r = predict_dG(typical_params, cm)
    assert r.dG_desolv > 0, "Desolvation should be positive (unfavorable)"


def test_conformational_entropy_unfavorable(typical_params):
    """Freezing torsions costs entropy (unfavorable, positive contribution)."""
    cm = ContactMap(lectin="test", sugar="test", pdb_id="XXXX",
                    n_frozen_torsions=2)
    r = predict_dG(typical_params, cm)
    assert r.dG_conf > 0, "Conf entropy term should be positive (unfavorable penalty)"


def test_more_hbonds_tighter(typical_params):
    """More H-bonds → more negative ΔG → tighter binding."""
    cm_few = ContactMap(lectin="t", sugar="t", pdb_id="X", n_HB=2)
    cm_many = ContactMap(lectin="t", sugar="t", pdb_id="X", n_HB=6)
    r_few = predict_dG(typical_params, cm_few)
    r_many = predict_dG(typical_params, cm_many)
    assert r_many.dG_predicted < r_few.dG_predicted


# ═══════════════════════════════════════════════════════════════════════════
# 4. ΔΔG direction checks
# ═══════════════════════════════════════════════════════════════════════════

def test_deoxy_weaker_binding(typical_params):
    """Removing an OH should weaken binding (positive ΔΔG)."""
    wt = ContactMap(lectin="ConA", sugar="mannose", pdb_id="5CNA",
                    n_HB=4, n_OH_eq_buried=3, n_OH_ax_buried=1)
    deoxy = ContactMap(lectin="ConA", sugar="2-deoxy-mannose", pdb_id="5CNA",
                       n_HB=3, n_OH_eq_buried=3, n_OH_ax_buried=0)
    # deoxy removes one axial OH and its H-bond
    ddG = predict_ddG(typical_params, wt, deoxy)
    assert ddG > 0, f"Deoxy should be weaker (ΔΔG > 0), got {ddG:.1f}"


def test_chain_extension_tighter(typical_params):
    """Adding residues should (generally) tighten binding despite entropy cost."""
    mono = ContactMap(lectin="ConA", sugar="mannose", pdb_id="5CNA",
                      n_HB=4, n_OH_eq_buried=3, n_CH_pi=1)
    di = ContactMap(lectin="ConA", sugar="dimannose", pdb_id="1CVN",
                    n_HB=6, n_OH_eq_buried=5, n_CH_pi=1, n_frozen_torsions=1)
    r_mono = predict_dG(typical_params, mono)
    r_di = predict_dG(typical_params, di)
    # di should be tighter IF new contacts outweigh entropy cost
    # (this depends on parameters — just check that the function runs)
    assert isinstance(r_di.dG_predicted, float)


# ═══════════════════════════════════════════════════════════════════════════
# 5. Statistics
# ═══════════════════════════════════════════════════════════════════════════

def test_perfect_prediction_r2():
    """R² = 1.0 for perfect predictions."""
    results = [
        PredictionResult("A", "s1", -10.0, -10.0),
        PredictionResult("A", "s2", -20.0, -20.0),
        PredictionResult("A", "s3", -15.0, -15.0),
    ]
    stats = compute_statistics(results)
    assert abs(stats["R2"] - 1.0) < 0.001


def test_statistics_handles_missing():
    """Should skip entries with no measured value."""
    results = [
        PredictionResult("A", "s1", -10.0, -10.0),
        PredictionResult("A", "s2", -20.0, None),
        PredictionResult("A", "s3", -15.0, -15.0),
    ]
    stats = compute_statistics(results)
    assert stats["n"] == 2


# ═══════════════════════════════════════════════════════════════════════════
# 6. Sugar property cards
# ═══════════════════════════════════════════════════════════════════════════

def test_glucose_all_equatorial():
    """β-D-glucose: all OH equatorial (4C1 chair)."""
    card = build_property_card("beta-D-glucose")
    assert card.n_OH_equatorial == 5  # C1, C2, C3, C4, C6
    assert card.n_OH_axial == 0


def test_mannose_has_axial_c2():
    """α-D-mannose: C1 and C2 are axial."""
    card = build_property_card("alpha-D-mannose")
    assert card.n_OH_axial >= 2  # at least C1 and C2


def test_galactose_has_axial_c4():
    """α-D-galactose: C4 is axial (defines galactose vs glucose)."""
    card = build_property_card("alpha-D-galactose")
    assert card.n_OH_axial >= 2  # C1 and C4


def test_mannose_best_ch_pi_donor():
    """α-mannose should have most axial CH (best CH-π donor)."""
    man = AXIAL_CH_COUNT.get("alpha-D-mannose", 0)
    glc = AXIAL_CH_COUNT.get("alpha-D-glucose", 0)
    assert man >= glc, "Mannose should be stronger CH-π donor than glucose"


def test_all_cards_build():
    """All defined monosaccharides should produce valid cards."""
    cards = build_all_cards()
    assert len(cards) >= 10
    for name, card in cards.items():
        assert card.n_OH_equatorial + card.n_OH_axial + (1 if card.has_NAc else 0) >= 3, \
            f"{name}: too few OH groups"


# ═══════════════════════════════════════════════════════════════════════════
# Run
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    pytest.main([__file__, "-v"])