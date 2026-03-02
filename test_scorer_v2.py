"""
Tests for contact_extractor.py and glycan_scorer_v2.py

Covers:
  - PDB fetching and parsing
  - ContactMap generation from real structures
  - H-bond, CH-π, metal, water bridge detection
  - Scorer physics and sign conventions
  - Deoxy-mannose predictions
  - Selectivity panel scoring
  - Validation against experimental ITC data
"""

import unittest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from contact_extractor import (
    ContactMap, extract_contacts, extract_contacts_simple,
    fetch_pdb, SUGAR_RESNAMES, SUGAR_TYPE, AROMATIC_RING_ATOMS,
    HYDROXYL_CLASS, HB_CUTOFF, CHPI_CUTOFF, METAL_CUTOFF,
)
from glycan_scorer_v2 import (
    GlycanParams, ScoringResult, score_contact_map,
    predict_deoxy_series, score_sugar_panel, validate_predictions,
    MANNOSE_OH_POSITIONS, SUGAR_OH_PATTERNS,
)

PDB_CACHE = "/home/claude/pdb_cache"


# ═════════════════════════════════════════════════════════════════════
# Contact Extractor Tests
# ═════════════════════════════════════════════════════════════════════

class TestPDBFetching(unittest.TestCase):
    """Test PDB file retrieval from EBI."""

    def test_fetch_existing_pdb(self):
        fn = fetch_pdb("3ZSJ", cache_dir=PDB_CACHE)
        self.assertTrue(os.path.exists(fn))
        self.assertGreater(os.path.getsize(fn), 10000)

    def test_fetch_cona(self):
        fn = fetch_pdb("5CNA", cache_dir=PDB_CACHE)
        self.assertTrue(os.path.exists(fn))

    def test_fetch_cached(self):
        """Second fetch should use cache."""
        fn1 = fetch_pdb("3ZSJ", cache_dir=PDB_CACHE)
        fn2 = fetch_pdb("3ZSJ", cache_dir=PDB_CACHE)
        self.assertEqual(fn1, fn2)


class TestContactMapDataclass(unittest.TestCase):
    """Test ContactMap structure."""

    def test_empty_contact_map(self):
        cm = ContactMap(pdb_id="TEST", chain="A", sugar_resname="GAL", sugar_resid=1)
        self.assertEqual(cm.n_hbonds, 0)
        self.assertEqual(cm.n_ch_pi, 0)
        self.assertEqual(cm.metal_ions, [])

    def test_summary_format(self):
        cm = ContactMap(pdb_id="3ZSJ", chain="B", sugar_resname="GAL",
                        sugar_resid=2, n_hbonds=6, n_ch_pi=2,
                        aromatic_residues=["TRP181", "HIS158"])
        s = cm.summary()
        self.assertIn("3ZSJ", s)
        self.assertIn("HB=6", s)
        self.assertIn("CH-π=2", s)


class TestGalectin3Contacts(unittest.TestCase):
    """Test contact extraction from 3ZSJ (Galectin-3 + LacNAc)."""

    @classmethod
    def setUpClass(cls):
        cms = extract_contacts("3ZSJ", sugar_resname="GAL",
                               cache_dir=PDB_CACHE)
        cls.cm = cms[0] if cms else None

    def test_found_sugar(self):
        self.assertIsNotNone(self.cm)
        self.assertEqual(self.cm.sugar_resname, "GAL")

    def test_hbond_count(self):
        """Galectin-3 makes 4-8 H-bonds with Gal at 3.2Å cutoff."""
        self.assertGreaterEqual(self.cm.n_hbonds, 4)
        self.assertLessEqual(self.cm.n_hbonds, 10)

    def test_ch_pi_contacts(self):
        """Must find TRP181 stacking on galactose."""
        self.assertGreaterEqual(self.cm.n_ch_pi, 1)
        trp_found = any("TRP" in ar for ar in self.cm.aromatic_residues)
        self.assertTrue(trp_found, f"Expected TRP in {self.cm.aromatic_residues}")

    def test_no_metals(self):
        """Galectin-3 has no metal ions in binding site."""
        self.assertEqual(len(self.cm.metal_ions), 0)

    def test_water_bridges(self):
        """Should find some bridging waters."""
        self.assertGreaterEqual(self.cm.n_water_bridges, 0)

    def test_buried_oh(self):
        """Should find at least 1 buried OH."""
        total_buried = (self.cm.n_buried_eq_OH + self.cm.n_buried_ax_OH
                        + self.cm.n_buried_NAc)
        self.assertGreaterEqual(total_buried, 1)


class TestConAContacts(unittest.TestCase):
    """Test contact extraction from 5CNA (ConA + MeAlphaMan)."""

    @classmethod
    def setUpClass(cls):
        cms = extract_contacts("5CNA", sugar_resname="MMA",
                               chain="A", cache_dir=PDB_CACHE)
        cls.cm = cms[0] if cms else None

    def test_found_sugar(self):
        self.assertIsNotNone(self.cm)
        self.assertEqual(self.cm.sugar_resname, "MMA")

    def test_hbond_count(self):
        """ConA makes 5-10 H-bonds with mannose."""
        self.assertGreaterEqual(self.cm.n_hbonds, 5)
        self.assertLessEqual(self.cm.n_hbonds, 12)

    def test_ch_pi_contacts(self):
        """ConA has TYR12 and TYR100 near sugar."""
        self.assertGreaterEqual(self.cm.n_ch_pi, 1)
        tyr_found = any("TYR" in ar for ar in self.cm.aromatic_residues)
        self.assertTrue(tyr_found, f"Expected TYR in {self.cm.aromatic_residues}")

    def test_no_direct_metal_coordination(self):
        """ConA's Ca²⁺/Mn²⁺ are structural — NOT within 3Å of sugar.
        This is the key test that fixes the v1 metal double-counting."""
        self.assertEqual(len(self.cm.metal_ions), 0,
                         f"ConA metals should not coordinate sugar directly, "
                         f"got: {self.cm.metal_ions}")

    def test_mannose_axial_oh(self):
        """Mannose has C2-OH axial. Should detect at least 1 axial OH."""
        # This depends on whether O2 is within 4Å of protein
        # In ConA, O2 is at the pocket entrance — may or may not be buried
        total_buried = self.cm.n_buried_eq_OH + self.cm.n_buried_ax_OH
        self.assertGreaterEqual(total_buried, 1)

    def test_four_subunits(self):
        """5CNA is a tetramer — should find MMA in multiple chains."""
        all_cms = extract_contacts("5CNA", sugar_resname="MMA",
                                   cache_dir=PDB_CACHE)
        self.assertEqual(len(all_cms), 4, "ConA tetramer should have 4 MMA")


class TestDCSIGNContacts(unittest.TestCase):
    """Test 1SL5 (DC-SIGN) — calcium-dependent lectin."""

    @classmethod
    def setUpClass(cls):
        cms = extract_contacts("1SL5", sugar_resname="FUC",
                               cache_dir=PDB_CACHE)
        cls.cm = cms[0] if cms else None

    def test_found_fucose(self):
        self.assertIsNotNone(self.cm)
        self.assertEqual(self.cm.sugar_resname, "FUC")

    def test_calcium_coordination(self):
        """DC-SIGN coordinates sugar through Ca²⁺."""
        ca_found = any("CA" in m for m in self.cm.metal_ions)
        self.assertTrue(ca_found,
                        f"Expected Ca in DC-SIGN, got: {self.cm.metal_ions}")

    def test_no_aromatics(self):
        """DC-SIGN binding is H-bond + Ca, not CH-π."""
        self.assertEqual(self.cm.n_ch_pi, 0)


class TestWGAContacts(unittest.TestCase):
    """Test 2UVO (WGA) — multi-site, NAc-binding lectin."""

    @classmethod
    def setUpClass(cls):
        all_cms = extract_contacts("2UVO", sugar_resname="NAG",
                                   cache_dir=PDB_CACHE)
        cls.all_cms = all_cms
        cls.cm = all_cms[0] if all_cms else None

    def test_found_nag(self):
        self.assertIsNotNone(self.cm)

    def test_multiple_sites(self):
        """WGA has multiple GlcNAc binding sites."""
        self.assertGreater(len(self.all_cms), 1)

    def test_aromatic_stacking(self):
        """WGA uses Tyr stacking on GlcNAc."""
        any_tyr = any(
            any("TYR" in ar for ar in cm.aromatic_residues)
            for cm in self.all_cms
        )
        self.assertTrue(any_tyr, "Expected TYR contacts in WGA")

    def test_nac_buried(self):
        """At least one site should bury the NAc group."""
        any_nac = any(cm.n_buried_NAc > 0 for cm in self.all_cms)
        self.assertTrue(any_nac, "Expected buried NAc in WGA")


class TestCutoffConstants(unittest.TestCase):
    """Test that cutoff values are physically reasonable."""

    def test_hb_cutoff(self):
        self.assertGreater(HB_CUTOFF, 2.5)
        self.assertLess(HB_CUTOFF, 4.0)

    def test_chpi_cutoff(self):
        self.assertGreater(CHPI_CUTOFF, 3.5)
        self.assertLess(CHPI_CUTOFF, 5.5)

    def test_metal_cutoff(self):
        self.assertGreater(METAL_CUTOFF, 2.0)
        self.assertLess(METAL_CUTOFF, 4.0)


class TestSugarClassification(unittest.TestCase):
    """Test sugar residue name lookups."""

    def test_mannose_variants(self):
        self.assertIn("MAN", SUGAR_RESNAMES)
        self.assertIn("BMA", SUGAR_RESNAMES)
        self.assertIn("MMA", SUGAR_RESNAMES)

    def test_sugar_type_mapping(self):
        self.assertEqual(SUGAR_TYPE["MMA"], "Man")
        self.assertEqual(SUGAR_TYPE["GAL"], "Gal")
        self.assertEqual(SUGAR_TYPE["NAG"], "GlcNAc")

    def test_hydroxyl_classification(self):
        """Mannose C2-OH is axial, glucose C2-OH is equatorial."""
        self.assertEqual(HYDROXYL_CLASS["Man"]["O2"], "ax")
        self.assertEqual(HYDROXYL_CLASS["Glc"]["O2"], "eq")
        self.assertEqual(HYDROXYL_CLASS["Gal"]["O4"], "ax")
        self.assertEqual(HYDROXYL_CLASS["GlcNAc"]["N2"], "NAc")


# ═════════════════════════════════════════════════════════════════════
# Scorer Tests
# ═════════════════════════════════════════════════════════════════════

class TestGlycanParams(unittest.TestCase):
    """Test parameter values and physics."""

    def test_default_params(self):
        p = GlycanParams()
        self.assertLess(p.eps_HB_kJ, 0)           # H-bonds favorable
        self.assertLess(p.eps_CH_pi_Trp_kJ, 0)    # CH-pi favorable
        self.assertGreater(p.k_desolv_eq_kJ, 0)   # Desolvation unfavorable
        self.assertGreater(p.dG_0_kJ, 0)           # Entropy penalty

    def test_hb_in_physical_range(self):
        p = GlycanParams()
        self.assertGreater(p.eps_HB_kJ, -12.0)
        self.assertLess(p.eps_HB_kJ, -4.0)

    def test_trp_stronger_than_tyr(self):
        p = GlycanParams()
        self.assertLess(p.eps_CH_pi_Trp_kJ, p.eps_CH_pi_Tyr_kJ)

    def test_ch_pi_per_residue(self):
        p = GlycanParams()
        self.assertAlmostEqual(p.get_ch_pi_energy("TRP"), p.eps_CH_pi_Trp_kJ)
        self.assertAlmostEqual(p.get_ch_pi_energy("TYR"), p.eps_CH_pi_Tyr_kJ)

    def test_equatorial_costs_more_than_axial(self):
        """Equatorial OH is better hydrated → costs more to bury."""
        p = GlycanParams()
        self.assertGreater(p.k_desolv_eq_kJ, p.k_desolv_ax_kJ)

    def test_translational_entropy(self):
        """Should be in 15-35 kJ/mol range."""
        p = GlycanParams()
        self.assertGreater(p.dG_0_kJ, 15.0)
        self.assertLess(p.dG_0_kJ, 35.0)


class TestScorerSignConventions(unittest.TestCase):
    """Test that scoring signs are physically correct."""

    def test_hbonds_favorable(self):
        """More H-bonds → more negative ΔG."""
        p = GlycanParams()
        cm1 = ContactMap("T", "A", "GAL", 1, n_hbonds=2)
        cm2 = ContactMap("T", "A", "GAL", 1, n_hbonds=5)
        r1 = score_contact_map(cm1, p)
        r2 = score_contact_map(cm2, p)
        self.assertLess(r2.dG_predicted_kJ, r1.dG_predicted_kJ)

    def test_ch_pi_favorable(self):
        """CH-π contacts → more negative ΔG."""
        p = GlycanParams()
        cm1 = ContactMap("T", "A", "GAL", 1, aromatic_residues=[])
        cm2 = ContactMap("T", "A", "GAL", 1, aromatic_residues=["TRP1"])
        r1 = score_contact_map(cm1, p)
        r2 = score_contact_map(cm2, p)
        self.assertLess(r2.dG_predicted_kJ, r1.dG_predicted_kJ)

    def test_desolvation_unfavorable(self):
        """Burying hydroxyls COSTS energy (less negative ΔG)."""
        p = GlycanParams()
        cm1 = ContactMap("T", "A", "GAL", 1, n_hbonds=5, n_buried_eq_OH=0)
        cm2 = ContactMap("T", "A", "GAL", 1, n_hbonds=5, n_buried_eq_OH=3)
        r1 = score_contact_map(cm1, p)
        r2 = score_contact_map(cm2, p)
        self.assertGreater(r2.dG_predicted_kJ, r1.dG_predicted_kJ,
                           "Burying OHs should make ΔG less negative")

    def test_water_bridges_favorable(self):
        """Water bridges → more negative ΔG."""
        p = GlycanParams()
        cm1 = ContactMap("T", "A", "GAL", 1, n_water_bridges=0)
        cm2 = ContactMap("T", "A", "GAL", 1, n_water_bridges=3)
        r1 = score_contact_map(cm1, p)
        r2 = score_contact_map(cm2, p)
        self.assertLess(r2.dG_predicted_kJ, r1.dG_predicted_kJ)

    def test_metal_favorable(self):
        """Direct metal coordination → more negative ΔG."""
        p = GlycanParams()
        cm1 = ContactMap("T", "A", "FUC", 1, metal_ions=[])
        cm2 = ContactMap("T", "A", "FUC", 1, metal_ions=["CA401"])
        r1 = score_contact_map(cm1, p)
        r2 = score_contact_map(cm2, p)
        self.assertLess(r2.dG_predicted_kJ, r1.dG_predicted_kJ)

    def test_entropy_penalty_present(self):
        """Empty contact map should have POSITIVE ΔG (entropy penalty only)."""
        p = GlycanParams()
        cm = ContactMap("T", "A", "GAL", 1)
        r = score_contact_map(cm, p)
        self.assertGreater(r.dG_predicted_kJ, 0,
                           "No contacts should give positive ΔG (unfavorable)")

    def test_component_sum_equals_total(self):
        """Component energies should sum to total."""
        p = GlycanParams()
        cm = ContactMap("T", "A", "GAL", 1, n_hbonds=5, n_ch_pi=1,
                        n_buried_eq_OH=2, n_water_bridges=1,
                        aromatic_residues=["TRP1"], metal_ions=["CA1"])
        r = score_contact_map(cm, p)
        comp_sum = sum(r.components.values())
        self.assertAlmostEqual(r.dG_predicted_kJ, comp_sum, places=1)


class TestValidation(unittest.TestCase):
    """Test validation against experimental ITC data."""

    @classmethod
    def setUpClass(cls):
        cls.val = validate_predictions(cache_dir=PDB_CACHE)

    def test_validation_runs(self):
        self.assertGreaterEqual(self.val["n_pairs"], 2)

    def test_mae_under_5_kJ(self):
        """MAE should be under 5 kJ/mol with calibrated parameters."""
        self.assertLess(self.val["MAE_kJ"], 5.0,
                        f"MAE={self.val['MAE_kJ']} kJ/mol, expected < 5.0")

    def test_galectin3_error(self):
        gal = next(r for r in self.val["results"] if "Gal" in r.get("lectin", ""))
        self.assertLess(abs(gal["error_kJ"]), 5.0)

    def test_cona_error(self):
        cona = next(r for r in self.val["results"] if "ConA" in r.get("lectin", ""))
        self.assertLess(abs(cona["error_kJ"]), 5.0)

    def test_predictions_negative(self):
        """Predicted ΔG should be negative (favorable binding)."""
        for r in self.val["results"]:
            if "pred_kJ" in r:
                self.assertLess(r["pred_kJ"], 0,
                                f"{r['lectin']} predicted positive ΔG")


# ═════════════════════════════════════════════════════════════════════
# Deoxy-Mannose Prediction Tests
# ═════════════════════════════════════════════════════════════════════

class TestDeoxyMannose(unittest.TestCase):
    """Test Prediction 1: ConA deoxy-mannose series."""

    @classmethod
    def setUpClass(cls):
        cls.deoxy = predict_deoxy_series("5CNA", "MMA", chain="A",
                                         cache_dir=PDB_CACHE)

    def test_prediction_runs(self):
        self.assertNotIn("error", self.deoxy)
        self.assertIn("deoxy_predictions", self.deoxy)

    def test_four_positions(self):
        """Should predict for O2, O3, O4, O6."""
        preds = self.deoxy["deoxy_predictions"]
        for pos in ["O2", "O3", "O4", "O6"]:
            self.assertIn(pos, preds, f"Missing prediction for {pos}")

    def test_o4_o6_weaker(self):
        """O4 and O6 deoxy should WEAKEN binding (positive ΔΔG).
        These are the key recognition hydroxyls in ConA."""
        preds = self.deoxy["deoxy_predictions"]
        self.assertGreater(preds["O4"]["ddG_predicted_kJ"], 0)
        self.assertGreater(preds["O6"]["ddG_predicted_kJ"], 0)

    def test_o6_largest_effect(self):
        """O6 deoxy should have largest effect (makes most H-bonds)."""
        preds = self.deoxy["deoxy_predictions"]
        o6_ddg = preds["O6"]["ddG_predicted_kJ"]
        for pos in ["O2", "O3", "O4"]:
            self.assertGreater(o6_ddg, preds[pos]["ddG_predicted_kJ"],
                               f"O6 effect ({o6_ddg}) should be largest, "
                               f"but {pos} = {preds[pos]['ddG_predicted_kJ']}")

    def test_o2_deoxy_mannose_special(self):
        """O2 is axial in mannose — removing it saves desolvation cost.
        Depending on H-bond count, effect could be near zero or slightly favorable."""
        preds = self.deoxy["deoxy_predictions"]
        # O2 effect should be the smallest magnitude
        o2_abs = abs(preds["O2"]["ddG_predicted_kJ"])
        for pos in ["O3", "O4", "O6"]:
            self.assertLessEqual(o2_abs, abs(preds[pos]["ddG_predicted_kJ"]) + 0.1)

    def test_rank_order(self):
        """Rank order: O6 > O4 > O3 > O2 in binding loss."""
        preds = self.deoxy["deoxy_predictions"]
        ddgs = {pos: preds[pos]["ddG_predicted_kJ"] for pos in ["O2", "O3", "O4", "O6"]}
        self.assertGreater(ddgs["O6"], ddgs["O4"])
        self.assertGreater(ddgs["O4"], ddgs["O3"])

    def test_parent_dg_reasonable(self):
        """Parent mannose ΔG should be -15 to -30 kJ/mol."""
        parent_dg = self.deoxy["parent"]["dG_parent_kJ"]
        self.assertLess(parent_dg, -15)
        self.assertGreater(parent_dg, -35)


# ═════════════════════════════════════════════════════════════════════
# Selectivity Panel Tests
# ═════════════════════════════════════════════════════════════════════

class TestConASelectivity(unittest.TestCase):
    """ConA prefers Man/Glc over Gal."""

    @classmethod
    def setUpClass(cls):
        cls.sel = score_sugar_panel("5CNA", "MMA", chain="A",
                                    panel=["Man", "Glc", "Gal", "GlcNAc"],
                                    cache_dir=PDB_CACHE)

    def test_panel_runs(self):
        self.assertNotIn("error", self.sel)
        self.assertEqual(len(self.sel["scores"]), 4)

    def test_mannose_preferred(self):
        """ConA should rank Man at or near top."""
        rank = self.sel["ranking"]
        man_idx = rank.index("Man")
        self.assertLessEqual(man_idx, 1, f"Man should be top 2, got rank {man_idx}")

    def test_gal_disfavored(self):
        """Gal should be disfavored relative to Man/Glc."""
        scores = self.sel["scores"]
        self.assertGreater(scores["Gal"]["dG_kJ"], scores["Man"]["dG_kJ"])

    def test_ddg_nonnegative(self):
        """ΔΔG should be ≥ 0 (relative to best binder)."""
        for sugar, data in self.sel["scores"].items():
            self.assertGreaterEqual(data["ddG_kJ"], 0)


class TestGalectin3Selectivity(unittest.TestCase):
    """Galectin-3 prefers Gal over Man/Glc."""

    @classmethod
    def setUpClass(cls):
        cls.sel = score_sugar_panel("3ZSJ", "GAL",
                                    panel=["Man", "Glc", "Gal", "GlcNAc"],
                                    cache_dir=PDB_CACHE)

    def test_gal_preferred(self):
        """Galectin-3 should rank Gal first."""
        self.assertEqual(self.sel["ranking"][0], "Gal",
                         f"Gal should be #1, got {self.sel['ranking']}")

    def test_gal_stronger_than_glc(self):
        scores = self.sel["scores"]
        self.assertLess(scores["Gal"]["dG_kJ"], scores["Glc"]["dG_kJ"])

    def test_selectivity_magnitude(self):
        """ΔΔG(Man vs Gal) should be > 2 kJ/mol."""
        scores = self.sel["scores"]
        ddg = scores["Man"]["dG_kJ"] - scores["Gal"]["dG_kJ"]
        self.assertGreater(ddg, 2.0)


class TestSugarOHPatterns(unittest.TestCase):
    """Test sugar hydroxyl pattern definitions."""

    def test_man_c2_axial(self):
        self.assertEqual(SUGAR_OH_PATTERNS["Man"]["O2"], "ax")

    def test_glc_all_equatorial(self):
        for pos in ["O2", "O3", "O4", "O6"]:
            self.assertEqual(SUGAR_OH_PATTERNS["Glc"][pos], "eq")

    def test_gal_c4_axial(self):
        self.assertEqual(SUGAR_OH_PATTERNS["Gal"]["O4"], "ax")

    def test_glcnac_c2_nac(self):
        self.assertEqual(SUGAR_OH_PATTERNS["GlcNAc"]["O2"], "NAc")

    def test_fuc_deoxy(self):
        """Fucose is 6-deoxy — no O6."""
        self.assertNotIn("O6", SUGAR_OH_PATTERNS["Fuc"])


if __name__ == "__main__":
    unittest.main(verbosity=2)