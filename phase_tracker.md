# Phase Tracker — Glycan Scoring Capstone
**Last updated:** 2026-02-25

---

## Phase 0: Literature + Data Compilation

### Structural Data (PDB + Contacts)
- [x] PDB structures downloaded: 5CNA, 1CVN, 2UVO, 2PEL, 1A3K, 1LZB
- [x] BioPython contact analyzer with geometry-verified CH-π
- [x] contact_maps_auto.csv: 35 binding sites across 6 structures
- [x] Trp62 CH-π confirmed in 1LZB (3.77 Å, 5.3° from normal)
- [x] Trp181 CH-π confirmed in 1A3K (3.68 Å, 5.7° from normal)
- [x] TYR12 confirmed in 5CNA + 1CVN (ConA conserved)
- [x] TYR66+TYR73 tandem pair in 2UVO (WGA)
- [x] TYR125 in 2PEL (PNA)

### Database Integration
- [x] GlyCosmos master lectin list parsed (5,099 lectins)
- [x] CarboGrove CFG array data extracted (27,277 entries, 6 project lectins)
- [x] Project lectin cross-reference built (GL→UniProt→LfDB→PDB→Predictions)
- [ ] **LfDB Ka values download** — separate from CarboGrove; quantitative Ka data

### Answer Key Construction
- [x] Answer key templates created (5 CSVs in data/answer_keys/)
- [x] ConA monosaccharide Ka values populated (MeαMan, MeαGlc, MeαGal, D-Man)
- [x] ConA oligosaccharide partial (mono, tri, penta populated; disaccharides gap)
- [x] WGA ITC values populated (GlcNAc, chitobiose, chitotriose — Bains 1992)
- [ ] ConA deoxy series ΔΔG: **NEEDS Chervenak 1995 Table**
- [ ] Lysozyme mutant ΔΔG: **NEEDS Kumagai 1993 Table**
- [ ] PNA monosaccharide Ka: **NEEDS Swaminathan 1998 or LfDB**
- [ ] Galectin-3 Ka: **NEEDS Seetharaman 1998 or LfDB**
- [ ] SBA Ka: **NEEDS Dam & Brewer 2002 Table or LfDB**
- [ ] Davis receptor Ka panel: **NEEDS Davis 2012 Nature Chem**

### Reference Corrections (documented in DATA_GAPS_AND_CORRECTIONS.md)
- [x] Schwarz citation corrected: J. Solution Chem. 26:471 (NOT Thermochim. Acta 284)
- [x] Schwarz scope documented: glucose deoxy C1,C2,C3,C6 only; no C4-deoxy, no mannose deoxy, no galactose
- [x] Lysozyme mutant ref corrected: Kumagai 1993 (NOT Muraki 1997)
- [x] Fallback strategy documented for mannose deoxy gap

### Parameter Source Data
- [ ] **Schwarz ΔH_sol table values**: paper paywalled, need institutional access
- [ ] **Jasra ΔH_sol table values**: paper paywalled
- [ ] **Laughrey CH-π ε values**: paper paywalled
- [ ] **GLYCAM06 torsion profiles**: need QM data from Kirschner 2008

---

## Summary: What's Done vs What Needs Humans

### Computationally complete:
- All code modules (glycan_scorer, sugar_properties, parameter_fitting, run_predictions)
- All data templates (8 CSVs, schwarz/jasra/laughrey/etc)
- PDB download + contact analysis (6 structures, 35 sites, geometry-verified CH-π)
- GlyCosmos/CarboGrove integration (27,277 CFG array entries)
- Answer key templates with 7/15 cross-lectin entries populated
- Test suite: 20/20 passing
- Reference corrections documented

### Blocked on human paper retrieval:
1. Schwarz 1997 J. Sol. Chem. — dissolution table values (k_desolv_eq, k_desolv_ax)
2. Jasra 1982 — polyol dissolution (k_desolv_NAc)
3. Chervenak 1995 — ConA deoxy ITC (Prediction 1 answer key)
4. Kumagai 1993 — lysozyme mutant Ka (Prediction 3 answer key)
5. Davis 2012 — synthetic receptor Ka panel (Prediction 7 answer key)
6. Laughrey 2008 — CH-π per-contact ε values
7. Kirschner 2008 / GLYCAM06 — QM torsion potential surfaces

### Recommended next action:
Students retrieve papers 1-7 via Northeastern library access.
LfDB CSV download from GlyCosmos may fill PNA/Gal-3/SBA Ka gaps.
Once paper tables are in hand, Phase 1 parameter determination can begin immediately.

---

## Test Status
- **20/20 tests passing** (test_glycan_scorer.py)
- Covers: sign conventions, scoring function, batch scoring, statistics, deoxy ΔΔG, epimer selectivity

## Repository Files
- **~45 files** in glycan-scoring/
- Core code: 4 Python modules
- Data: 8 template CSVs + 5 answer key CSVs + 6 PDB files + contact maps
- Scripts: fetch_data.py, contact_analyzer.py
- Docs: README.md, PHASE_TRACKER.md, DATA_GAPS_AND_CORRECTIONS.md, data/README.md