# Glycan Scoring — Phase Tracker

**Target:** Poster at 9th Canadian Glycomics Symposium, May 12–14, 2026, Canmore AB
**Timeline:** ~10 weeks from March 1, 2026

---

## Phase 0: Literature + Data Compilation (Weeks 1–2)

**Status: IN PROGRESS — database pulls complete, papers still needed**

### Database pulls (DONE 2026-03-02)
- [x] LfDB Ka values from GlyCosmos (418 entries, 6 lectins)
  - WGA: 88 Ka values (GlcNAc oligomers up to Ka=1358 × 10⁴ M⁻¹)
  - PNA: 40 Ka values (Gal-containing, up to Ka=64.3 × 10⁴ M⁻¹)
  - Galectin-3: 132 Ka values (full + CRD, up to Ka=591 × 10⁴ M⁻¹)
  - SBA: 96 Ka values (GalNAc-specific, up to Ka=92.4 × 10⁴ M⁻¹)
  - Conarva (ConA-like): 62 Ka values (Man/Glc, up to Ka=150 × 10⁴ M⁻¹)
  - **NOTE: ConA (LfDB0170) has NO FAC data in LfDB**
- [x] GlyTouCan IDs resolved to IUPAC names (250K glycan lookup loaded)
- [x] CarboGrove CFG glycan array data (26,863 entries for 5 project lectins)
- [x] GlyCosmos lectins cross-reference (5,099 lectins)
- [x] ChEMBL ConA activities (12 entries, all IC50 for glycoconjugates — not useful)
- [ ] ~~BindingDB lectin binding~~ (domain blocked: 403 Forbidden)
- [ ] ~~ProCarbDB mutagenesis~~ (server down: procarbdb.bioinfo.se)
- [x] PDB structures downloaded (6 files: 5CNA, 1CVN, 2UVO, 2PEL, 3GAL→fixed to 1A3K, 1LZB)
- [x] Contact maps auto-extracted (BioPython-based, all 6 structures)

### Answer keys populated from databases
- [x] Prediction 2: Cross-lectin selectivity — 12 glycans × 5 lectins (LfDB Ka)
- [x] Prediction 4: Chain enhancement — WGA GlcNAc series (n=2→5) + PNA Gal series (LfDB Ka)
- [x] Prediction 6: Selectivity matrix — 86 glycans × 6 lectins with IUPAC names (LfDB Ka)

### Paper retrieval (STILL NEEDED — no database substitute)
- [x] Schwarz 1996 (Biochem. J. 316:123) — ConA/pea/lentil deoxy binding Kb + ΔH
  - 16 derivatives × 3 lectins: Kb, ΔG, ΔH, TΔS at 10°C and 25°C
  - Critical: C3, C4, C6 deoxy = NO BINDING (defines pharmacophore)
  - ΔΔH values extracted for parameter fitting (Table 4)
  - NOTE: This is the BINDING paper. Dissolution paper (J. Solution Chem.) still needed
- [ ] Jasra 1982 (J. Solution Chem. 11:325) — polyol ΔH_sol
- [x] Schwarz 1996 (Biochem. J. 316:123) — **Prediction 1 answer key** (ConA deoxy Kb+ΔH)
  - [ ] Chervenak 1995 (Biochemistry 34:5685) — ConA oligosaccharide ITC (still useful, lower priority)
- [ ] Kumagai 1993 (Eur J Biochem 212:151) — **Prediction 3 answer key** (lysozyme mutants)
- [ ] Dam & Brewer 2002 (Chem. Rev. 102:387) — ConA monosaccharide Ka + Prediction 5
- [ ] Davis 2012 (Nature Chem. 4:548) — **Prediction 7 answer key** (synthetic receptor)
- [ ] Laughrey 2008 — CH-π per-contact ε values
- [ ] Kirschner 2008 / GLYCAM06 — QM torsion surfaces

### RISK CHECK
- [ ] **CONFIRM** Schwarz has systematic deoxy coverage
  - If not: locate Galema, Banipal, or Kozak alternatives
  - **DO NOT proceed past Phase 0 without a confirmed desolvation source**

### LfDB data note
LfDB Ka values are for PA-labeled oligosaccharides (disaccharides+), measured by
frontal affinity chromatography at 25°C. These are NOT monosaccharide Ka values.
For monosaccharide Ka (ConA + Glc/Man/Gal), Dam & Brewer 2002 is required.
Ka units: 10⁴ M⁻¹ (Hirabayashi et al. 2002 Biochim Biophys Acta 1572:232)

---

## Phase 1: Parameter Determination (Week 3)

**Status: BLOCKED — waiting for Schwarz/Jasra/Laughrey papers**

### G1: Desolvation
- [ ] Extract ΔΔH per position from Schwarz data → schwarz_dissolution.csv
- [ ] Compute k_desolv_eq, k_desolv_ax, k_desolv_NAc
- [ ] **CHECK:** k_desolv_eq > k_desolv_ax

### G2: Conformational entropy
- [ ] Tabulate QM torsion profiles from GLYCAM06
- [ ] Compute TΔS_freeze for each linkage type
- [ ] Set eps_conf = mean across linkages

### G3: CH-π
- [ ] Collect Ka values from Laughrey/Diederich/Asensio
- [ ] Subtract locked HG terms
- [ ] Compute eps_CH_pi = residual / n_contacts

---

## Test Status
- **20/20 tests passing** (test_glycan_scorer.py, from previous session)
- Covers: sign conventions, scoring function, batch scoring, statistics, deoxy ΔΔG, epimer selectivity

## Repository Files (as of 2026-03-02)
### Code (4 modules)
- glycan_scorer.py — scoring function + dataclasses
- sugar_properties.py — Phase G0 hydroxyl maps for 16 monosaccharides
- parameter_fitting.py — Phase G1-G3 parameter extraction
- run_predictions.py — all 7 predictions with figure generators

### Data
- data/lfdb/ — 6 per-lectin Ka CSVs + resolved names + raw LfDB dump (3,246 entries)
- data/answer_keys/ — Predictions 2, 4, 6 populated from LfDB
- data/carbogrove_project_lectins.csv — 26,863 CFG array binding intensities
- data/glycosmos_lectins_list.csv — 5,099 lectins cross-reference
- data/chembl_cona_activities.csv — ChEMBL ConA (IC50 only, limited)
- data/contact_maps/ — PDB-derived contact counts

### Scripts
- scripts/fetch_lfdb.py — re-runnable LfDB download + processing
- scripts/fetch_data.py — PDB + contact analysis
- scripts/contact_analyzer.py — BioPython contact extractor
