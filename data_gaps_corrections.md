# Data Gaps & Reference Corrections
## Glycan Scoring Capstone Project — Phase 0 Status

**Last updated:** 2026-02-25

---

## CRITICAL REFERENCE CORRECTION

### Schwarz dissolution calorimetry paper

**Spec cites:** Schwarz FP (1996) Thermochim. Acta 284:129
**Actual paper:** Schwarz FP (1997) J. Solution Chem. 26:471-485
- "Enthalpy of solution of carbohydrates using a modified DSC"
- DOI: 10.1007/BF00972993

**Actual compound coverage (from abstract):**

| Compound class | Positions covered | Notes |
|---|---|---|
| Deoxy-glucose | C1, C2, C3, C6 | **NO C4-deoxy** (only C4-fluoro) |
| Fluoro-glucose | C1, C2, C3, C4, C6 | F replaces OH — different physics |
| Mannose | Parent ONLY | **No mannose deoxy series** |
| Galactose | NOT covered | — |
| Methoxy derivatives | αMeOGlu, αPheOGlu, 3MeOGlu | — |
| Temperatures | 15.1, 25.0, 35.0, 45.1°C | — |

**Impact on Phase G1:**
- Glucose deoxy series sufficient for equatorial OH calibration (C2, C3, C6 all eq)
- C1 axial in α-glucose gives one axial data point
- **FALLBACK for mannose deoxy:** use glucose ΔΔH values where same stereochemistry applies (physically justified: per-OH hydration depends on eq/ax orientation, not ring identity)
- C4-deoxy gap: 4-fluoro data available but F≠H substitution changes physics. Alternative: use 4-deoxy-galactose from Banipal (if available) or treat as interpolation.

### Lysozyme mutant reference

**Spec cites:** Muraki 1997
**Actual key papers:**
- Kumagai et al. (1993) Eur. J. Biochem. 212:151 — binding constants for W62 mutants
- Maenaka & Kumagai (1995) J. Mol. Biol. 247:281 — crystal structures of mutants
- García-Hernández & Hernández-Arana (2002) — ITC thermodynamics (ΔCp = -83 cal/(mol·K))

---

## ANSWER KEY DATA — POPULATED vs GAPS

### ✓ POPULATED (from web-accessible sources)

**ConA monosaccharides (Prediction 2):**
- MeαMan: Ka = 8,240 M⁻¹, ΔG = -5.35, ΔH = -8.4 kcal/mol (Chervenak 1995 / Dam & Brewer 2002)
- MeαGlc: Ka = 2,100 M⁻¹, ΔG = -4.53 kcal/mol
- D-mannose: Ka ≈ 8,200 M⁻¹
- MeαGal: no detectable binding
- GlcNAc: Ka ≈ 920 M⁻¹ (from ConA general specificity)
- ManNAc: Ka weak

**ConA oligosaccharides (Prediction 4):**
- MeαMan: Ka = 8,240 M⁻¹
- Trimannoside core: Ka = 353,000 M⁻¹ (Naismith 1998)
- Pentasaccharide: Ka = 1.41×10⁶ M⁻¹

**WGA binding (Prediction 6):**
- GlcNAc: Ka = 400 M⁻¹, ΔH = -6.1 kcal/mol at 26°C (Bains et al. 1992)
- (GlcNAc)₂: Ka = 5,300 M⁻¹, ΔH = -15.6 kcal/mol
- (GlcNAc)₃: Ka = 11,100 M⁻¹, ΔH = -19.4 kcal/mol
- Enthalpically driven (TΔS always negative)

### ✗ GAPS (require paper retrieval)

| Answer key | What's missing | Paper needed | Priority |
|---|---|---|---|
| ConA deoxy (Pred 1) | All ΔΔG values | Chervenak & Toone 1995 Table | HIGH |
| ConA oligosaccharides | Man-α1,3-Man, Man-α1,6-Man | Dam & Brewer 2002 Table 2 | MED |
| Lysozyme mutants (Pred 3) | All mutant Ka/ΔΔG | Kumagai 1993 Table | HIGH |
| PNA (Pred 6) | Gal Ka, Glc Ka | Swaminathan 1998 | MED |
| Galectin-3 (Pred 6) | LacNAc Ka, Gal Ka | Seetharaman 1998 | MED |
| SBA (Pred 6) | GalNAc Ka, Gal Ka | Dam & Brewer 2002 Table | MED |
| Davis receptor (Pred 7) | Structure + Ka panel | Davis 2012 Nature Chem 4:548 | HIGH |
| Dissolution calor. | Actual ΔH_sol table values | Schwarz 1997 J. Sol. Chem. | HIGH |
| Dissolution calor. | Jasra sugar/polyol ΔH_sol | Jasra 1982 J. Sol. Chem. 11:325 | HIGH |
| CH-π calibration | ε values per contact | Laughrey 2008 | MED |

---

## STRATEGY FOR REMAINING GAPS

1. **Papers behind paywall** (Schwarz, Jasra, Chervenak, Kumagai, Davis):
   Students or PI retrieve via institutional access. These are the Phase 0 deliverables.

2. **CarboGrove CFG array data** (27,277 entries already extracted):
   Use for Prediction 6 VALIDATION only (rank-ordering), NOT calibration.
   Semi-quantitative RFU, not Ka.

3. **LfDB Ka values** (frontal affinity chromatography):
   Download from https://glycosmos.org/lfdb
   Quantitative Ka for lectin-monosaccharide pairs. May fill PNA/Gal-3/SBA gaps.

4. **Fallback strategy for Schwarz gaps:**
   If mannose deoxy dissolution data truly don't exist, use glucose ΔΔH by stereochemistry.
   Document assumption explicitly. This is the RISK CHECK flagged in Phase 0.