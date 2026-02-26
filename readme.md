# Predicting Glycan Recognition Selectivity from Pure Chemistry

**Can we predict which sugars bind which receptors — and how tightly — using parameters measured from pure chemistry, without training on any protein-sugar binding data?**

## The Approach

Lectins discriminate between monosaccharides differing by a single hydroxyl group. We measure:

1. **Desolvation cost** per hydroxyl from published dissolution calorimetry (Schwarz 1996)
2. **CH-π stacking energy** from published synthetic host-guest data (Laughrey 2008)
3. **Conformational entropy** of glycosidic linkages from QM calculations (GLYCAM06)

Then test whether these three numbers — derived entirely from chemistry — predict lectin sugar preferences using published ITC binding data as the answer key.

## The Scoring Function

```
ΔG_pred(S,R) = ε_HB × n_HB           H-bonds (favorable)
             + ε_CH-π × n_CH-π        CH-π stacking (favorable)
             − k_desolv × n_OH_buried  desolvation penalty
             − ε_conf × n_frozen_tor   conformational entropy cost
             + ΔG_0                    offset
```

6 parameters. All from pure chemistry. ~50 lines of Python.

## Project Structure

```
glycan-scoring/
├── glycan_scorer.py          Core scoring function + dataclasses
├── sugar_properties.py       Phase G0: sugar property cards
├── parameter_fitting.py      Phases G1-G3: parameter determination
├── run_predictions.py        Predictions 1-7 + figure generation
├── data/
│   ├── schwarz_dissolution.csv      Phase G1 source (desolvation)
│   ├── jasra_polyols.csv            Phase G1 source (NAc)
│   ├── laughrey_ch_pi.csv           Phase G3 source (CH-π)
│   ├── chervenak_conA_deoxy.csv     Answer key: Prediction 1
│   ├── dam_brewer_compilation.csv   Answer keys: Predictions 2, 4, 6
│   ├── davis_synthetic_lectin.csv   Answer key: Prediction 7
│   ├── muraki_mutants.csv           Answer key: Prediction 3
│   └── contact_maps/                PDB-derived contact counts
├── predictions/                     (prediction output files)
├── figures/                         (generated poster figures)
├── tests/
│   └── test_glycan_scorer.py
├── PHASE_TRACKER.md                 Phase-by-phase checklist
└── requirements.txt
```

## Quick Start

```bash
pip install -r requirements.txt
python -m pytest tests/ -v          # run sanity checks
python glycan_scorer.py             # demo with placeholder params
python sugar_properties.py          # print sugar property cards
```

## Timeline

| Phase | Weeks | What |
|-------|-------|------|
| 0 | 1–2 | Literature compilation, PDB contact extraction |
| 1 | 3 | Parameter determination from pure chemistry |
| 2 | 3–4 | Predictions 1 + 4 (ConA deoxy + chain extension) |
| 3 | 5 | Prediction 3 (CH-π mutant validation) |
| 4 | 6–7 | Predictions 2, 6, 7 (epimers, cross-lectin, Davis) |
| 5 | 8–9 | Full validation scatter plot + R² |
| 6 | 10 | Poster assembly |

## Target

Poster at 9th Canadian Glycomics Symposium, May 12–14, 2026, Canmore AB.

## Authors

Hannah Sanford-Crane, PhD — MAAD Scientist Technologies Inc.
Busayo Gbenga-Ajose — Northeastern University, Toronto
Likhita Kombathula — Northeastern University, Toronto

## What Makes This Publishable

Parameters from sugar dissolution calorimetry (no receptor), synthetic host-guest binding (no protein), and QM calculations (no experiment) correctly predict which sugars bind which lectins and which synthetic receptors — including **opposite selectivity preferences in different binding geometries**.