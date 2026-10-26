# Changelog

## 1.1.0 (2026-10-26)

- Evidence levels recalibrated with ClinGen's PP3/BP4 method (Pejaver et al. 2022): prior 0.0441, local likelihood ratio, thresholds at the one-sided 95% bootstrap bound, calibrated per unit (missense; other classes; with AlphaMissense). Replaces the Tavtigian OddsPath binning at the dataset base rate.
- New `acmg_code` (missense only; PP3 capped at Strong). Old label kept as `evidence_v1`. `lr` is now the local LR.
- `mechanism_pred` coverage follows the new levels (moderate or stronger).
- Scores identical to 1.0.2.

## 1.0.2 (2026-10-26)

- 2,198 missense variants that had only a consequence-class score now carry a protein-level score (protein change placed on UniProt by aligning the ClinVar transcript's reported residues).
- 94 missense variants that had been scored at the wrong residue (AlphaMissense numbering coinciding with ours) rescored at the correct one.
- `protein_variant` for these rows is in UniProt numbering. All other rows identical to 1.0.1.

## 1.0.1 (2026-10-25)

- Activity and mechanism layers completed for 18,530 missense variants.

## 1.0.0 (2026-10-21)

- First release.
