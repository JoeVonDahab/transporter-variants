# TripleVar transporter variant release

349,056 ClinVar variants in 783 transporter, ion-channel, ABC and ATPase genes, scored for pathogenicity and calibrated to ACMG evidence strength.

## Files

| File | What it is | Licence |
|---|---|---|
| `transporter_variants_routed.csv.gz` | Each variant served by the better track: `plus` where AlphaMissense exists, `native` elsewhere. Highest accuracy. | Inherits CC BY-NC-SA 4.0 from AlphaMissense |
| `transporter_variants_native.csv.gz` | Our model only, every variant, no third-party predictor. | Freely redistributable |

## Columns

| Column | Meaning |
|---|---|
| `gene`, `family` | HGNC symbol and transporter family |
| `variation_id`, `rsid`, `hgvs` | ClinVar and dbSNP identifiers |
| `consequence`, `protein_variant` | Variant class and protein change where defined |
| `clinsig`, `review_status` | The ClinVar assertion at download, for reference only |
| `in_training` | TRUE if the model saw this variant. **Exclude these rows for any performance estimate.** |
| `track` | `plus` (with AlphaMissense) or `native` (ours alone) |
| `has_alphamissense` | Whether AlphaMissense genuinely covers this variant |
| `score` | Pathogenicity score, 0 to 1, from the named track |
| `lr` | Local likelihood ratio at this score (ClinGen/Pejaver 2022 method, see `CALIBRATION_CHECK.md`) |
| `evidence` | Evidence strength from the ClinGen PP3/BP4 calibration: pathogenic very strong / strong / moderate / supporting, indeterminate, benign supporting / moderate / strong |
| `acmg_code` | The matching ACMG/AMP code for missense variants (PP3_Supporting ... PP3_Strong, BP4_Supporting ... BP4_Strong); blank for indeterminate and for non-missense variants |
| `evidence_v1` | The 1.0.x evidence label (Tavtigian OddsPath bins at the dataset base rate), kept for comparison |

## How to use it

Use `evidence` (or `acmg_code` for missense), not `score`. A 0.5 cut-off over-calls badly: the model's training set is 45% pathogenic while a gene panel is mostly not. From 1.1.0 the levels follow ClinGen's calibration for computational evidence (Pejaver et al. 2022): a 0.0441 prior, a local likelihood ratio, and thresholds set at the one-sided 95% bootstrap bound, so a level is assigned only where even the pessimistic estimate reaches it. Maximum levels (pathogenic / benign): missense, ours alone pathogenic moderate / benign moderate; missense with AlphaMissense pathogenic moderate / benign strong; other classes pathogenic strong / benign strong. PP3 is capped at Strong in `acmg_code`, as ClinGen recommends for computational evidence; do not combine PP3 with PVS1 for predicted loss-of-function variants.

## What it is not

- `in_training` rows are not held out and must not be used to estimate accuracy.
- The `plus` track imputes AlphaMissense where it is absent, so trust it only where `has_alphamissense` is TRUE; otherwise read the `native` track.
- Our model alone loses to AlphaMissense on human missense variants (0.9189 against 0.9455); the routed track wins (0.9553) because it combines them. Across all classes ours reaches 0.9853, and AlphaMissense scores nothing outside missense.
- ClinVar assertions are partly informed by computational evidence, so some circularity cannot be excluded from any benchmark built on them.

## Version 1.0.1

The activity and mechanism layers are completed for missense variants whose ESM-2 embeddings were in the training-pool feature files rather than the published-variant file, including every variant in SLC19A2 and SLC19A3. Activity now covers 176,153 missense variants. Pathogenicity scores, likelihood ratios and evidence bins are unchanged, and every activity or mechanism value present in 1.0.0 is unchanged. Missense variants whose ClinVar reference residue does not match the UniProt sequence still carry no activity or mechanism.

## Version 1.0.2

2,198 missense variants that carried only a consequence-class score now carry a protein-level (missense-branch) score: multi-nucleotide substitutions and missense variants on non-canonical transcripts that AlphaMissense does not map, and variants whose AlphaMissense numbering did not fit the UniProt sequence. Their ClinVar protein change was placed on the UniProt sequence by aligning the ClinVar transcript's reported residues; `protein_variant` for these rows is now in UniProt numbering. A further 94 missense variants in the same genes had been scored at the wrong residue (AlphaMissense's residue matched the UniProt sequence by coincidence) and are rescored at the correct one. Evidence bins use the 1.0.1 calibration unchanged, and every other row is identical to 1.0.1. 654 missense variants could not be placed unambiguously (mostly residues absent from the UniProt isoform) and keep their class-level score.

## Version 1.1.0

Evidence levels are recalibrated with ClinGen's PP3/BP4 method (Pejaver et al. 2022) in place of the 1.0.x binning, which on held-out genes had over-stated its top level ("pathogenic very strong" observed at strong) and under-stated its middle one ("indeterminate" observed at supporting). New column `acmg_code`; the old label is kept as `evidence_v1`; `lr` is now the local likelihood ratio. Scores are unchanged from 1.0.2. `mechanism_pred` coverage follows the new levels (moderate or stronger). See `CALIBRATION_CHECK.md` and `CHANGELOG.md`.
