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
| `lr` | Likelihood ratio of this score bin against the base rate |
| `evidence` | ACMG evidence strength from `lr`, at the Tavtigian 2018 thresholds |

## How to use it

Use `evidence`, not `score`. A 0.5 cut-off over-calls badly: the model's training set is 45% pathogenic while a gene panel is mostly not, so a raw probability is miscalibrated for that setting. The likelihood ratio is what transfers, and the bins map onto the ACMG combining rules directly.

## What it is not

- `in_training` rows are not held out and must not be used to estimate accuracy.
- The `plus` track imputes AlphaMissense where it is absent, so trust it only where `has_alphamissense` is TRUE; otherwise read the `native` track.
- Our model alone loses to AlphaMissense on human missense variants (0.9163 against 0.9455); the routed track wins (0.9538) because it combines them. On every other variant class AlphaMissense offers nothing and ours reaches 0.9851.
- ClinVar assertions are partly informed by computational evidence, so some circularity cannot be excluded from any benchmark built on them.

## Version 1.0.1

The activity and mechanism layers are completed for missense variants whose ESM-2 embeddings were in the training-pool feature files rather than the published-variant file, including every variant in SLC19A2 and SLC19A3. Activity now covers 176,153 missense variants. Pathogenicity scores, likelihood ratios and evidence bins are unchanged, and every activity or mechanism value present in 1.0.0 is unchanged. Missense variants whose ClinVar reference residue does not match the UniProt sequence still carry no activity or mechanism.
