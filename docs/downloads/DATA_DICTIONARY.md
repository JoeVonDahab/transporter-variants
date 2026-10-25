# Data dictionary

349,056 rows, one per ClinVar variant, in 783 genes.

| Column | Type | Meaning |
|---|---|---|
| `gene` | string | HGNC symbol |
| `family` | string | SLC, ion channel, ABC, P-type ATPase, connexin/pannexin, V-type ATPase, aquaporin or CFTR |
| `variation_id` | integer | ClinVar VariationID |
| `rsid` | string | dbSNP identifier where one exists |
| `hgvs` | string | HGVS expression as given by ClinVar |
| `consequence` | string | missense, synonymous, intronic, splice_region, splice_canonical, frameshift, nonsense, inframe_indel, utr3, utr5, start_loss, stop_loss, other |
| `protein_variant` | string | protein change, where the consequence defines one |
| `clinsig` | string | the ClinVar assertion at download; reference only, not a target |
| `review_status` | string | ClinVar review status for that assertion |
| `in_training` | boolean | **TRUE means the model saw this variant. Exclude these rows from any performance estimate.** |
| `track` | string | `plus` where AlphaMissense covers the variant, else `native` |
| `has_alphamissense` | boolean | whether AlphaMissense genuinely covers this variant |
| `score` | float 0-1 | pathogenicity score from the named track |
| `lr` | float | likelihood ratio of this score bin against the base rate |
| `evidence` | string | ACMG strength derived from `lr`: pathogenic very strong / strong / moderate / supporting, indeterminate, benign supporting / moderate / strong |
| `v2_native` | float | our model alone, for every variant |
| `v2_plus` | float | our model with AlphaMissense; **imputes AlphaMissense where absent, so read it only where `has_alphamissense` is TRUE** |
| `alphamissense` | float | AlphaMissense score where available |

## Evidence thresholds

Likelihood ratios are binned at the odds of pathogenicity behind the ACMG combining rules (Tavtigian et al. 2018): 350 very strong, 18.7 strong, 4.33 moderate, 2.08 supporting, and the reciprocals for benign evidence.

## Observed purity of each bin

Measured on gene-disjoint halves of the held-out labelled set, so no gene calibrated the bin its own variants fall in. See `CALIBRATION_CHECK.md` for the full table.