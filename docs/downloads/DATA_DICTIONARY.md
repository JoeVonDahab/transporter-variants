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

## Activity and mechanism layers

| Column | Meaning |
|---|---|
| `activity_pred` | predicted transport or enzymatic activity, 0 (abolished) to 1 (wild-type-like), on the rank scale the deep mutational scans were normalised to. Defined for missense variants with ESM-2 features. |
| `mechanism_pred` | among variants called confidently pathogenic, the predicted probability that loss of function is **activity-specific** (protein reaches its compartment but cannot work) rather than **expression-driven** (never arrives). A regression output, so a few values fall slightly outside 0 to 1. Blank elsewhere, because the question is meaningless for a variant that does not abolish function. |
| `layers_in_sample` | TRUE if the variant's protein carries a deep mutational scan that trained these two models. Their predictions are in-sample for those proteins. |

### How much to trust them

Activity transfers to a protein with no measurements at Spearman 0.50 alone and 0.54 combined with AlphaMissense, against 0.52 for AlphaMissense alone; it is a graded estimate, not a call.

Mechanism transfers at AUC 0.78 across proteins and replicates at 0.73 on assay pairs curated by ProteinGym, including 0.70 on a fungal enzyme absent from our data. But **per-variant calls on genes with no measurements have not been shown to be reliable**: on SLC19A2 and SLC19A3 two independent models agreed at 0.009 where they agree at 0.536 on proteins with known answers. Treat this column as a hypothesis generator, not as evidence about any individual variant.

The two layers were fitted on deep mutational scans, not on ClinVar labels, so `in_training` (which refers to the pathogenicity model) does not make them in-sample.
