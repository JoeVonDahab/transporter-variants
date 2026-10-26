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
| `lr` | float | local likelihood ratio at this score: (P/LP density / B/LB density) in a sliding score window, ClinGen method (Pejaver et al. 2022); prior-free |
| `evidence` | string | evidence strength from the ClinGen PP3/BP4 calibration: pathogenic very strong / strong / moderate / supporting, indeterminate, benign supporting / moderate / strong |
| `acmg_code` | string | ACMG/AMP code for missense variants: PP3_Supporting, PP3_Moderate, PP3_Strong (also for `evidence` = pathogenic very strong; PP3 is capped at Strong), BP4_Supporting, BP4_Moderate, BP4_Strong; blank for indeterminate and for every non-missense variant |
| `evidence_v1` | string | the 1.0.x label (Tavtigian OddsPath bins at the dataset base rate), kept for comparison; superseded |
| `v2_native` | float | our model alone, for every variant |
| `v2_plus` | float | our model with AlphaMissense; **imputes AlphaMissense where absent, so read it only where `has_alphamissense` is TRUE** |
| `alphamissense` | float | AlphaMissense score where available |

## Evidence thresholds

ClinGen's calibration for computational evidence (Pejaver et al. 2022): prior probability of pathogenicity 0.0441; local likelihood ratio estimated in a sliding window around each score (at least 100 labelled variants and 3% of the unit's variants); thresholds where the one-sided 95% bootstrap bound reaches LR 2.406 (supporting), 5.790 (moderate), 33.53 (strong), 1,124 (very strong), and the reciprocals for benign. A level the bound never reaches is not assigned. Calibrated separately for our score on missense, our score on other classes, and our score with AlphaMissense where AlphaMissense exists, on labelled variants the model never trained on.

| Unit | benign strong | benign moderate | benign supporting | pathogenic supporting | pathogenic moderate | pathogenic strong | pathogenic very strong |
|---|---|---|---|---|---|---|---|
| native_missense | - | <= 0.1747 | <= 0.2469 | >= 0.9289 | >= 0.9724 | - | - |
| native_other | <= 0.0045 | <= 0.0203 | <= 0.4184 | >= 0.5677 | >= 0.5677 | >= 0.9441 | - |
| plus | <= 0.0094 | <= 0.0920 | <= 0.3499 | >= 0.9525 | >= 0.9660 | - | - |

Units: `native_missense` = `v2_native` on missense (native track; routed rows without AlphaMissense); `native_other` = `v2_native` on every other class; `plus` = `v2_plus` where `has_alphamissense` is TRUE (routed track).

## Observed evidence of each level

Measured on gene-disjoint halves of the held-out labelled set, so no gene calibrated the level its own variants fall in. See `CALIBRATION_CHECK.md` for the full table, including the observed likelihood ratio of every level.

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
