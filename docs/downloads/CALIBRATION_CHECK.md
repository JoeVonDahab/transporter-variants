# Do the evidence levels mean what they say? (release 1.1.0)

Evidence levels follow ClinGen's PP3/BP4 calibration (Pejaver et al. 2022): prior 0.0441, local likelihood ratio over a sliding score window, thresholds at the one-sided 95% bound of 10,000 bootstrap resamples. Calibrated per unit (our score on missense; our score on other classes; our score with AlphaMissense where it exists) on labelled out-of-sample variants only.

## Score thresholds

| Unit | serves | labelled (P/LP) | benign strong | benign moderate | benign supporting | pathogenic supporting | pathogenic moderate | pathogenic strong | pathogenic very strong | maximum levels |
|---|---|---|---|---|---|---|---|---|---|---|
| native_missense | missense, native track; missense without AlphaMissense, routed track | 2,011 (1,496) | not reached | <= 0.1747 | <= 0.2469 | >= 0.9289 | >= 0.9724 | not reached | not reached | pathogenic moderate; benign moderate |
| native_other | every other class, both tracks | 6,468 (2,985) | <= 0.0045 | <= 0.0203 | <= 0.4184 | >= 0.5677 | >= 0.5677 | >= 0.9441 | not reached | pathogenic strong; benign strong |
| plus | routed track where AlphaMissense exists | 1,998 (1,487) | <= 0.0094 | <= 0.0920 | <= 0.3499 | >= 0.9525 | >= 0.9660 | not reached | not reached | pathogenic moderate; benign strong |

Scores between the innermost benign and pathogenic thresholds are indeterminate.

## Gene-disjoint validation

Calibration refitted on half the genes and read off the other half, both directions (112's split), so no gene calibrates the level its own variants fall in. 8,479 labelled out-of-sample variants. Observed LR is class-conditional ((P/LP in level / all P/LP) / (B/LB in level / all B/LB)), so it does not depend on the dataset's mix; the posterior applies it to the 0.0441 prior. A flag marks a level whose observed LR falls outside its range; 'over-claims' means weaker evidence than the label says.

### routed track, ClinGen calibration (1.1.0)

| Level | n | observed P/LP | observed LR [95% CI] | posterior at prior 0.0441 | flag |
|---|---|---|---|---|---|
| benign strong | 849 | 0.8% | 0.00742 [0.00215, 0.0131] | 0.000 |  |
| benign moderate | 2,684 | 2.8% | 0.026 [0.0204, 0.0318] | 0.001 | below range: conservative |
| benign supporting | 498 | 15.3% | 0.161 [0.123, 0.199] | 0.007 | below range: conservative |
| indeterminate | 479 | 85.4% | 5.21 [4.15, 6.88] | 0.194 | above range: informative (pathogenic) |
| pathogenic supporting | 408 | 96.1% | 21.9 [14.1, 40.3] | 0.502 | above range: conservative |
| pathogenic moderate | 1,246 | 97.6% | 36.2 [26.3, 54.7] | 0.625 | above range: conservative |
| pathogenic strong | 2,315 | 99.6% | 206 [123, 505] | 0.905 |  |

### routed track, for comparison: the 1.0.x binning on the same halves, against its own claimed ranges

| Level | n | observed P/LP | observed LR [95% CI] | posterior at prior 0.0441 | flag |
|---|---|---|---|---|---|
| benign strong | 3,387 | 1.6% | 0.0145 [0.0109, 0.0183] | 0.001 |  |
| benign moderate | 335 | 17.0% | 0.183 [0.135, 0.236] | 0.008 |  |
| benign supporting | 318 | 30.8% | 0.397 [0.308, 0.498] | 0.018 |  |
| indeterminate | 232 | 82.3% | 4.16 [3.08, 6.08] | 0.161 | above range: informative (pathogenic) |
| pathogenic supporting | 220 | 80.0% | 3.57 [2.63, 5.13] | 0.141 |  |
| pathogenic moderate | 734 | 95.2% | 17.8 [12.9, 26.4] | 0.451 |  |
| pathogenic strong | 2,091 | 99.5% | 169 [103, 370] | 0.886 |  |
| pathogenic very strong | 1,135 | 99.2% | 112 [65.8, 259] | 0.837 | BELOW range: over-claims (CI excludes range) |
| unscored | 27 | 0.0% | 0 [0, 0] | 0.000 |  |

### native track, ClinGen calibration (1.1.0)

| Level | n | observed P/LP | observed LR [95% CI] | posterior at prior 0.0441 | flag |
|---|---|---|---|---|---|
| benign strong | 610 | 0.3% | 0.00293 [0, 0.00765] | 0.000 |  |
| benign moderate | 2,850 | 3.1% | 0.0284 [0.0224, 0.0342] | 0.001 | below range: conservative |
| benign supporting | 552 | 18.5% | 0.202 [0.159, 0.247] | 0.009 |  |
| indeterminate | 622 | 81.0% | 3.81 [3.16, 4.63] | 0.150 | above range: informative (pathogenic) |
| pathogenic supporting | 561 | 96.4% | 24.1 [16.3, 42.8] | 0.527 | above range: conservative |
| pathogenic moderate | 969 | 96.9% | 27.9 [20.3, 42.9] | 0.563 |  |
| pathogenic strong | 2,315 | 99.6% | 206 [123, 505] | 0.905 |  |

### native track, for comparison: the 1.0.x binning on the same halves, against its own claimed ranges

| Level | n | observed P/LP | observed LR [95% CI] | posterior at prior 0.0441 | flag |
|---|---|---|---|---|---|
| benign strong | 3,198 | 1.6% | 0.0145 [0.0106, 0.0185] | 0.001 |  |
| benign moderate | 281 | 12.5% | 0.127 [0.086, 0.176] | 0.006 |  |
| benign supporting | 798 | 49.7% | 0.883 [0.769, 1.01] | 0.039 | ABOVE range: over-claims (CI excludes range) |
| pathogenic supporting | 499 | 73.3% | 2.46 [2.03, 3.02] | 0.102 |  |
| pathogenic strong | 2,234 | 99.0% | 85.8 [59.7, 141] | 0.798 |  |
| pathogenic very strong | 1,442 | 98.5% | 60.4 [41.7, 103] | 0.736 | BELOW range: over-claims (CI excludes range) |
| unscored | 27 | 0.0% | 0 [0, 0] | 0.000 |  |

### Verdict

- routed: levels over-claiming 0 (0 variants) with the new calibration, against 1 (1,135) with 1.0.x.
- native: levels over-claiming 0 (0 variants) with the new calibration, against 2 (2,240) with 1.0.x.

## Why missense and other classes are calibrated apart

The same v2_native value carries very different evidence in the two branches:

| v2_native score | local LR, missense | local LR, other classes | local LR, pooled |
|---|---|---|---|
| 0.05 | 0.0531 | 0.169 | 0.174 |
| 0.10 | 0.0877 | 0.153 | 0.149 |
| 0.20 | 0.156 | 0.17 | 0.463 |
| 0.30 | 0.318 | 0.212 | 0.892 |
| 0.50 | 0.834 | 2.58 | 2.27 |
| 0.70 | 1.2 | 13.5 | 3.3 |
| 0.80 | 2.44 | 17.8 | 5.64 |
| 0.90 | 5.98 | 18 | 15 |
| 0.95 | 6.36 | 116 | 23.7 |
| 0.99 | 13.9 | 117 | 88.7 |

Pooling them would let the consequence-class branch (whose scores cluster at a few values per class) set the thresholds for missense variants. PP3/BP4 is in any case a missense criterion, so `acmg_code` is written for missense only.

## Missense only (the rows that carry `acmg_code`)

Same gene-disjoint procedure and split, 10,000 bootstrap resamples per fit; levels of the 2,011 labelled out-of-sample missense variants, with observed LR computed against the missense totals.

### routed track, missense

| Level | n | observed P/LP | observed LR [95% CI] | posterior at prior 0.0441 | flag |
|---|---|---|---|---|---|
| benign strong (BP4_Strong) | 239 | 2.1% | 0.00736 [0.00149, 0.0145] | 0.000 |  |
| benign moderate (BP4_Moderate) | 195 | 19.0% | 0.0806 [0.0557, 0.111] | 0.004 |  |
| benign supporting (BP4_Supporting) | 48 | 39.6% | 0.226 [0.12, 0.391] | 0.010 |  |
| indeterminate (no code) | 463 | 85.5% | 2.03 [1.62, 2.63] | 0.086 |  |
| pathogenic supporting (PP3_Supporting) | 404 | 96.0% | 8.35 [5.39, 15.1] | 0.278 | above range: conservative |
| pathogenic moderate (PP3_Moderate) | 662 | 98.3% | 20.4 [12.3, 44.9] | 0.485 |  |

Over-claiming levels: 0 (0 variants).

### native track, missense

| Level | n | observed P/LP | observed LR [95% CI] | posterior at prior 0.0441 | flag |
|---|---|---|---|---|---|
| benign moderate (BP4_Moderate) | 361 | 13.6% | 0.0541 [0.0402, 0.0705] | 0.002 |  |
| benign supporting (BP4_Supporting) | 102 | 44.1% | 0.272 [0.182, 0.398] | 0.012 |  |
| indeterminate (no code) | 606 | 81.0% | 1.47 [1.24, 1.77] | 0.064 |  |
| pathogenic supporting (PP3_Supporting) | 557 | 96.4% | 9.24 [6.34, 15.6] | 0.299 | above range: conservative |
| pathogenic moderate (PP3_Moderate) | 385 | 97.1% | 11.7 [7.06, 25.3] | 0.351 |  |

Over-claiming levels: 0 (0 variants).

