# Do the evidence bins mean what they say?

Calibration refitted on half the genes and read off the other half, so no gene calibrates the bin its own variants fall in. 8,479 labelled out-of-sample variants across 444 genes; the base rate is 52.8% pathogenic.

| Evidence bin | n | observed P/LP | vs base rate |
|---|---|---|---|
| benign strong | 3,392 | 1.8% | 0.03x |
| benign moderate | 319 | 16.9% | 0.32x |
| benign supporting | 315 | 31.1% | 0.59x |
| indeterminate | 247 | 76.1% | 1.44x |
| pathogenic supporting | 223 | 80.7% | 1.53x |
| pathogenic moderate | 734 | 95.2% | 1.80x |
| pathogenic strong | 2,089 | 99.5% | 1.88x |
| pathogenic very strong | 1,133 | 99.2% | 1.88x |
| unscored | 27 | 0.0% | 0.00x |

## Reading

The observed rate rises across the bins apart from one inversion of 0.3 percentage points (pathogenic strong to pathogenic very strong), which is saturation at the top of the scale rather than miscalibration: both bins are above 99% pathogenic and the ordering between them carries no practical information.
The pathogenic very strong bin is 99.2% pathogenic over 1,133 variants.
The pathogenic strong bin is 99.5% pathogenic over 2,089 variants.
The benign strong bin is 1.8% pathogenic over 3,392 variants.

These are read on genes that did not contribute to the calibration, so they are what a user should expect on a gene the model has not been tuned for.