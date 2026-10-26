# Evidence distribution, 1.0.2 -> 1.1.0

## transporter_variants_routed

| Level | 1.0.2 (`evidence_v1`) | 1.1.0 (`evidence`) |
|---|---|---|
| benign strong | 168,714 | 98,524 |
| benign moderate | 68,915 | 91,966 |
| benign supporting | 0 | 41,391 |
| indeterminate | 27,412 | 53,006 |
| pathogenic supporting | 11,986 | 4,735 |
| pathogenic moderate | 25,312 | 32,836 |
| pathogenic strong | 36,061 | 26,598 |
| pathogenic very strong | 10,656 | 0 |

`acmg_code`: (blank) 222,755, BP4_Moderate 45,336, PP3_Moderate 28,081, BP4_Supporting 25,903, BP4_Strong 22,246, PP3_Supporting 4,735

Rows: 1.0.2 level; columns: 1.1.0 level.

| 1.0.2 \ 1.1.0 | benign strong | benign moderate | benign supporting | indeterminate | pathogenic supporting | pathogenic moderate | pathogenic strong | pathogenic very strong |
|---|---|---|---|---|---|---|---|---|
| benign strong | 98,524 | 70,190 | 0 | 0 | 0 | 0 | 0 | 0 |
| benign moderate | 0 | 21,776 | 41,391 | 5,748 | 0 | 0 | 0 | 0 |
| indeterminate | 0 | 0 | 0 | 27,412 | 0 | 0 | 0 | 0 |
| pathogenic supporting | 0 | 0 | 0 | 11,155 | 831 | 0 | 0 | 0 |
| pathogenic moderate | 0 | 0 | 0 | 8,691 | 3,904 | 11,957 | 760 | 0 |
| pathogenic strong | 0 | 0 | 0 | 0 | 0 | 20,868 | 15,193 | 0 |
| pathogenic very strong | 0 | 0 | 0 | 0 | 0 | 11 | 10,645 | 0 |

## transporter_variants_native

| Level | 1.0.2 (`evidence_v1`) | 1.1.0 (`evidence`) |
|---|---|---|
| benign strong | 117,270 | 76,278 |
| benign moderate | 121,618 | 113,666 |
| benign supporting | 0 | 26,574 |
| indeterminate | 0 | 74,459 |
| pathogenic supporting | 0 | 14,996 |
| pathogenic moderate | 76,385 | 16,485 |
| pathogenic strong | 22,310 | 26,598 |
| pathogenic very strong | 11,473 | 0 |

`acmg_code`: (blank) 244,208, BP4_Moderate 67,036, PP3_Supporting 14,996, PP3_Moderate 11,730, BP4_Supporting 11,086

Rows: 1.0.2 level; columns: 1.1.0 level.

| 1.0.2 \ 1.1.0 | benign strong | benign moderate | benign supporting | indeterminate | pathogenic supporting | pathogenic moderate | pathogenic strong | pathogenic very strong |
|---|---|---|---|---|---|---|---|---|
| benign strong | 76,278 | 40,992 | 0 | 0 | 0 | 0 | 0 | 0 |
| benign moderate | 0 | 72,674 | 26,574 | 22,370 | 0 | 0 | 0 | 0 |
| pathogenic moderate | 0 | 0 | 0 | 52,089 | 14,996 | 8,540 | 760 | 0 |
| pathogenic strong | 0 | 0 | 0 | 0 | 0 | 7,117 | 15,193 | 0 |
| pathogenic very strong | 0 | 0 | 0 | 0 | 0 | 828 | 10,645 | 0 |

## Mechanism coverage

- transporter_variants_routed: mechanism refit reproduces 5,430 kept values (max |diff| 4.2e-07); filled 0 newly confident rows, blanked 12,303 no longer confident; now 28,081
- transporter_variants_native: mechanism refit reproduces 3,868 kept values (max |diff| 3.7e-07); filled 0 newly confident rows, blanked 66,793 no longer confident; now 11,730
