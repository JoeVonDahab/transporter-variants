# Transporter Variants

A static website for browsing calibrated pathogenicity predictions for all 349,056 ClinVar variants in 783 human transporter, ion-channel, ABC-transporter and ATPase genes.

**Site:** https://joevondahab.github.io/transporter-variants/

Each variant carries an ACMG evidence strength derived from a likelihood ratio, validated on gene-disjoint halves of the held-out labelled set. Gene pages show a sortable, filterable variant table and a protein map: the chain drawn across the membrane from UniProt topology, plus the AlphaFold structure, both coloured by the regional average of the predictions.

## Two tracks, two licences

| Track | Contents | Licence |
|---|---|---|
| native | Our model alone, every variant, no third-party predictions | CC BY 4.0 |
| routed | Our model with AlphaMissense where AlphaMissense covers the variant, our model elsewhere | CC BY-NC-SA 4.0, inherited from AlphaMissense |

The tracks are kept separate throughout: separate files in `docs/data/native/` and `docs/data/routed/`, separate downloads in `docs/downloads/`. Do not merge them.

## Using the data correctly

- Read `evidence`, not `score`. A 0.5 cut-off over-calls, because the training set was 45% pathogenic.
- Rows with `in_training` TRUE were seen during model fitting. Exclude them from any accuracy estimate.
- On the routed track, read `v2_plus` only where `has_alphamissense` is TRUE.
- `mechanism_pred` is exploratory: a hypothesis generator, not a per-variant call.

## Layout

```
docs/                 the site, served by GitHub Pages
  index.html          landing page
  gene.html           gene page (?g=SYMBOL, optional &t=routed, &q=, &pos=)
  assets/             CSS, JavaScript, vendored 3Dmol.js 2.5.4 (BSD-3-Clause)
  data/index.json     gene index, landing-page numbers, shared dictionaries
  data/native/        one JSON file per gene, native track
  data/routed/        one JSON file per gene, routed track
  data/protein/       one JSON file per gene: sequence and UniProt features
  data/lookup/        rsID and ClinVar VariationID to gene, in 100 shards each
  downloads/          the release files, unchanged, with their licences and checksums
tools/build_data.py   rebuilds docs/data and docs/downloads from the release
```

There is no build step for the site itself and no server: everything is static files.

## Rebuilding the data

```bash
python -I tools/build_data.py --release ~/projects/TripleVar/outputs/transporters/release
```

The script reads the two release files, fetches UniProt features once (cached in `tools/cache/`, which is not committed), checks the release checksums, and writes `docs/data/`. AlphaFold models are not stored here; gene pages fetch them from AlphaFold DB in the browser.

## Licences

- Site code (HTML, CSS, JavaScript, `tools/`): MIT, see `LICENSE`.
- Native track and the native-track data files under `docs/data/`: CC BY 4.0.
- Routed track and `docs/data/routed/`: CC BY-NC-SA 4.0.
- UniProtKB annotations and AlphaFold DB models: CC BY 4.0, from their providers.

## Citation

Abo-Dahab Y, Chun J, Giacomini KM. *Calibrated pathogenicity predictions for 783 human transporter and ion-channel genes.* Version 1.0.0, 2026. See `docs/downloads/CITATION.cff`. If you use the routed track, also cite AlphaMissense (Cheng et al., *Science* 2023).

Predictions are for research use and are not a clinical classification.
