# Changelog

## deepSTRAPP 1.1.1

- [`convert_BSM_to_simmap()`](https://maeldore.github.io/deepSTRAPP/reference/convert_BSM_to_simmap.md)
  and
  [`convert_BSMs_to_simmaps()`](https://maeldore.github.io/deepSTRAPP/reference/convert_BSM_to_simmap.md)
  now repair paths to phylo and tip ranges objects used in the example
  `eel_biogeo_data` dataset.

## deepSTRAPP 1.1.0

CRAN release: 2026-09-28

- Handle uncertainty in trait estimates. See the ‘uncertainty_strategy’
  argument in run_deepSTRAPP\_\*() functions.
- Add functions to load results from external BAMM analyses. See the
  dedicated vignette/tutorial “import_external_analyses”.

### Accounting for uncertainty in ancestral trait estimates

- STRAPP tests can now be run across a posterior sample of trait
  histories rather than a single reconstruction. The
  `uncertainty_strategy` argument of
  [`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md),
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  and
  [`compute_STRAPP_test_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/compute_STRAPP_test_for_focal_time.md)
  selects how trait and rate uncertainty are combined: `"rates_only"`
  (BAMM posterior samples with a single trait reconstruction),
  `"paired"` (each stochastic map paired with one BAMM sample), or
  `"full"` (all stochastic maps crossed with all BAMM samples).
- Stochastic maps can be supplied directly through the new `contMaps`
  and `simmaps` arguments, or simulated from posterior probabilities
  with `nb_simulations`.
- `trait_maps_vs_BAMM_samples_list` records, and lets you impose, which
  stochastic map is paired with which BAMM sample, so that the same
  pairing is used at every time step of a trajectory.
- See the new vignette `handle_uncertainty`.

### Importing results from external analyses

- New functions to bring analyses run outside deepSTRAPP into the
  workflow:
  [`build_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/build_BAMM_object.md),
  [`subset_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/subset_BAMM_object.md)
  and
  [`prune_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/prune_BAMM_object.md)
  for BAMM output;
  [`convert_BSM_to_simmap()`](https://maeldore.github.io/deepSTRAPP/reference/convert_BSM_to_simmap.md)
  and
  [`convert_BSMs_to_simmaps()`](https://maeldore.github.io/deepSTRAPP/reference/convert_BSM_to_simmap.md)
  for BioGeoBEARS biogeographic stochastic maps;
  [`convert_contsimmap_to_contMaps()`](https://maeldore.github.io/deepSTRAPP/reference/convert_contsimmap_to_contMaps.md)
  for `contsimmap` output; and
  [`convert_simmaps_to_densityMaps()`](https://maeldore.github.io/deepSTRAPP/reference/convert_simmaps_to_densityMaps.md)
  to summarise stochastic maps as posterior densities.
- [`aggregate_contMaps()`](https://maeldore.github.io/deepSTRAPP/reference/aggregate_contMaps.md)
  summarises a list of continuous stochastic maps into a single mean or
  median `contMap`.
- See the new vignette `import_external_analyses`.

### Consistency of the statistical method across time steps

- [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  now selects the statistical method once, before any test is run, from
  all states/ranges described in the complete trait mapping, and applies
  it at every time step. The method previously depended on the
  states/ranges still present at each time step, so a trajectory could
  mix Kruskal-Wallis and Mann-Whitney U p-values on a single curve when
  a state/range was absent from the deeper time steps. This is not the
  case anymore, and all p-values across time-steps are prodcued by the
  same type of test.
- The states/ranges actually observed are now reported:
  `$states_observed` and `$nb_states_observed` per time step, and
  `$states_observed_overall`, `$states_observed_per_time_steps` and
  `$nb_states_observed_per_time_steps` in the output of
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).

### Other additions

- New
  [`extract_all_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_all_trait_values_for_focal_time.md)
  to extract trait data from every stochastic map at a given time, and
  [`extract_trait_data_melted_df_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_trait_data_melted_df_for_focal_time.md)
  to return it in long dataframe format.
- New
  [`cut_contMaps_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_contMaps_for_focal_time.md),
  [`cut_simmap_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_simmap_for_focal_time.md)
  and
  [`cut_simmaps_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_simmaps_for_focal_time.md)
  to cut lists of stochastic maps at a focal time.
- [`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
  and
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  gain `return_updated_Maps` to return the mappings cut at each focal
  time, and
  [`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
  gains `extract_trait_data_melted_df` to return the underlying trait
  data in long dataframe format.

## deepSTRAPP 1.0.0

CRAN release: 2026-01-19

- First release on CRAN
