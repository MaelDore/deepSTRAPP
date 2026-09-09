# Package index

## Prepare data for deepSTRAPP

All-in-one functions to prepare trait and diversification data to use as
inputs for a deepSTRAPP run.

- [`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md)
  : Map trait evolution on a time-calibrated phylogeny
- [`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md)
  : Run a full BAMM (Bayesian Analysis of Macroevolutionary Mixtures)
  workflow

## Run deepSTRAPP worflow

Core functions to run deepSTRAPP.

- [`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
  : Run deepSTRAPP to test for a relationship between diversification
  rates and trait data at a given focal time
- [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  : Run deepSTRAPP to test for a relationship between diversification
  rates and trait data over multiple time steps
- [`compute_STRAPP_test_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/compute_STRAPP_test_for_focal_time.md)
  : Compute STRAPP to test for a relationship between diversification
  rates and trait data

## Plot deepSTRAPP outputs

Plot results from deepSTRAPP runs.

- [`plot_STRAPP_pvalues_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_STRAPP_pvalues_over_time.md)
  : Plot evolution of p-values of STRAPP tests over time
- [`plot_histogram_STRAPP_test_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_histogram_STRAPP_test_for_focal_time.md)
  : Plot histogram of STRAPP test statistics to assess results
- [`plot_histograms_STRAPP_tests_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_histograms_STRAPP_tests_over_time.md)
  : Plot multiple histograms of STRAPP test statistics over time-steps
- [`plot_rates_through_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_rates_through_time.md)
  : Plot evolution of diversification rates in relation to trait values
  over time
- [`plot_rates_vs_trait_data_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_rates_vs_trait_data_for_focal_time.md)
  : Plot rates vs. trait data for a given focal time
- [`plot_rates_vs_trait_data_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_rates_vs_trait_data_over_time.md)
  : Plot mean rates vs. trait data over time-steps

## Plot mapped phylogenies

Plot phylogenies with mapped trait or diversification rates evolution on
branches.

- [`plot_contMap()`](https://maeldore.github.io/deepSTRAPP/reference/plot_contMap.md)
  : Plot continuous trait evolution on the tree
- [`plot_densityMaps_overlay()`](https://maeldore.github.io/deepSTRAPP/reference/plot_densityMaps_overlay.md)
  : Plot posterior probabilities of states/ranges on phylogeny from
  densityMaps
- [`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
  : Plot diversification rates and regime shifts from BAMM on phylogeny
- [`plot_traits_vs_rates_on_phylogeny_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_traits_vs_rates_on_phylogeny_for_focal_time.md)
  : Plot trait/range evolution vs. diversification rates and regime
  shifts on phylogeny
- [`plot_traits_vs_rates_on_phylogeny_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_traits_vs_rates_on_phylogeny_over_time.md)
  : Plot multiple mapped phylogenies of trait/range evolution vs.
  diversification rates and regime shifts over time-steps

## Import results from external analyses

Format results from external analyses as inputs for a deepSTRAPP run.

- [`convert_contsimmap_to_contMaps()`](https://maeldore.github.io/deepSTRAPP/reference/convert_contsimmap_to_contMaps.md)
  : Convert a contsimmap object into a list of contMaps
- [`aggregate_contMaps()`](https://maeldore.github.io/deepSTRAPP/reference/aggregate_contMaps.md)
  : Aggregate a list of contMaps into a unique mean/median contMap
- [`convert_BSM_to_simmap()`](https://maeldore.github.io/deepSTRAPP/reference/convert_BSM_to_simmap.md)
  [`convert_BSMs_to_simmaps()`](https://maeldore.github.io/deepSTRAPP/reference/convert_BSM_to_simmap.md)
  : Convert Biogeographic Stochastic Map (BSM) to phytools SIMMAP
  stochastic map (SM) format
- [`convert_simmaps_to_densityMaps()`](https://maeldore.github.io/deepSTRAPP/reference/convert_simmaps_to_densityMaps.md)
  : Convert simmaps into densityMaps suitable for deepSTRAPP
- [`build_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/build_BAMM_object.md)
  : Build a BAMM object for a deepSTRAPP run
- [`subset_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/subset_BAMM_object.md)
  : Subset a BAMM object before a deepSTRAPP run
- [`prune_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/prune_BAMM_object.md)
  : Prune a BAMM object to a subset of tips

## Model selection

Select best models based on AICc.

- [`select_best_trait_model_from_geiger()`](https://maeldore.github.io/deepSTRAPP/reference/select_best_trait_model_from_geiger.md)
  : Compare trait evolutionary model fits with AICc and Akaike's weights
- [`select_best_model_from_BioGeoBEARS()`](https://maeldore.github.io/deepSTRAPP/reference/select_best_model_from_BioGeoBEARS.md)
  : Compare model fits with AICc and Akaike's weights

## Extract data from mapped phylogenies

Extract trait data from contMap(s), densityMaps, simmaps.

- [`extract_most_likely_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_most_likely_trait_values_for_focal_time.md)
  : Extract most likely trait data mapped on a phylogeny at a given time
  in the past
- [`extract_all_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_all_trait_values_for_focal_time.md)
  : Extract all trait data from stochastic maps at a given time in the
  past
- [`extract_trait_data_melted_df_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_trait_data_melted_df_for_focal_time.md)
  : Extract trait data from a trait_data object

## Extract data from BAMM objects

Extract diversification rates and regimes from BAMM objects.

- [`extract_diversification_data_melted_df_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_diversification_data_melted_df_for_focal_time.md)
  : Extract diversification data from a BAMM_object
- [`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md)
  : Update diversification rates/regimes mapped on a phylogeny up to a
  given time in the past

## Cut phylogenies

Cut different types of mapped phylogenies to a given focal-time.

- [`cut_phylo_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_phylo_for_focal_time.md)
  : Cut the phylogeny for a given time in the past
- [`cut_contMap_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_contMap_for_focal_time.md)
  : Cut the phylogeny and continuous trait mapping for a given focal
  time in the past
- [`cut_contMaps_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_contMaps_for_focal_time.md)
  : Cut the phylogenies and continuous trait mappings of a list of
  contMaps for a given focal time in the past
- [`cut_densityMap_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_densityMap_for_focal_time.md)
  : Cut the phylogeny and posterior probability mapping of a categorical
  trait for a given focal time in the past
- [`cut_densityMaps_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_densityMaps_for_focal_time.md)
  : Cut phylogenies and posterior probability mapping of each state for
  a given focal time in the past
- [`cut_simmap_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_simmap_for_focal_time.md)
  : Cut the phylogeny and categorical trait/range mapping for a given
  focal time in the past
- [`cut_simmaps_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_simmaps_for_focal_time.md)
  : Cut the phylogenies and categorical trait/range mappings of a list
  of simmaps for a given focal time in the past

## Datasets

Datasets used to run examples and display expected results.

### Phylogenies

- [`Ponerinae_tree`](https://maeldore.github.io/deepSTRAPP/reference/Ponerinae_tree.md)
  : Dataset providing the extensive time-calibrated phylogeny of extant
  ponerine ants
- [`Ponerinae_tree_old_calib`](https://maeldore.github.io/deepSTRAPP/reference/Ponerinae_tree_old_calib.md)
  : Dataset providing the extensive time-calibrated phylogeny of extant
  ponerine ants using an old calibration for illustrative purposes

### Trait datasets

- [`mammals`](https://maeldore.github.io/deepSTRAPP/reference/mammals.md)
  : Phylogeny and body mass data for extant and extinct mammal
  families/genera from Slater, 2013
- [`Ponerinae_trait_tip_data`](https://maeldore.github.io/deepSTRAPP/reference/Ponerinae_trait_tip_data.md)
  : Dataset providing fake trait data for extant ponerine ants for
  illustrative purposes
- [`Ponerinae_binary_range_table`](https://maeldore.github.io/deepSTRAPP/reference/Ponerinae_binary_range_table.md)
  : Dataset providing biogeographic range data for extant ponerine ants

### Trait evolution data

- [`Ponerinae_trait_cont_tip_data_10My`](https://maeldore.github.io/deepSTRAPP/reference/Ponerinae_trait_cont_tip_data_10My.md)
  : Data summarizing the evolution of a fake continuous trait in
  Ponerinae ants extracted for 10 Mya
- [`eel_cat_3lvl_data`](https://maeldore.github.io/deepSTRAPP/reference/eel_cat_3lvl_data.md)
  : Data summarizing the evolution of feeding habits in eels using a
  3-level factor as categorical trait
- [`Ponerinae_cat_2lvl_data_old_calib`](https://maeldore.github.io/deepSTRAPP/reference/Ponerinae_cat_2lvl_data_old_calib.md)
  : Data summarizing the evolution of fake size data in Ponerinae ants
  using a 2-level factor as categorical trait
- [`Ponerinae_cat_3lvl_data_old_calib`](https://maeldore.github.io/deepSTRAPP/reference/Ponerinae_cat_3lvl_data_old_calib.md)
  : Data summarizing the evolution of fake habitat data in Ponerinae
  ants using a 3-level factor as categorical trait
- [`eel_biogeo_data`](https://maeldore.github.io/deepSTRAPP/reference/eel_biogeo_data.md)
  : Data summarizing the evolution of geographic ranges in eels
- [`Ponerinae_biogeo_data_old_calib`](https://maeldore.github.io/deepSTRAPP/reference/Ponerinae_biogeo_data_old_calib.md)
  : Data summarizing the evolution of geographic ranges in Ponerinae
  ants using an old ill-calibrated phylogeny for illustrative purposes

### BAMM data

- [`BAMM_template_diversification`](https://maeldore.github.io/deepSTRAPP/reference/BAMM_template_diversification.md)
  : Template file for BAMM diversification analyses
- [`whale_BAMM_object`](https://maeldore.github.io/deepSTRAPP/reference/whale_BAMM_object.md)
  : Dataset summarizing 1000 posterior samples of BAMM for extant whales
