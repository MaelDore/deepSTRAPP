# Plot multiple mapped phylogenies of trait/range evolution vs. diversification rates and regime shifts over time-steps

Plot two mapped phylogenies with evolutionary data with branches cut off
for each time-step. Main input = output of a deepSTRAPP run over time
using
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)).

- Left facet: plot the evolution of trait data/geographic ranges on the
  left time-calibrated phylogeny.

  - For continuous data: Each branch is colored according to the
    estimated value of the traits.

  - For categorical and biogeographic data: Each branch is colored
    according to the posterior probability of being in a given
    state/range. Colors for each state/range are overlaid using
    transparency to produce a single plot for all states/ranges.

- Right facet: plot the evolution of diversification rates and location
  of regime shifts estimated from a BAMM (Bayesian Analysis of
  Macroevolutionary Mixtures). Each branch is colored according to the
  estimated rates of speciation, extinction, or net diversification
  stored in an object of class `"bammdata"`. Rates can vary along time,
  thus colors evolve along individual branches.

This function is a wrapper of
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
used to plot mapped phylogenies for a unique given `focal_time`.

Returns one plot per focal time in `$time_steps`. If a PDF file path is
provided in `PDF_file_path`, the plots will be saved directly in a PDF
file, with one page per focal time in `$time_steps`.

## Usage

``` r
plot_traits_vs_rates_on_phylogeny_over_time(
  deepSTRAPP_outputs,
  color_scale = NULL,
  colors_per_levels = NULL,
  add_ACE_pies = TRUE,
  rate_type = "net_diversification",
  keep_initial_colorbreaks = FALSE,
  add_regime_shifts = TRUE,
  configuration_type = "MAP",
  sample_index = 1,
  regimes_fill = "grey",
  regimes_size = 1,
  regimes_pch = 21,
  regimes_border_col = "black",
  regimes_border_width = 1,
  ...,
  cex_pies = 0.5,
  adjust_size_to_prob = TRUE,
  display_plot = TRUE,
  PDF_file_path = NULL
)
```

## Arguments

- deepSTRAPP_outputs:

  List of elements generated with
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md),
  that summarize the results of a STRAPP run over multiple time-steps in
  the past.

- color_scale:

  Character vector. List of colors to use to build the color scale with
  [`grDevices::colorRampPalette()`](https://rdrr.io/r/grDevices/colorRamp.html)
  showing the evolution of a continuous trait. From lowest values to
  highest values. (For continuous trait data only)

- colors_per_levels:

  Named character string. To set the colors to use to map each
  state/range posterior probabilities. Names = states/ranges; values =
  colors. If `NULL` (default), the color scale provided in `densityMaps`
  will be used. (For categorical and biogeographic data only)

- add_ACE_pies:

  Logical. Whether to add pies of posterior probabilities of
  states/ranges at internal nodes on the mapped phylogeny. Default =
  `TRUE`.

- rate_type:

  A character string specifying the type of diversification rates to
  plot. Must be one of 'speciation', 'extinction' or
  'net_diversification' (default).

- keep_initial_colorbreaks:

  Logical. Whether to keep the same color breaks as used for the most
  recent focal time. Typically, the current time (t = 0). Default =
  `FALSE`.

- add_regime_shifts:

  Logical. Whether to add the location of regime shifts on the phylogeny
  (Step 2). Default is `TRUE`.

- configuration_type:

  A character string specifying how to select the location of regime
  shifts across posterior samples.

  - `configuration_type = "MAP"`: Use the average locations recorded in
    posterior samples with the Maximum A Posteriori probability (MAP)
    configuration. This regime shift configuration is the most frequent
    configuration among the posterior samples (See
    [`BAMMtools::getBestShiftConfiguration()`](https://rdrr.io/pkg/BAMMtools/man/getBestShiftConfiguration.html)).
    This is the default option.

  - `configuration_type = "MSC"`: Use the average locations recorded in
    posterior samples with the Maximum Shift Credibility (MSC)
    configuration. This regime shift configuration has the highest
    product of marginal probabilities across branches (See
    [`BAMMtools::maximumShiftCredibility()`](https://rdrr.io/pkg/BAMMtools/man/maximumShiftCredibility.html)).

  - `configuration_type = "index"`: Use the configuration of a unique
    posterior sample whose index is provided in `sample_index`.

- sample_index:

  Integer. Index of the posterior samples to use to plot the location of
  regime shifts. Used only if `configuration_type = index`. Default =
  `1`.

- regimes_fill:

  Character string. Set the color of the background of the symbols
  showing the location of regime shifts. Equivalent to the `bg` argument
  in
  [`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html).
  Default is `"grey"`.

- regimes_size:

  Numeric. Set the size of the symbols showing the location of regime
  shifts. Equivalent to the `cex` argument in
  [`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html).
  Default is `1`.

- regimes_pch:

  Integer. Set the shape of the symbols showing the location of regime
  shifts. Equivalent to the `pch` argument in
  [`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html).
  Default is `21`.

- regimes_border_col:

  Character string. Set the color of the border of the symbols showing
  the location of regime shifts. Equivalent to the `col` argument in
  [`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html).
  Default is `"black"`.

- regimes_border_width:

  Numeric. Set the width of the border of the symbols showing the
  location of regime shifts. Equivalent to the `lwd` argument in
  [`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html).
  Default is `1`.

- ...:

  Additional graphical arguments to pass down to
  [`phytools::plot.contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html),
  [`phytools::plot.simmap()`](https://rdrr.io/pkg/phytools/man/plotSimmap.html),
  [`plot_densityMaps_overlay()`](https://maeldore.github.io/deepSTRAPP/reference/plot_densityMaps_overlay.md),
  [`BAMMtools::plot.bammdata()`](https://rdrr.io/pkg/BAMMtools/man/plot.html),
  [`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html),
  and [`par()`](https://rdrr.io/r/graphics/par.html). If provided,
  `par.reset` is ignored with a message. See
  [`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
  for details.

- cex_pies:

  Numeric. To adjust the size of the ACE pies. Default = `0.5`. Provide
  the full argument name to avoid ambiguity with the `cex` graphical
  argument.

- adjust_size_to_prob:

  Logical. Whether to scale the size of the symbols showing the location
  of regime shifts according to the marginal shift probability of the
  shift happening on each location/branch. This will only work if there
  is an `$MSP_tree` element summarizing the marginal shift probabilities
  across branches in the `BAMM_object`. Default is `TRUE`. Provide the
  full argument name to avoid ambiguity with the `adj` graphical
  argument.

- display_plot:

  Logical. Whether to display the plot generated in the R console.
  Default is `TRUE`.

- PDF_file_path:

  Character string. If provided, the plot will be saved in a PDF file
  following the path provided here. The path must end with ".pdf".

## Value

If `display_plot = TRUE`, the function displays the successive plots
with two facets in the R console:

- (Left) A time-calibrated phylogeny displaying the evolution of
  trait/biogeographic data.

- (Right) A time-calibrated phylogeny displaying diversification rates
  and regime shifts. A plot per time-steps is generated and displayed
  successively in the console.

If `PDF_file_path` is provided, the plots will be exported in a unique
PDF file with one page per time-step.

## Details

The main input `deepSTRAPP_outputs` is the typical output of
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).
It provides information on results of a STRAPP run performed over
multiple time-steps, and can also encompass updated phylogenies with
mapped trait evolution and diversification rates and regime shifts if
appropriate arguments are set.

- `return_updated_Maps` must be set to `TRUE` so that the updated
  version of mapped phylogeny (`contMap`/`densityMaps`) for the given
  `$time_steps` are returned among the outputs under
  `$updated_Maps_over_time`. The update of the `contMap`/`densityMaps`
  consists of cutting off branches and mappings that are younger than
  the successive `$time_steps`. `contMap`/`densityMaps` must have been
  provided among the inputs to
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  in the first place.

- `return_updated_BAMM_objects` must be set to `TRUE` so that the list
  of `updated_BAMM_object` with phylogenies and mapped diversification
  rates cut-off at the successive `$time_steps` are returned among the
  outputs under `$updated_BAMM_objects_over_time`.

`$MAP_BAMM_object` and `$MSC_BAMM_object` elements are required in the
elements of `$updated_BAMM_objects_over_time` to plot regime shift
locations following the "MAP" or "MSC" `configuration_type`
respectively. A `$MSP_tree` element is required to scale the size of the
symbols showing the location of regime shifts according to marginal
shift probabilities. (If `adjust_size_to_prob = TRUE`).

For plotting a single time-step (i.e., `$focal_time`) at once, see
[`plot_traits_vs_rates_on_phylogeny_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_traits_vs_rates_on_phylogeny_for_focal_time.md).

## See also

[`plot_contMap()`](https://maeldore.github.io/deepSTRAPP/reference/plot_contMap.md)
[`phytools::plot.densityMap()`](https://rdrr.io/pkg/phytools/man/densityMap.html)
[`plot_densityMaps_overlay()`](https://maeldore.github.io/deepSTRAPP/reference/plot_densityMaps_overlay.md)
[`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)

Function in deepSTRAPP needed to produce the `deepSTRAPP_outputs` as
input:
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
Function in deepSTRAPP to plot a single time-step at once:
[`plot_traits_vs_rates_on_phylogeny_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_traits_vs_rates_on_phylogeny_for_focal_time.md)

## Author

Maël Doré

## Examples

``` r
if (deepSTRAPP::is_dev_version())
{
 # ----- Example 1: Continuous trait ----- #

 ## Load data

 # Load trait df
 data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
 # Load phylogeny with old calibration
 data(Ponerinae_tree_old_calib, package = "deepSTRAPP")

 # Load the BAMM_object summarizing 1000 posterior samples of BAMM
 data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 ## Prepare trait data

 # Extract continuous trait data as a named vector
 Ponerinae_cont_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cont_tip_data,
                                     nm = Ponerinae_trait_tip_data$Taxa)

 # Select a color scheme from lowest to highest values
 color_scale = c("darkgreen", "limegreen", "orange", "red")

  # Get Ancestral Character Estimates based on a Brownian Motion model
 # To obtain values at internal nodes
 Ponerinae_ACE <- phytools::fastAnc(tree = Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data)

  # (May take several minutes to run)
 # Run a Stochastic Mapping based on a Brownian Motion model
 # to interpolate values along branches and obtain a "contMap" object
 Ponerinae_contMap <- phytools::contMap(Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data,
                                        res = 100, # Number of time steps
                                        plot = FALSE)
 # Plot contMap = stochastic mapping of continuous trait
 plot_contMap(contMap = Ponerinae_contMap,
              color_scale = color_scale)

 ## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
 # nb_time_steps <- 5
 time_step_duration <- 5
 time_range <- c(0, 40)

 ## Run deepSTRAPP on net diversification rates
 Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- run_deepSTRAPP_over_time(
    contMap = Ponerinae_contMap,
    ace = Ponerinae_ACE,
    tip_data = Ponerinae_cont_tip_data,
    trait_data_type = "continuous",
    BAMM_object = Ponerinae_BAMM_object_old_calib,
    # nb_time_steps = nb_time_steps,
    time_range = time_range,
    time_step_duration = time_step_duration,
    uncertainty_strategy = "rates_only",
    return_perm_data = TRUE,
    extract_trait_data_melted_df = TRUE,
    extract_diversification_data_melted_df = TRUE,
    return_STRAPP_results = TRUE,
    return_updated_Maps = TRUE,
    return_updated_BAMM_objects = TRUE,
    verbose = TRUE,
    verbose_extended = TRUE) 

 ## Load directly trait data output
 Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- readRDS(system.file("extdata",
     "Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", package = "deepSTRAPP"))
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 ## Plot updated contMap vs. updated diversification rates
 plot_traits_vs_rates_on_phylogeny_over_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
    # To use the same color breaks as for t = 0 My in all BAMM plots
    keep_initial_colorbreaks = TRUE,
    color_scale = c("limegreen", "orange", "red"), # Adjust color scale on contMap
    legend = TRUE, labels = TRUE, # Show legend and label on BAMM plot
    cex = 0.7) # Adjust label size on contMap
    # PDF_file_path = "Updated_maps_cont_old_calib_0_40My.pdf")

 # ----- Example 2: Biogeographic data ----- #

 ## Load data

 # Load phylogeny
 data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
 # Load trait df
 data(Ponerinae_binary_range_table, package = "deepSTRAPP")

 # Load the BAMM_object summarizing 1000 posterior samples of BAMM
 data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 ## Prepare range data for Old World vs. New World

 # No overlap in ranges
 table(Ponerinae_binary_range_table$Old_World, Ponerinae_binary_range_table$New_World)

 Ponerinae_NO_tip_data <- stats::setNames(object = Ponerinae_binary_range_table$Old_World,
                                      nm = Ponerinae_binary_range_table$Taxa)
 Ponerinae_NO_tip_data <- as.character(Ponerinae_NO_tip_data)
 Ponerinae_NO_tip_data[Ponerinae_NO_tip_data == "TRUE"] <- "O" # O = Old World
 Ponerinae_NO_tip_data[Ponerinae_NO_tip_data == "FALSE"] <- "N" # N = New World
 names(Ponerinae_NO_tip_data) <- Ponerinae_binary_range_table$Taxa
 table(Ponerinae_NO_tip_data)

 colors_per_levels <- c("mediumpurple2", "peachpuff2")
 names(colors_per_levels) <- c("N", "O")

  # (May take several minutes to run)

 # Load ape and BioGeoBEARS to use models internally
 library(ape)
 library(BioGeoBEARS)

 ## Run evolutionary models
 Ponerinae_biogeo_data <- prepare_trait_data(
    tip_data = Ponerinae_NO_tip_data,
    trait_data_type = "biogeographic",
    phylo = Ponerinae_tree_old_calib,
    evolutionary_models = "DEC+J", # Default = "DEC" for biogeographic
    BioGeoBEARS_directory_path = tempdir(), # Ex: "./BioGeoBEARS_directory/"
    keep_BioGeoBEARS_files = FALSE,
    prefix_for_files = "Ponerinae_old_calib",
    max_range_size = 2,
    split_multi_area_ranges = TRUE, # Set to TRUE to use only unique areas "A" and "B"
    nb_simulations = 100, # Reduce to save time (Default = '1000')
    colors_per_levels = colors_per_levels,
    return_model_selection_df = TRUE,
    verbose = TRUE) 

 # Load directly output
 data(Ponerinae_biogeo_data_old_calib, package = "deepSTRAPP")

 ## Set for time steps of 5 My. Will generate deepSTRAPP workflows from 0 to 40 Mya.
 time_range <- c(0, 40)
 time_step_duration <- 5

  # (May take several minutes to run)
 ## Run deepSTRAPP on net diversification rates for time-steps = 0 to 40 Mya.
 Ponerinae_deepSTRAPP_biogeo_old_calib_0_40 <- run_deepSTRAPP_over_time(
    densityMaps = Ponerinae_biogeo_data_old_calib$densityMaps,
    ace = Ponerinae_biogeo_data_old_calib$ace,
    tip_data = Ponerinae_NO_tip_data,
    nb_simulations = 100,
    trait_data_type = "biogeographic",
    BAMM_object = Ponerinae_BAMM_object_old_calib,
    time_range = time_range,
    time_step_duration = time_step_duration,
    uncertainty_strategy = "paired",
    seed = 1234, # Set seed for reproducibility
    alpha = 0.10, # Select a generous level of significance for the sake of the example
    rate_type = "net_diversification",
    return_perm_data = TRUE,
    extract_trait_data_melted_df = TRUE,
    extract_diversification_data_melted_df = TRUE,
    return_STRAPP_results = TRUE,
    return_updated_Maps = TRUE,
    return_updated_BAMM_objects = TRUE,
    verbose = TRUE,
    verbose_extended = TRUE) 

 ## Load directly output
 Ponerinae_deepSTRAPP_biogeo_old_calib_0_40 <- readRDS(system.file("extdata",
     "Ponerinae_deepSTRAPP_biogeo_old_calib_0_40.rds", package = "deepSTRAPP"))
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 ## Explore output
 str(Ponerinae_deepSTRAPP_biogeo_old_calib_0_40, max.level =  1)

 ## Plot updated contMap vs. updated diversification rates
 plot_traits_vs_rates_on_phylogeny_over_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_biogeo_old_calib_0_40,
    # Adjust colors on densityMaps
    colors_per_levels = c("N" = "dodgerblue2", "O" = "orange"),
    legend = TRUE, labels = TRUE, # Show legend and label on BAMM plot
    cex_pies = 0.2, # Adjust size of ACE pies on densityMaps
    cex = 0.7) # Adjust label size on contMap
    # PDF_file_path = "Updated_maps_biogeo_old_calib_0_40My.pdf")
}
```
