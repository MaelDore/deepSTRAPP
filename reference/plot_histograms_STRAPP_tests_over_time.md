# Plot multiple histograms of STRAPP test statistics over time-steps

Plot a histogram of the distribution of the test statistics obtained
from a deepSTRAPP workflow carried out for each focal time in
`$time_steps`. Main input = output of a deepSTRAPP run over time using
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)).

Returns one histogram for overall tests for each focal time in
`$time_steps`. If `plot_posthoc_tests = TRUE`, it will return one
faceted plot with a histogram per post hoc test for each focal time in
`$time_steps`.

If a PDF file path is provided in `PDF_file_path`, the plots will be
saved directly in a PDF file, with one page per focal time in
`$time_steps`.

## Usage

``` r
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs,
  display_plots = TRUE,
  plot_posthoc_tests = FALSE,
  PDF_file_path = NULL
)
```

## Arguments

- deepSTRAPP_outputs:

  List of elements generated with
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md),
  that summarize the results of multiple deepSTRAPP runs across
  `$time_steps`. It needs to include the `$STRAPP_results_over_time`
  element with `$perm_data_df` obtained when setting both
  `return_STRAPP_results = TRUE` and `return_perm_data = TRUE`.

- display_plots:

  Logical. Whether to display the histograms generated in the R console.
  Default is `TRUE`.

- plot_posthoc_tests:

  Logical. For multinomial data only. Whether to plot the histograms for
  the overall Kruskal-Wallis test across all states
  (`plot_posthoc_tests = FALSE`), or plot the histograms for all the
  pairwise post hoc Dunn's tests across pairs of states
  (`plot_posthoc_tests = TRUE`). Time-steps at which the data does not
  yield more than two states/ranges will show a warning and generate no
  plot. Default is `FALSE`.

- PDF_file_path:

  Character string. If provided, the plots will be saved in a unique PDF
  file following the path provided here. The path must end with ".pdf".
  Each page of the PDF corresponds to a focal time in `$time_steps`.

## Value

By default, the function returns a list of sub-lists of classes `gg` and
`ggplot` ordered as in `$time_steps`. Each sub-list corresponds to a
ggplot for a given `focal_time` `i` that can be displayed on the console
with `print(output[[i]])`. If `display_plots = TRUE`, the histograms are
being displayed on the console one by one while generated.

If using multinomial data and set `plot_posthoc_tests = TRUE`, the
function will return a list of sub-lists of objects ordered as in
`$time_steps`. Each sub-list is a list of the ggplots associated with
pairwise post hoc tests carried out for a given `focal_time`. For a
given `focal_time` i, to plot each histogram j individually, use
`print(output_list[[i]][[j]])`. To plot all histograms of a given
`focal_time` `i` at once in a multifaceted plot, as displayed
sequentially on the console if `display_plots = TRUE`, use
`cowplot::plot_grid(plotlist = output_list[[i]])`.

Each plot also displays summary statistics for the STRAPP test
associated with the data displayed.

- QX% = The quantile of null statistic distribution at the significance
  threshold used to define test significance. This is the value found on
  the red dashed line. The test will be considered significant (i.e.,
  the null hypothesis is rejected) if this value is higher than zero
  (the black dashed line).

- P-value = The p-value of the STRAPP test which corresponds to the
  proportion of cases in which the statistic was lower than expected
  under the null hypothesis (i.e., the proportion of the histogram found
  below / on the left-side of the black dashed line).

- N = The number of test statistics in the null distribution.

If a `PDF_file_path` is provided, the function will also generate a PDF
file of the plots with one page per `$time_steps`. For post hoc tests,
this will save the multifaceted plots.

## Details

The main input `deepSTRAPP_outputs` is the typical output of
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).
It provides information on results of STRAPP tests performed over
multiple time-steps.

Histograms are built based on the distribution of the test statistics.
Such distributions are recorded in the outputs of a deepSTRAPP run
carried out with
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
when `return_STRAPP_results = TRUE` AND `return_perm_data = TRUE`. The
`$STRAPP_results_over_time` objects provided within the input are lists
that must contain a `$perm_data_df` element that summarizes test
statistics computed across posterior samples.

For multinomial data (categorical or biogeographic data with more than 2
states), it is possible to plot the histograms of post hoc pairwise
tests. Set `plot_posthoc_tests = TRUE` to generate histograms for all
the pairwise post hoc Dunn's tests across pairs of states. To achieve
this, the `$STRAPP_results_over_time` objects must contain a
`$posthoc_pairwise_tests$perm_data_array` element that summarizes test
statistics computed across posterior samples for all pairwise post hoc
tests. This is obtained from
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
when setting `return_STRAPP_results = TRUE` to return the STRAPP
results, `posthoc_pairwise_tests = TRUE` to carry out post hoc tests,
and `return_perm_data = TRUE` to record distributions of test
statistics. Time-steps for which the data do not yield more than two
states/ranges will show a warning and generate no plot.

## See also

Associated functions in deepSTRAPP:
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
[`plot_histogram_STRAPP_test_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_histogram_STRAPP_test_for_focal_time.md)

## Author

Maël Doré

## Examples

``` r
if (deepSTRAPP::is_dev_version())
{
 # ----- Example 1: Continuous trait ----- #

 # Load fake trait df
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
 Ponerinae_contMap <- phytools::contMap(tree = Ponerinae_tree_old_calib,
                                        x = Ponerinae_cont_tip_data,
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

 ## Plot histograms of STRAPP overall test results
 # Tests are Spearman's rank correlation tests

 # Plot all histograms
 histogram_ggplots <- plot_histograms_STRAPP_tests_over_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
    display_plots = TRUE,
    # PDF_file_path = "./plot_STRAPP_histogram_overall_test.pdf",
    plot_posthoc_tests = FALSE)

 # Print histogram for time step 1 = 0 My
 print(histogram_ggplots[[1]])
 # Adjust aesthetics of plot for time step 1 a posteriori
 histogram_ggplot_adj <- histogram_ggplots[[1]] +
    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
 print(histogram_ggplot_adj)

 # ----- Example 2: Categorical data ----- #

 ## Load data

 # Load trait df
 data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
 # Load phylogeny
 data(Ponerinae_tree_old_calib, package = "deepSTRAPP")

 # Load the BAMM_object summarizing 1000 posterior samples of BAMM
 data(Ponerinae_BAMM_object_old_calib, package = "deepSTRAPP")
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 ## Prepare trait data

 # Extract categorical data with 3-levels
 Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
                                         nm = Ponerinae_trait_tip_data$Taxa)
 table(Ponerinae_cat_3lvl_tip_data)

 # Select color scheme for states
 colors_per_states <- c("forestgreen", "sienna", "goldenrod")
 names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")

  # (May take several minutes to run)
 ## Produce densityMaps using stochastic character mapping based on an ARD Mk model
 Ponerinae_cat_3lvl_data_old_calib <- prepare_trait_data(
    tip_data = Ponerinae_cat_3lvl_tip_data,
    phylo = Ponerinae_tree_old_calib,
    trait_data_type = "categorical",
    colors_per_levels = colors_per_states,
    evolutionary_models = "ARD",
    nb_simulations = 100,
    return_best_model_fit = TRUE,
    return_model_selection_df = TRUE,
    plot_map = FALSE) 

 # Load directly trait data output
 data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")

 ## Set for time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 40 Mya.
 # nb_time_steps <- 5
 time_step_duration <- 5
 time_range <- c(0, 40)

  # (May take several minutes to run)
 ## Run deepSTRAPP on net diversification rates across time-steps.
 Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- run_deepSTRAPP_over_time(
    densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
    ace = Ponerinae_cat_3lvl_data_old_calib$ace,
    tip_data = Ponerinae_cat_3lvl_tip_data,
    nb_simulations = 100,
    trait_data_type = "categorical",
    BAMM_object = Ponerinae_BAMM_object_old_calib,
    # nb_time_steps = nb_time_steps,
    time_range = time_range,
    time_step_duration = time_step_duration,
    rate_type = "net_diversification",
    uncertainty_strategy = "paired",
    seed = 1234, # Set for reproducibility
    alpha = 0.10, # Select a generous level of significance for the sake of the example
    posthoc_pairwise_tests = TRUE,
    return_perm_data = TRUE,
    extract_trait_data_melted_df = TRUE,
    extract_diversification_data_melted_df = TRUE,
    return_STRAPP_results = TRUE,
    return_updated_Maps = TRUE,
    return_updated_BAMM_objects = TRUE,
    verbose = TRUE,
    verbose_extended = TRUE) 

 ## Load directly deepSTRAPP output
 Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- readRDS(system.file("extdata",
    "Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rds", package = "deepSTRAPP"))
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 ## Explore output
 str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, max.level = 1)

 ## Plot histograms of STRAPP overall test results #
 # Tests are Kruskall-Wallis H tests when more than two states/ranges are present.
 # Tests are Mann-Whitney-Wilcoxon rank-sum tests when only two states/ranges are present.

 histogram_ggplots <- plot_histograms_STRAPP_tests_over_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
    display_plots = TRUE,
    # PDF_file_path = "./plot_STRAPP_histograms_overall_tests.pdf",
    plot_posthoc_tests = FALSE)

 # Print histogram for time step 1 = 0 My
 print(histogram_ggplots[[1]])
 # Adjust aesthetics of plot for time step 1 a posteriori
 histogram_ggplot_adj <- histogram_ggplots[[1]] +
     ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
 print(histogram_ggplot_adj)

 ## Plot histograms of STRAPP post hoc test results ------ #
 # Tests are Dunn's multiple comparison pairwise post hoc tests possible
 # only when more than two states/ranges are present.

 histograms_ggplots_list <- plot_histograms_STRAPP_tests_over_time(
     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
     display_plots = TRUE,
     # PDF_file_path = "./plot_STRAPP_histograms_posthoc_tests.pdf",
     plot_posthoc_tests = TRUE)

 # Print all histograms for time step 1 (= 0 My) one by one
 print(histograms_ggplots_list[[1]])
 # Plot all histograms for time step 1 (= 0 My) on one faceted plot
 cowplot::plot_grid(plotlist = histograms_ggplots_list[[1]])
}
```
