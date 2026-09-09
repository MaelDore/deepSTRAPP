# Plot histogram of STRAPP test statistics to assess results

Plot a histogram of the distribution of the test statistics obtained
from a STRAPP test carried out for a unique `focal_time`. (See
[`compute_STRAPP_test_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/compute_STRAPP_test_for_focal_time.md)
and
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)).

Returns a single histogram for overall tests. If
`plot_posthoc_tests = TRUE`, it will return a faceted plot with a
histogram per post hoc test.

If a PDF file path is provided in `PDF_file_path`, the plot will be
saved directly in a PDF file.

## Usage

``` r
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs,
  focal_time = NULL,
  display_plot = TRUE,
  plot_posthoc_tests = FALSE,
  PDF_file_path = NULL
)
```

## Arguments

- deepSTRAPP_outputs:

  List of elements generated with
  [`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md),
  that summarize the results of a STRAPP test for a specific time in the
  past (i.e. the `focal_time`). `deepSTRAPP_outputs` can also be
  extracted from the output of
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  that run the whole deepSTRAPP workflow and store results inside
  `$STRAPP_results_over_time`.

- focal_time:

  Numeric. (Optional) If `deepSTRAPP_outputs` comprises results over
  multiple time-steps (i.e., output of
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md),
  this is the time of the STRAPP test targeted for plotting.

- display_plot:

  Logical. Whether to display the histogram(s) generated in the R
  console. Default is `TRUE`.

- plot_posthoc_tests:

  Logical. For multinomial data only. Whether to plot the histogram for
  the overall Kruskal-Wallis test across all states
  (`plot_posthoc_tests = FALSE`), or plot the histograms for all the
  pairwise post hoc Dunn's tests across pairs of states
  (`plot_posthoc_tests = TRUE`). Default is `FALSE`.

- PDF_file_path:

  Character string. If provided, the plot will be saved in a PDF file
  following the path provided here. The path must end with ".pdf".

## Value

By default, the function returns a list of classes `gg` and `ggplot`.
This object is a ggplot that can be displayed on the console with
`print(output)`. It corresponds to the histogram being displayed on the
console when the function is run, if `display_plot = TRUE`, and can be
further modify for aesthetics using the ggplot2 grammar.

If using multinomial data and setting `plot_posthoc_tests = TRUE`, the
function will return a list of objects. Each object is the ggplot
associated with a pairwise post hoc test. To plot each histogram `i`
individually, use `print(output_list[[i]])`. To plot all histograms at
once in a multifaceted plot, as displayed on the console if
`display_plot = TRUE`, use `cowplot::plot_grid(plotlist = output_list)`.

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
file of the plot. For post hoc tests, this will save the multifaceted
plot.

## Details

The main input `deepSTRAPP_outputs` is the typical output of
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md).
It provides information on results of a STRAPP test performed at a given
`focal_time`.

Histograms are built based on the distribution of the test statistics.
Such distributions are recorded in the outputs of a deepSTRAPP run
carried out with
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
when `return_perm_data = TRUE` so that the distributions of test stats
computed across posterior samples are returned among the outputs under
`$STRAPP_results$perm_data_df`.

For multinomial data (categorical or biogeographic data with more than 2
states), it is possible to plot the histograms of post hoc pairwise
tests. Set `plot_posthoc_tests = TRUE` to generate histograms for all
the pairwise post hoc Dunn's tests across pairs of states. To achieve
this, the `deepSTRAPP_outputs` input object must contain a
`$STRAPP_results$posthoc_pairwise_tests$perm_data_array` element that
summarizes test statistics computed across posterior samples for all
pairwise post hoc tests. This is obtained from
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
when setting both `posthoc_pairwise_tests = TRUE` to carry out post hoc
tests, and `return_perm_data = TRUE` to record distributions of test
statistics.

Alternatively, the main input `deepSTRAPP_outputs` can be the output of
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md),
providing results of STRAPP tests over multiple time-steps. In this
case, you must provide a `focal_time` to select the unique time-step
used for plotting.

- `return_perm_data` must be set to `TRUE` so that the permutation data
  used to compute the tests are returned among the outputs under
  `$STRAPP_results_over_time[[i]]$perm_data_df`.

- `posthoc_pairwise_tests` must be set to `TRUE` so that the permutation
  data used to perform the post hoc tests are also returned among the
  outputs under
  `$STRAPP_results_over_time[[i]]$posthoc_pairwise_tests$perm_data_array`.

For plotting all time-steps at once, see
[`plot_histograms_STRAPP_tests_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_histograms_STRAPP_tests_over_time.md).

## See also

Associated functions in deepSTRAPP:
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
[`plot_histograms_STRAPP_tests_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_histograms_STRAPP_tests_over_time.md)

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

 ## Set focal time to 10 Mya
 focal_time <- 10

 ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.

 Ponerinae_deepSTRAPP_cont_old_calib_10My <- run_deepSTRAPP_for_focal_time(
   contMap = Ponerinae_contMap,
   ace = Ponerinae_ACE,
   tip_data = Ponerinae_cont_tip_data,
   trait_data_type = "continuous",
   BAMM_object = Ponerinae_BAMM_object_old_calib,
   focal_time = focal_time,
   rate_type = "net_diversification",
   uncertainty_strategy = "rates_only",
   return_perm_data = TRUE)

 ## Explore output
 str(Ponerinae_deepSTRAPP_cont_old_calib_10My, max.level = 1)

 # ----- Plot histogram of STRAPP overall test results from run_deepSTRAPP_for_focal_time() ----- #

 histogram_ggplot <- plot_histogram_STRAPP_test_for_focal_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_10My,
    display_plot = TRUE,
    # PDF_file_path = "./plot_STRAPP_histogram_10My.pdf",
    plot_posthoc_tests = FALSE)
 # Adjust aesthetics a posteriori
 histogram_ggplot_adj <- histogram_ggplot +
    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
 print(histogram_ggplot_adj) 

 # ----- Plot histogram of STRAPP overall test results from run_deepSTRAPP_over_time() ----- #

 ## Load directly outputs from run_deepSTRAPP_over_time()
 Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- readRDS(system.file("extdata",
    "Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", package = "deepSTRAPP"))
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 # Select focal_time = 10My
 focal_time <- 10

 ## Plot histogram for overall test
 histogram_ggplot <- plot_histogram_STRAPP_test_for_focal_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
    focal_time = focal_time,
    display_plot = TRUE,
    # PDF_file_path = "./plot_STRAPP_histogram_10My.pdf",
    plot_posthoc_tests = FALSE)
 # Adjust aesthetics a posteriori
 histogram_ggplot_adj <- histogram_ggplot +
   ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
 print(histogram_ggplot_adj)

 # ----- Example 2: Categorical trait ----- #

 ## Load data

 # Load phylogeny
 data(Ponerinae_tree_old_calib, package = "deepSTRAPP")
 # Load trait df
 data(Ponerinae_trait_tip_data, package = "deepSTRAPP")

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
    colors_per_states = colors_per_states,
    evolutionary_models = "ARD", # Use default ARD model
    nb_simulations = 100, # Reduce number of simulations to save time
    seed = 1234, # Seet seed for reproducibility
    return_best_model_fit = TRUE,
    return_model_selection_df = TRUE,
    plot_map = FALSE) 

 # Load directly output
 data(Ponerinae_cat_3lvl_data_old_calib, package = "deepSTRAPP")

 ## Set focal time to 10 Mya
 focal_time <- 10

  # (May take several minutes to run)
 ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.

 Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My <- run_deepSTRAPP_for_focal_time(
    densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
    nb_simulations = 100,
    ace = Ponerinae_cat_3lvl_data_old_calib$ace,
    tip_data = Ponerinae_cat_3lvl_tip_data,
    trait_data_type = "categorical",
    BAMM_object = Ponerinae_BAMM_object_old_calib,
    focal_time = focal_time,
    rate_type = "net_diversification",
    uncertainty_strategy = "paired",
    posthoc_pairwise_tests = TRUE,
    return_perm_data = TRUE,
    return_updated_Maps = TRUE,
    return_updated_BAMM_object = TRUE)

 ## Explore output
 str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My, max.level = 1)

 # ------ Plot histogram of STRAPP overall test results ------ #

 histogram_ggplot <- plot_histogram_STRAPP_test_for_focal_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My,
    display_plot = TRUE,
    # PDF_file_path = "./plot_STRAPP_histogram_overall_test.pdf",
    plot_posthoc_tests = FALSE)
 # Adjust aesthetics a posteriori
 histogram_ggplot_adj <- histogram_ggplot +
    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
 print(histogram_ggplot_adj)

 # ------ Plot histograms of STRAPP post hoc test results ------ #

 histograms_ggplot_list <- plot_histogram_STRAPP_test_for_focal_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My,
    display_plot = TRUE,
    # PDF_file_path = "./plot_STRAPP_histograms_posthoc_tests.pdf",
    plot_posthoc_tests = TRUE)
 # Plot all histograms one by one
 print(histograms_ggplot_list)
 # Plot all histograms on one faceted plot
 cowplot::plot_grid(plotlist = histograms_ggplot_list) 
}
```
