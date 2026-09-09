# Plot rates vs. trait data for a given focal time

Plot rates vs. trait data of branches as extracted for a given focal
time. Data are extracted from the output of a deepSTRAPP run carried out
with
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
or
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)).

Returns a single plot showing mean rates vs. trait data extracted across
all branches for a given focal time. If the trait data are 'continuous',
the plot is a scatter plot. If the trait data are 'categorical' or
'biogeographic', the plot is a boxplot.

If a PDF file path is provided in `PDF_file_path`, the plot will be
saved directly in a PDF file.

## Usage

``` r
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs,
  focal_time = NULL,
  rate_type = "net_diversification",
  select_trait_levels = "all",
  color_scale = NULL,
  colors_per_levels = NULL,
  display_plot = TRUE,
  PDF_file_path = NULL,
  return_mean_data_per_branches_df = FALSE
)
```

## Arguments

- deepSTRAPP_outputs:

  List of elements generated with
  [`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md),
  that summarize the results of a STRAPP test for a specific time in the
  past (i.e. the `focal_time`). The list needs to include two
  data.frames: `$trait_data_df_over_time` and
  `$diversification_data_df_over_time` by setting
  `extract_trait_data_melted_df = TRUE` and
  `extract_diversification_data_melted_df = TRUE`. `deepSTRAPP_outputs`
  can also be extracted from the output of
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  that runs the whole deepSTRAPP workflow over multiple time-steps.

- focal_time:

  Numeric. (Optional) If `deepSTRAPP_outputs` comprises results over
  multiple time-steps (i.e., output of
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md),
  this is the time of the STRAPP test targeted for plotting.

- rate_type:

  A character string specifying the type of diversification rates to
  plot. Must be one of 'speciation', 'extinction' or
  'net_diversification' (default). Even if the `deepSTRAPP_outputs`
  object was generated with
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  for testing another type of rates, the object will contain data for
  all types of rates.

- select_trait_levels:

  (Vector of) character string. Only for categorical and biogeographic
  trait data. To provide a list of a subset of states/ranges to plot.
  Names must match the ones found in the `deepSTRAPP_outputs`. Default
  is `all` which means all states/ranges will be plotted.

- color_scale:

  Character vector. List of colors to use to build the color scale with
  [`grDevices::colorRampPalette()`](https://rdrr.io/r/grDevices/colorRamp.html)
  to display the points. Color scale from lowest values to highest rate
  values. Only for continuous data. Default = `NULL` will use the
  'Spectral' color palette in
  [`RColorBrewer::brewer.pal()`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).

- colors_per_levels:

  Named character string. To set the colors to use to plot data points
  and box for each state/range. Names = states/ranges; values = colors.
  If `NULL` (default), the default ggplot2 color palette
  ([`scales::hue_pal()`](https://scales.r-lib.org/reference/pal_hue.html))
  will be used. Only for categorical and biogeographic data.

- display_plot:

  Logical. Whether to display the plot generated in the R console.
  Default is `TRUE`.

- PDF_file_path:

  Character string. If provided, the plot will be saved in a PDF file
  following the path provided here. The path must end with ".pdf".

- return_mean_data_per_branches_df:

  Logical. Whether to include in the output the data.frame of mean rates
  per trait values/states/ranges of all branches computed at the focal
  time and used for the plot. Default is `FALSE`.

## Value

The function returns a list with at least one element.

- `rates_vs_trait_ggplot` An object of classes `gg` and `ggplot`. This
  is a ggplot that can be displayed on the console with
  `print(output$rates_vs_trait_ggplot)`. It corresponds to the plot
  being displayed on the console when the function is run, if
  `display_plot = TRUE`, and can be further modified for aesthetics
  using the ggplot2 grammar.

If the trait data are 'continuous', the plot is a scatter plot showing
how mean diversification rates vary with mean trait values across
branches. If the trait data are 'categorical' or 'biogeographic', the
plot is a boxplot showing mean diversification rates by states/ranges
per branch.

Each plot also displays summary statistics for the STRAPP test
associated with the data displayed:

- An observed statistic computed across the mean traits/ranges and rates
  values shown on the plot. This is not the statistic of the STRAPP test
  itself, which is conducted across stochastic maps X BAMM posterior
  samples (i.e., it is a distribution, not a unique value).

- The quantile of null statistic distribution at the significance
  threshold used to define test significance. The test will be
  considered significant (i.e., the null hypothesis is rejected) if this
  value is higher than zero.

- The p-value of the associated STRAPP test.

Optional summary data.frame:

- `mean_data_per_branches_df` A data.frame with four columns providing
  the `$mean_trait_values` (continuous) / `$states` (categorical) /
  `$ranges` (biogeographic) and `$mean_rates` computed along branches
  (`$tip_ID`) at `focal_time`. Rates/Traits are averaged for each
  stochastic maps X BAMM posterior samples used for testing, in
  accordance with the `uncertainty_strategy` selected when running
  deepSTRAPP. This is the raw data used to draw the plot. Included if
  `return_mean_data_per_branches_df = TRUE`.

If a `PDF_file_path` is provided, the function will also generate a PDF
file of the plot.

## Details

The main input `deepSTRAPP_outputs` is the typical output of
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md).
It provides information on results of a STRAPP test performed at a given
`focal_time`.

Plots are built based on both trait data and diversification data as
extracted for the given `focal_time`. Such data are recorded in the
outputs of a deepSTRAPP run carried out with
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
when `extract_trait_data_melted_df = TRUE` for trait data, and
`extract_diversification_data_melted_df = TRUE` for diversification
data. Please ensure to select those arguments when running deepSTRAPP.

Alternatively, the main input `deepSTRAPP_outputs` can be the output of
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md),
providing results of STRAPP tests over multiple time-steps. In this
case, you must provide a `focal_time` to select the unique time-step
used for plotting.

- `extract_trait_data_melted_df` must be set to `TRUE` so that trait
  data are returned in a melted data.frame among the outputs under
  `$trait_data_df_over_time`.

- `extract_diversification_data_melted_df` must be set to `TRUE` so that
  the diversification rates are returned in a melted data.frame among
  the outputs under `$diversification_data_df_over_time`.

- `return_STRAPP_results` must be set to `TRUE` so that the STRAPP
  results objects recording the summary statistics for the STRAPP tests
  are returned among the outputs under `$STRAPP_results_over_time`.

For plotting all time-steps at once, see
[`plot_rates_vs_trait_data_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_rates_vs_trait_data_over_time.md).

## See also

Associated functions in deepSTRAPP:
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
[`plot_rates_vs_trait_data_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_rates_vs_trait_data_over_time.md)

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
 Ponerinae_contMap <- phytools::contMap(Ponerinae_tree_old_calib, x = Ponerinae_cont_tip_data,
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
    return_perm_data = TRUE,
    # Need to be set to TRUE to save diversification data
    extract_diversification_data_melted_df = TRUE,
    # Need to be set to TRUE to save trait data
    extract_trait_data_melted_df = TRUE,
    return_updated_BAMM_object = TRUE)

 ## Explore output
 str(Ponerinae_deepSTRAPP_cont_old_calib_10My, max.level = 1)

 # ----- Plot scatterplot of rates vs. trait values from run_deepSTRAPP_for_focal_time() ----- #

 # Get plot
 rates_vs_trait_output <- plot_rates_vs_trait_data_for_focal_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_10My,
    color_scale = c("grey80", "orange"),
    display_plot = TRUE,
    # PDF_file_path = "./plot_rates_vs_trait_10My.pdf"
    return_mean_data_per_branches_df = TRUE)
 # Adjust aesthetics a posteriori
 rates_vs_trait_ggplot_adj <- rates_vs_trait_output$rates_vs_trait_ggplot +
    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
 print(rates_vs_trait_ggplot_adj)

 # Explore melted data.frame of mean rates and trait values extracted for the given focal time.
 head(rates_vs_trait_output$mean_data_per_branches_df) 

 # ----- Plot scatterplot of rates vs. trait values from run_deepSTRAPP_over_time() ----- #

 ## Load directly outputs from run_deepSTRAPP_over_time()
 Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- readRDS(system.file("extdata",
    "Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", package = "deepSTRAPP"))
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 # Select focal_time = 10My
 focal_time <- 10

 # Get plot
 rates_vs_trait_output <- plot_rates_vs_trait_data_for_focal_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
    focal_time = focal_time,
    color_scale = c("grey80", "purple"),
    display_plot = TRUE)
    # PDF_file_path = "./plot_rates_vs_trait_10My.pdf"

 # Adjust aesthetics a posteriori
 rates_vs_trait_ggplot_adj <- rates_vs_trait_output$rates_vs_trait_ggplot +
     ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
 print(rates_vs_trait_ggplot_adj)


 # ----- Example 2: Categorical trait ----- #

 ## Load data

 # Load phylogeny
 data(Ponerinae_tree, package = "deepSTRAPP")
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

 ## Run deepSTRAPP on net diversification rates for focal time = 10 Mya.

 Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My <- run_deepSTRAPP_for_focal_time(
    densityMaps = Ponerinae_cat_3lvl_data_old_calib$densityMaps,
    ace = Ponerinae_cat_3lvl_data_old_calib$ace,
    tip_data = Ponerinae_cat_3lvl_tip_data,
    nb_simulations = 100,
    trait_data_type = "categorical",
    BAMM_object = Ponerinae_BAMM_object_old_calib,
    focal_time = focal_time,
    rate_type = "net_diversification",
    uncertainty_strategy = "paired",
    posthoc_pairwise_tests = TRUE,
    return_perm_data = TRUE,
    # Need to be set to TRUE to save diversification data
    extract_diversification_data_melted_df = TRUE,
    # Need to be set to TRUE to save trait data
    extract_trait_data_melted_df = TRUE,
    return_updated_BAMM_object = TRUE)

 ## Explore output
 str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My, max.level = 1)

 ## Plot rates vs. states
 rates_vs_trait_output <- plot_rates_vs_trait_data_for_focal_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_10My,
    focal_time = 10,
    select_trait_levels = c("arboreal", "terricolous"), # Select only two levels
    colors_per_levels = colors_per_states[c("arboreal", "terricolous")], # Adjust colors
    display_plot = TRUE,
    # PDF_file_path = "./plot_rates_vs_trait_10My.pdf",
    return_mean_data_per_branches_df = TRUE)

 # Adjust aesthetics a posteriori
 rates_vs_trait_ggplot_adj <- rates_vs_trait_output$rates_vs_trait_ggplot +
    ggplot2::theme(plot.title = ggplot2::element_text(color = "red", size = 15))
 print(rates_vs_trait_ggplot_adj)

 # Explore melted data.frame of mean rates and states extracted for the given focal time.
 head(rates_vs_trait_output$mean_data_per_branches_df) 
 }
```
