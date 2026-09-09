# Plot mean rates vs. trait data over time-steps

Plot mean rates vs. trait data of branches as extracted for all
time-steps. Data are extracted from the output of a deepSTRAPP run
carried out with
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md))
over multiple time-steps.

For each time-step, returns a plot showing rates vs. trait data. If the
trait data are 'continuous', the plot is a scatter plot. If the trait
data are 'categorical' or 'biogeographic', the plot is a boxplot.

If a PDF file path is provided in `PDF_file_path`, all plots will be
saved directly in a unique PDF file, with one page per plot/time-step.

## Usage

``` r
plot_rates_vs_trait_data_over_time(
  deepSTRAPP_outputs,
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
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  that runs the whole deepSTRAPP workflow over multiple time-steps. The
  list needs to include two data.frames: `$trait_data_df_over_time` and
  `$diversification_data_df_over_time` by setting
  `extract_trait_data_melted_df = TRUE` and
  `extract_diversification_data_melted_df = TRUE`.

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
  per trait values/states/ranges of all branches computed over all
  time-steps and used for the plots. Default is `FALSE`.

## Value

The function returns a list with at least one element.

- `rates_vs_trait_ggplots` A list of objects of classes `gg` and
  `ggplot` ordered as in `$time_steps`. Each element corresponds to a
  ggplot for a given `focal_time`. They can be displayed on the console
  with `print(output$rates_vs_trait_ggplots[[i]])`. They correspond to
  the plots being displayed on the console one by one when the function
  is run, if `display_plot = TRUE`, and can be further modified for
  aesthetics using the ggplot2 grammar.

If the trait data are 'continuous', the plots are scatter plots showing
how mean diversification rates vary with mean trait values across
branches. If the trait data are 'categorical' or 'biogeographic', the
plots are boxplots showing mean diversification rates per state/ranges
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
  the `$mean_trait_value` (continuous) / `$states` (categorical) /
  `$ranges` (biogeographic) and `$mean_rates` computed along branches
  (`$tip_ID`) at the different `focal_time`. Rates/Traits are averaged
  for each stochastic maps X BAMM posterior samples used for testing, in
  accordance with the `uncertainty_strategy` selected when running
  deepSTRAPP. This is the raw data used to draw each plot for each
  `focal_time`. Included if `return_mean_data_per_branches_df = TRUE`.

If a `PDF_file_path` is provided, the function will also generate a
unique PDF file with one plot/page per `$time_steps`.

## Details

The main input `deepSTRAPP_outputs` is the typical output of
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).
It provides information on results of a STRAPP test performed over
multiple time-steps.

Plots are built based on both trait data and diversification data as
extracted for the time-steps. Such data are recorded in the outputs of a
deepSTRAPP run carried out with
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).

- `extract_trait_data_melted_df` must be set to `TRUE` so that trait
  data are returned in a melted data.frame among the outputs under
  `$trait_data_df_over_time`.

- `extract_diversification_data_melted_df` must be set to `TRUE` so that
  the diversification rates are returned among the outputs under
  `$diversification_data_df_over_time`.

- `return_STRAPP_results` must be set to `TRUE` so that the STRAPP
  results objects recording the summary statistics for the STRAPP tests
  are returned among the outputs under `$STRAPP_results_over_time`.

For plotting a single `focal_time`, see
[`plot_rates_vs_trait_data_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_rates_vs_trait_data_for_focal_time.md).

## See also

Associated functions in deepSTRAPP:
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
[`plot_rates_vs_trait_data_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/plot_rates_vs_trait_data_for_focal_time.md)

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
    # Need to be set to TRUE to save trait data
    extract_trait_data_melted_df = TRUE,
    # Need to be set to TRUE to save diversification data
    extract_diversification_data_melted_df = TRUE,
    return_STRAPP_results = TRUE,
    return_updated_BAMM_object = TRUE,
    verbose = TRUE,
    verbose_extended = TRUE) 

 ## Load directly trait data output
 Ponerinae_deepSTRAPP_cont_old_calib_0_40 <- readRDS(system.file("extdata",
    "Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", package = "deepSTRAPP"))
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 # Explore output
 str(Ponerinae_deepSTRAPP_cont_old_calib_0_40, max.level = 1)

 ## Plot for all time-steps
 rates_vs_trait_outputs <- plot_rates_vs_trait_data_over_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cont_old_calib_0_40,
    color_scale = c("grey80", "purple"), # Adjust color scale
    display_plot = TRUE,
    # PDF_file_path = "./plot_rates_vs_trait_0_40My.pdf",
    return_mean_data_per_branches_df = TRUE)

 ## Print plot for time step 3 = 10 My
 print(rates_vs_trait_outputs$rates_vs_trait_ggplots[[3]])

 ## Explore melted data.frame of rates and trait data
 head(rates_vs_trait_outputs$mean_data_per_branches_df)

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
    # Need to be set to TRUE to save trait data
    extract_trait_data_melted_df = TRUE,
    # Need to be set to TRUE to save diversification data
    extract_diversification_data_melted_df = TRUE,
    return_STRAPP_results = TRUE,
    return_updated_BAMM_object = TRUE,
    verbose = TRUE,
    verbose_extended = TRUE) 

 ## Load directly deepSTRAPP output
 Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- readRDS(system.file("extdata",
    "Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rds", package = "deepSTRAPP"))
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 # Explore output
 str(Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40, max.level = 1)

 # Adjust color scheme
 colors_per_states <- c("orange", "dodgerblue", "red")
 names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")

 ## Plot for all time-steps
 rates_vs_trait_outputs <- plot_rates_vs_trait_data_over_time(
    deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
    colors_per_levels = colors_per_states, # Adjust color scheme
    display_plot = TRUE,
    # PDF_file_path = "./plot_rates_vs_trait_0_40My.pdf",
    return_mean_data_per_branches_df = TRUE)

 ## Print plot for time step 3 = 10 My
 print(rates_vs_trait_outputs$rates_vs_trait_ggplots[[3]])

 ## Explore melted data.frame of rates and states
 head(rates_vs_trait_outputs$mean_data_per_branches_df)
}
```
