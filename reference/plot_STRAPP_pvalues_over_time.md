# Plot evolution of p-values of STRAPP tests over time

Plot the evolution of the p-values of STRAPP tests carried out across
multiple `time_steps`, obtained from
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).

By default, return a plot with a single line for p-values of overall
tests. If `plot_posthoc_tests = TRUE`, it will return a plot with
multiple lines, one per pair in post hoc tests (only for multinomial
data, with more than two states).

If a PDF file path is provided in `PDF_file_path`, the plot will be
saved directly in a PDF file.

## Usage

``` r
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs,
  time_range = NULL,
  pvalues_max = NULL,
  alpha = 0.05,
  display_plot = TRUE,
  plot_significant_time_frame = TRUE,
  plot_posthoc_tests = FALSE,
  select_posthoc_pairs = "all",
  plot_adjusted_pvalues = FALSE,
  PDF_file_path = NULL
)
```

## Arguments

- deepSTRAPP_outputs:

  List of elements generated with
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md),
  that summarize the results of multiple deepSTRAPP runs across
  `$time_steps`.

- time_range:

  Numeric vector of length 2. Time boundaries used for the X-axis of the
  plot. If `NULL` (the default), the range of data provided in
  `deepSTRAPP_outputs` will be used.

- pvalues_max:

  Numeric. Set the max boundary used for the Y-axis of the plot. If
  `NULL` (the default), the maximum p-value provided in
  `deepSTRAPP_outputs` will be used.

- alpha:

  Numeric. Significance level to display as a red dashed line on the
  plot. If set to `NULL`, no line will be added. Default is `0.05`.

- display_plot:

  Logical. Whether to display the plot generated in the R console.
  Default is `TRUE`.

- plot_significant_time_frame:

  Logical. Whether to display a green band over the time frame that
  yields significant results according to the chosen alpha level.
  Default is `TRUE`.

- plot_posthoc_tests:

  Logical. For multinomial data only. Whether to plot the p-values for
  the overall Kruskal-Wallis test across all states
  (`plot_posthoc_tests = FALSE`), or plot the p-values for the pairwise
  post hoc Dunn's test across pairs of states
  (`plot_posthoc_tests = TRUE`). Default is `FALSE`. This is only
  possible if `deepSTRAPP_outputs` contains the
  `$pvalues_summary_df_for_posthoc_pairwise_tests` element returned by
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
  when `posthoc_pairwise_tests = TRUE`.

- select_posthoc_pairs:

  Character vector used to specify the pairs to include in the plot.
  Names of pairs must match the pairs found in
  `deepSTRAPP_outputs$pvalues_summary_df_for_posthoc_pairwise_tests$pair`.
  Default is "all" to include all pairs.

- plot_adjusted_pvalues:

  Logical. Whether to display the p-values adjusted for multiple testing
  rather than the raw p-values. See argument 'p.adjust_method' in
  [`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
  or
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).
  Default is `FALSE`.

- PDF_file_path:

  Character string. If provided, the plot will be saved in a PDF file
  following the path provided here. The path must end with ".pdf".

## Value

The function returns a list of classes `gg` and `ggplot`. This object is
a ggplot that can be displayed on the console with `print(output)`. It
corresponds to the plot being displayed on the console when the function
is run, if `display_plot = TRUE`, and can be further modified for
aesthetics using the ggplot2 grammar.

If a `PDF_file_path` is provided, the function will also generate a PDF
file of the plot.

## Details

Plots are built based on the p-values recorded in summary_df provided by
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).

For overall tests, those p-values are found in `$pvalues_summary_df`.

For multinomial data (categorical or biogeographic data with more than 2
states), it is possible to plot p-values of post hoc pairwise tests. Set
`plot_posthoc_tests = TRUE` to generate plots for the pairwise post hoc
Dunn's test across pairs of states. To achieve this, the
`deepSTRAPP_outputs` input object must contain a
`$pvalues_summary_df_for_posthoc_pairwise_tests` element that summarizes
p-values computed across pairs of states for all post hoc tests. This is
obtained from
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)
when setting `posthoc_pairwise_tests = TRUE` to carry out post hoc
tests.

## See also

[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)

## Author

Maël Doré

## Examples

``` r
if (deepSTRAPP::is_dev_version())
{
  ## Load results of run_deepSTRAPP_over_time() for categorical data with 3-levels
  Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40 <- readRDS(system.file("extdata",
     "Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40.rds", package = "deepSTRAPP"))
  ## This dataset is only available in development versions installed from GitHub.
  # It is not available in CRAN versions.
  # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

  ## Plot results of overall Kruskal-Wallis / Mann-Whitney-Wilcoxon tests across all time-steps
  plot_overall <- plot_STRAPP_pvalues_over_time(
     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
     alpha = 0.10,
     time_range = c(0, 30), # Adjust time range if needed
     display_plot = FALSE)

  # Adjust aesthetics a posteriori
  plot_overall <- plot_overall +
     ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16))

  print(plot_overall)

  ## Plot results of post hoc pairwise Dunn's tests between selected pairs of states
  plot_posthoc <- plot_STRAPP_pvalues_over_time(
     deepSTRAPP_outputs = Ponerinae_deepSTRAPP_cat_3lvl_old_calib_0_40,
     alpha = 0.10,
     plot_posthoc_tests = TRUE,
     # PDF_file_path = "./pvalues_over_time.pdf",
     select_posthoc_pairs = c("arboreal != subterranean",
                              "arboreal != terricolous"),
     display_plot = FALSE)

  # Adjust aesthetics a posteriori
  plot_posthoc <- plot_posthoc +
     ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size  = 14),
        legend.position.inside = c(0.25, 0.25))

  print(plot_posthoc)
}
```
