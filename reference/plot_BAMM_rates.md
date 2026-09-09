# Plot diversification rates and regime shifts from BAMM on phylogeny

Plot on a time-calibrated phylogeny the evolution of diversification
rates and the location of regime shifts estimated from a BAMM (Bayesian
Analysis of Macroevolutionary Mixtures). Each branch is colored
according to the estimated rates of speciation, extinction, or net
diversification stored in an object of class `"bammdata"`. Rates can
vary along time, thus colors evolve along individual branches.

This function is a wrapper of original functions from the R package
`{BAMMtools}`:

- Step 1: Use
  [`BAMMtools::plot.bammdata()`](https://rdrr.io/pkg/BAMMtools/man/plot.html)
  to map rates on the phylogeny.

- Step 2: Add the location of regime shifts with
  [`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html)
  (if `add_regime_shifts = TRUE`).

## Usage

``` r
plot_BAMM_rates(
  BAMM_object,
  rate_type = "net_diversification",
  method = "phylogram",
  configuration_type = "MAP",
  sample_index = 1,
  regimes_fill = "grey",
  regimes_size = 1,
  regimes_pch = 21,
  regimes_border_col = "black",
  regimes_border_width = 1,
  ...,
  add_regime_shifts = TRUE,
  adjust_size_to_prob = TRUE,
  display_plot = TRUE,
  PDF_file_path = NULL
)
```

## Arguments

- BAMM_object:

  Object of class `"bammdata"`, typically generated with
  [`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md),
  that contains a phylogenetic tree and associated diversification rate
  mapping across selected posterior samples. It works also for
  `BAMM_object` updated for a specific `focal_time` using
  [`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md),
  or the deepSTRAPP workflow with
  [run_deepSTRAPP_for_focal_time](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
  and
  [`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).

- rate_type:

  A character string specifying the type of diversification rates to
  plot. Must be one of 'speciation', 'extinction' or
  'net_diversification' (default).

- method:

  A character string indicating the method for plotting the phylogenetic
  tree.

  - `method = "phylogram"` (default) plots the phylogenetic tree using
    rectangular coordinates.

  - `method = "polar"` plots the phylogenetic tree using polar
    coordinates.

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
  [`BAMMtools::plot.bammdata()`](https://rdrr.io/pkg/BAMMtools/man/plot.html),
  [`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html),
  and [`par()`](https://rdrr.io/r/graphics/par.html). Among them,
  `par.reset` is ignored, with a message: it would make
  [`BAMMtools::plot.bammdata()`](https://rdrr.io/pkg/BAMMtools/man/plot.html)
  restore a full par() snapshot, discarding the coordinate system of the
  phylogeny and rewinding the panel counter of a multi-panel layout.
  Graphical parameters are instead restored by `plot_BAMM_rates()` when
  it exits: every parameter the call actually changed is put back,
  except those describing the plot itself, which are kept so that the
  phylogeny can be annotated afterwards with for instance
  [`graphics::abline()`](https://rdrr.io/r/graphics/abline.html) or
  [`graphics::title()`](https://rdrr.io/r/graphics/title.html), as in
  the examples.

- add_regime_shifts:

  Logical. Whether to add the location of regime shifts on the phylogeny
  (Step 2). Default is `TRUE`. Provide the full argument name to avoid
  ambiguity with the `add` graphical argument.

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

The function returns (invisibly) a list with three elements similarly to
[`BAMMtools::plot.bammdata()`](https://rdrr.io/pkg/BAMMtools/man/plot.html).

- `$coords`: A matrix of plot coordinates. Rows correspond to branches.
  Columns 1-2 are starting (x,y) coordinates of each branch and columns
  3-4 are ending (x,y) coordinates of each branch. If method = "polar" a
  fifth column gives the angle (in radians) of each branch.

- `$colorbreaks`: A vector of percentiles used to group
  macroevolutionary rates into color bins.

- `$colordens`: A matrix of the kernel density estimates (column 2) of
  evolutionary rates (column 1) and the color (column 3) corresponding
  to each rate value.

## Details

The main input `BAMM_object` is the typical output of
[`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md).
It provides information on rates and regime shifts across the posterior
samples of a BAMM.

`$MAP_BAMM_object` and `$MSC_BAMM_object` elements are required to plot
regime shift locations following the "MAP" or "MSC" `configuration_type`
respectively. A `$MSP_tree` element is required to scale the size of the
symbols showing the location of regime shifts according to marginal
shift probabilities. (If `adjust_size_to_prob = TRUE`).

The default option to display regime shift is to use the average
locations from the posterior samples with the Maximum A Posteriori
probability (MAP) configuration. However, sometimes, multiple
configurations have similarly high frequency in the posterior samples
(See
[`BAMMtools::credibleShiftSet()`](https://rdrr.io/pkg/BAMMtools/man/credibleShiftSet.html)).
An alternative is to use the average locations from posterior samples
with the Maximum Shift Credibility (MSC) configuration instead. This
regime shift configuration has the highest product of marginal
probabilities across branches where a shift is estimated. It may differ
from the MAP configuration. (See
[`BAMMtools::maximumShiftCredibility()`](https://rdrr.io/pkg/BAMMtools/man/maximumShiftCredibility.html)).

## See also

Initial functions in BAMMtools:
[`BAMMtools::plot.bammdata()`](https://rdrr.io/pkg/BAMMtools/man/plot.html)
[`BAMMtools::addBAMMshifts()`](https://rdrr.io/pkg/BAMMtools/man/addBAMMshifts.html)

Associated functions in deepSTRAPP:
[`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md)
[`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md)
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)

## Author

Maël Doré

Original functions by Mike Grundler & Pascal Title in R package
`{BAMMtools}`.

## Examples

``` r
# Load BAMM output
data(whale_BAMM_object, package = "deepSTRAPP")

mfrow_ini <- par()$mfrow
par(mfrow = c(1,2))

## Plot overall mean rates with MAP configuration for regime shifts
# (rates are averaged across all posterior samples)
plot_BAMM_rates(whale_BAMM_object, # Use the full object to plot mean overall rates
                add_regime_shifts = TRUE,
                configuration_type = "MAP",
                regimes_size = 3, bg = "black")
title("Overall mean - MAP shift config", col.main = "white")

## Plot overall mean rates with MSC configuration for regime shifts
# (rates are averaged across all posterior samples)
plot_BAMM_rates(whale_BAMM_object, add_regime_shifts = TRUE,
                configuration_type = "MSC",
                regimes_size = 3, bg = "black")
title("Overall mean - MSC shift config", col.main = "white")


par(mfrow_ini)
#> Warning: argument 1 does not name a graphical parameter
#> NULL


mfrow_ini <- par()$mfrow
par(mfrow = c(1,2))

## Plot mean MAP rates with regime shifts
# (rates are averaged only across MAP samples)
plot_BAMM_rates(whale_BAMM_object$MAP_BAMM_object, # Use the MAP object to plot mean MAP rates
                add_regime_shifts = TRUE,
                configuration_type = "index",
                # Set to index to use the regime shift location from MAP configuration
                regimes_size = 3)
title("MAP mean - MAP shift config")

## Plot mean MSC rates with regime shifts
# (rates averaged only across MSC samples)
plot_BAMM_rates(whale_BAMM_object$MSC_BAMM_object, # Use the MSC object to plot mean MSC rates
                add_regime_shifts = TRUE,
                configuration_type = "index",
                # Set to index to use the regime shift data from MSC configuration
                regimes_size = 3)
title("MSC mean - MSC shift config")


par(mfrow_ini)
#> Warning: argument 1 does not name a graphical parameter
#> NULL

```
