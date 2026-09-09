# Aggregate a list of contMaps into a unique mean/median contMap

Aggregate a list of contMaps into a unique mean/median contMap as
produced by the
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html)
function.

## Usage

``` r
aggregate_contMaps(
  contMaps,
  fun = "mean",
  display_plot = TRUE,
  ...,
  color_scale = NULL,
  verbose = TRUE
)
```

## Arguments

- contMaps:

  List of objects of class `"contMap"` that represent independent
  simulations of the evolution of a continuous trait (i.e., continuous
  stochastic maps).

- fun:

  Character string. Select the aggregating function. Available options
  are `"mean"` and `"median"`. Default = `"mean"`.

- display_plot:

  Logical. Whether to plot the resulting aggregated contMap. Default =
  `TRUE`.

- ...:

  List of named arguments. Additional arguments to be passed down to
  [`deepSTRAPP::plot_contMap()`](https://maeldore.github.io/deepSTRAPP/reference/plot_contMap.md).

- color_scale:

  Character vector. List of colors to use to build the color scale with
  [`grDevices::colorRampPalette()`](https://rdrr.io/r/grDevices/colorRamp.html)
  showing the evolution of the continuous trait. From lowest values to
  highest values. Default (`color_scale = NULL`) is using the color
  palette recorded in the `contMap$cols` item. If none was provided, the
  [`rainbow()`](https://rdrr.io/r/grDevices/palettes.html) palette is
  used. Provide the full argument name to avoid ambiguity with the `col`
  graphical argument.

- verbose:

  Logical. Whether to display progress every 100 edges.

## Value

Returns a unique `"contMap"` object representing the average trait
evolutionary history recorded across all continuous stochastic maps
provided as input in `contMaps`.

The resulting aggregated `"contMap"` is a list with three items:

- `$tree` A list of classes `"simmap"` and `"phylo"`. The mapped
  phylogeny including:

  - `$maps` A list of named numeric vectors. Provides the mapping of
    aggregated trait values along each edge (scaled between 0 and 1000).

  - `$mapped.edge` A numeric matrix. Provides the evolutionary time
    spent across trait values (columns) along the edges (rows).

- `$cols` A named character vector mapping values to colors used to
  display trait evolution along the branches.

- `$lims` A numeric vector. Minimum and maximum aggregated trait values
  recorded.

## Details

The function is primarily designed to average trait values across
multiple continuous stochastic maps produced from the same simulation,
typically with `contsimmap::make.contsimmap()` and then converted to a
list of contMaps (phytools format) with the
[`convert_contsimmap_to_contMaps()`](https://maeldore.github.io/deepSTRAPP/reference/convert_contsimmap_to_contMaps.md)
function.

The result is a unique `contMap` which should be consistent with the
`contMap` produced by
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html)
when interpolating maximum likelihood estimates of ancestral trait
values between nodes.

## References

For continuous stochastic mapping: Martin, B. S., & Weber, M. G. (2026).
Stochastic character mapping of continuous traits on phylogenies.
Systematic Biology, syag031.
[doi:10.1093/sysbio/syag031](https://doi.org/10.1093/sysbio/syag031) .

## See also

`contsimmap::make.contsimmap()`
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html)
[`plot_contMap()`](https://maeldore.github.io/deepSTRAPP/reference/plot_contMap.md)

## Author

Maël Doré

## Examples

``` r
if (deepSTRAPP::is_dev_version())
{
 ## The R package 'contsimmap' is needed for this example to work.
 # Please install it manually from: https://github.com/bstaggmartin/contsimmap.

 # ----- Prepare data ----- #

 # Load eel phylogeny and tip data from the R package phytools
 # Source: Collar et al., 2014; DOI: 10.1038/ncomms6505

 data("eel.tree", package = "phytools")
 data("eel.data", package = "phytools")

 # Extract body size
 eel_tip_data <- stats::setNames(eel.data$Max_TL_cm,
                                 rownames(eel.data))

 # ----- Run continuous stochastic mapping ----- #

  # (May take several seconds to run)
 eel_contsimmap <- contsimmap::make.contsimmap(
    tree = eel.tree,
    trait.data = eel_tip_data ,
    nsim = 100, res = 100,
    verbose = TRUE)

 dim(eel_contsimmap) # 121 edges, 1 trait, 100 simulations

 # ----- Convert to a list of contMaps ----- #

 eel_contMaps <- convert_contsimmap_to_contMaps(
    contsimmap = eel_contsimmap, verbose = TRUE)

 # ----- Explore contMaps ----- #

 # Plot contMap from simulation n°1
 plot_contMap(eel_contMaps[[1]])
 # Plot contMap from simulation n°50
 plot_contMap(eel_contMaps[[50]])

 # ----- Aggregate all simulations ----- #

 eel_aggregated_contMap <- aggregate_contMaps(
    contMaps = eel_contMaps,
    verbose = TRUE, display_plot = FALSE)

 plot_contMap(eel_aggregated_contMap)

 # ----- Compare to interpolated contMap from phytools ----- #

 # Produce interpolated contMap with phytools
 eel_contMap <- phytools::contMap(tree = eel.tree, x = eel_tip_data,
                                  res = 100, # Number of time steps
                                  plot = FALSE)

 # Plot the interpolated contMap mapping Maximum Likelihood estimates
 # of ancestral trait values interpolated between nodes
 plot_contMap(eel_contMap)
 # Plot the aggregated contMap representing average trait values
 # from true continuous stochastic mapping simulations
 plot_contMap(eel_aggregated_contMap)
 # Both should be visually similar
 
}
```
