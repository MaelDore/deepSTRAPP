# Convert a contsimmap object into a list of contMaps

Convert a contsimmap object generated with
`contsimmap::make.contsimmap()` into a list of `contMaps` as produced by
the
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html)
function.

## Usage

``` r
convert_contsimmap_to_contMaps(contsimmap, verbose = FALSE)
```

## Arguments

- contsimmap:

  List of class `"contsimmap"`, typically generated with
  `contsimmap::make.contsimmap()`, that summarizes a set of continuous
  stochastic maps. Each map represents a simulation of the evolution of
  a continuous trait compatible with the associated evolutionary model
  fit on observed trait data.

- verbose:

  Logical. Whether to display conversion progress every 10 maps.

## Value

Returns a list of `"contMap"` objects. Each `"contMap"` represents a
unique evolutionary history simulation conditioned on the observed trait
data and model fit (i.e., a continuous stochastic map).

A `"contMap"` typically contains:

- `$tree` A list of classes `"simmap"` and `"phylo"`. The mapped
  phylogeny including:

  - `$maps` A list of named numeric vectors. Provides the mapping of
    trait values along each edge (scaled between 0 and 1000).

  - `$mapped.edge` A numeric matrix. Provides the evolutionary time
    spent across trait values (columns) along the edges (rows).

- `$cols` A named vector mapping values to colors used to display trait
  evolution along the branches.

- `$lims` A numeric vector. Minimum and maximum trait values recorded
  during the simulation.

## Details

A contsimmap object is typically produced by the
`contsimmap::make.contsimmap()` function. It stores results of a
continuous stochastic mapping simulation that represent independent
evolutionary history of a continuous trait mapped along a
time-calibrated phylogeny, and compatible with the associated
evolutionary model fit on observed trait data. The method is described
in Martin & Weber, 2026.

A typical `"contsimmap"` object stores simulated trait values across
maps in a 3D array, recording values for all time-points used to
represent the continuous evolution. While `"contsimmap"` objects can
record evolution of multivariate traits, deepSTRAPP currently works only
with simulation of univariate trait evolution. See the
`contsimmap::make.contsimmap()` help section for a detailed description
of the structure of the object.

A `"contMap"` object typically produced by the
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html)
function also records the evolution of a continuous trait alongside
branches of a time-calibrated phylogeny. Trait values are stored in
`$maps` as list of named vectors recording trait values (as names scaled
between 0 and 1000) and length of the associated branch segment between
two time-points (as values). Therefore, a unique `"contsimmap"` object
is converted into a list of `"contMaps"`, with each item being a
`"contMap"` that represents a unique evolutionary history simulation
conditioned on the observed trait data and model fit.

## References

For continuous stochastic mapping: Martin, B. S., & Weber, M. G. (2026).
Stochastic character mapping of continuous traits on phylogenies.
Systematic Biology, syag031.
[doi:10.1093/sysbio/syag031](https://doi.org/10.1093/sysbio/syag031) .

## See also

`contsimmap::make.contsimmap()`
[`phytools::contMap()`](https://rdrr.io/pkg/phytools/man/contMap.html)
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md)

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
    trait.data = eel_tip_data,
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
 
}
```
