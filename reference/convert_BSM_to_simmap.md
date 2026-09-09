# Convert Biogeographic Stochastic Map (BSM) to phytools SIMMAP stochastic map (SM) format

These functions convert a Biogeographic Stochastic Map (BSM) output from
BioGeoBEARS into a `simmap` object from R package `{phytools}` (See
[`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html)).

They require a model fit with `BioGeoBEARS::bears_optim_run()` and the
output of a Biogeographic Stochastic Mapping performed with
`BioGeoBEARS::runBSM()` to produce `simmap` objects as phylogenies with
the associated mapping of range evolution along branches across
simulations.

- `convert_BSM_to_simmap()`: Produce one `simmap` for the required
  simulation (index of the simulation provided with `sim_index`).

- `convert_BSMs_to_simmaps()`: Produce all `simmap` objects for all
  simulations stored in a unique `multiSimmap` object.

Initial functions in R package BioGeoBEARS by Nicholas J. Matzke:

- `BioGeoBEARS::BSM_to_phytools_SM()`

- `BioGeoBEARS::BSMs_to_phytools_SMs()`

## Usage

``` r
convert_BSM_to_simmap(model_fit, phylo, BSM_output, sim_index)

convert_BSMs_to_simmaps(model_fit, phylo, BSM_output)
```

## Arguments

- model_fit:

  A BioGeoBEARS results object, produced by ML inference via
  `BioGeoBEARS::bears_optim_run()`.

- phylo:

  Time-calibrated phylogeny used in the BioGeoBEARS analyses to produce
  the historical biogeographic inference and run the Biogeographic
  Stochastic Mapping. Object of class `"phylo"` as defined in `{ape}`.

- BSM_output:

  A list with two objects, a cladogenetic events table and an anagenetic
  events table, as the result of Biogeographic Stochastic Mapping
  conducted with `BioGeoBEARS::runBSM()`.

- sim_index:

  Integer. Index of the biogeographic simulation targeted to produce the
  `simmap` with `convert_BSM_to_simmap()`.

## Value

The `convert_BSM_to_simmap()` function returns a list with two elements:

- `$simmap` A unique `simmap` for a given biogeographic simulation as an
  object of classes `c("simmap", "phylo")`. This is a modified `{ape}`
  tree with additional elements to report range mapping, model
  parameters and likelihood.

  - `$maps` A list of named numeric vectors. Provides the mapping of
    ranges along each remaining edge. Names are the ranges. Values are
    residence times in each state across segments

  - `$mapped.edge` A numeric matrix. Provides the evolutionary time
    spent across ranges (columns) along the edges (rows). row.names()
    are the node ID at the rootward and tipward ends of each edge.

  - `$Q` Numerical matrix. The transition rates across ranges calculated
    from the ML parameter estimates of the model.

  - `$logL` Numeric. The log-likelihood of the data under the ML model.

- `$residence_times` Data.frame with two rows. Summarizes the residence
  time spent in each range along all branches, in (raw) evolutionary
  time (i.e., branch lengths), and in percentage (perc).

The `convert_BSMs_to_simmaps()` function loops around the
`convert_BSM_to_simmap()` function to aggregate all `simmaps` from all
biogeographic simulations in a unique list of classes
`c("multiSimmap", "multiPhylo")`.

- Each element is the `$simmap` of a biogeographic simulation obtained
  with `convert_BSM_to_simmap()`.

- `$residence_times` summary data.frames are not preserved.

## Details

These functions are slight adaptations of original functions from the R
Package `BioGeoBEARS` by N. Matzke.

Initial functions: `BioGeoBEARS::BSM_to_phytools_SM()`
`BioGeoBEARS::BSMs_to_phytools_SMs()`

Changes:

- Solves issue with differences in ranges allowed across time-strata.

- Requires directly the output of `BioGeoBEARS::runBSM()` instead of
  separated cladogenetic and anagenetic event tables.

- Update the documentation.

## Notes on BioGeoBEARS

The R package `BioGeoBEARS` is needed for this function to work with
biogeographic data. Please install it manually from:
<https://github.com/nmatzke/BioGeoBEARS>.

## Notes on using the resulting simmap object in phytools (adapted from Nicholas J. Matzke)

The phytools functions, like
[`phytools::countSimmap()`](https://rdrr.io/pkg/phytools/man/countSimmap.html),
will only count the anagenetic events (range transitions occurring along
branches) as they were written assuming purely anagenetic models.

It remains possible to extract cladogenetic events (range transitions
occurring at speciation) by comparing the last-state-below-a-node with
the descendant-pairs-above-a-node. However, it is recommended to use the
built-in functions from BioGeoBEARS to summarize the biogeographic
history based on the tables of cladogenetic and anagenetic events
obtained from `BioGeoBEARS::runBSM()`. `simmap` objects should primarily
be considered as a tool for visualization.

Associated functions in R package `BioGeoBEARS`:

- `BioGeoBEARS::simulate_source_areas_ana_clado()`: To select randomly a
  unique area source for transition from a multi-area state to a single
  area.

- `BioGeoBEARS::get_dmat_times_from_res()`: To generate matrices of
  range expansion from source area to destination area.

- `BioGeoBEARS::count_ana_clado_events()`: To count the number and type
  of events from BSM tables.

- `BioGeoBEARS::hist_event_counts()`: To plot histograms of event counts
  across BSM tables.

Please note carefully that area-to-area dispersal events are not
identical with the state transitions. For example, a state can be a
geographic range with multiple areas, but under the logic of DEC-type
models, a range-expansion event like ABC-\>ABCD actually means that a
dispersal happened from some specific area (A, B, or C) to the new area.
BSMs track this area-to-area sourcing in their cladogenetic and
anagenetic event tables, at least if
`BioGeoBEARS::simulate_source_areas_ana_clado()` has been run on the
output of `BioGeoBEARS::runBSM()`.

## References

For BioGeoBEARS: Matzke, Nicholas J. (2018). BioGeoBEARS: BioGeography
with Bayesian (and likelihood) Evolutionary Analysis with R Scripts.
version 1.1.1, published on GitHub on November 6, 2018.
[doi:10.5281/zenodo.1478250](https://doi.org/10.5281/zenodo.1478250) .
Website: <http://phylo.wikidot.com/biogeobears>.

## See also

[`phytools::countSimmap()`](https://rdrr.io/pkg/phytools/man/countSimmap.html)
[`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html)
`BioGeoBEARS::simulate_source_areas_ana_clado()`
`BioGeoBEARS::get_dmat_times_from_res()`
`BioGeoBEARS::count_ana_clado_events()`
`BioGeoBEARS::hist_event_counts()`

## Author

Nicholas J. Matzke. Contact: <matzke@berkeley.edu>

Changes by Maël Doré (see Details)

## Examples

``` r
if (deepSTRAPP::is_dev_version())
{
 ## Run only if you have R package 'BioGeoBEARS' installed.
 # Please install it manually from: https://github.com/nmatzke/BioGeoBEARS")

 ## Load phylogeny and tip data
 library(phytools)
 data(eel.tree)

  # (May take several minutes to run)
 ## Load directly output of prepare_trait_data() run on biogeographic data
 data(eel_biogeo_data, package = "deepSTRAPP")

 ## Convert BSM output into a unique simmap, including residence times
 simmap_1 <- convert_BSM_to_simmap(model_fit = eel_biogeo_data$best_model_fit,
                                    phylo = eel.tree,
                                    BSM_output = eel_biogeo_data$BSM_output,
                                    sim_index = 1)
 # Explore output
 str(simmap_1, max.level = 1)
 # Print residence times in each range
 simmap_1$residence_times
 # Plot simmap
 plot(simmap_1$simmap)

 ## Convert BSM output into all simmaps in a multiSimmap/multiPhylo object
 all_simmaps <- convert_BSMs_to_simmaps(model_fit = eel_biogeo_data$best_model_fit,
                                        phylo = eel.tree,
                                        BSM_output = eel_biogeo_data$BSM_output)
 # Explore output
 str(all_simmaps, max.level = 1)
 # Plot simmap n°1
 plot(all_simmaps[[1]]) 
}
```
