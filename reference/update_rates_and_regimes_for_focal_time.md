# Update diversification rates/regimes mapped on a phylogeny up to a given time in the past

Updates an object of class `"bammdata"` to obtain the diversification
rates/regimes found along branches at a specific time in the past (i.e.
the `focal_time`). Optionally, the function can update the object to
display a mapped phylogeny such that branches overlapping the
`focal_time` are shortened to the `focal_time`.

## Usage

``` r
update_rates_and_regimes_for_focal_time(
  BAMM_object,
  focal_time,
  update_rates = TRUE,
  update_regimes = TRUE,
  update_tree = FALSE,
  update_plot = FALSE,
  update_all_elements = FALSE,
  keep_tip_labels = TRUE,
  verbose = TRUE
)
```

## Arguments

- BAMM_object:

  Object of class `"bammdata"`, typically generated with
  [`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md),
  that contains a phylogenetic tree and associated diversification rate
  mapping across selected posterior samples. The phylogenetic tree must
  be rooted and fully resolved/dichotomous, but it does not need to be
  ultrametric (it can include fossils).

- focal_time:

  Numeric. The time, in terms of time distance from the present, at
  which the tree and rate mapping must be cut. It must be smaller than
  the root age of the phylogeny.

- update_rates:

  Logical. Specify whether diversification rates stored in `$tipLambda`
  (speciation) and `$tipMu` (extinction) must be updated to summarize
  rates found along branches at `focal_time`. Default is `TRUE`.

- update_regimes:

  Logical. Specify whether diversification regimes stored in
  `$tipStates` must be updated to summarize regimes found along branches
  at `focal_time`. Default is `TRUE`.

- update_tree:

  Logical. Specify whether to update the phylogeny such that all
  branches that are younger than the `focal_time` are cut-off. Default
  is `FALSE`.

- update_plot:

  Logical. Specify whether to update the phylogeny AND the elements used
  by
  [`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
  to plot diversification rates on the phylogeny such that all branches
  that are younger than the `focal_time` are cut-off. Default is
  `FALSE`. If set to `TRUE`, it will override the `update_tree`
  parameter and update the phylogeny.

- update_all_elements:

  Logical. Specify whether to update all the elements in the object,
  including rates/regimes/phylogeny/elements for
  [`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)/all
  other elements. Default is `FALSE`. If set to `TRUE`, it will override
  other `update_*` parameters and update all elements.

- keep_tip_labels:

  Logical. Should terminal branches with a single descendant tip retain
  their initial `tip.label` on the updated phylogeny and diversification
  rate mapping? Default is `TRUE`. If set to `FALSE`, the tipward node
  ID will be used as label for all tips.

- verbose:

  Logical. Should progression be displayed? A message will be printed
  for every batch of 100 BAMM posterior samples updated. Default is
  `TRUE`.

## Value

The function returns a list as an updated `BAMM_object` of class
`"bammdata"`.

Phylogeny-related elements used to plot a phylogeny with
[`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html):

- `$edge` Integer matrix. Defines the tree topology by providing
  rootward and tipward node ID of each edge.

- `$Nnode` Integer. Number of internal nodes.

- `$tip.label` Character vector. Labels of all tips, including fossils
  older than `focal_time` if present.

  - If `keep_tip_labels = TRUE`, cut-off branches with a single
    descendant tip retain their initial `tip.label`.

  - If `keep_tip_labels = FALSE`, all cut-off branches are labeled using
    their tipward node ID.

- `$edge.length` Numeric vector. Length of edges/branches.

- `$node.label` Character vector. Labels of all internal nodes. (Present
  only if present in the initial `BAMM_object`)

BAMM internal elements used for tree exploration:

- `$begin` Numeric vector. Absolute time since root of edge/branch start
  (rootward).

- `$end` Numeric vector. Absolute time since root of edge/branch end
  (tipward).

- `$downseq` Integer vector. Order of node visits when using a pre-order
  tree traversal.

- `$lastvisit` ID of the last node visited when starting from the node
  in the corresponding position in `$downseq`.

BAMM elements summarizing diversification data:

- `$numberEvents` Integer vector. Number of events/macroevolutionary
  regimes (k+1) recorded in each posterior configuration. k = number of
  shifts.

- `$eventData` List of data.frames. One per posterior sample. Records
  shift events and macroevolutionary regimes parameters. 1st line =
  Background root regime.

- `$eventVectors` List of integer vectors. One per posterior sample.
  Record regime ID per branch.

- `$tipStates` List of named integer vectors. One per posterior sample.
  Record regime ID per tip present at `focal_time`. Updated if
  `update_regimes = TRUE`.

- `$tipLambda` List of named numeric vectors. One per posterior sample.
  Record speciation rates per tip present at `focal_time`. Updated if
  `update_rates = TRUE`.

- `$tipMu` List of named numeric vectors. One per posterior sample.
  Record extinction rates per tip present at `focal_time`. Updated if
  `update_rates = TRUE`.

- `$eventBranchSegs` List of numeric matrices. One per posterior sample.
  Record regime ID per segment of branches.

- `$meanTipLambda` Named numeric vector. Mean tip speciation rates
  across all posterior configurations of tips present at `focal_time`
  (does not include older fossils).

- `$meanTipMu` Named numeric vector. Mean tip extinction rates across
  all posterior configurations of tips present at `focal_time` (does not
  include older fossils).

- `$type` Character string. Set the type of data modeled with BAMM.
  Should be "diversification".

Additional elements providing key information for downstream analyses:

- `$expectedNumberOfShifts` Integer. The expected number of regime
  shifts used to set the prior in BAMM.

- `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating
  information from the Phylogeny-related elements above, except
  `$MSP_tree$edge.length` is recording the Marginal Shift Probability of
  each branch (i.e., the probability of a regime shift to occur along
  each branch) whose origin is older than `focal_time`.

- `$MAP_indices` Integer vector. The indices of the Maximum A Posteriori
  probability (MAP) configurations among the posterior samples.

- `$MAP_BAMM_object`. List of 18 elements of class `"bammdata"`
  recording the mean rates and regime shift locations found across the
  Maximum A Posteriori probability (MAP) configuration. All BAMM
  elements summarizing diversification data hold a single entry
  describing this the mean diversification history, updated for the
  `focal_time`.

- `$MSC_indices` Integer vector. The indices of the Maximum Shift
  Credibility (MSC) configurations among the posterior samples.

- `$MSC_BAMM_object` List of 18 elements of class `"bammdata"` recording
  the mean rates and regime shift locations found across the Maximum
  Shift Credibility (MSC) configurations. All BAMM elements summarizing
  diversification data hold a single entry describing this mean
  diversification history, updated for the `focal_time`.

New elements added to provide update information:

- `$root_age` Integer. Stores the age of the root of the tree.

- `$nodes_ID_df` Data.frame with two columns. Provides the conversion
  from the `new_node_ID` to the `initial_node_ID`. Each row is a node.

- `$initial_nodes_ID` Character vector. Provides the initial ID of
  internal nodes. Used to plot internal node IDs as labels with
  [`ape::nodelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html).

- `$edges_ID_df` Data.frame with two columns. Provides the conversion
  from the `new_edge_ID` to the `initial_edge_ID`. Each row is an
  edge/branch.

- `$initial_edges_ID` Character vector. Provides the initial ID of
  edges/branches. Used to plot edge/branch IDs as labels with
  [`ape::edgelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html).

- `$dtrates` List of three elements.

  - 1/ `$dtrates$tau` Numeric. Resolution factor describing the fraction
    of each segment length used in
    [`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
    compared to the full depth of the initial tree (i.e., the root_age)

  - 2/ `$dtrates$rates` List of two numeric vectors. Speciation and
    extinction rates along segments used by
    [`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md).

  - 3/ `$dtrates$tmat` Numeric matrix. Start and end times of segments
    in term of distance to the root.

- `$initial_colorbreaks` List of three numeric vectors. Rate values of
  the percentiles delimiting the bins for mapping rates to colors with
  [`BAMMtools::plot.bammdata()`](https://rdrr.io/pkg/BAMMtools/man/plot.html).
  Each element provides values for different type of rates
  (`$speciation`, `$extinction`, `$net_diversification`).

- `$focal_time` Integer. The time, in terms of time distance from the
  present, at which the rates/regimes were extracted and the tree was
  possibly cut.

## Details

The object of class `"bammdata"` (`BAMM_object`) is cut-off at a
specific time in the past (i.e. the `focal_time`) and the current
diversification rate values of the overlapping edges/branches are
extracted.

—– Update diversification rate data —–

If `update_rates = TRUE`, diversification rates of branches overlapping
with `focal_time` will be updated. Each cut-off branch forms a new tip
present at `focal_time` with updated associated diversification rate
data. Fossils older than `focal_time` do not yield any data.

- `$tipLambda` contains speciation rates.

- `$tipMu` contains extinction rates.

If `update_regimes = TRUE`, diversification regimes of branches
overlapping with `focal_time` will be updated. Each cut-off branch forms
a new tip present at `focal_time` with updated associated
diversification regime ID found in `$tipStates`. Fossils older than
`focal_time` do not yield any data.

—– Update the phylogeny —–

If `update_tree = TRUE`, elements defining the phylogeny, as in a
`"phylo"` object will be updated such that all branches that are younger
than the `focal_time` are cut-off:

- `$edge` defines the tree topology.

- `$Nnode` defines the number of internal nodes.

- `$tip.label` provides the labels of all tips, including fossils older
  than `focal_time` if present.

- `$edge.length` provides length of all branches.

- `$node.label` provides the labels of all internal nodes. (Optional)

—– Update the plot from
[`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
—–

If `update_plot = TRUE`, elements used to plot diversification rates on
the phylogeny using
[`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
will be updated such that all branches that are younger than the
`focal_time` are cut-off:

- `$begin` provides absolute time since root of edge/branch start
  (rootward).

- `$end` provides absolute time since root of edge/branch end (tipward).

- `$eventVectors` provides regime membership per branch in each
  posterior sample configuration.

- `$eventBranchSegs` provides regime membership per segment of branches
  in each posterior sample configuration.

- `$dtrates` provides mean speciation and extinction rates along
  segments of branches, and resolution fraction (tau) describing the
  fraction of each segment length compared to the full depth of the
  initial tree (i.e., the root_age).

## See also

[`cut_phylo_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_phylo_for_focal_time.md)
[`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)

## Author

Maël Doré

## Examples

``` r
# ----- Example 1: Extant whales (87 taxa) ----- #

## Load the BAMM_object summarizing 1000 posterior samples of BAMM
data(whale_BAMM_object, package = "deepSTRAPP")

## Set focal-time to 5 My
focal_time = 5

 # (May take several minutes to run)
## Update the BAMM object
whale_BAMM_object_5My <- update_rates_and_regimes_for_focal_time(
   BAMM_object = whale_BAMM_object,
   focal_time = 5,
   update_rates = TRUE, update_regimes = TRUE,
   update_tree = TRUE, update_plot = TRUE,
   update_all_elements = TRUE,
   keep_tip_labels = TRUE,
   verbose = TRUE)
#> 2026-09-09 06:38:30.831165 - Tip regimes/rates updated for BAMM posterior sample n°100/1000
#> 2026-09-09 06:38:30.861592 - Tip regimes/rates updated for BAMM posterior sample n°200/1000
#> 2026-09-09 06:38:30.892336 - Tip regimes/rates updated for BAMM posterior sample n°300/1000
#> 2026-09-09 06:38:30.953043 - Tip regimes/rates updated for BAMM posterior sample n°400/1000
#> 2026-09-09 06:38:30.980129 - Tip regimes/rates updated for BAMM posterior sample n°500/1000
#> 2026-09-09 06:38:31.021999 - Tip regimes/rates updated for BAMM posterior sample n°600/1000
#> 2026-09-09 06:38:31.047008 - Tip regimes/rates updated for BAMM posterior sample n°700/1000
#> 2026-09-09 06:38:31.070095 - Tip regimes/rates updated for BAMM posterior sample n°800/1000
#> 2026-09-09 06:38:31.093395 - Tip regimes/rates updated for BAMM posterior sample n°900/1000
#> 2026-09-09 06:38:31.118679 - Tip regimes/rates updated for BAMM posterior sample n°1000/1000

# Add "phylo" class to be compatible with phytools::nodeHeights()
class(whale_BAMM_object) <- unique(c(class(whale_BAMM_object), "phylo"))
root_age <- max(phytools::nodeHeights(whale_BAMM_object)[,2])
# Remove temporary "phylo" class
class(whale_BAMM_object) <- setdiff(class(whale_BAMM_object), "phylo")

## Plot initial BAMM_object for t = 0 My
plot_BAMM_rates(whale_BAMM_object, add_regime_shifts = TRUE,
                labels = TRUE, legend = TRUE, cex = 0.5)
abline(v = root_age - focal_time,
      col = "red", lty = 2, lwd = 2)


## Plot updated BAMM_object for t = 5 My
plot_BAMM_rates(whale_BAMM_object_5My, add_regime_shifts = TRUE,
                labels = TRUE, legend = TRUE, cex = 0.8) 


# ----- Example 2: Extant Ponerinae (1,534 taxa) ----- #

if (deepSTRAPP::is_dev_version())
{
 ## Load the BAMM_object summarizing 1000 posterior samples of BAMM
 data(Ponerinae_BAMM_object, package = "deepSTRAPP")
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 ## Set focal-time to 10 My
 focal_time = 10

  # (May take several minutes to run)
 ## Update the BAMM object
 Ponerinae_BAMM_object_10My <- update_rates_and_regimes_for_focal_time(
    BAMM_object = Ponerinae_BAMM_object,
    focal_time = focal_time,
    update_rates = TRUE, update_regimes = TRUE,
    update_tree = TRUE, update_plot = TRUE,
    update_all_elements = TRUE,
    keep_tip_labels = TRUE,
    verbose = TRUE) 

 ## Load results to save time
 data(Ponerinae_BAMM_object_10My, package = "deepSTRAPP")
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 ## Extract root age
 # Add temporarily the "phylo" class to be compatible with phytools::nodeHeights()
 class(Ponerinae_BAMM_object) <- unique(c(class(Ponerinae_BAMM_object), "phylo"))
 root_age <- max(phytools::nodeHeights(Ponerinae_BAMM_object)[,2])
 # Remove temporary "phylo" class
 class(Ponerinae_BAMM_object) <- setdiff(class(Ponerinae_BAMM_object), "phylo")

 ## Plot diversification rates on the initial tree
 plot_BAMM_rates(Ponerinae_BAMM_object,
                 legend = TRUE, labels = FALSE)
 abline(v = root_age - focal_time,
        col = "red", lty = 2, lwd = 2)

 ## Plot diversification rates and regime shifts on the updated tree (cut-off for 10 My)
 # Keep the initial color scheme
 plot_BAMM_rates(Ponerinae_BAMM_object_10My, legend = TRUE, labels = FALSE,
                 colorbreaks = Ponerinae_BAMM_object_10My$initial_colorbreaks$net_diversification)

 # Use a new color scheme mapped on the new distribution of rates
 plot_BAMM_rates(Ponerinae_BAMM_object_10My, legend = TRUE, labels = FALSE)
}
```
