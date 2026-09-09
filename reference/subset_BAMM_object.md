# Subset a BAMM object before a deepSTRAPP run

Subset a BAMM object of class `bammdata` while keeping additional
information for deepSTRAPP.

The `BAMM_object` is typically generated directly with
[`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md)
or from external BAMM output files with
[`build_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/build_BAMM_object.md).

This is a wrapper of the original
[`BAMMtools::subsetEventData()`](https://rdrr.io/pkg/BAMMtools/man/subsetEventData.html)
function that additionally preserves information on:

- the Marginal Shift Probability (MSP) = the probability of a regime
  shift to occur along each branch.

- the Maximum A Posteriori probability (MAP) configurations among the
  posterior samples = the configurations of regime shifts that were
  sampled most frequently (See
  [`BAMMtools::getBestShiftConfiguration()`](https://rdrr.io/pkg/BAMMtools/man/getBestShiftConfiguration.html)).

- the Maximum Shift Credibility (MSC) configurations among the posterior
  samples = the configurations of regime shift location with the highest
  product of marginal probabilities across branches (See
  [`BAMMtools::maximumShiftCredibility()`](https://rdrr.io/pkg/BAMMtools/man/maximumShiftCredibility.html))
  that are stored in a `BAMM_object` when built with deepSTRAPP
  functions.

Those additional elements are used by
[`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
to display regime shift probabilities and locations.

## Usage

``` r
subset_BAMM_object(
  BAMM_object,
  nb_posterior_samples = NULL,
  sample_indices = NULL,
  seed = NULL
)
```

## Arguments

- BAMM_object:

  Object of class `"bammdata"`, typically generated with
  [`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md)
  or
  [`build_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/build_BAMM_object.md),
  that contains a phylogenetic tree and associated diversification rate
  mapping across selected posterior samples.

- nb_posterior_samples:

  Integer. Number of posterior samples to extract from `BAMM_object`.
  Default = `NULL`.

- sample_indices:

  Integer or vector of integers. Indices of the posterior samples to
  extract. Default = `NULL`.

- seed:

  Integer. Set the seed to ensure reproducibility when `sample_indices`
  is not provided and posterior samples are drawn randomly. Default is
  `NULL` (a random seed is used).

## Value

The function returns a `BAMM_object` of class `"bammdata"` which is a
list with at least 23 elements.

Phylogeny-related elements used to plot a phylogeny with
[`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html):

- `$edge` Integer matrix. Defines the tree topology by providing
  rootward and tipward node ID of each edge.

- `$Nnode` Integer. Number of internal nodes.

- `$tip.label` Character vector. Labels of all tips.

- `$edge.length` Numeric vector. Length of edges/branches.

- `$node.label` Character vector. Labels of all internal nodes. (Present
  only if present in the initial `phylo`)

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
  Record regime ID per tip.

- `$tipLambda` List of named numeric vectors. One per posterior sample.
  Record speciation rates per tip.

- `$tipMu` List of named numeric vectors. One per posterior sample.
  Record extinction rates per tip.

- `$eventBranchSegs` List of numeric matrices. One per posterior sample.
  Record regime ID per segment of branches.

- `$meanTipLambda` Named numeric vector. Mean tip speciation rates
  across all posterior configurations of tips.

- `$meanTipMu` Named numeric vector. Mean tip extinction rates across
  all posterior configurations of tips.

- `$type` Character string. Set the type of data modeled with BAMM.
  Should be "diversification".

Additional elements providing key information for downstream analyses:

- `$expectedNumberOfShifts` Integer. The expected number of regime
  shifts used to set the prior in BAMM.

- `$MSP_tree` Object of class `phylo`. List of 4 elements duplicating
  information from the Phylogeny-related elements above, except
  `$MSP_tree$edge.length` is recording the Marginal Shift Probability of
  each branch (i.e., the probability of a regime shift to occur along
  each branch)

- `$MAP_indices` Integer vector. The indices of the Maximum A Posteriori
  probability (MAP) configurations among the posterior samples.

- `$MAP_BAMM_object` List of 18 elements of class `"bammdata"` recording
  the mean rates and regime shift locations found across the Maximum A
  Posteriori probability (MAP) configurations. All BAMM elements
  summarizing diversification data hold a single entry describing this
  mean diversification history.

- `$MSC_indices` Integer vector. The indices of the Maximum Shift
  Credibility (MSC) configurations among the posterior samples.

- `$MSC_BAMM_object` List of 18 elements of class `"bammdata"` recording
  the mean rates and regime shift locations found across the Maximum
  Shift Credibility (MSC) configurations. All BAMM elements summarizing
  diversification data hold a single entry describing this mean
  diversification history.

## References

For BAMM: Rabosky, D. L. (2014). Automatic detection of key innovations,
rate shifts, and diversity-dependence on phylogenetic trees. PloS one,
9(2), e89543.
[doi:10.1371/journal.pone.0089543](https://doi.org/10.1371/journal.pone.0089543)
. Website: <http://bamm-project.org/>.

For `{BAMMtools}`: Rabosky, D. L., Grundler, M., Anderson, C., Title,
P., Shi, J. J., Brown, J. W., ... & Larson, J. G. (2014). BAMM tools: an
R package for the analysis of evolutionary dynamics on phylogenetic
trees. Methods in Ecology and Evolution, 5(7), 701-707.
[doi:10.1111/2041-210X.12199](https://doi.org/10.1111/2041-210X.12199)

## See also

Initial function in BAMMtools:
[`BAMMtools::subsetEventData()`](https://rdrr.io/pkg/BAMMtools/man/subsetEventData.html)

Associated functions in deepSTRAPP:
[`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md)
[`build_BAMM_object()`](https://maeldore.github.io/deepSTRAPP/reference/build_BAMM_object.md)

For a guided tutorial, see this vignette:
[`vignette("model_diversification_dynamics", package = "deepSTRAPP")`](https://maeldore.github.io/deepSTRAPP/articles/model_diversification_dynamics.md)

## Author

Maël Doré

## Examples

``` r
if (deepSTRAPP::is_dev_version())
{
 ## Load BAMM object
 # data(Ponerinae_BAMM_object_old_calib)
 # This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

 # Check structure of BAMM_object
 str(Ponerinae_BAMM_object_old_calib, 1)
 # Check current number of BAMM posterior samples
 length(Ponerinae_BAMM_object_old_calib$eventData)
 # We have initially 1000 posterior samples in the updated BAMM object

 ## Subset BAMM_object
 BAMM_object_subset <- subset_BAMM_object(
    BAMM_object = Ponerinae_BAMM_object_old_calib,
    nb_posterior_samples = 100,
    seed = 1234)

 # Check structure of the updated BAMM_object
 str(BAMM_object_subset, 1)
 # Check updated number of BAMM posterior samples
 length(BAMM_object_subset$eventData)
 # We have now 100 posterior samples in the updated BAMM object
}
```
