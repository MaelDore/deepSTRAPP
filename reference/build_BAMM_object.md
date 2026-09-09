# Build a BAMM object for a deepSTRAPP run

Build a BAMM object of class `bammdata` based on the output file of a
BAMM run that contains a phylogenetic tree and associated
diversification rates mapped along branches across BAMM posterior
samples.

The `BAMM_object` output is typically used as input to run deepSTRAPP
with
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
or
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).

This is a wrapper of the original
[`BAMMtools::getEventData()`](https://rdrr.io/pkg/BAMMtools/man/getEventData.html)
function that additionally provides information on:

- the Marginal Shift Probability (MSP) = the probability of a regime
  shift to occur along each branch.

- the Maximum A Posteriori probability (MAP) configurations among the
  posterior samples = the configurations of regime shifts that were
  sampled most frequently (See
  [`BAMMtools::getBestShiftConfiguration()`](https://rdrr.io/pkg/BAMMtools/man/getBestShiftConfiguration.html)).

- the Maximum Shift Credibility (MSC) configurations among the posterior
  samples = the configurations of regime shift location with the highest
  product of marginal probabilities across branches (See
  [`BAMMtools::maximumShiftCredibility()`](https://rdrr.io/pkg/BAMMtools/man/maximumShiftCredibility.html)).

Those additional elements are used by
[`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
to display regime shift probabilities and locations.

This function is meant to enable users to inject into the deepSTRAPP
framework the results of their own BAMM analyses. Alternatively, a full
BAMM analysis starting from a time-calibrated phylogeny alone can be
carried out within deepSTRAPP with
[`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md).

## Usage

``` r
build_BAMM_object(
  phylo,
  eventdata,
  burn_in = 0.25,
  nb_posterior_samples = NULL,
  seed = NULL,
  expectedNumberOfShifts,
  MAP_odds_ratio_threshold = 5,
  verbose = FALSE
)
```

## Arguments

- phylo:

  Object of class `"phylo"` as defined in R package `{ape}`.
  Time-calibrated phylogeny that was used to produce the BAMM run. The
  phylogeny must be rooted and fully resolved.

- eventdata:

  Character string specifying the path to a BAMM event-data file.
  Alternatively, an object of class data.frame that includes the event
  data from a BAMM run.

- burn_in:

  Numeric. Proportion of posterior samples removed from the BAMM output
  to ensure that the remaining samples were drawn once the equilibrium
  distribution was reached. Default is `0.25`

- nb_posterior_samples:

  Integer. Number of posterior samples to extract, after removing the
  burn-in, in the final `BAMM_object` to use for downstream analyses. If
  set to `NULL` (default), all samples remaining after removing the
  burn-in will be kept.

- seed:

  Integer. Set the seed to ensure reproducibility when drawing random
  posterior samples. Default is `NULL` (a random seed is used).

- expectedNumberOfShifts:

  Integer. The expected number of regime shifts set during the BAMM run
  as a hyperparameter controlling the exponential prior distribution
  used to modulate reversible jumps across model configurations in the
  rjMCMC run. This is needed to compute priors for regime shift along
  branches.

- MAP_odds_ratio_threshold:

  Numeric. Controls the definition of 'core-shifts' used to distinguish
  across configurations when fetching the MAP samples. Shifts that have
  an odds ratio of marginal posterior probability / prior lower than
  `MAP_odds_ratio_threshold` are ignored. See
  [`BAMMtools::getBestShiftConfiguration()`](https://rdrr.io/pkg/BAMMtools/man/getBestShiftConfiguration.html).
  Default = `5`.

- verbose:

  Logical. Whether to display progress in the console. Default =
  `FALSE`.

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

## Note on Bayesian Analysis of Macroevolutionary Mixtures (BAMM)

BAMM is a model of diversification for time-calibrated phylogenies that
explores complex diversification dynamics by allowing multiple regime
shifts across clades without a priori hypotheses on the location of such
shifts. It uses reversible jump Markov chain Monte Carlo (rjMCMC) to
automatically explore a vast range of models with different speciation
and extinction rates, and different number and location of regime
shifts.

BAMM is one option among others for modeling diversification on
phylogenies. You may wish to explore alternative models such as LSBDS
model in RevBayes (Höhna et al., 2016), the MTBD model (Barido-Sottani
et al., 2020), or the ClaDS2 model (Maliet et al., 2019) for your own
data. However, you will need Bayesian models that infer regime shifts to
be able to perform STRAPP tests (Rabosky & Huang, 2016). Additionally,
you need to format the model output such as in `BAMM_object`, so it can
be used in a deepSTRAPP workflow.

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

Initial functions in BAMMtools:
[`BAMMtools::getEventData()`](https://rdrr.io/pkg/BAMMtools/man/getEventData.html)
[`BAMMtools::getBestShiftConfiguration()`](https://rdrr.io/pkg/BAMMtools/man/getBestShiftConfiguration.html)
[`BAMMtools::maximumShiftCredibility()`](https://rdrr.io/pkg/BAMMtools/man/maximumShiftCredibility.html)

Associated functions in deepSTRAPP:
[`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md)
[`plot_BAMM_rates()`](https://maeldore.github.io/deepSTRAPP/reference/plot_BAMM_rates.md)
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md)

For a guided tutorial, see this vignette:
[`vignette("model_diversification_dynamics", package = "deepSTRAPP")`](https://maeldore.github.io/deepSTRAPP/articles/model_diversification_dynamics.md)

## Author

Maël Doré

## Examples

``` r
# The key output from a BAMM is the 'event_data.txt' file
# It can be loaded in R to use as input for deepSTRAPP

library(phytools)
#> Loading required package: ape
#> Loading required package: maps
data(whale.tree)

if (FALSE) { # \dontrun{
## The 'whale_event_data.txt' file used for example here is not provided within deepSTRAPP
BAMM_object <- build_BAMM_object(
   phylo = whale.tree,
   eventdata = "./BAMM_outputs/whale_event_data.txt",
   burn_in = 0.25, # Remove 25% as burn-in
   nb_posterior_samples = 1000, # Retain 1000 samples
   expectedNumberOfShifts = 1,
   verbose = TRUE)
str(BAMM_object, 1)
} # }
```
