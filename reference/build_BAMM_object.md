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

# The 'whale_event_data.txt' file used for example here is provided within deepSTRAPP in the 'extdata' directory

library(phytools)
#> Loading required package: ape
#> Loading required package: maps
data(whale.tree)

BAMM_object <- build_BAMM_object(
   phylo = whale.tree,
   eventdata = system.file("extdata", "whale_event_data.txt", package = "deepSTRAPP"),
   burn_in = 0.25, # Remove 25% as burn-in
   nb_posterior_samples = 1000, # Retain 1000 samples
   expectedNumberOfShifts = 1,
   verbose = TRUE)
#> Reading event datafile:  /home/runner/.cache/R/renv/library/deepSTRAPP-c9e56a57/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu/deepSTRAPP/extdata/whale_event_data.txt 
#>      ...........
#> Read a total of 2000 samples from posterior
#> 
#> Discarded as burnin: GENERATIONS <  24950
#> Analyzing  1501  samples from posterior
#> 
#> Setting recursive sequence on tree...
#> 
#> Done with recursive sequence
#> 
#> Start preprocessing MRCA pairs....
#> Done preprocessing MRCA pairs....
#> Processing event:  1 
#> Processing event:  2 
#> Processing event:  3 
#> Processing event:  4 
#> Processing event:  5 
#> Processing event:  6 
#> Processing event:  7 
#> Processing event:  8 
#> Processing event:  9 
#> Processing event:  10 
#> Processing event:  11 
#> Processing event:  12 
#> Processing event:  13 
#> Processing event:  14 
#> Processing event:  15 
#> Processing event:  16 
#> Processing event:  17 
#> Processing event:  18 
#> Processing event:  19 
#> Processing event:  20 
#> Processing event:  21 
#> Processing event:  22 
#> Processing event:  23 
#> Processing event:  24 
#> Processing event:  25 
#> Processing event:  26 
#> Processing event:  27 
#> Processing event:  28 
#> Processing event:  29 
#> Processing event:  30 
#> Processing event:  31 
#> Processing event:  32 
#> Processing event:  33 
#> Processing event:  34 
#> Processing event:  35 
#> Processing event:  36 
#> Processing event:  37 
#> Processing event:  38 
#> Processing event:  39 
#> Processing event:  40 
#> Processing event:  41 
#> Processing event:  42 
#> Processing event:  43 
#> Processing event:  44 
#> Processing event:  45 
#> Processing event:  46 
#> Processing event:  47 
#> Processing event:  48 
#> Processing event:  49 
#> Processing event:  50 
#> Processing event:  51 
#> Processing event:  52 
#> Processing event:  53 
#> Processing event:  54 
#> Processing event:  55 
#> Processing event:  56 
#> Processing event:  57 
#> Processing event:  58 
#> Processing event:  59 
#> Processing event:  60 
#> Processing event:  61 
#> Processing event:  62 
#> Processing event:  63 
#> Processing event:  64 
#> Processing event:  65 
#> Processing event:  66 
#> Processing event:  67 
#> Processing event:  68 
#> Processing event:  69 
#> Processing event:  70 
#> Processing event:  71 
#> Processing event:  72 
#> Processing event:  73 
#> Processing event:  74 
#> Processing event:  75 
#> Processing event:  76 
#> Processing event:  77 
#> Processing event:  78 
#> Processing event:  79 
#> Processing event:  80 
#> Processing event:  81 
#> Processing event:  82 
#> Processing event:  83 
#> Processing event:  84 
#> Processing event:  85 
#> Processing event:  86 
#> Processing event:  87 
#> Processing event:  88 
#> Processing event:  89 
#> Processing event:  90 
#> Processing event:  91 
#> Processing event:  92 
#> Processing event:  93 
#> Processing event:  94 
#> Processing event:  95 
#> Processing event:  96 
#> Processing event:  97 
#> Processing event:  98 
#> Processing event:  99 
#> Processing event:  100 
#> Processing event:  101 
#> Processing event:  102 
#> Processing event:  103 
#> Processing event:  104 
#> Processing event:  105 
#> Processing event:  106 
#> Processing event:  107 
#> Processing event:  108 
#> Processing event:  109 
#> Processing event:  110 
#> Processing event:  111 
#> Processing event:  112 
#> Processing event:  113 
#> Processing event:  114 
#> Processing event:  115 
#> Processing event:  116 
#> Processing event:  117 
#> Processing event:  118 
#> Processing event:  119 
#> Processing event:  120 
#> Processing event:  121 
#> Processing event:  122 
#> Processing event:  123 
#> Processing event:  124 
#> Processing event:  125 
#> Processing event:  126 
#> Processing event:  127 
#> Processing event:  128 
#> Processing event:  129 
#> Processing event:  130 
#> Processing event:  131 
#> Processing event:  132 
#> Processing event:  133 
#> Processing event:  134 
#> Processing event:  135 
#> Processing event:  136 
#> Processing event:  137 
#> Processing event:  138 
#> Processing event:  139 
#> Processing event:  140 
#> Processing event:  141 
#> Processing event:  142 
#> Processing event:  143 
#> Processing event:  144 
#> Processing event:  145 
#> Processing event:  146 
#> Processing event:  147 
#> Processing event:  148 
#> Processing event:  149 
#> Processing event:  150 
#> Processing event:  151 
#> Processing event:  152 
#> Processing event:  153 
#> Processing event:  154 
#> Processing event:  155 
#> Processing event:  156 
#> Processing event:  157 
#> Processing event:  158 
#> Processing event:  159 
#> Processing event:  160 
#> Processing event:  161 
#> Processing event:  162 
#> Processing event:  163 
#> Processing event:  164 
#> Processing event:  165 
#> Processing event:  166 
#> Processing event:  167 
#> Processing event:  168 
#> Processing event:  169 
#> Processing event:  170 
#> Processing event:  171 
#> Processing event:  172 
#> Processing event:  173 
#> Processing event:  174 
#> Processing event:  175 
#> Processing event:  176 
#> Processing event:  177 
#> Processing event:  178 
#> Processing event:  179 
#> Processing event:  180 
#> Processing event:  181 
#> Processing event:  182 
#> Processing event:  183 
#> Processing event:  184 
#> Processing event:  185 
#> Processing event:  186 
#> Processing event:  187 
#> Processing event:  188 
#> Processing event:  189 
#> Processing event:  190 
#> Processing event:  191 
#> Processing event:  192 
#> Processing event:  193 
#> Processing event:  194 
#> Processing event:  195 
#> Processing event:  196 
#> Processing event:  197 
#> Processing event:  198 
#> Processing event:  199 
#> Processing event:  200 
#> Processing event:  201 
#> Processing event:  202 
#> Processing event:  203 
#> Processing event:  204 
#> Processing event:  205 
#> Processing event:  206 
#> Processing event:  207 
#> Processing event:  208 
#> Processing event:  209 
#> Processing event:  210 
#> Processing event:  211 
#> Processing event:  212 
#> Processing event:  213 
#> Processing event:  214 
#> Processing event:  215 
#> Processing event:  216 
#> Processing event:  217 
#> Processing event:  218 
#> Processing event:  219 
#> Processing event:  220 
#> Processing event:  221 
#> Processing event:  222 
#> Processing event:  223 
#> Processing event:  224 
#> Processing event:  225 
#> Processing event:  226 
#> Processing event:  227 
#> Processing event:  228 
#> Processing event:  229 
#> Processing event:  230 
#> Processing event:  231 
#> Processing event:  232 
#> Processing event:  233 
#> Processing event:  234 
#> Processing event:  235 
#> Processing event:  236 
#> Processing event:  237 
#> Processing event:  238 
#> Processing event:  239 
#> Processing event:  240 
#> Processing event:  241 
#> Processing event:  242 
#> Processing event:  243 
#> Processing event:  244 
#> Processing event:  245 
#> Processing event:  246 
#> Processing event:  247 
#> Processing event:  248 
#> Processing event:  249 
#> Processing event:  250 
#> Processing event:  251 
#> Processing event:  252 
#> Processing event:  253 
#> Processing event:  254 
#> Processing event:  255 
#> Processing event:  256 
#> Processing event:  257 
#> Processing event:  258 
#> Processing event:  259 
#> Processing event:  260 
#> Processing event:  261 
#> Processing event:  262 
#> Processing event:  263 
#> Processing event:  264 
#> Processing event:  265 
#> Processing event:  266 
#> Processing event:  267 
#> Processing event:  268 
#> Processing event:  269 
#> Processing event:  270 
#> Processing event:  271 
#> Processing event:  272 
#> Processing event:  273 
#> Processing event:  274 
#> Processing event:  275 
#> Processing event:  276 
#> Processing event:  277 
#> Processing event:  278 
#> Processing event:  279 
#> Processing event:  280 
#> Processing event:  281 
#> Processing event:  282 
#> Processing event:  283 
#> Processing event:  284 
#> Processing event:  285 
#> Processing event:  286 
#> Processing event:  287 
#> Processing event:  288 
#> Processing event:  289 
#> Processing event:  290 
#> Processing event:  291 
#> Processing event:  292 
#> Processing event:  293 
#> Processing event:  294 
#> Processing event:  295 
#> Processing event:  296 
#> Processing event:  297 
#> Processing event:  298 
#> Processing event:  299 
#> Processing event:  300 
#> Processing event:  301 
#> Processing event:  302 
#> Processing event:  303 
#> Processing event:  304 
#> Processing event:  305 
#> Processing event:  306 
#> Processing event:  307 
#> Processing event:  308 
#> Processing event:  309 
#> Processing event:  310 
#> Processing event:  311 
#> Processing event:  312 
#> Processing event:  313 
#> Processing event:  314 
#> Processing event:  315 
#> Processing event:  316 
#> Processing event:  317 
#> Processing event:  318 
#> Processing event:  319 
#> Processing event:  320 
#> Processing event:  321 
#> Processing event:  322 
#> Processing event:  323 
#> Processing event:  324 
#> Processing event:  325 
#> Processing event:  326 
#> Processing event:  327 
#> Processing event:  328 
#> Processing event:  329 
#> Processing event:  330 
#> Processing event:  331 
#> Processing event:  332 
#> Processing event:  333 
#> Processing event:  334 
#> Processing event:  335 
#> Processing event:  336 
#> Processing event:  337 
#> Processing event:  338 
#> Processing event:  339 
#> Processing event:  340 
#> Processing event:  341 
#> Processing event:  342 
#> Processing event:  343 
#> Processing event:  344 
#> Processing event:  345 
#> Processing event:  346 
#> Processing event:  347 
#> Processing event:  348 
#> Processing event:  349 
#> Processing event:  350 
#> Processing event:  351 
#> Processing event:  352 
#> Processing event:  353 
#> Processing event:  354 
#> Processing event:  355 
#> Processing event:  356 
#> Processing event:  357 
#> Processing event:  358 
#> Processing event:  359 
#> Processing event:  360 
#> Processing event:  361 
#> Processing event:  362 
#> Processing event:  363 
#> Processing event:  364 
#> Processing event:  365 
#> Processing event:  366 
#> Processing event:  367 
#> Processing event:  368 
#> Processing event:  369 
#> Processing event:  370 
#> Processing event:  371 
#> Processing event:  372 
#> Processing event:  373 
#> Processing event:  374 
#> Processing event:  375 
#> Processing event:  376 
#> Processing event:  377 
#> Processing event:  378 
#> Processing event:  379 
#> Processing event:  380 
#> Processing event:  381 
#> Processing event:  382 
#> Processing event:  383 
#> Processing event:  384 
#> Processing event:  385 
#> Processing event:  386 
#> Processing event:  387 
#> Processing event:  388 
#> Processing event:  389 
#> Processing event:  390 
#> Processing event:  391 
#> Processing event:  392 
#> Processing event:  393 
#> Processing event:  394 
#> Processing event:  395 
#> Processing event:  396 
#> Processing event:  397 
#> Processing event:  398 
#> Processing event:  399 
#> Processing event:  400 
#> Processing event:  401 
#> Processing event:  402 
#> Processing event:  403 
#> Processing event:  404 
#> Processing event:  405 
#> Processing event:  406 
#> Processing event:  407 
#> Processing event:  408 
#> Processing event:  409 
#> Processing event:  410 
#> Processing event:  411 
#> Processing event:  412 
#> Processing event:  413 
#> Processing event:  414 
#> Processing event:  415 
#> Processing event:  416 
#> Processing event:  417 
#> Processing event:  418 
#> Processing event:  419 
#> Processing event:  420 
#> Processing event:  421 
#> Processing event:  422 
#> Processing event:  423 
#> Processing event:  424 
#> Processing event:  425 
#> Processing event:  426 
#> Processing event:  427 
#> Processing event:  428 
#> Processing event:  429 
#> Processing event:  430 
#> Processing event:  431 
#> Processing event:  432 
#> Processing event:  433 
#> Processing event:  434 
#> Processing event:  435 
#> Processing event:  436 
#> Processing event:  437 
#> Processing event:  438 
#> Processing event:  439 
#> Processing event:  440 
#> Processing event:  441 
#> Processing event:  442 
#> Processing event:  443 
#> Processing event:  444 
#> Processing event:  445 
#> Processing event:  446 
#> Processing event:  447 
#> Processing event:  448 
#> Processing event:  449 
#> Processing event:  450 
#> Processing event:  451 
#> Processing event:  452 
#> Processing event:  453 
#> Processing event:  454 
#> Processing event:  455 
#> Processing event:  456 
#> Processing event:  457 
#> Processing event:  458 
#> Processing event:  459 
#> Processing event:  460 
#> Processing event:  461 
#> Processing event:  462 
#> Processing event:  463 
#> Processing event:  464 
#> Processing event:  465 
#> Processing event:  466 
#> Processing event:  467 
#> Processing event:  468 
#> Processing event:  469 
#> Processing event:  470 
#> Processing event:  471 
#> Processing event:  472 
#> Processing event:  473 
#> Processing event:  474 
#> Processing event:  475 
#> Processing event:  476 
#> Processing event:  477 
#> Processing event:  478 
#> Processing event:  479 
#> Processing event:  480 
#> Processing event:  481 
#> Processing event:  482 
#> Processing event:  483 
#> Processing event:  484 
#> Processing event:  485 
#> Processing event:  486 
#> Processing event:  487 
#> Processing event:  488 
#> Processing event:  489 
#> Processing event:  490 
#> Processing event:  491 
#> Processing event:  492 
#> Processing event:  493 
#> Processing event:  494 
#> Processing event:  495 
#> Processing event:  496 
#> Processing event:  497 
#> Processing event:  498 
#> Processing event:  499 
#> Processing event:  500 
#> Processing event:  501 
#> Processing event:  502 
#> Processing event:  503 
#> Processing event:  504 
#> Processing event:  505 
#> Processing event:  506 
#> Processing event:  507 
#> Processing event:  508 
#> Processing event:  509 
#> Processing event:  510 
#> Processing event:  511 
#> Processing event:  512 
#> Processing event:  513 
#> Processing event:  514 
#> Processing event:  515 
#> Processing event:  516 
#> Processing event:  517 
#> Processing event:  518 
#> Processing event:  519 
#> Processing event:  520 
#> Processing event:  521 
#> Processing event:  522 
#> Processing event:  523 
#> Processing event:  524 
#> Processing event:  525 
#> Processing event:  526 
#> Processing event:  527 
#> Processing event:  528 
#> Processing event:  529 
#> Processing event:  530 
#> Processing event:  531 
#> Processing event:  532 
#> Processing event:  533 
#> Processing event:  534 
#> Processing event:  535 
#> Processing event:  536 
#> Processing event:  537 
#> Processing event:  538 
#> Processing event:  539 
#> Processing event:  540 
#> Processing event:  541 
#> Processing event:  542 
#> Processing event:  543 
#> Processing event:  544 
#> Processing event:  545 
#> Processing event:  546 
#> Processing event:  547 
#> Processing event:  548 
#> Processing event:  549 
#> Processing event:  550 
#> Processing event:  551 
#> Processing event:  552 
#> Processing event:  553 
#> Processing event:  554 
#> Processing event:  555 
#> Processing event:  556 
#> Processing event:  557 
#> Processing event:  558 
#> Processing event:  559 
#> Processing event:  560 
#> Processing event:  561 
#> Processing event:  562 
#> Processing event:  563 
#> Processing event:  564 
#> Processing event:  565 
#> Processing event:  566 
#> Processing event:  567 
#> Processing event:  568 
#> Processing event:  569 
#> Processing event:  570 
#> Processing event:  571 
#> Processing event:  572 
#> Processing event:  573 
#> Processing event:  574 
#> Processing event:  575 
#> Processing event:  576 
#> Processing event:  577 
#> Processing event:  578 
#> Processing event:  579 
#> Processing event:  580 
#> Processing event:  581 
#> Processing event:  582 
#> Processing event:  583 
#> Processing event:  584 
#> Processing event:  585 
#> Processing event:  586 
#> Processing event:  587 
#> Processing event:  588 
#> Processing event:  589 
#> Processing event:  590 
#> Processing event:  591 
#> Processing event:  592 
#> Processing event:  593 
#> Processing event:  594 
#> Processing event:  595 
#> Processing event:  596 
#> Processing event:  597 
#> Processing event:  598 
#> Processing event:  599 
#> Processing event:  600 
#> Processing event:  601 
#> Processing event:  602 
#> Processing event:  603 
#> Processing event:  604 
#> Processing event:  605 
#> Processing event:  606 
#> Processing event:  607 
#> Processing event:  608 
#> Processing event:  609 
#> Processing event:  610 
#> Processing event:  611 
#> Processing event:  612 
#> Processing event:  613 
#> Processing event:  614 
#> Processing event:  615 
#> Processing event:  616 
#> Processing event:  617 
#> Processing event:  618 
#> Processing event:  619 
#> Processing event:  620 
#> Processing event:  621 
#> Processing event:  622 
#> Processing event:  623 
#> Processing event:  624 
#> Processing event:  625 
#> Processing event:  626 
#> Processing event:  627 
#> Processing event:  628 
#> Processing event:  629 
#> Processing event:  630 
#> Processing event:  631 
#> Processing event:  632 
#> Processing event:  633 
#> Processing event:  634 
#> Processing event:  635 
#> Processing event:  636 
#> Processing event:  637 
#> Processing event:  638 
#> Processing event:  639 
#> Processing event:  640 
#> Processing event:  641 
#> Processing event:  642 
#> Processing event:  643 
#> Processing event:  644 
#> Processing event:  645 
#> Processing event:  646 
#> Processing event:  647 
#> Processing event:  648 
#> Processing event:  649 
#> Processing event:  650 
#> Processing event:  651 
#> Processing event:  652 
#> Processing event:  653 
#> Processing event:  654 
#> Processing event:  655 
#> Processing event:  656 
#> Processing event:  657 
#> Processing event:  658 
#> Processing event:  659 
#> Processing event:  660 
#> Processing event:  661 
#> Processing event:  662 
#> Processing event:  663 
#> Processing event:  664 
#> Processing event:  665 
#> Processing event:  666 
#> Processing event:  667 
#> Processing event:  668 
#> Processing event:  669 
#> Processing event:  670 
#> Processing event:  671 
#> Processing event:  672 
#> Processing event:  673 
#> Processing event:  674 
#> Processing event:  675 
#> Processing event:  676 
#> Processing event:  677 
#> Processing event:  678 
#> Processing event:  679 
#> Processing event:  680 
#> Processing event:  681 
#> Processing event:  682 
#> Processing event:  683 
#> Processing event:  684 
#> Processing event:  685 
#> Processing event:  686 
#> Processing event:  687 
#> Processing event:  688 
#> Processing event:  689 
#> Processing event:  690 
#> Processing event:  691 
#> Processing event:  692 
#> Processing event:  693 
#> Processing event:  694 
#> Processing event:  695 
#> Processing event:  696 
#> Processing event:  697 
#> Processing event:  698 
#> Processing event:  699 
#> Processing event:  700 
#> Processing event:  701 
#> Processing event:  702 
#> Processing event:  703 
#> Processing event:  704 
#> Processing event:  705 
#> Processing event:  706 
#> Processing event:  707 
#> Processing event:  708 
#> Processing event:  709 
#> Processing event:  710 
#> Processing event:  711 
#> Processing event:  712 
#> Processing event:  713 
#> Processing event:  714 
#> Processing event:  715 
#> Processing event:  716 
#> Processing event:  717 
#> Processing event:  718 
#> Processing event:  719 
#> Processing event:  720 
#> Processing event:  721 
#> Processing event:  722 
#> Processing event:  723 
#> Processing event:  724 
#> Processing event:  725 
#> Processing event:  726 
#> Processing event:  727 
#> Processing event:  728 
#> Processing event:  729 
#> Processing event:  730 
#> Processing event:  731 
#> Processing event:  732 
#> Processing event:  733 
#> Processing event:  734 
#> Processing event:  735 
#> Processing event:  736 
#> Processing event:  737 
#> Processing event:  738 
#> Processing event:  739 
#> Processing event:  740 
#> Processing event:  741 
#> Processing event:  742 
#> Processing event:  743 
#> Processing event:  744 
#> Processing event:  745 
#> Processing event:  746 
#> Processing event:  747 
#> Processing event:  748 
#> Processing event:  749 
#> Processing event:  750 
#> Processing event:  751 
#> Processing event:  752 
#> Processing event:  753 
#> Processing event:  754 
#> Processing event:  755 
#> Processing event:  756 
#> Processing event:  757 
#> Processing event:  758 
#> Processing event:  759 
#> Processing event:  760 
#> Processing event:  761 
#> Processing event:  762 
#> Processing event:  763 
#> Processing event:  764 
#> Processing event:  765 
#> Processing event:  766 
#> Processing event:  767 
#> Processing event:  768 
#> Processing event:  769 
#> Processing event:  770 
#> Processing event:  771 
#> Processing event:  772 
#> Processing event:  773 
#> Processing event:  774 
#> Processing event:  775 
#> Processing event:  776 
#> Processing event:  777 
#> Processing event:  778 
#> Processing event:  779 
#> Processing event:  780 
#> Processing event:  781 
#> Processing event:  782 
#> Processing event:  783 
#> Processing event:  784 
#> Processing event:  785 
#> Processing event:  786 
#> Processing event:  787 
#> Processing event:  788 
#> Processing event:  789 
#> Processing event:  790 
#> Processing event:  791 
#> Processing event:  792 
#> Processing event:  793 
#> Processing event:  794 
#> Processing event:  795 
#> Processing event:  796 
#> Processing event:  797 
#> Processing event:  798 
#> Processing event:  799 
#> Processing event:  800 
#> Processing event:  801 
#> Processing event:  802 
#> Processing event:  803 
#> Processing event:  804 
#> Processing event:  805 
#> Processing event:  806 
#> Processing event:  807 
#> Processing event:  808 
#> Processing event:  809 
#> Processing event:  810 
#> Processing event:  811 
#> Processing event:  812 
#> Processing event:  813 
#> Processing event:  814 
#> Processing event:  815 
#> Processing event:  816 
#> Processing event:  817 
#> Processing event:  818 
#> Processing event:  819 
#> Processing event:  820 
#> Processing event:  821 
#> Processing event:  822 
#> Processing event:  823 
#> Processing event:  824 
#> Processing event:  825 
#> Processing event:  826 
#> Processing event:  827 
#> Processing event:  828 
#> Processing event:  829 
#> Processing event:  830 
#> Processing event:  831 
#> Processing event:  832 
#> Processing event:  833 
#> Processing event:  834 
#> Processing event:  835 
#> Processing event:  836 
#> Processing event:  837 
#> Processing event:  838 
#> Processing event:  839 
#> Processing event:  840 
#> Processing event:  841 
#> Processing event:  842 
#> Processing event:  843 
#> Processing event:  844 
#> Processing event:  845 
#> Processing event:  846 
#> Processing event:  847 
#> Processing event:  848 
#> Processing event:  849 
#> Processing event:  850 
#> Processing event:  851 
#> Processing event:  852 
#> Processing event:  853 
#> Processing event:  854 
#> Processing event:  855 
#> Processing event:  856 
#> Processing event:  857 
#> Processing event:  858 
#> Processing event:  859 
#> Processing event:  860 
#> Processing event:  861 
#> Processing event:  862 
#> Processing event:  863 
#> Processing event:  864 
#> Processing event:  865 
#> Processing event:  866 
#> Processing event:  867 
#> Processing event:  868 
#> Processing event:  869 
#> Processing event:  870 
#> Processing event:  871 
#> Processing event:  872 
#> Processing event:  873 
#> Processing event:  874 
#> Processing event:  875 
#> Processing event:  876 
#> Processing event:  877 
#> Processing event:  878 
#> Processing event:  879 
#> Processing event:  880 
#> Processing event:  881 
#> Processing event:  882 
#> Processing event:  883 
#> Processing event:  884 
#> Processing event:  885 
#> Processing event:  886 
#> Processing event:  887 
#> Processing event:  888 
#> Processing event:  889 
#> Processing event:  890 
#> Processing event:  891 
#> Processing event:  892 
#> Processing event:  893 
#> Processing event:  894 
#> Processing event:  895 
#> Processing event:  896 
#> Processing event:  897 
#> Processing event:  898 
#> Processing event:  899 
#> Processing event:  900 
#> Processing event:  901 
#> Processing event:  902 
#> Processing event:  903 
#> Processing event:  904 
#> Processing event:  905 
#> Processing event:  906 
#> Processing event:  907 
#> Processing event:  908 
#> Processing event:  909 
#> Processing event:  910 
#> Processing event:  911 
#> Processing event:  912 
#> Processing event:  913 
#> Processing event:  914 
#> Processing event:  915 
#> Processing event:  916 
#> Processing event:  917 
#> Processing event:  918 
#> Processing event:  919 
#> Processing event:  920 
#> Processing event:  921 
#> Processing event:  922 
#> Processing event:  923 
#> Processing event:  924 
#> Processing event:  925 
#> Processing event:  926 
#> Processing event:  927 
#> Processing event:  928 
#> Processing event:  929 
#> Processing event:  930 
#> Processing event:  931 
#> Processing event:  932 
#> Processing event:  933 
#> Processing event:  934 
#> Processing event:  935 
#> Processing event:  936 
#> Processing event:  937 
#> Processing event:  938 
#> Processing event:  939 
#> Processing event:  940 
#> Processing event:  941 
#> Processing event:  942 
#> Processing event:  943 
#> Processing event:  944 
#> Processing event:  945 
#> Processing event:  946 
#> Processing event:  947 
#> Processing event:  948 
#> Processing event:  949 
#> Processing event:  950 
#> Processing event:  951 
#> Processing event:  952 
#> Processing event:  953 
#> Processing event:  954 
#> Processing event:  955 
#> Processing event:  956 
#> Processing event:  957 
#> Processing event:  958 
#> Processing event:  959 
#> Processing event:  960 
#> Processing event:  961 
#> Processing event:  962 
#> Processing event:  963 
#> Processing event:  964 
#> Processing event:  965 
#> Processing event:  966 
#> Processing event:  967 
#> Processing event:  968 
#> Processing event:  969 
#> Processing event:  970 
#> Processing event:  971 
#> Processing event:  972 
#> Processing event:  973 
#> Processing event:  974 
#> Processing event:  975 
#> Processing event:  976 
#> Processing event:  977 
#> Processing event:  978 
#> Processing event:  979 
#> Processing event:  980 
#> Processing event:  981 
#> Processing event:  982 
#> Processing event:  983 
#> Processing event:  984 
#> Processing event:  985 
#> Processing event:  986 
#> Processing event:  987 
#> Processing event:  988 
#> Processing event:  989 
#> Processing event:  990 
#> Processing event:  991 
#> Processing event:  992 
#> Processing event:  993 
#> Processing event:  994 
#> Processing event:  995 
#> Processing event:  996 
#> Processing event:  997 
#> Processing event:  998 
#> Processing event:  999 
#> Processing event:  1000 
#> Processing event:  1001 
#> Processing event:  1002 
#> Processing event:  1003 
#> Processing event:  1004 
#> Processing event:  1005 
#> Processing event:  1006 
#> Processing event:  1007 
#> Processing event:  1008 
#> Processing event:  1009 
#> Processing event:  1010 
#> Processing event:  1011 
#> Processing event:  1012 
#> Processing event:  1013 
#> Processing event:  1014 
#> Processing event:  1015 
#> Processing event:  1016 
#> Processing event:  1017 
#> Processing event:  1018 
#> Processing event:  1019 
#> Processing event:  1020 
#> Processing event:  1021 
#> Processing event:  1022 
#> Processing event:  1023 
#> Processing event:  1024 
#> Processing event:  1025 
#> Processing event:  1026 
#> Processing event:  1027 
#> Processing event:  1028 
#> Processing event:  1029 
#> Processing event:  1030 
#> Processing event:  1031 
#> Processing event:  1032 
#> Processing event:  1033 
#> Processing event:  1034 
#> Processing event:  1035 
#> Processing event:  1036 
#> Processing event:  1037 
#> Processing event:  1038 
#> Processing event:  1039 
#> Processing event:  1040 
#> Processing event:  1041 
#> Processing event:  1042 
#> Processing event:  1043 
#> Processing event:  1044 
#> Processing event:  1045 
#> Processing event:  1046 
#> Processing event:  1047 
#> Processing event:  1048 
#> Processing event:  1049 
#> Processing event:  1050 
#> Processing event:  1051 
#> Processing event:  1052 
#> Processing event:  1053 
#> Processing event:  1054 
#> Processing event:  1055 
#> Processing event:  1056 
#> Processing event:  1057 
#> Processing event:  1058 
#> Processing event:  1059 
#> Processing event:  1060 
#> Processing event:  1061 
#> Processing event:  1062 
#> Processing event:  1063 
#> Processing event:  1064 
#> Processing event:  1065 
#> Processing event:  1066 
#> Processing event:  1067 
#> Processing event:  1068 
#> Processing event:  1069 
#> Processing event:  1070 
#> Processing event:  1071 
#> Processing event:  1072 
#> Processing event:  1073 
#> Processing event:  1074 
#> Processing event:  1075 
#> Processing event:  1076 
#> Processing event:  1077 
#> Processing event:  1078 
#> Processing event:  1079 
#> Processing event:  1080 
#> Processing event:  1081 
#> Processing event:  1082 
#> Processing event:  1083 
#> Processing event:  1084 
#> Processing event:  1085 
#> Processing event:  1086 
#> Processing event:  1087 
#> Processing event:  1088 
#> Processing event:  1089 
#> Processing event:  1090 
#> Processing event:  1091 
#> Processing event:  1092 
#> Processing event:  1093 
#> Processing event:  1094 
#> Processing event:  1095 
#> Processing event:  1096 
#> Processing event:  1097 
#> Processing event:  1098 
#> Processing event:  1099 
#> Processing event:  1100 
#> Processing event:  1101 
#> Processing event:  1102 
#> Processing event:  1103 
#> Processing event:  1104 
#> Processing event:  1105 
#> Processing event:  1106 
#> Processing event:  1107 
#> Processing event:  1108 
#> Processing event:  1109 
#> Processing event:  1110 
#> Processing event:  1111 
#> Processing event:  1112 
#> Processing event:  1113 
#> Processing event:  1114 
#> Processing event:  1115 
#> Processing event:  1116 
#> Processing event:  1117 
#> Processing event:  1118 
#> Processing event:  1119 
#> Processing event:  1120 
#> Processing event:  1121 
#> Processing event:  1122 
#> Processing event:  1123 
#> Processing event:  1124 
#> Processing event:  1125 
#> Processing event:  1126 
#> Processing event:  1127 
#> Processing event:  1128 
#> Processing event:  1129 
#> Processing event:  1130 
#> Processing event:  1131 
#> Processing event:  1132 
#> Processing event:  1133 
#> Processing event:  1134 
#> Processing event:  1135 
#> Processing event:  1136 
#> Processing event:  1137 
#> Processing event:  1138 
#> Processing event:  1139 
#> Processing event:  1140 
#> Processing event:  1141 
#> Processing event:  1142 
#> Processing event:  1143 
#> Processing event:  1144 
#> Processing event:  1145 
#> Processing event:  1146 
#> Processing event:  1147 
#> Processing event:  1148 
#> Processing event:  1149 
#> Processing event:  1150 
#> Processing event:  1151 
#> Processing event:  1152 
#> Processing event:  1153 
#> Processing event:  1154 
#> Processing event:  1155 
#> Processing event:  1156 
#> Processing event:  1157 
#> Processing event:  1158 
#> Processing event:  1159 
#> Processing event:  1160 
#> Processing event:  1161 
#> Processing event:  1162 
#> Processing event:  1163 
#> Processing event:  1164 
#> Processing event:  1165 
#> Processing event:  1166 
#> Processing event:  1167 
#> Processing event:  1168 
#> Processing event:  1169 
#> Processing event:  1170 
#> Processing event:  1171 
#> Processing event:  1172 
#> Processing event:  1173 
#> Processing event:  1174 
#> Processing event:  1175 
#> Processing event:  1176 
#> Processing event:  1177 
#> Processing event:  1178 
#> Processing event:  1179 
#> Processing event:  1180 
#> Processing event:  1181 
#> Processing event:  1182 
#> Processing event:  1183 
#> Processing event:  1184 
#> Processing event:  1185 
#> Processing event:  1186 
#> Processing event:  1187 
#> Processing event:  1188 
#> Processing event:  1189 
#> Processing event:  1190 
#> Processing event:  1191 
#> Processing event:  1192 
#> Processing event:  1193 
#> Processing event:  1194 
#> Processing event:  1195 
#> Processing event:  1196 
#> Processing event:  1197 
#> Processing event:  1198 
#> Processing event:  1199 
#> Processing event:  1200 
#> Processing event:  1201 
#> Processing event:  1202 
#> Processing event:  1203 
#> Processing event:  1204 
#> Processing event:  1205 
#> Processing event:  1206 
#> Processing event:  1207 
#> Processing event:  1208 
#> Processing event:  1209 
#> Processing event:  1210 
#> Processing event:  1211 
#> Processing event:  1212 
#> Processing event:  1213 
#> Processing event:  1214 
#> Processing event:  1215 
#> Processing event:  1216 
#> Processing event:  1217 
#> Processing event:  1218 
#> Processing event:  1219 
#> Processing event:  1220 
#> Processing event:  1221 
#> Processing event:  1222 
#> Processing event:  1223 
#> Processing event:  1224 
#> Processing event:  1225 
#> Processing event:  1226 
#> Processing event:  1227 
#> Processing event:  1228 
#> Processing event:  1229 
#> Processing event:  1230 
#> Processing event:  1231 
#> Processing event:  1232 
#> Processing event:  1233 
#> Processing event:  1234 
#> Processing event:  1235 
#> Processing event:  1236 
#> Processing event:  1237 
#> Processing event:  1238 
#> Processing event:  1239 
#> Processing event:  1240 
#> Processing event:  1241 
#> Processing event:  1242 
#> Processing event:  1243 
#> Processing event:  1244 
#> Processing event:  1245 
#> Processing event:  1246 
#> Processing event:  1247 
#> Processing event:  1248 
#> Processing event:  1249 
#> Processing event:  1250 
#> Processing event:  1251 
#> Processing event:  1252 
#> Processing event:  1253 
#> Processing event:  1254 
#> Processing event:  1255 
#> Processing event:  1256 
#> Processing event:  1257 
#> Processing event:  1258 
#> Processing event:  1259 
#> Processing event:  1260 
#> Processing event:  1261 
#> Processing event:  1262 
#> Processing event:  1263 
#> Processing event:  1264 
#> Processing event:  1265 
#> Processing event:  1266 
#> Processing event:  1267 
#> Processing event:  1268 
#> Processing event:  1269 
#> Processing event:  1270 
#> Processing event:  1271 
#> Processing event:  1272 
#> Processing event:  1273 
#> Processing event:  1274 
#> Processing event:  1275 
#> Processing event:  1276 
#> Processing event:  1277 
#> Processing event:  1278 
#> Processing event:  1279 
#> Processing event:  1280 
#> Processing event:  1281 
#> Processing event:  1282 
#> Processing event:  1283 
#> Processing event:  1284 
#> Processing event:  1285 
#> Processing event:  1286 
#> Processing event:  1287 
#> Processing event:  1288 
#> Processing event:  1289 
#> Processing event:  1290 
#> Processing event:  1291 
#> Processing event:  1292 
#> Processing event:  1293 
#> Processing event:  1294 
#> Processing event:  1295 
#> Processing event:  1296 
#> Processing event:  1297 
#> Processing event:  1298 
#> Processing event:  1299 
#> Processing event:  1300 
#> Processing event:  1301 
#> Processing event:  1302 
#> Processing event:  1303 
#> Processing event:  1304 
#> Processing event:  1305 
#> Processing event:  1306 
#> Processing event:  1307 
#> Processing event:  1308 
#> Processing event:  1309 
#> Processing event:  1310 
#> Processing event:  1311 
#> Processing event:  1312 
#> Processing event:  1313 
#> Processing event:  1314 
#> Processing event:  1315 
#> Processing event:  1316 
#> Processing event:  1317 
#> Processing event:  1318 
#> Processing event:  1319 
#> Processing event:  1320 
#> Processing event:  1321 
#> Processing event:  1322 
#> Processing event:  1323 
#> Processing event:  1324 
#> Processing event:  1325 
#> Processing event:  1326 
#> Processing event:  1327 
#> Processing event:  1328 
#> Processing event:  1329 
#> Processing event:  1330 
#> Processing event:  1331 
#> Processing event:  1332 
#> Processing event:  1333 
#> Processing event:  1334 
#> Processing event:  1335 
#> Processing event:  1336 
#> Processing event:  1337 
#> Processing event:  1338 
#> Processing event:  1339 
#> Processing event:  1340 
#> Processing event:  1341 
#> Processing event:  1342 
#> Processing event:  1343 
#> Processing event:  1344 
#> Processing event:  1345 
#> Processing event:  1346 
#> Processing event:  1347 
#> Processing event:  1348 
#> Processing event:  1349 
#> Processing event:  1350 
#> Processing event:  1351 
#> Processing event:  1352 
#> Processing event:  1353 
#> Processing event:  1354 
#> Processing event:  1355 
#> Processing event:  1356 
#> Processing event:  1357 
#> Processing event:  1358 
#> Processing event:  1359 
#> Processing event:  1360 
#> Processing event:  1361 
#> Processing event:  1362 
#> Processing event:  1363 
#> Processing event:  1364 
#> Processing event:  1365 
#> Processing event:  1366 
#> Processing event:  1367 
#> Processing event:  1368 
#> Processing event:  1369 
#> Processing event:  1370 
#> Processing event:  1371 
#> Processing event:  1372 
#> Processing event:  1373 
#> Processing event:  1374 
#> Processing event:  1375 
#> Processing event:  1376 
#> Processing event:  1377 
#> Processing event:  1378 
#> Processing event:  1379 
#> Processing event:  1380 
#> Processing event:  1381 
#> Processing event:  1382 
#> Processing event:  1383 
#> Processing event:  1384 
#> Processing event:  1385 
#> Processing event:  1386 
#> Processing event:  1387 
#> Processing event:  1388 
#> Processing event:  1389 
#> Processing event:  1390 
#> Processing event:  1391 
#> Processing event:  1392 
#> Processing event:  1393 
#> Processing event:  1394 
#> Processing event:  1395 
#> Processing event:  1396 
#> Processing event:  1397 
#> Processing event:  1398 
#> Processing event:  1399 
#> Processing event:  1400 
#> Processing event:  1401 
#> Processing event:  1402 
#> Processing event:  1403 
#> Processing event:  1404 
#> Processing event:  1405 
#> Processing event:  1406 
#> Processing event:  1407 
#> Processing event:  1408 
#> Processing event:  1409 
#> Processing event:  1410 
#> Processing event:  1411 
#> Processing event:  1412 
#> Processing event:  1413 
#> Processing event:  1414 
#> Processing event:  1415 
#> Processing event:  1416 
#> Processing event:  1417 
#> Processing event:  1418 
#> Processing event:  1419 
#> Processing event:  1420 
#> Processing event:  1421 
#> Processing event:  1422 
#> Processing event:  1423 
#> Processing event:  1424 
#> Processing event:  1425 
#> Processing event:  1426 
#> Processing event:  1427 
#> Processing event:  1428 
#> Processing event:  1429 
#> Processing event:  1430 
#> Processing event:  1431 
#> Processing event:  1432 
#> Processing event:  1433 
#> Processing event:  1434 
#> Processing event:  1435 
#> Processing event:  1436 
#> Processing event:  1437 
#> Processing event:  1438 
#> Processing event:  1439 
#> Processing event:  1440 
#> Processing event:  1441 
#> Processing event:  1442 
#> Processing event:  1443 
#> Processing event:  1444 
#> Processing event:  1445 
#> Processing event:  1446 
#> Processing event:  1447 
#> Processing event:  1448 
#> Processing event:  1449 
#> Processing event:  1450 
#> Processing event:  1451 
#> Processing event:  1452 
#> Processing event:  1453 
#> Processing event:  1454 
#> Processing event:  1455 
#> Processing event:  1456 
#> Processing event:  1457 
#> Processing event:  1458 
#> Processing event:  1459 
#> Processing event:  1460 
#> Processing event:  1461 
#> Processing event:  1462 
#> Processing event:  1463 
#> Processing event:  1464 
#> Processing event:  1465 
#> Processing event:  1466 
#> Processing event:  1467 
#> Processing event:  1468 
#> Processing event:  1469 
#> Processing event:  1470 
#> Processing event:  1471 
#> Processing event:  1472 
#> Processing event:  1473 
#> Processing event:  1474 
#> Processing event:  1475 
#> Processing event:  1476 
#> Processing event:  1477 
#> Processing event:  1478 
#> Processing event:  1479 
#> Processing event:  1480 
#> Processing event:  1481 
#> Processing event:  1482 
#> Processing event:  1483 
#> Processing event:  1484 
#> Processing event:  1485 
#> Processing event:  1486 
#> Processing event:  1487 
#> Processing event:  1488 
#> Processing event:  1489 
#> Processing event:  1490 
#> Processing event:  1491 
#> Processing event:  1492 
#> Processing event:  1493 
#> Processing event:  1494 
#> Processing event:  1495 
#> Processing event:  1496 
#> Processing event:  1497 
#> Processing event:  1498 
#> Processing event:  1499 
#> Processing event:  1500 
#> Processing event:  1501 
#> Processing event data from data.frame
#> 
#> Discarded as burnin: GENERATIONS <  0
#> Analyzing  1  samples from posterior
#> 
#> Setting recursive sequence on tree...
#> 
#> Done with recursive sequence
#> 
#> Processing event data from data.frame
#> 
#> Discarded as burnin: GENERATIONS <  0
#> Analyzing  1  samples from posterior
#> 
#> Setting recursive sequence on tree...
#> 
#> Done with recursive sequence
#> 
str(BAMM_object, 1)
#> List of 24
#>  $ edge                  : int [1:172, 1:2] 88 89 90 90 91 91 92 92 89 93 ...
#>  $ Nnode                 : int 86
#>  $ tip.label             : chr [1:87] "Balaena_mysticetus" "Eubalaena_australis" "Eubalaena_glacialis" "Eubalaena_japonica" ...
#>  $ edge.length           : num [1:172] 7.59 19.26 8.58 7.3 1.28 ...
#>  $ begin                 : num [1:172] 0 7.59 26.85 26.85 34.15 ...
#>  $ end                   : num [1:172] 7.59 26.85 35.42 34.15 35.42 ...
#>  $ downseq               : int [1:173] 88 89 90 1 91 2 92 3 4 93 ...
#>  $ lastvisit             : int [1:173] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ numberEvents          : int [1:1000] 2 2 2 2 2 2 2 2 2 2 ...
#>  $ eventData             :List of 1000
#>  $ eventVectors          :List of 1000
#>  $ tipStates             :List of 1000
#>  $ tipLambda             :List of 1000
#>  $ tipMu                 :List of 1000
#>  $ eventBranchSegs       :List of 1000
#>  $ meanTipLambda         : Named num [1:87] 0.0618 0.0721 0.0738 0.0738 0.0593 ...
#>   ..- attr(*, "names")= chr [1:87] "Balaena_mysticetus" "Eubalaena_australis" "Eubalaena_glacialis" "Eubalaena_japonica" ...
#>  $ meanTipMu             : Named num [1:87] 0.00253 0.00938 0.01219 0.01219 0 ...
#>   ..- attr(*, "names")= chr [1:87] "Balaena_mysticetus" "Eubalaena_australis" "Eubalaena_glacialis" "Eubalaena_japonica" ...
#>  $ type                  : chr "diversification"
#>  $ expectedNumberOfShifts: num 1
#>  $ MSP_tree              :List of 4
#>   ..- attr(*, "class")= chr "phylo"
#>   ..- attr(*, "order")= chr "cladewise"
#>  $ MAP_indices           : int [1:500] 2 6 9 10 11 12 14 15 20 21 ...
#>  $ MAP_BAMM_object       :List of 18
#>   ..- attr(*, "class")= chr "bammdata"
#>   ..- attr(*, "order")= chr "cladewise"
#>  $ MSC_indices           : int [1:452] 2 6 9 10 11 12 14 15 20 21 ...
#>  $ MSC_BAMM_object       :List of 18
#>   ..- attr(*, "class")= chr "bammdata"
#>   ..- attr(*, "order")= chr "cladewise"
#>  - attr(*, "class")= chr "bammdata"
#>  - attr(*, "order")= chr "cladewise"
```
