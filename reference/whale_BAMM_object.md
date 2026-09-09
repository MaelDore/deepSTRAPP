# Dataset summarizing 1000 posterior samples of BAMM for extant whales

An object of class `"bammdata"` containing information on
diversification dynamics of extant whales (Cetacea order) modeled with
BAMM.

Source: Steeman, M. E., M. B. Hebsgaard, R. E. Fordyce, S. Y. W. Ho, D.
L. Rabosky, R. Nielsen, C. Rahbek, H. Glenner, M. V. Sorensen, and E.
Willerslev (2009) Radiation of extant cetaceans driven by restructuring
of the oceans. Systematic Biology, 58, 573-585.

## Usage

``` r
data(whale_BAMM_object)
```

## Format

A list with 24 elements.

## Details

An object of class `"bammdata"` containing information on
diversification dynamics of extant whales (Cetacea order) modeled with
BAMM.

Phylogeny-related elements used to plot a phylogeny with
[`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html):

- `$edge` Integer matrix. Defines the tree topology by providing
  rootward and tipward node ID of each edge.

- `$Nnode` Integer. Number of internal nodes.

- `$tip.label` Character vector. Labels of all tips.

- `$edge.length` Numeric vector. Length of edges/branches.

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

- `$tipStates` List of integer vectors. One per posterior sample. Record
  regime ID per tip.

- `$tipLambda` List of numeric vectors. One per posterior sample. Record
  speciation rates per tip.

- `$tipMu` List of numeric vectors. One per posterior sample. Record
  extinction rates per tip.

- `$eventBranchSegs` List of numeric matrices. One per posterior sample.
  Record regime ID per segment of branches.

- `$meanTipLambda` Numeric vector. Mean tip speciation rates across all
  posterior configurations of tips.

- `$meanTipMu` Numeric vector. Mean tip extinction rates across all
  posterior configurations of tips.

- `$type` Character string. Set the type of data modeled with BAMM.
  Here, type = "diversification".

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

Steeman, M. E., M. B. Hebsgaard, R. E. Fordyce, S. Y. W. Ho, D. L.
Rabosky, R. Nielsen, C. Rahbek, H. Glenner, M. V. Sorensen, and E.
Willerslev (2009) Radiation of extant cetaceans driven by restructuring
of the oceans. Systematic Biology, 58, 573-585.

## See also

BAMM software website: <http://bamm-project.org/>
