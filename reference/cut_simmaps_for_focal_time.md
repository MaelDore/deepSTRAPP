# Cut the phylogenies and categorical trait/range mappings of a list of simmaps for a given focal time in the past

Cuts off all the branches of the phylogeny which are younger than a
specific time in the past (i.e. the `focal_time`). Branches overlapping
the `focal_time` are shortened to the `focal_time`. Likewise, remove
continuous trait mapping for the cut off branches by updating the
`$maps` and `$mapped.edge` elements.

## Usage

``` r
cut_simmaps_for_focal_time(simmaps, focal_time, keep_tip_labels = TRUE)
```

## Arguments

- simmaps:

  List of objects with the classes `"phylo"` and `"simmap"`, typically
  generated with
  [`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html),
  that contains a phylogenetic tree and associated categorical
  trait/range mapping. The phylogenetic tree must be rooted and fully
  resolved/dichotomous, but it does not need to be ultrametric (it can
  include fossils).

- focal_time:

  Numeric. The time, in terms of time distance from the present, for
  which the tree and mapping must be cut. It must be smaller than the
  root age of the phylogeny.

- keep_tip_labels:

  Logical. Specify whether terminal branches with a single descendant
  tip must retain their initial `tip.label`. Default is `TRUE`.

## Value

The function returns an updated list of objects as cut/updated simmaps
of classes `"phylo"` and `"simmap"`.

Each simmap object represents an updated stochastic map cut/updated
simmap as a list of at least six elements:

Initial updated elements:

- `$edge` A numeric matrix with two columns listing the updated ID of
  rootward and tipward nodes of the remaining edges

- `$edge.length` A numeric vector. Updated length of remaining edges.

- `$Nnode` Integer. Number of remaining internal nodes.

- `$tip.label` Character vector. Labels of the branches cut at the
  `focal_time`.

  - If `keep_tip_labels = TRUE`, the cut terminal branches are labeled
    with the tip.label of their unique descendant tip.

  - If `keep_tip_labels = FALSE`, the cut terminal branches are labeled
    with the node ID of their unique descendant tip.

- `$maps` An updated list of named numeric vectors. Provides the mapping
  of trait values along each remaining edge.

- `$mapped.edge` An updated matrix. Provides the evolutionary time spent
  across trait values (columns) along the remaining edges (rows).

Additional elements:

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

## Details

The phylogenetic trees are cut for a specific time in the past (i.e. the
`focal_time`).

When a branch with a single descendant tip is cut and
`keep_tip_labels = TRUE`, the leaf left is labeled with the tip.label of
the unique descendant tip.

When a branch with a single descendant tip is cut and
`keep_tip_labels = FALSE`, the leaf left is labeled with the node ID of
the unique descendant tip.

In all cases, when a branch with multiple descendant tips (i.e., a
clade) is cut, the leaf left is labeled with the node ID of the MRCA of
the cut-off clade.

The categorical trait/range mappings are updated accordingly by removing
mapping associated with the cut off branches.

## See also

[`cut_phylo_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/cut_phylo_for_focal_time.md)
[`extract_most_likely_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_most_likely_trait_values_for_focal_time.md)
[`extract_all_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_all_trait_values_for_focal_time.md)

For a guided tutorial, see this vignette:
[`vignette("cut_phylogenies", package = "deepSTRAPP")`](https://maeldore.github.io/deepSTRAPP/articles/cut_phylogenies.md)

## Author

Maël Doré

## Examples

``` r
# ----- Prepare data ----- #

## Load mammals phylogeny and data from the R package motmot, and implemented in deepSTRAPP
# Data source: Slater, 2013; DOI: 10.1111/2041-210X.12084
data("mammals", package = "deepSTRAPP")

# Obtain mammal tree
mammals_tree <- mammals$mammal.phy
# Convert mass data into categories
mammals_mass <- setNames(object = mammals$mammal.mass$mean,
                         nm = row.names(mammals$mammal.mass))[mammals_tree$tip.label]
mammals_data <- mammals_mass
mammals_data[seq_along(mammals_data)] <- "small"
mammals_data[mammals_mass > 5] <- "medium"
mammals_data[mammals_mass > 10] <- "large"
table(mammals_data)
#> mammals_data
#>  large medium  small 
#>     36     83     92 

 # (May take several seconds to run)
# Produce densityMaps using stochastic character mapping based on an equal-rates (ER) Mk model
mammals_cat_data <- prepare_trait_data(tip_data = mammals_data, phylo = mammals_tree,
                                       trait_data_type = "categorical",
                                       evolutionary_models = "ER",
                                       nb_simulations = 100,
                                       return_simmaps = TRUE,
                                       plot_map = FALSE)
#> 
#> 2026-09-09 06:30:42.970124 - Fit 1 evolutionary model(s): ER.
#> 
#> ------ ER model ------ 
#> 
#> GEIGER-fitted comparative model of discrete data
#>  fitted Q matrix:
#>                   large       medium        small
#>     large  -0.006199429  0.003099714  0.003099714
#>     medium  0.003099714 -0.006199429  0.003099714
#>     small   0.003099714  0.003099714 -0.006199429
#> 
#>  model summary:
#>  log-likelihood = -152.749035
#>  AIC = 307.498071
#>  AICc = 307.517210
#>  free parameters = 1
#> 
#> Convergence diagnostics:
#>  optimization iterations = 100
#>  failed iterations = 0
#>  number of iterations with same best fit = 100
#>  frequency of best fit = 1.000
#> 
#>  object summary:
#>  'lik' -- likelihood function
#>  'bnd' -- bounds for likelihood search
#>  'res' -- optimization iteration summary
#>  'opt' -- maximum likelihood parameter estimates
#> 
#> 2026-09-09 06:30:44.711403 - Compare model fits.
#> 
#>    model     logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -152.749 1 307.4981 307.5172          0            100    1
#> 2026-09-09 06:30:44.713203 - Run simulations for stochastic mapping.
#> 
#> make.simmap is sampling character histories conditioned on
#> the transition matrix
#> 
#> Q =
#>               large       medium        small
#> large  -0.006199429  0.003099714  0.003099714
#> medium  0.003099714 -0.006199429  0.003099714
#> small   0.003099714  0.003099714 -0.006199429
#> (specified by the user);
#> and (mean) root node prior probabilities
#> pi =
#>     large    medium     small 
#> 0.3333333 0.3333333 0.3333333 
#> Done.
#> 2026-09-09 06:31:01.580484 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-09-09 06:31:02.183175 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-09-09 06:31:03.654649 - Posterior probability computed for edge n°100/420
#> 2026-09-09 06:31:04.676962 - Posterior probability computed for edge n°200/420
#> 2026-09-09 06:31:05.904813 - Posterior probability computed for edge n°300/420
#> 2026-09-09 06:31:07.164445 - Posterior probability computed for edge n°400/420
#> 2026-09-09 06:31:07.637747 - Posterior probabilities computed for State = large - n°1/3
#> 2026-09-09 06:31:09.113945 - Posterior probability computed for edge n°100/420
#> 2026-09-09 06:31:10.130366 - Posterior probability computed for edge n°200/420
#> 2026-09-09 06:31:11.360672 - Posterior probability computed for edge n°300/420
#> 2026-09-09 06:31:12.741893 - Posterior probability computed for edge n°400/420
#> 2026-09-09 06:31:13.127865 - Posterior probabilities computed for State = medium - n°2/3
#> 2026-09-09 06:31:14.610414 - Posterior probability computed for edge n°100/420
#> 2026-09-09 06:31:15.615315 - Posterior probability computed for edge n°200/420
#> 2026-09-09 06:31:16.907034 - Posterior probability computed for edge n°300/420
#> 2026-09-09 06:31:18.363908 - Posterior probability computed for edge n°400/420
#> 2026-09-09 06:31:18.72267 - Posterior probabilities computed for State = small - n°3/3

# Extract the simulated evolutionary histories (simmaps)
mammals_simmaps <- mammals_cat_data$simmaps

# Set focal time to 80Mya
focal_time <- 80

# ----- Example 1: keep_tip_labels = TRUE ----- #

# Cut the simmap to 80 Mya while keeping tip.label
# on terminal branches with a unique descending tip.
updated_simmaps <- cut_simmaps_for_focal_time(simmaps = mammals_simmaps,
                                              focal_time = focal_time,
                                              keep_tip_labels = TRUE)

# Plot node labels on initial stochastic map n°1 with cut-off
plot(mammals_simmaps[[1]], lwd = 2, fsize = 0.5,
     colors = c(large = "black", medium = "blue", small = "dodgerblue"))
ape::nodelabels(cex = 0.5)
abline(v = max(phytools::nodeHeights(mammals_simmaps[[1]])[,2]) - focal_time,
       col = "red", lty = 2, lwd = 2)


# Plot initial node labels on cut stochastic map n°1
plot(updated_simmaps[[1]], fsize = 0.8,
     colors = c(large = "black", medium = "blue", small = "dodgerblue"))
ape::nodelabels(cex = 0.8, text = updated_simmaps[[1]]$initial_nodes_ID)


# ----- Example 2: keep_tip_labels = FALSE ----- #

# Cut the simmap to 80 Mya while NOT keeping tip.label.
updated_simmaps <- cut_simmaps_for_focal_time(simmaps = mammals_simmaps,
                                              focal_time = focal_time,
                                              keep_tip_labels = FALSE)

# Plot node labels on initial stochastic map n°1 with cut-off
plot(mammals_simmaps[[1]], lwd = 2, fsize = 0.5,
     colors = c(large = "black", medium = "blue", small = "dodgerblue"))
ape::nodelabels(cex = 0.5)
abline(v = max(phytools::nodeHeights(mammals_simmaps[[1]])[,2]) - focal_time,
       col = "red", lty = 2, lwd = 2)


# Plot initial node labels on cut stochastic map n°1
plot(updated_simmaps[[1]], fsize = 0.8,
     colors = c(large = "black", medium = "blue", small = "dodgerblue"))
ape::nodelabels(cex = 0.8, text = updated_simmaps[[1]]$initial_nodes_ID)


```
