# Convert simmaps into densityMaps suitable for deepSTRAPP

Convert `simmaps` representing independent simulations of evolutionary
history into a `densityMaps` object suitable to use as input for a
deepSTRAPP run.

`simmaps` are mapping discrete character/geographic evolutionary history
as transitions in character states/geographic ranges across nodes and
branches of a phylogeny.

`densityMaps` are recording the posterior probability/frequencies of
being in a given state/range along branches. In deepSTRAPP format,
`densityMaps` is a list of objects of class `"densityMap"` where each
object corresponds to the frequency of presence/absence of a
state/range.

## Usage

``` r
convert_simmaps_to_densityMaps(
  simmaps,
  colors_per_levels = NULL,
  tol = 1e-05,
  verbose = TRUE
)
```

## Arguments

- simmaps:

  List of objects of class `"simmap"`, typically generated with
  [`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md)
  or
  [`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html),
  that represent discrete character/geographic evolutionary history
  (i.e., transitions in character states/geographic ranges) mapped along
  branches.

- colors_per_levels:

  Named character string. To set the colors to use to map each
  state/range posterior probabilities. Names = states/ranges; values =
  colors. If `NULL` (default), the
  [`rainbow()`](https://rdrr.io/r/grDevices/palettes.html) color scale
  will be used.

- tol:

  Positive numeric. To set the tolerance used to match node ages and
  time steps (i.e., consider them equal). Default = 1e-5.

- verbose:

  Logical. Whether to display progress every 100 edges. Default =
  `TRUE`.

## Value

The function returns a `densityMaps` as a list of objects of class
`"densityMap"`, where each object is mapping the frequency of
presence/absence of a state/range.

The number of objects depends on the number of states/ranges observed
across the `simmaps` provided for conversion. Each `densityMap` is named
as "Density_map_X" with X being each of the state/range recorded across
the `simmaps`.

Each `densityMap` is a list of three elements:

- `$tree` List of classes `"simmap"` and `"phylo"` that contains the
  phylogeny in [ape](https://rdrr.io/pkg/ape/man/ape-package.html)
  format and the mapping of states/ranges frequencies along branches in
  `$tree$maps`.

- `$cols` Named character strings. Colors mapped to the 0 to 1000 scale
  used to record frequencies in `$tree$maps`.

- `$states` Character string with two values. First entry is the absence
  of the state/range recorded as "Not X". Second entry is the presence
  of the state/range recorded as "X", X being the state/range name.

## Details

The function is a wrapper of the original
[`phytools::densityMap()`](https://rdrr.io/pkg/phytools/man/densityMap.html)
by Liam Revell. However, it does not produce a single `densityMap` for
binary states, but rather can handle `simmaps` mapping any number of
states/ranges, and produce a `densityMaps` list that summarizes the
frequency of presence/absence of each state/range in subsequent
`densityMap` objects.

Both `simmaps` and `densityMaps` objects can be used as inputs for a
deepSTRAPP run with
[`run_deepSTRAPP_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_for_focal_time.md)
or
[`run_deepSTRAPP_over_time()`](https://maeldore.github.io/deepSTRAPP/reference/run_deepSTRAPP_over_time.md).

`simmaps` retain the identity of each simulated history, allowing
deepSTRAPP to keep track of which simulation generated each set of trait
values (i.e., states/ranges). However, retaining all simulated histories
can require substantial RAM, particularly for large phylogenies or many
simulations.

Alternatively, `densityMaps` summarize the frequency of states/ranges
across simulations. They require substantially less RAM and can be used
to visualize the overall uncertainty in trait evolution with
[`plot_densityMaps_overlay()`](https://maeldore.github.io/deepSTRAPP/reference/plot_densityMaps_overlay.md).

The trade-off is that `densityMaps` discard the identity of individual
simulations and therefore cannot be used to track which simulated
history generated a given set of trait values.

## See also

[`phytools::densityMap()`](https://rdrr.io/pkg/phytools/man/densityMap.html)
[`plot_densityMaps_overlay()`](https://maeldore.github.io/deepSTRAPP/reference/plot_densityMaps_overlay.md)
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md)

## Author

Maël Doré

## Examples

``` r
## Load data

# Load trait df
data(Ponerinae_trait_tip_data, package = "deepSTRAPP")
# Load phylogeny
data(Ponerinae_tree_old_calib, package = "deepSTRAPP")

# Extract categorical data with 3-levels
Ponerinae_cat_3lvl_tip_data <- setNames(object = Ponerinae_trait_tip_data$fake_cat_3lvl_tip_data,
                                         nm = Ponerinae_trait_tip_data$Taxa)
table(Ponerinae_cat_3lvl_tip_data)
#> Ponerinae_cat_3lvl_tip_data
#>     arboreal subterranean  terricolous 
#>          193          953          388 

# Select color scheme for states
colors_per_states <- c("forestgreen", "sienna", "goldenrod")
names(colors_per_states) <- c("arboreal", "subterranean", "terricolous")

 # (May take several minutes to run)
## Produce densityMaps using stochastic character mapping based on an ER Mk model
Ponerinae_cat_3lvl_data_old_calib <- prepare_trait_data(
    tip_data = Ponerinae_cat_3lvl_tip_data,
    phylo = Ponerinae_tree_old_calib,
    trait_data_type = "categorical",
    colors_per_levels = colors_per_states,
    evolutionary_models = "ER", # Use ER model
    run_stochastic_maps = TRUE,
    nb_simulations = 100, # Reduce number of simulations to save time
    seed = 1234, # Set seed for reproducibility
    return_simmaps = TRUE, # Return simmaps in the output
    plot_map = FALSE)
#> Warning: Entries in 'tip_data' were reordered to match 'phylo$tip.label.
#> 
#> 2026-09-09 04:26:48.152448 - Fit 1 evolutionary model(s): ER.
#> 
#> ------ ER model ------ 
#> 
#> GEIGER-fitted comparative model of discrete data
#>  fitted Q matrix:
#>                       arboreal  subterranean   terricolous
#>     arboreal     -0.0002508455  0.0001254228  0.0001254228
#>     subterranean  0.0001254228 -0.0002508455  0.0001254228
#>     terricolous   0.0001254228  0.0001254228 -0.0002508455
#> 
#>  model summary:
#>  log-likelihood = -45.665085
#>  AIC = 93.330171
#>  AICc = 93.332782
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
#> 2026-09-09 04:26:52.862327 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-09-09 04:26:52.863791 - Run simulations for stochastic mapping.
#> 
#> make.simmap is sampling character histories conditioned on
#> the transition matrix
#> 
#> Q =
#>                   arboreal  subterranean   terricolous
#> arboreal     -0.0002508455  0.0001254228  0.0001254228
#> subterranean  0.0001254228 -0.0002508455  0.0001254228
#> terricolous   0.0001254228  0.0001254228 -0.0002508455
#> (specified by the user);
#> and (mean) root node prior probabilities
#> pi =
#>     arboreal subterranean  terricolous 
#>    0.3333333    0.3333333    0.3333333 
#> Done.
#> 2026-09-09 04:29:14.476297 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-09-09 04:29:17.962216 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-09-09 04:29:21.796884 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 04:29:22.887474 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 04:29:24.183585 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 04:29:26.557785 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 04:29:27.448889 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 04:29:28.297863 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 04:29:28.73746 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 04:29:29.076267 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 04:29:29.457887 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 04:29:29.857433 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 04:29:30.244592 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 04:29:30.951986 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 04:29:31.641279 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 04:29:32.715061 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 04:29:33.410419 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 04:29:34.635818 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 04:29:35.814981 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 04:29:36.683394 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 04:29:37.350288 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 04:29:37.908486 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 04:29:38.656933 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 04:29:39.534876 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 04:29:40.62518 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 04:29:41.373658 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 04:29:42.07622 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 04:29:42.937977 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 04:29:44.194309 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 04:29:46.130086 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 04:29:48.362587 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 04:29:50.290279 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 04:29:51.481278 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-09 04:29:54.777461 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 04:29:55.713662 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 04:29:56.840185 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 04:29:59.015032 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 04:29:59.794653 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 04:30:00.119693 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 04:30:00.543216 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 04:30:00.8904 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 04:30:01.278094 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 04:30:01.682271 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 04:30:02.099037 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 04:30:02.825659 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 04:30:03.535182 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 04:30:04.673845 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 04:30:05.387501 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 04:30:06.676788 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 04:30:07.883584 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 04:30:08.759679 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 04:30:09.450129 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 04:30:10.034208 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 04:30:11.354097 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 04:30:12.214719 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 04:30:13.272962 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 04:30:14.004575 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 04:30:14.694465 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 04:30:15.551386 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 04:30:16.766581 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 04:30:18.664914 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 04:30:20.883614 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 04:30:22.774385 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 04:30:23.965978 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-09 04:30:27.318983 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 04:30:28.253654 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 04:30:29.38147 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 04:30:31.441 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 04:30:32.270909 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 04:30:32.597608 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 04:30:33.002862 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 04:30:33.345993 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 04:30:33.732288 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 04:30:34.153967 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 04:30:34.550367 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 04:30:35.263811 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 04:30:35.983101 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 04:30:37.090518 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 04:30:37.806147 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 04:30:39.060283 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 04:30:40.298467 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 04:30:41.178459 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 04:30:41.880349 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 04:30:42.459701 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 04:30:43.225005 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 04:30:44.124917 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 04:30:45.248963 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 04:30:46.027084 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 04:30:46.746383 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 04:30:47.639438 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 04:30:48.911504 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 04:30:51.451405 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 04:30:53.651424 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 04:30:55.54285 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 04:30:56.712492 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-09-09 04:31:00.121429 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 04:31:01.05538 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 04:31:02.194092 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 04:31:04.245768 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 04:31:05.058773 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 04:31:05.378841 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 04:31:05.780712 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 04:31:06.124419 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 04:31:06.526803 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 04:31:06.927899 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 04:31:07.318288 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 04:31:08.042275 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 04:31:08.737472 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 04:31:09.841002 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 04:31:10.553288 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 04:31:11.811564 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 04:31:13.002347 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 04:31:13.877224 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 04:31:14.56646 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 04:31:15.148361 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 04:31:15.927269 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 04:31:16.835731 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 04:31:17.952231 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 04:31:18.718192 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 04:31:19.463234 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 04:31:20.358078 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 04:31:21.619105 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 04:31:23.628742 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 04:31:25.926912 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 04:31:27.918798 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 04:31:29.721688 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-09 04:31:32.985353 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 04:31:33.93558 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 04:31:35.056073 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 04:31:37.232116 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 04:31:38.006142 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 04:31:38.325315 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 04:31:38.745271 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 04:31:39.083781 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 04:31:39.465607 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 04:31:39.880438 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 04:31:40.271204 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 04:31:40.973691 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 04:31:41.681031 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 04:31:42.79022 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 04:31:43.512988 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 04:31:44.75311 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 04:31:45.944374 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 04:31:46.806015 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 04:31:47.50415 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 04:31:48.080862 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 04:31:48.841028 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 04:31:49.737963 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 04:31:50.861817 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 04:31:51.643462 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 04:31:52.363903 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 04:31:53.25662 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 04:31:54.533273 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 04:31:56.52481 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 04:31:58.8132 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 04:32:00.804686 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 04:32:02.052824 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-09 04:32:05.626183 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 04:32:07.30364 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 04:32:08.356117 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 04:32:10.286388 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 04:32:11.026866 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 04:32:11.333244 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 04:32:11.720615 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 04:32:12.046918 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 04:32:12.415906 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 04:32:12.799484 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 04:32:13.171692 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 04:32:13.870986 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 04:32:14.553705 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 04:32:15.577615 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 04:32:16.251347 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 04:32:17.434464 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 04:32:18.542378 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 04:32:19.351045 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 04:32:19.995279 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 04:32:20.542105 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 04:32:21.269897 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 04:32:22.115719 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 04:32:23.152823 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 04:32:23.874856 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 04:32:24.579886 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 04:32:25.404637 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 04:32:26.58552 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 04:32:28.453143 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 04:32:30.571316 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 04:32:32.425801 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 04:32:33.559677 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
