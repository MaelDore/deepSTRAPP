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
#> 2026-10-14 09:32:12.743018 - Fit 1 evolutionary model(s): ER.
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
#> 2026-10-14 09:32:16.119423 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-10-14 09:32:16.120503 - Run simulations for stochastic mapping.
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
#> 2026-10-14 09:33:55.379045 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-10-14 09:33:58.284866 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-10-14 09:34:01.205909 - Posterior probability computed for edge n°100/3066
#> 2026-10-14 09:34:01.969097 - Posterior probability computed for edge n°200/3066
#> 2026-10-14 09:34:02.898726 - Posterior probability computed for edge n°300/3066
#> 2026-10-14 09:34:04.596053 - Posterior probability computed for edge n°400/3066
#> 2026-10-14 09:34:05.251755 - Posterior probability computed for edge n°500/3066
#> 2026-10-14 09:34:05.510314 - Posterior probability computed for edge n°600/3066
#> 2026-10-14 09:34:05.843054 - Posterior probability computed for edge n°700/3066
#> 2026-10-14 09:34:06.131072 - Posterior probability computed for edge n°800/3066
#> 2026-10-14 09:34:06.451724 - Posterior probability computed for edge n°900/3066
#> 2026-10-14 09:34:06.787059 - Posterior probability computed for edge n°1000/3066
#> 2026-10-14 09:34:07.176524 - Posterior probability computed for edge n°1100/3066
#> 2026-10-14 09:34:07.77373 - Posterior probability computed for edge n°1200/3066
#> 2026-10-14 09:34:08.364979 - Posterior probability computed for edge n°1300/3066
#> 2026-10-14 09:34:09.286379 - Posterior probability computed for edge n°1400/3066
#> 2026-10-14 09:34:09.885012 - Posterior probability computed for edge n°1500/3066
#> 2026-10-14 09:34:11.487279 - Posterior probability computed for edge n°1600/3066
#> 2026-10-14 09:34:12.392553 - Posterior probability computed for edge n°1700/3066
#> 2026-10-14 09:34:13.043503 - Posterior probability computed for edge n°1800/3066
#> 2026-10-14 09:34:13.556942 - Posterior probability computed for edge n°1900/3066
#> 2026-10-14 09:34:13.989729 - Posterior probability computed for edge n°2000/3066
#> 2026-10-14 09:34:14.560154 - Posterior probability computed for edge n°2100/3066
#> 2026-10-14 09:34:15.236942 - Posterior probability computed for edge n°2200/3066
#> 2026-10-14 09:34:16.073797 - Posterior probability computed for edge n°2300/3066
#> 2026-10-14 09:34:16.667761 - Posterior probability computed for edge n°2400/3066
#> 2026-10-14 09:34:17.208293 - Posterior probability computed for edge n°2500/3066
#> 2026-10-14 09:34:17.880172 - Posterior probability computed for edge n°2600/3066
#> 2026-10-14 09:34:18.824326 - Posterior probability computed for edge n°2700/3066
#> 2026-10-14 09:34:20.323188 - Posterior probability computed for edge n°2800/3066
#> 2026-10-14 09:34:22.048604 - Posterior probability computed for edge n°2900/3066
#> 2026-10-14 09:34:23.542835 - Posterior probability computed for edge n°3000/3066
#> 2026-10-14 09:34:24.465178 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-14 09:34:27.280362 - Posterior probability computed for edge n°100/3066
#> 2026-10-14 09:34:27.999605 - Posterior probability computed for edge n°200/3066
#> 2026-10-14 09:34:28.875256 - Posterior probability computed for edge n°300/3066
#> 2026-10-14 09:34:30.47772 - Posterior probability computed for edge n°400/3066
#> 2026-10-14 09:34:31.093813 - Posterior probability computed for edge n°500/3066
#> 2026-10-14 09:34:31.423174 - Posterior probability computed for edge n°600/3066
#> 2026-10-14 09:34:31.737157 - Posterior probability computed for edge n°700/3066
#> 2026-10-14 09:34:32.005949 - Posterior probability computed for edge n°800/3066
#> 2026-10-14 09:34:32.300922 - Posterior probability computed for edge n°900/3066
#> 2026-10-14 09:34:32.621135 - Posterior probability computed for edge n°1000/3066
#> 2026-10-14 09:34:32.923793 - Posterior probability computed for edge n°1100/3066
#> 2026-10-14 09:34:33.475824 - Posterior probability computed for edge n°1200/3066
#> 2026-10-14 09:34:34.025535 - Posterior probability computed for edge n°1300/3066
#> 2026-10-14 09:34:34.909625 - Posterior probability computed for edge n°1400/3066
#> 2026-10-14 09:34:35.464347 - Posterior probability computed for edge n°1500/3066
#> 2026-10-14 09:34:36.448643 - Posterior probability computed for edge n°1600/3066
#> 2026-10-14 09:34:37.444452 - Posterior probability computed for edge n°1700/3066
#> 2026-10-14 09:34:38.114321 - Posterior probability computed for edge n°1800/3066
#> 2026-10-14 09:34:38.655432 - Posterior probability computed for edge n°1900/3066
#> 2026-10-14 09:34:39.721774 - Posterior probability computed for edge n°2000/3066
#> 2026-10-14 09:34:40.260249 - Posterior probability computed for edge n°2100/3066
#> 2026-10-14 09:34:40.891402 - Posterior probability computed for edge n°2200/3066
#> 2026-10-14 09:34:41.673644 - Posterior probability computed for edge n°2300/3066
#> 2026-10-14 09:34:42.212469 - Posterior probability computed for edge n°2400/3066
#> 2026-10-14 09:34:42.718271 - Posterior probability computed for edge n°2500/3066
#> 2026-10-14 09:34:43.343691 - Posterior probability computed for edge n°2600/3066
#> 2026-10-14 09:34:44.244076 - Posterior probability computed for edge n°2700/3066
#> 2026-10-14 09:34:45.686539 - Posterior probability computed for edge n°2800/3066
#> 2026-10-14 09:34:47.347718 - Posterior probability computed for edge n°2900/3066
#> 2026-10-14 09:34:48.783202 - Posterior probability computed for edge n°3000/3066
#> 2026-10-14 09:34:49.676357 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-14 09:34:52.362213 - Posterior probability computed for edge n°100/3066
#> 2026-10-14 09:34:53.097248 - Posterior probability computed for edge n°200/3066
#> 2026-10-14 09:34:54.149503 - Posterior probability computed for edge n°300/3066
#> 2026-10-14 09:34:55.794851 - Posterior probability computed for edge n°400/3066
#> 2026-10-14 09:34:56.420547 - Posterior probability computed for edge n°500/3066
#> 2026-10-14 09:34:56.671275 - Posterior probability computed for edge n°600/3066
#> 2026-10-14 09:34:56.996306 - Posterior probability computed for edge n°700/3066
#> 2026-10-14 09:34:57.265392 - Posterior probability computed for edge n°800/3066
#> 2026-10-14 09:34:57.575674 - Posterior probability computed for edge n°900/3066
#> 2026-10-14 09:34:57.89921 - Posterior probability computed for edge n°1000/3066
#> 2026-10-14 09:34:58.215416 - Posterior probability computed for edge n°1100/3066
#> 2026-10-14 09:34:58.781304 - Posterior probability computed for edge n°1200/3066
#> 2026-10-14 09:34:59.342098 - Posterior probability computed for edge n°1300/3066
#> 2026-10-14 09:35:00.225997 - Posterior probability computed for edge n°1400/3066
#> 2026-10-14 09:35:00.799806 - Posterior probability computed for edge n°1500/3066
#> 2026-10-14 09:35:01.812055 - Posterior probability computed for edge n°1600/3066
#> 2026-10-14 09:35:02.767725 - Posterior probability computed for edge n°1700/3066
#> 2026-10-14 09:35:03.462709 - Posterior probability computed for edge n°1800/3066
#> 2026-10-14 09:35:04.012632 - Posterior probability computed for edge n°1900/3066
#> 2026-10-14 09:35:04.475704 - Posterior probability computed for edge n°2000/3066
#> 2026-10-14 09:35:05.089402 - Posterior probability computed for edge n°2100/3066
#> 2026-10-14 09:35:06.461365 - Posterior probability computed for edge n°2200/3066
#> 2026-10-14 09:35:07.249303 - Posterior probability computed for edge n°2300/3066
#> 2026-10-14 09:35:07.79305 - Posterior probability computed for edge n°2400/3066
#> 2026-10-14 09:35:08.306252 - Posterior probability computed for edge n°2500/3066
#> 2026-10-14 09:35:08.93797 - Posterior probability computed for edge n°2600/3066
#> 2026-10-14 09:35:09.829211 - Posterior probability computed for edge n°2700/3066
#> 2026-10-14 09:35:11.258035 - Posterior probability computed for edge n°2800/3066
#> 2026-10-14 09:35:12.921192 - Posterior probability computed for edge n°2900/3066
#> 2026-10-14 09:35:14.351736 - Posterior probability computed for edge n°3000/3066
#> 2026-10-14 09:35:15.230403 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-10-14 09:35:18.087643 - Posterior probability computed for edge n°100/3066
#> 2026-10-14 09:35:18.831851 - Posterior probability computed for edge n°200/3066
#> 2026-10-14 09:35:19.727155 - Posterior probability computed for edge n°300/3066
#> 2026-10-14 09:35:21.496208 - Posterior probability computed for edge n°400/3066
#> 2026-10-14 09:35:22.118289 - Posterior probability computed for edge n°500/3066
#> 2026-10-14 09:35:22.367923 - Posterior probability computed for edge n°600/3066
#> 2026-10-14 09:35:22.694874 - Posterior probability computed for edge n°700/3066
#> 2026-10-14 09:35:22.961654 - Posterior probability computed for edge n°800/3066
#> 2026-10-14 09:35:23.272506 - Posterior probability computed for edge n°900/3066
#> 2026-10-14 09:35:23.589352 - Posterior probability computed for edge n°1000/3066
#> 2026-10-14 09:35:23.907398 - Posterior probability computed for edge n°1100/3066
#> 2026-10-14 09:35:24.477104 - Posterior probability computed for edge n°1200/3066
#> 2026-10-14 09:35:25.034159 - Posterior probability computed for edge n°1300/3066
#> 2026-10-14 09:35:25.903461 - Posterior probability computed for edge n°1400/3066
#> 2026-10-14 09:35:26.470951 - Posterior probability computed for edge n°1500/3066
#> 2026-10-14 09:35:27.468278 - Posterior probability computed for edge n°1600/3066
#> 2026-10-14 09:35:28.424815 - Posterior probability computed for edge n°1700/3066
#> 2026-10-14 09:35:29.115826 - Posterior probability computed for edge n°1800/3066
#> 2026-10-14 09:35:29.667339 - Posterior probability computed for edge n°1900/3066
#> 2026-10-14 09:35:30.174415 - Posterior probability computed for edge n°2000/3066
#> 2026-10-14 09:35:30.78778 - Posterior probability computed for edge n°2100/3066
#> 2026-10-14 09:35:31.501616 - Posterior probability computed for edge n°2200/3066
#> 2026-10-14 09:35:33.079419 - Posterior probability computed for edge n°2300/3066
#> 2026-10-14 09:35:33.618281 - Posterior probability computed for edge n°2400/3066
#> 2026-10-14 09:35:34.132652 - Posterior probability computed for edge n°2500/3066
#> 2026-10-14 09:35:34.762205 - Posterior probability computed for edge n°2600/3066
#> 2026-10-14 09:35:35.658407 - Posterior probability computed for edge n°2700/3066
#> 2026-10-14 09:35:37.069591 - Posterior probability computed for edge n°2800/3066
#> 2026-10-14 09:35:38.714072 - Posterior probability computed for edge n°2900/3066
#> 2026-10-14 09:35:40.125077 - Posterior probability computed for edge n°3000/3066
#> 2026-10-14 09:35:41.002036 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-14 09:35:43.71613 - Posterior probability computed for edge n°100/3066
#> 2026-10-14 09:35:44.449143 - Posterior probability computed for edge n°200/3066
#> 2026-10-14 09:35:45.330651 - Posterior probability computed for edge n°300/3066
#> 2026-10-14 09:35:47.029802 - Posterior probability computed for edge n°400/3066
#> 2026-10-14 09:35:47.653261 - Posterior probability computed for edge n°500/3066
#> 2026-10-14 09:35:47.912792 - Posterior probability computed for edge n°600/3066
#> 2026-10-14 09:35:48.425578 - Posterior probability computed for edge n°700/3066
#> 2026-10-14 09:35:48.696566 - Posterior probability computed for edge n°800/3066
#> 2026-10-14 09:35:49.003556 - Posterior probability computed for edge n°900/3066
#> 2026-10-14 09:35:49.322144 - Posterior probability computed for edge n°1000/3066
#> 2026-10-14 09:35:49.63177 - Posterior probability computed for edge n°1100/3066
#> 2026-10-14 09:35:50.197696 - Posterior probability computed for edge n°1200/3066
#> 2026-10-14 09:35:50.756949 - Posterior probability computed for edge n°1300/3066
#> 2026-10-14 09:35:51.644808 - Posterior probability computed for edge n°1400/3066
#> 2026-10-14 09:35:52.212576 - Posterior probability computed for edge n°1500/3066
#> 2026-10-14 09:35:53.216005 - Posterior probability computed for edge n°1600/3066
#> 2026-10-14 09:35:54.174058 - Posterior probability computed for edge n°1700/3066
#> 2026-10-14 09:35:54.869596 - Posterior probability computed for edge n°1800/3066
#> 2026-10-14 09:35:55.412789 - Posterior probability computed for edge n°1900/3066
#> 2026-10-14 09:35:55.879725 - Posterior probability computed for edge n°2000/3066
#> 2026-10-14 09:35:56.4877 - Posterior probability computed for edge n°2100/3066
#> 2026-10-14 09:35:57.204534 - Posterior probability computed for edge n°2200/3066
#> 2026-10-14 09:35:58.091121 - Posterior probability computed for edge n°2300/3066
#> 2026-10-14 09:35:59.401249 - Posterior probability computed for edge n°2400/3066
#> 2026-10-14 09:35:59.905078 - Posterior probability computed for edge n°2500/3066
#> 2026-10-14 09:36:00.538251 - Posterior probability computed for edge n°2600/3066
#> 2026-10-14 09:36:01.433508 - Posterior probability computed for edge n°2700/3066
#> 2026-10-14 09:36:02.901352 - Posterior probability computed for edge n°2800/3066
#> 2026-10-14 09:36:04.538791 - Posterior probability computed for edge n°2900/3066
#> 2026-10-14 09:36:05.966778 - Posterior probability computed for edge n°3000/3066
#> 2026-10-14 09:36:06.848875 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-14 09:36:09.666068 - Posterior probability computed for edge n°100/3066
#> 2026-10-14 09:36:10.398012 - Posterior probability computed for edge n°200/3066
#> 2026-10-14 09:36:11.278005 - Posterior probability computed for edge n°300/3066
#> 2026-10-14 09:36:12.920509 - Posterior probability computed for edge n°400/3066
#> 2026-10-14 09:36:13.541301 - Posterior probability computed for edge n°500/3066
#> 2026-10-14 09:36:13.794703 - Posterior probability computed for edge n°600/3066
#> 2026-10-14 09:36:14.110514 - Posterior probability computed for edge n°700/3066
#> 2026-10-14 09:36:14.383418 - Posterior probability computed for edge n°800/3066
#> 2026-10-14 09:36:14.686037 - Posterior probability computed for edge n°900/3066
#> 2026-10-14 09:36:15.009219 - Posterior probability computed for edge n°1000/3066
#> 2026-10-14 09:36:15.312973 - Posterior probability computed for edge n°1100/3066
#> 2026-10-14 09:36:15.875123 - Posterior probability computed for edge n°1200/3066
#> 2026-10-14 09:36:16.431911 - Posterior probability computed for edge n°1300/3066
#> 2026-10-14 09:36:17.319686 - Posterior probability computed for edge n°1400/3066
#> 2026-10-14 09:36:17.938524 - Posterior probability computed for edge n°1500/3066
#> 2026-10-14 09:36:18.929448 - Posterior probability computed for edge n°1600/3066
#> 2026-10-14 09:36:19.877537 - Posterior probability computed for edge n°1700/3066
#> 2026-10-14 09:36:20.576449 - Posterior probability computed for edge n°1800/3066
#> 2026-10-14 09:36:21.132234 - Posterior probability computed for edge n°1900/3066
#> 2026-10-14 09:36:21.585455 - Posterior probability computed for edge n°2000/3066
#> 2026-10-14 09:36:22.190557 - Posterior probability computed for edge n°2100/3066
#> 2026-10-14 09:36:22.902203 - Posterior probability computed for edge n°2200/3066
#> 2026-10-14 09:36:23.792364 - Posterior probability computed for edge n°2300/3066
#> 2026-10-14 09:36:24.504874 - Posterior probability computed for edge n°2400/3066
#> 2026-10-14 09:36:25.820559 - Posterior probability computed for edge n°2500/3066
#> 2026-10-14 09:36:26.450634 - Posterior probability computed for edge n°2600/3066
#> 2026-10-14 09:36:27.343193 - Posterior probability computed for edge n°2700/3066
#> 2026-10-14 09:36:28.92723 - Posterior probability computed for edge n°2800/3066
#> 2026-10-14 09:36:30.570099 - Posterior probability computed for edge n°2900/3066
#> 2026-10-14 09:36:32.085786 - Posterior probability computed for edge n°3000/3066
#> 2026-10-14 09:36:32.956843 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
