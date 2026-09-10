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
#> 2026-09-10 00:35:32.404216 - Fit 1 evolutionary model(s): ER.
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
#> 2026-09-10 00:35:36.028277 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-09-10 00:35:36.030121 - Run simulations for stochastic mapping.
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
#> 2026-09-10 00:37:57.613342 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-09-10 00:38:01.089564 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-09-10 00:38:04.889829 - Posterior probability computed for edge n°100/3066
#> 2026-09-10 00:38:05.971196 - Posterior probability computed for edge n°200/3066
#> 2026-09-10 00:38:07.25668 - Posterior probability computed for edge n°300/3066
#> 2026-09-10 00:38:09.590083 - Posterior probability computed for edge n°400/3066
#> 2026-09-10 00:38:10.477428 - Posterior probability computed for edge n°500/3066
#> 2026-09-10 00:38:11.323729 - Posterior probability computed for edge n°600/3066
#> 2026-09-10 00:38:11.759014 - Posterior probability computed for edge n°700/3066
#> 2026-09-10 00:38:12.095956 - Posterior probability computed for edge n°800/3066
#> 2026-09-10 00:38:12.47503 - Posterior probability computed for edge n°900/3066
#> 2026-09-10 00:38:12.872049 - Posterior probability computed for edge n°1000/3066
#> 2026-09-10 00:38:13.256326 - Posterior probability computed for edge n°1100/3066
#> 2026-09-10 00:38:13.953329 - Posterior probability computed for edge n°1200/3066
#> 2026-09-10 00:38:14.637217 - Posterior probability computed for edge n°1300/3066
#> 2026-09-10 00:38:15.704181 - Posterior probability computed for edge n°1400/3066
#> 2026-09-10 00:38:16.39425 - Posterior probability computed for edge n°1500/3066
#> 2026-09-10 00:38:17.603899 - Posterior probability computed for edge n°1600/3066
#> 2026-09-10 00:38:18.753122 - Posterior probability computed for edge n°1700/3066
#> 2026-09-10 00:38:19.606362 - Posterior probability computed for edge n°1800/3066
#> 2026-09-10 00:38:20.270977 - Posterior probability computed for edge n°1900/3066
#> 2026-09-10 00:38:20.827039 - Posterior probability computed for edge n°2000/3066
#> 2026-09-10 00:38:21.570869 - Posterior probability computed for edge n°2100/3066
#> 2026-09-10 00:38:22.442428 - Posterior probability computed for edge n°2200/3066
#> 2026-09-10 00:38:23.516606 - Posterior probability computed for edge n°2300/3066
#> 2026-09-10 00:38:24.264513 - Posterior probability computed for edge n°2400/3066
#> 2026-09-10 00:38:24.958955 - Posterior probability computed for edge n°2500/3066
#> 2026-09-10 00:38:25.814873 - Posterior probability computed for edge n°2600/3066
#> 2026-09-10 00:38:27.052874 - Posterior probability computed for edge n°2700/3066
#> 2026-09-10 00:38:28.967897 - Posterior probability computed for edge n°2800/3066
#> 2026-09-10 00:38:31.181312 - Posterior probability computed for edge n°2900/3066
#> 2026-09-10 00:38:33.088799 - Posterior probability computed for edge n°3000/3066
#> 2026-09-10 00:38:34.273292 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-10 00:38:37.707699 - Posterior probability computed for edge n°100/3066
#> 2026-09-10 00:38:38.653539 - Posterior probability computed for edge n°200/3066
#> 2026-09-10 00:38:39.780736 - Posterior probability computed for edge n°300/3066
#> 2026-09-10 00:38:41.850015 - Posterior probability computed for edge n°400/3066
#> 2026-09-10 00:38:42.632053 - Posterior probability computed for edge n°500/3066
#> 2026-09-10 00:38:42.956854 - Posterior probability computed for edge n°600/3066
#> 2026-09-10 00:38:43.379858 - Posterior probability computed for edge n°700/3066
#> 2026-09-10 00:38:43.735764 - Posterior probability computed for edge n°800/3066
#> 2026-09-10 00:38:44.127071 - Posterior probability computed for edge n°900/3066
#> 2026-09-10 00:38:44.532097 - Posterior probability computed for edge n°1000/3066
#> 2026-09-10 00:38:44.940948 - Posterior probability computed for edge n°1100/3066
#> 2026-09-10 00:38:45.657911 - Posterior probability computed for edge n°1200/3066
#> 2026-09-10 00:38:46.356775 - Posterior probability computed for edge n°1300/3066
#> 2026-09-10 00:38:48.043514 - Posterior probability computed for edge n°1400/3066
#> 2026-09-10 00:38:48.726671 - Posterior probability computed for edge n°1500/3066
#> 2026-09-10 00:38:49.913431 - Posterior probability computed for edge n°1600/3066
#> 2026-09-10 00:38:51.048846 - Posterior probability computed for edge n°1700/3066
#> 2026-09-10 00:38:51.886003 - Posterior probability computed for edge n°1800/3066
#> 2026-09-10 00:38:52.55194 - Posterior probability computed for edge n°1900/3066
#> 2026-09-10 00:38:53.092643 - Posterior probability computed for edge n°2000/3066
#> 2026-09-10 00:38:53.822004 - Posterior probability computed for edge n°2100/3066
#> 2026-09-10 00:38:54.685317 - Posterior probability computed for edge n°2200/3066
#> 2026-09-10 00:38:55.739187 - Posterior probability computed for edge n°2300/3066
#> 2026-09-10 00:38:56.467889 - Posterior probability computed for edge n°2400/3066
#> 2026-09-10 00:38:57.158397 - Posterior probability computed for edge n°2500/3066
#> 2026-09-10 00:38:57.998566 - Posterior probability computed for edge n°2600/3066
#> 2026-09-10 00:38:59.240282 - Posterior probability computed for edge n°2700/3066
#> 2026-09-10 00:39:01.125381 - Posterior probability computed for edge n°2800/3066
#> 2026-09-10 00:39:03.32447 - Posterior probability computed for edge n°2900/3066
#> 2026-09-10 00:39:05.201608 - Posterior probability computed for edge n°3000/3066
#> 2026-09-10 00:39:06.383957 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-10 00:39:09.75582 - Posterior probability computed for edge n°100/3066
#> 2026-09-10 00:39:10.693317 - Posterior probability computed for edge n°200/3066
#> 2026-09-10 00:39:11.816098 - Posterior probability computed for edge n°300/3066
#> 2026-09-10 00:39:13.865507 - Posterior probability computed for edge n°400/3066
#> 2026-09-10 00:39:14.659656 - Posterior probability computed for edge n°500/3066
#> 2026-09-10 00:39:14.981942 - Posterior probability computed for edge n°600/3066
#> 2026-09-10 00:39:15.386738 - Posterior probability computed for edge n°700/3066
#> 2026-09-10 00:39:15.733644 - Posterior probability computed for edge n°800/3066
#> 2026-09-10 00:39:16.133304 - Posterior probability computed for edge n°900/3066
#> 2026-09-10 00:39:16.538617 - Posterior probability computed for edge n°1000/3066
#> 2026-09-10 00:39:16.93547 - Posterior probability computed for edge n°1100/3066
#> 2026-09-10 00:39:17.661861 - Posterior probability computed for edge n°1200/3066
#> 2026-09-10 00:39:18.365605 - Posterior probability computed for edge n°1300/3066
#> 2026-09-10 00:39:19.471517 - Posterior probability computed for edge n°1400/3066
#> 2026-09-10 00:39:20.188733 - Posterior probability computed for edge n°1500/3066
#> 2026-09-10 00:39:21.431964 - Posterior probability computed for edge n°1600/3066
#> 2026-09-10 00:39:22.638991 - Posterior probability computed for edge n°1700/3066
#> 2026-09-10 00:39:23.492668 - Posterior probability computed for edge n°1800/3066
#> 2026-09-10 00:39:24.196171 - Posterior probability computed for edge n°1900/3066
#> 2026-09-10 00:39:24.774378 - Posterior probability computed for edge n°2000/3066
#> 2026-09-10 00:39:25.534874 - Posterior probability computed for edge n°2100/3066
#> 2026-09-10 00:39:26.430098 - Posterior probability computed for edge n°2200/3066
#> 2026-09-10 00:39:27.548947 - Posterior probability computed for edge n°2300/3066
#> 2026-09-10 00:39:28.907837 - Posterior probability computed for edge n°2400/3066
#> 2026-09-10 00:39:29.577322 - Posterior probability computed for edge n°2500/3066
#> 2026-09-10 00:39:30.427543 - Posterior probability computed for edge n°2600/3066
#> 2026-09-10 00:39:31.613026 - Posterior probability computed for edge n°2700/3066
#> 2026-09-10 00:39:33.503672 - Posterior probability computed for edge n°2800/3066
#> 2026-09-10 00:39:35.678696 - Posterior probability computed for edge n°2900/3066
#> 2026-09-10 00:39:37.542893 - Posterior probability computed for edge n°3000/3066
#> 2026-09-10 00:39:38.702503 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-09-10 00:39:42.009539 - Posterior probability computed for edge n°100/3066
#> 2026-09-10 00:39:43.065369 - Posterior probability computed for edge n°200/3066
#> 2026-09-10 00:39:44.185271 - Posterior probability computed for edge n°300/3066
#> 2026-09-10 00:39:46.222823 - Posterior probability computed for edge n°400/3066
#> 2026-09-10 00:39:47.00883 - Posterior probability computed for edge n°500/3066
#> 2026-09-10 00:39:47.328977 - Posterior probability computed for edge n°600/3066
#> 2026-09-10 00:39:47.732417 - Posterior probability computed for edge n°700/3066
#> 2026-09-10 00:39:48.071262 - Posterior probability computed for edge n°800/3066
#> 2026-09-10 00:39:48.472074 - Posterior probability computed for edge n°900/3066
#> 2026-09-10 00:39:48.872874 - Posterior probability computed for edge n°1000/3066
#> 2026-09-10 00:39:49.264376 - Posterior probability computed for edge n°1100/3066
#> 2026-09-10 00:39:49.984101 - Posterior probability computed for edge n°1200/3066
#> 2026-09-10 00:39:50.67762 - Posterior probability computed for edge n°1300/3066
#> 2026-09-10 00:39:51.772237 - Posterior probability computed for edge n°1400/3066
#> 2026-09-10 00:39:52.481392 - Posterior probability computed for edge n°1500/3066
#> 2026-09-10 00:39:53.734979 - Posterior probability computed for edge n°1600/3066
#> 2026-09-10 00:39:54.920786 - Posterior probability computed for edge n°1700/3066
#> 2026-09-10 00:39:55.779661 - Posterior probability computed for edge n°1800/3066
#> 2026-09-10 00:39:56.458525 - Posterior probability computed for edge n°1900/3066
#> 2026-09-10 00:39:57.033724 - Posterior probability computed for edge n°2000/3066
#> 2026-09-10 00:39:57.803146 - Posterior probability computed for edge n°2100/3066
#> 2026-09-10 00:39:58.693433 - Posterior probability computed for edge n°2200/3066
#> 2026-09-10 00:39:59.785347 - Posterior probability computed for edge n°2300/3066
#> 2026-09-10 00:40:00.559826 - Posterior probability computed for edge n°2400/3066
#> 2026-09-10 00:40:01.282479 - Posterior probability computed for edge n°2500/3066
#> 2026-09-10 00:40:02.167802 - Posterior probability computed for edge n°2600/3066
#> 2026-09-10 00:40:03.413043 - Posterior probability computed for edge n°2700/3066
#> 2026-09-10 00:40:05.38793 - Posterior probability computed for edge n°2800/3066
#> 2026-09-10 00:40:08.211403 - Posterior probability computed for edge n°2900/3066
#> 2026-09-10 00:40:10.067758 - Posterior probability computed for edge n°3000/3066
#> 2026-09-10 00:40:11.226863 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-10 00:40:14.524412 - Posterior probability computed for edge n°100/3066
#> 2026-09-10 00:40:15.444052 - Posterior probability computed for edge n°200/3066
#> 2026-09-10 00:40:16.55874 - Posterior probability computed for edge n°300/3066
#> 2026-09-10 00:40:18.5914 - Posterior probability computed for edge n°400/3066
#> 2026-09-10 00:40:19.376056 - Posterior probability computed for edge n°500/3066
#> 2026-09-10 00:40:19.696376 - Posterior probability computed for edge n°600/3066
#> 2026-09-10 00:40:20.09701 - Posterior probability computed for edge n°700/3066
#> 2026-09-10 00:40:20.437541 - Posterior probability computed for edge n°800/3066
#> 2026-09-10 00:40:20.840109 - Posterior probability computed for edge n°900/3066
#> 2026-09-10 00:40:21.242058 - Posterior probability computed for edge n°1000/3066
#> 2026-09-10 00:40:21.63028 - Posterior probability computed for edge n°1100/3066
#> 2026-09-10 00:40:22.347276 - Posterior probability computed for edge n°1200/3066
#> 2026-09-10 00:40:23.040809 - Posterior probability computed for edge n°1300/3066
#> 2026-09-10 00:40:24.133823 - Posterior probability computed for edge n°1400/3066
#> 2026-09-10 00:40:24.838042 - Posterior probability computed for edge n°1500/3066
#> 2026-09-10 00:40:26.08982 - Posterior probability computed for edge n°1600/3066
#> 2026-09-10 00:40:27.275757 - Posterior probability computed for edge n°1700/3066
#> 2026-09-10 00:40:28.136888 - Posterior probability computed for edge n°1800/3066
#> 2026-09-10 00:40:28.818269 - Posterior probability computed for edge n°1900/3066
#> 2026-09-10 00:40:29.395713 - Posterior probability computed for edge n°2000/3066
#> 2026-09-10 00:40:30.162401 - Posterior probability computed for edge n°2100/3066
#> 2026-09-10 00:40:31.053741 - Posterior probability computed for edge n°2200/3066
#> 2026-09-10 00:40:32.161036 - Posterior probability computed for edge n°2300/3066
#> 2026-09-10 00:40:32.920436 - Posterior probability computed for edge n°2400/3066
#> 2026-09-10 00:40:33.64334 - Posterior probability computed for edge n°2500/3066
#> 2026-09-10 00:40:34.530105 - Posterior probability computed for edge n°2600/3066
#> 2026-09-10 00:40:35.775717 - Posterior probability computed for edge n°2700/3066
#> 2026-09-10 00:40:37.750166 - Posterior probability computed for edge n°2800/3066
#> 2026-09-10 00:40:40.031517 - Posterior probability computed for edge n°2900/3066
#> 2026-09-10 00:40:42.005479 - Posterior probability computed for edge n°3000/3066
#> 2026-09-10 00:40:43.239103 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-10 00:40:47.113584 - Posterior probability computed for edge n°100/3066
#> 2026-09-10 00:40:48.013371 - Posterior probability computed for edge n°200/3066
#> 2026-09-10 00:40:49.07414 - Posterior probability computed for edge n°300/3066
#> 2026-09-10 00:40:50.992778 - Posterior probability computed for edge n°400/3066
#> 2026-09-10 00:40:51.744332 - Posterior probability computed for edge n°500/3066
#> 2026-09-10 00:40:52.032688 - Posterior probability computed for edge n°600/3066
#> 2026-09-10 00:40:52.460542 - Posterior probability computed for edge n°700/3066
#> 2026-09-10 00:40:52.769917 - Posterior probability computed for edge n°800/3066
#> 2026-09-10 00:40:53.136965 - Posterior probability computed for edge n°900/3066
#> 2026-09-10 00:40:53.532655 - Posterior probability computed for edge n°1000/3066
#> 2026-09-10 00:40:53.903536 - Posterior probability computed for edge n°1100/3066
#> 2026-09-10 00:40:54.572379 - Posterior probability computed for edge n°1200/3066
#> 2026-09-10 00:40:55.218602 - Posterior probability computed for edge n°1300/3066
#> 2026-09-10 00:40:56.265713 - Posterior probability computed for edge n°1400/3066
#> 2026-09-10 00:40:56.93732 - Posterior probability computed for edge n°1500/3066
#> 2026-09-10 00:40:58.105627 - Posterior probability computed for edge n°1600/3066
#> 2026-09-10 00:40:59.220595 - Posterior probability computed for edge n°1700/3066
#> 2026-09-10 00:41:00.054478 - Posterior probability computed for edge n°1800/3066
#> 2026-09-10 00:41:00.70451 - Posterior probability computed for edge n°1900/3066
#> 2026-09-10 00:41:01.239969 - Posterior probability computed for edge n°2000/3066
#> 2026-09-10 00:41:01.960557 - Posterior probability computed for edge n°2100/3066
#> 2026-09-10 00:41:02.809628 - Posterior probability computed for edge n°2200/3066
#> 2026-09-10 00:41:03.852395 - Posterior probability computed for edge n°2300/3066
#> 2026-09-10 00:41:04.57941 - Posterior probability computed for edge n°2400/3066
#> 2026-09-10 00:41:05.265668 - Posterior probability computed for edge n°2500/3066
#> 2026-09-10 00:41:06.114105 - Posterior probability computed for edge n°2600/3066
#> 2026-09-10 00:41:07.290028 - Posterior probability computed for edge n°2700/3066
#> 2026-09-10 00:41:09.17427 - Posterior probability computed for edge n°2800/3066
#> 2026-09-10 00:41:11.332157 - Posterior probability computed for edge n°2900/3066
#> 2026-09-10 00:41:13.209884 - Posterior probability computed for edge n°3000/3066
#> 2026-09-10 00:41:14.376191 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
