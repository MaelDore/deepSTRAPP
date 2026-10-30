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
#> 2026-10-30 00:59:07.09644 - Fit 1 evolutionary model(s): ER.
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
#> 2026-10-30 00:59:10.009869 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-10-30 00:59:10.011661 - Run simulations for stochastic mapping.
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
#> 2026-10-30 01:01:27.660177 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-10-30 01:01:30.963099 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-10-30 01:01:34.825133 - Posterior probability computed for edge n°100/3066
#> 2026-10-30 01:01:35.777763 - Posterior probability computed for edge n°200/3066
#> 2026-10-30 01:01:36.942395 - Posterior probability computed for edge n°300/3066
#> 2026-10-30 01:01:39.051235 - Posterior probability computed for edge n°400/3066
#> 2026-10-30 01:01:39.857297 - Posterior probability computed for edge n°500/3066
#> 2026-10-30 01:01:40.189281 - Posterior probability computed for edge n°600/3066
#> 2026-10-30 01:01:40.603086 - Posterior probability computed for edge n°700/3066
#> 2026-10-30 01:01:40.956158 - Posterior probability computed for edge n°800/3066
#> 2026-10-30 01:01:41.357194 - Posterior probability computed for edge n°900/3066
#> 2026-10-30 01:01:41.771922 - Posterior probability computed for edge n°1000/3066
#> 2026-10-30 01:01:42.178664 - Posterior probability computed for edge n°1100/3066
#> 2026-10-30 01:01:42.912864 - Posterior probability computed for edge n°1200/3066
#> 2026-10-30 01:01:43.625137 - Posterior probability computed for edge n°1300/3066
#> 2026-10-30 01:01:44.753599 - Posterior probability computed for edge n°1400/3066
#> 2026-10-30 01:01:45.491881 - Posterior probability computed for edge n°1500/3066
#> 2026-10-30 01:01:46.779598 - Posterior probability computed for edge n°1600/3066
#> 2026-10-30 01:01:48.503873 - Posterior probability computed for edge n°1700/3066
#> 2026-10-30 01:01:49.338429 - Posterior probability computed for edge n°1800/3066
#> 2026-10-30 01:01:50.00014 - Posterior probability computed for edge n°1900/3066
#> 2026-10-30 01:01:50.565262 - Posterior probability computed for edge n°2000/3066
#> 2026-10-30 01:01:51.310625 - Posterior probability computed for edge n°2100/3066
#> 2026-10-30 01:01:52.17602 - Posterior probability computed for edge n°2200/3066
#> 2026-10-30 01:01:53.258046 - Posterior probability computed for edge n°2300/3066
#> 2026-10-30 01:01:53.99849 - Posterior probability computed for edge n°2400/3066
#> 2026-10-30 01:01:54.710799 - Posterior probability computed for edge n°2500/3066
#> 2026-10-30 01:01:55.572669 - Posterior probability computed for edge n°2600/3066
#> 2026-10-30 01:01:56.794622 - Posterior probability computed for edge n°2700/3066
#> 2026-10-30 01:01:58.717611 - Posterior probability computed for edge n°2800/3066
#> 2026-10-30 01:02:00.943834 - Posterior probability computed for edge n°2900/3066
#> 2026-10-30 01:02:02.862512 - Posterior probability computed for edge n°3000/3066
#> 2026-10-30 01:02:04.04508 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-30 01:02:07.437098 - Posterior probability computed for edge n°100/3066
#> 2026-10-30 01:02:08.366303 - Posterior probability computed for edge n°200/3066
#> 2026-10-30 01:02:09.48506 - Posterior probability computed for edge n°300/3066
#> 2026-10-30 01:02:11.525775 - Posterior probability computed for edge n°400/3066
#> 2026-10-30 01:02:12.304008 - Posterior probability computed for edge n°500/3066
#> 2026-10-30 01:02:12.622633 - Posterior probability computed for edge n°600/3066
#> 2026-10-30 01:02:13.032969 - Posterior probability computed for edge n°700/3066
#> 2026-10-30 01:02:13.371911 - Posterior probability computed for edge n°800/3066
#> 2026-10-30 01:02:13.76139 - Posterior probability computed for edge n°900/3066
#> 2026-10-30 01:02:14.16507 - Posterior probability computed for edge n°1000/3066
#> 2026-10-30 01:02:14.562225 - Posterior probability computed for edge n°1100/3066
#> 2026-10-30 01:02:15.26748 - Posterior probability computed for edge n°1200/3066
#> 2026-10-30 01:02:15.960176 - Posterior probability computed for edge n°1300/3066
#> 2026-10-30 01:02:17.054079 - Posterior probability computed for edge n°1400/3066
#> 2026-10-30 01:02:17.769352 - Posterior probability computed for edge n°1500/3066
#> 2026-10-30 01:02:19.004962 - Posterior probability computed for edge n°1600/3066
#> 2026-10-30 01:02:20.183485 - Posterior probability computed for edge n°1700/3066
#> 2026-10-30 01:02:21.040306 - Posterior probability computed for edge n°1800/3066
#> 2026-10-30 01:02:21.721962 - Posterior probability computed for edge n°1900/3066
#> 2026-10-30 01:02:22.868667 - Posterior probability computed for edge n°2000/3066
#> 2026-10-30 01:02:23.577313 - Posterior probability computed for edge n°2100/3066
#> 2026-10-30 01:02:24.409154 - Posterior probability computed for edge n°2200/3066
#> 2026-10-30 01:02:25.445851 - Posterior probability computed for edge n°2300/3066
#> 2026-10-30 01:02:26.157017 - Posterior probability computed for edge n°2400/3066
#> 2026-10-30 01:02:26.826291 - Posterior probability computed for edge n°2500/3066
#> 2026-10-30 01:02:27.654743 - Posterior probability computed for edge n°2600/3066
#> 2026-10-30 01:02:28.820386 - Posterior probability computed for edge n°2700/3066
#> 2026-10-30 01:02:30.701336 - Posterior probability computed for edge n°2800/3066
#> 2026-10-30 01:02:32.852329 - Posterior probability computed for edge n°2900/3066
#> 2026-10-30 01:02:34.718882 - Posterior probability computed for edge n°3000/3066
#> 2026-10-30 01:02:35.876304 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-30 01:02:39.091194 - Posterior probability computed for edge n°100/3066
#> 2026-10-30 01:02:40.018603 - Posterior probability computed for edge n°200/3066
#> 2026-10-30 01:02:41.260596 - Posterior probability computed for edge n°300/3066
#> 2026-10-30 01:02:43.301849 - Posterior probability computed for edge n°400/3066
#> 2026-10-30 01:02:44.081375 - Posterior probability computed for edge n°500/3066
#> 2026-10-30 01:02:44.400788 - Posterior probability computed for edge n°600/3066
#> 2026-10-30 01:02:44.811776 - Posterior probability computed for edge n°700/3066
#> 2026-10-30 01:02:45.150471 - Posterior probability computed for edge n°800/3066
#> 2026-10-30 01:02:45.545448 - Posterior probability computed for edge n°900/3066
#> 2026-10-30 01:02:45.947295 - Posterior probability computed for edge n°1000/3066
#> 2026-10-30 01:02:46.347671 - Posterior probability computed for edge n°1100/3066
#> 2026-10-30 01:02:47.061274 - Posterior probability computed for edge n°1200/3066
#> 2026-10-30 01:02:47.76277 - Posterior probability computed for edge n°1300/3066
#> 2026-10-30 01:02:48.851162 - Posterior probability computed for edge n°1400/3066
#> 2026-10-30 01:02:49.563916 - Posterior probability computed for edge n°1500/3066
#> 2026-10-30 01:02:50.808735 - Posterior probability computed for edge n°1600/3066
#> 2026-10-30 01:02:51.996631 - Posterior probability computed for edge n°1700/3066
#> 2026-10-30 01:02:52.861786 - Posterior probability computed for edge n°1800/3066
#> 2026-10-30 01:02:53.548386 - Posterior probability computed for edge n°1900/3066
#> 2026-10-30 01:02:54.122455 - Posterior probability computed for edge n°2000/3066
#> 2026-10-30 01:02:54.884246 - Posterior probability computed for edge n°2100/3066
#> 2026-10-30 01:02:56.350233 - Posterior probability computed for edge n°2200/3066
#> 2026-10-30 01:02:57.382193 - Posterior probability computed for edge n°2300/3066
#> 2026-10-30 01:02:58.09416 - Posterior probability computed for edge n°2400/3066
#> 2026-10-30 01:02:58.759367 - Posterior probability computed for edge n°2500/3066
#> 2026-10-30 01:02:59.588267 - Posterior probability computed for edge n°2600/3066
#> 2026-10-30 01:03:00.75437 - Posterior probability computed for edge n°2700/3066
#> 2026-10-30 01:03:02.616084 - Posterior probability computed for edge n°2800/3066
#> 2026-10-30 01:03:04.758669 - Posterior probability computed for edge n°2900/3066
#> 2026-10-30 01:03:06.623703 - Posterior probability computed for edge n°3000/3066
#> 2026-10-30 01:03:07.764261 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-10-30 01:03:11.131484 - Posterior probability computed for edge n°100/3066
#> 2026-10-30 01:03:12.061855 - Posterior probability computed for edge n°200/3066
#> 2026-10-30 01:03:13.171021 - Posterior probability computed for edge n°300/3066
#> 2026-10-30 01:03:15.304906 - Posterior probability computed for edge n°400/3066
#> 2026-10-30 01:03:16.08184 - Posterior probability computed for edge n°500/3066
#> 2026-10-30 01:03:16.404857 - Posterior probability computed for edge n°600/3066
#> 2026-10-30 01:03:16.811847 - Posterior probability computed for edge n°700/3066
#> 2026-10-30 01:03:17.149988 - Posterior probability computed for edge n°800/3066
#> 2026-10-30 01:03:17.543569 - Posterior probability computed for edge n°900/3066
#> 2026-10-30 01:03:17.949584 - Posterior probability computed for edge n°1000/3066
#> 2026-10-30 01:03:18.344223 - Posterior probability computed for edge n°1100/3066
#> 2026-10-30 01:03:19.054763 - Posterior probability computed for edge n°1200/3066
#> 2026-10-30 01:03:19.752166 - Posterior probability computed for edge n°1300/3066
#> 2026-10-30 01:03:20.840706 - Posterior probability computed for edge n°1400/3066
#> 2026-10-30 01:03:21.550508 - Posterior probability computed for edge n°1500/3066
#> 2026-10-30 01:03:22.79743 - Posterior probability computed for edge n°1600/3066
#> 2026-10-30 01:03:23.989151 - Posterior probability computed for edge n°1700/3066
#> 2026-10-30 01:03:24.848089 - Posterior probability computed for edge n°1800/3066
#> 2026-10-30 01:03:25.528157 - Posterior probability computed for edge n°1900/3066
#> 2026-10-30 01:03:26.108411 - Posterior probability computed for edge n°2000/3066
#> 2026-10-30 01:03:26.873154 - Posterior probability computed for edge n°2100/3066
#> 2026-10-30 01:03:27.763571 - Posterior probability computed for edge n°2200/3066
#> 2026-10-30 01:03:29.443031 - Posterior probability computed for edge n°2300/3066
#> 2026-10-30 01:03:30.157765 - Posterior probability computed for edge n°2400/3066
#> 2026-10-30 01:03:30.826401 - Posterior probability computed for edge n°2500/3066
#> 2026-10-30 01:03:31.660485 - Posterior probability computed for edge n°2600/3066
#> 2026-10-30 01:03:32.828217 - Posterior probability computed for edge n°2700/3066
#> 2026-10-30 01:03:34.677108 - Posterior probability computed for edge n°2800/3066
#> 2026-10-30 01:03:36.826545 - Posterior probability computed for edge n°2900/3066
#> 2026-10-30 01:03:38.680327 - Posterior probability computed for edge n°3000/3066
#> 2026-10-30 01:03:39.818385 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-30 01:03:43.043103 - Posterior probability computed for edge n°100/3066
#> 2026-10-30 01:03:43.975386 - Posterior probability computed for edge n°200/3066
#> 2026-10-30 01:03:45.094002 - Posterior probability computed for edge n°300/3066
#> 2026-10-30 01:03:47.244216 - Posterior probability computed for edge n°400/3066
#> 2026-10-30 01:03:48.023839 - Posterior probability computed for edge n°500/3066
#> 2026-10-30 01:03:48.349714 - Posterior probability computed for edge n°600/3066
#> 2026-10-30 01:03:48.759274 - Posterior probability computed for edge n°700/3066
#> 2026-10-30 01:03:49.104447 - Posterior probability computed for edge n°800/3066
#> 2026-10-30 01:03:49.493055 - Posterior probability computed for edge n°900/3066
#> 2026-10-30 01:03:49.902605 - Posterior probability computed for edge n°1000/3066
#> 2026-10-30 01:03:50.295924 - Posterior probability computed for edge n°1100/3066
#> 2026-10-30 01:03:51.011615 - Posterior probability computed for edge n°1200/3066
#> 2026-10-30 01:03:51.714774 - Posterior probability computed for edge n°1300/3066
#> 2026-10-30 01:03:52.814303 - Posterior probability computed for edge n°1400/3066
#> 2026-10-30 01:03:53.529759 - Posterior probability computed for edge n°1500/3066
#> 2026-10-30 01:03:54.782444 - Posterior probability computed for edge n°1600/3066
#> 2026-10-30 01:03:55.975009 - Posterior probability computed for edge n°1700/3066
#> 2026-10-30 01:03:56.84812 - Posterior probability computed for edge n°1800/3066
#> 2026-10-30 01:03:57.539402 - Posterior probability computed for edge n°1900/3066
#> 2026-10-30 01:03:58.116934 - Posterior probability computed for edge n°2000/3066
#> 2026-10-30 01:03:58.879132 - Posterior probability computed for edge n°2100/3066
#> 2026-10-30 01:03:59.778966 - Posterior probability computed for edge n°2200/3066
#> 2026-10-30 01:04:00.893431 - Posterior probability computed for edge n°2300/3066
#> 2026-10-30 01:04:02.247598 - Posterior probability computed for edge n°2400/3066
#> 2026-10-30 01:04:02.921202 - Posterior probability computed for edge n°2500/3066
#> 2026-10-30 01:04:03.756985 - Posterior probability computed for edge n°2600/3066
#> 2026-10-30 01:04:04.936764 - Posterior probability computed for edge n°2700/3066
#> 2026-10-30 01:04:06.804032 - Posterior probability computed for edge n°2800/3066
#> 2026-10-30 01:04:08.957166 - Posterior probability computed for edge n°2900/3066
#> 2026-10-30 01:04:10.809526 - Posterior probability computed for edge n°3000/3066
#> 2026-10-30 01:04:11.958112 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-30 01:04:15.27533 - Posterior probability computed for edge n°100/3066
#> 2026-10-30 01:04:16.202525 - Posterior probability computed for edge n°200/3066
#> 2026-10-30 01:04:17.307394 - Posterior probability computed for edge n°300/3066
#> 2026-10-30 01:04:19.373539 - Posterior probability computed for edge n°400/3066
#> 2026-10-30 01:04:20.16134 - Posterior probability computed for edge n°500/3066
#> 2026-10-30 01:04:20.482484 - Posterior probability computed for edge n°600/3066
#> 2026-10-30 01:04:20.890759 - Posterior probability computed for edge n°700/3066
#> 2026-10-30 01:04:21.234683 - Posterior probability computed for edge n°800/3066
#> 2026-10-30 01:04:21.622027 - Posterior probability computed for edge n°900/3066
#> 2026-10-30 01:04:22.037463 - Posterior probability computed for edge n°1000/3066
#> 2026-10-30 01:04:22.435152 - Posterior probability computed for edge n°1100/3066
#> 2026-10-30 01:04:23.152138 - Posterior probability computed for edge n°1200/3066
#> 2026-10-30 01:04:23.853238 - Posterior probability computed for edge n°1300/3066
#> 2026-10-30 01:04:24.946291 - Posterior probability computed for edge n°1400/3066
#> 2026-10-30 01:04:25.659451 - Posterior probability computed for edge n°1500/3066
#> 2026-10-30 01:04:26.906688 - Posterior probability computed for edge n°1600/3066
#> 2026-10-30 01:04:28.108743 - Posterior probability computed for edge n°1700/3066
#> 2026-10-30 01:04:28.976084 - Posterior probability computed for edge n°1800/3066
#> 2026-10-30 01:04:29.669441 - Posterior probability computed for edge n°1900/3066
#> 2026-10-30 01:04:30.25296 - Posterior probability computed for edge n°2000/3066
#> 2026-10-30 01:04:31.024027 - Posterior probability computed for edge n°2100/3066
#> 2026-10-30 01:04:31.93247 - Posterior probability computed for edge n°2200/3066
#> 2026-10-30 01:04:33.045087 - Posterior probability computed for edge n°2300/3066
#> 2026-10-30 01:04:34.406492 - Posterior probability computed for edge n°2400/3066
#> 2026-10-30 01:04:35.077412 - Posterior probability computed for edge n°2500/3066
#> 2026-10-30 01:04:35.912967 - Posterior probability computed for edge n°2600/3066
#> 2026-10-30 01:04:37.08554 - Posterior probability computed for edge n°2700/3066
#> 2026-10-30 01:04:38.942529 - Posterior probability computed for edge n°2800/3066
#> 2026-10-30 01:04:41.083087 - Posterior probability computed for edge n°2900/3066
#> 2026-10-30 01:04:42.939564 - Posterior probability computed for edge n°3000/3066
#> 2026-10-30 01:04:44.072361 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
