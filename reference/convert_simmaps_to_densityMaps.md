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
#> 2026-10-09 01:22:24.300663 - Fit 1 evolutionary model(s): ER.
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
#> 2026-10-09 01:22:27.663851 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-10-09 01:22:27.665336 - Run simulations for stochastic mapping.
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
#> 2026-10-09 01:24:10.877733 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-10-09 01:24:13.707058 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-10-09 01:24:17.099399 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 01:24:18.000796 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 01:24:19.073107 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 01:24:21.024339 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 01:24:22.354307 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 01:24:22.602641 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 01:24:22.925573 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 01:24:23.19122 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 01:24:23.491996 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 01:24:23.804392 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 01:24:24.105999 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 01:24:24.656672 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 01:24:25.197562 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 01:24:26.041481 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 01:24:26.589476 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 01:24:27.557199 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 01:24:28.46973 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 01:24:29.1222 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 01:24:29.652097 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 01:24:30.102845 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 01:24:30.690014 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 01:24:31.375961 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 01:24:32.228522 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 01:24:32.81636 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 01:24:33.369361 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 01:24:34.056722 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 01:24:35.041311 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 01:24:36.587824 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 01:24:38.357417 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 01:24:39.897165 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 01:24:40.83957 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-09 01:24:43.759486 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 01:24:44.515608 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 01:24:45.431295 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 01:24:47.247697 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 01:24:47.898537 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 01:24:48.160305 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 01:24:48.488682 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 01:24:48.763267 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 01:24:49.071006 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 01:24:49.407095 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 01:24:49.725343 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 01:24:50.296814 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 01:24:50.880438 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 01:24:51.762117 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 01:24:52.357035 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 01:24:53.368705 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 01:24:54.340697 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 01:24:55.048259 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 01:24:55.604129 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 01:24:56.078734 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 01:24:57.34479 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 01:24:58.017887 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 01:24:58.843202 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 01:24:59.416571 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 01:24:59.956709 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 01:25:00.619515 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 01:25:01.571146 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 01:25:03.075763 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 01:25:04.846189 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 01:25:06.384971 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 01:25:07.306089 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-09 01:25:10.142407 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 01:25:10.883811 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 01:25:11.780313 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 01:25:13.460911 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 01:25:14.132955 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 01:25:14.384507 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 01:25:14.698178 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 01:25:14.979117 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 01:25:15.281933 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 01:25:15.596567 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 01:25:15.90683 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 01:25:16.482446 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 01:25:17.032684 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 01:25:17.907941 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 01:25:18.470551 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 01:25:19.485314 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 01:25:20.438002 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 01:25:21.129348 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 01:25:21.669926 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 01:25:22.128194 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 01:25:22.741341 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 01:25:23.453932 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 01:25:24.327544 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 01:25:24.942002 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 01:25:25.503515 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 01:25:26.210757 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 01:25:27.224076 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 01:25:29.462692 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 01:25:31.159379 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 01:25:32.612407 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 01:25:33.514727 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-10-09 01:25:36.471152 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 01:25:37.219084 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 01:25:38.093106 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 01:25:39.774373 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 01:25:40.389856 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 01:25:40.641464 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 01:25:40.974414 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 01:25:41.246062 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 01:25:41.549669 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 01:25:41.866812 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 01:25:42.188617 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 01:25:42.750396 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 01:25:43.302578 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 01:25:44.181852 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 01:25:44.755667 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 01:25:45.743329 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 01:25:46.691415 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 01:25:47.383345 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 01:25:47.920382 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 01:25:48.376754 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 01:25:48.988693 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 01:25:49.70848 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 01:25:50.600078 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 01:25:51.1985 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 01:25:51.785063 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 01:25:52.493907 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 01:25:53.501732 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 01:25:55.09792 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 01:25:56.932403 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 01:25:58.523216 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 01:26:00.24534 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-09 01:26:02.979815 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 01:26:03.716773 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 01:26:04.604834 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 01:26:06.389665 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 01:26:07.014248 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 01:26:07.262605 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 01:26:07.577655 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 01:26:07.840182 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 01:26:08.140216 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 01:26:08.468204 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 01:26:08.772097 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 01:26:09.326453 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 01:26:09.89348 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 01:26:10.766176 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 01:26:11.324149 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 01:26:12.3124 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 01:26:13.265476 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 01:26:13.947056 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 01:26:14.482929 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 01:26:14.941125 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 01:26:15.527706 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 01:26:16.246214 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 01:26:17.107698 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 01:26:17.716883 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 01:26:18.292855 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 01:26:18.998974 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 01:26:20.003237 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 01:26:21.588343 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 01:26:23.45999 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 01:26:25.054276 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 01:26:26.050607 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-09 01:26:29.090814 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 01:26:30.663367 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 01:26:31.465121 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 01:26:32.938612 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 01:26:33.500502 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 01:26:33.730562 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 01:26:34.019752 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 01:26:34.263979 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 01:26:34.541032 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 01:26:34.829553 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 01:26:35.10867 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 01:26:35.617895 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 01:26:36.122101 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 01:26:36.908924 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 01:26:37.420463 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 01:26:38.326633 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 01:26:39.177459 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 01:26:39.812713 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 01:26:40.295926 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 01:26:40.719182 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 01:26:41.281183 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 01:26:41.93487 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 01:26:42.730131 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 01:26:43.273608 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 01:26:43.787056 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 01:26:44.414048 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 01:26:45.318838 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 01:26:46.747882 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 01:26:48.384824 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 01:26:49.817419 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 01:26:50.690392 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
