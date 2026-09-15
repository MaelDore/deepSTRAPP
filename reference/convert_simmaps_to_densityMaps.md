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
#> 2026-09-15 09:04:55.884721 - Fit 1 evolutionary model(s): ER.
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
#> 2026-09-15 09:04:59.52268 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-09-15 09:04:59.524026 - Run simulations for stochastic mapping.
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
#> 2026-09-15 09:06:43.570449 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-09-15 09:06:46.349753 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-09-15 09:06:49.616311 - Posterior probability computed for edge n°100/3066
#> 2026-09-15 09:06:50.514697 - Posterior probability computed for edge n°200/3066
#> 2026-09-15 09:06:51.58143 - Posterior probability computed for edge n°300/3066
#> 2026-09-15 09:06:53.528496 - Posterior probability computed for edge n°400/3066
#> 2026-09-15 09:06:54.2545 - Posterior probability computed for edge n°500/3066
#> 2026-09-15 09:06:55.073131 - Posterior probability computed for edge n°600/3066
#> 2026-09-15 09:06:55.402561 - Posterior probability computed for edge n°700/3066
#> 2026-09-15 09:06:55.662953 - Posterior probability computed for edge n°800/3066
#> 2026-09-15 09:06:55.957399 - Posterior probability computed for edge n°900/3066
#> 2026-09-15 09:06:56.266168 - Posterior probability computed for edge n°1000/3066
#> 2026-09-15 09:06:56.565645 - Posterior probability computed for edge n°1100/3066
#> 2026-09-15 09:06:57.116274 - Posterior probability computed for edge n°1200/3066
#> 2026-09-15 09:06:57.656368 - Posterior probability computed for edge n°1300/3066
#> 2026-09-15 09:06:58.491912 - Posterior probability computed for edge n°1400/3066
#> 2026-09-15 09:06:59.040744 - Posterior probability computed for edge n°1500/3066
#> 2026-09-15 09:06:59.990915 - Posterior probability computed for edge n°1600/3066
#> 2026-09-15 09:07:00.91003 - Posterior probability computed for edge n°1700/3066
#> 2026-09-15 09:07:01.564889 - Posterior probability computed for edge n°1800/3066
#> 2026-09-15 09:07:02.089968 - Posterior probability computed for edge n°1900/3066
#> 2026-09-15 09:07:02.515994 - Posterior probability computed for edge n°2000/3066
#> 2026-09-15 09:07:03.112483 - Posterior probability computed for edge n°2100/3066
#> 2026-09-15 09:07:03.793477 - Posterior probability computed for edge n°2200/3066
#> 2026-09-15 09:07:04.643475 - Posterior probability computed for edge n°2300/3066
#> 2026-09-15 09:07:05.231744 - Posterior probability computed for edge n°2400/3066
#> 2026-09-15 09:07:05.783226 - Posterior probability computed for edge n°2500/3066
#> 2026-09-15 09:07:06.453887 - Posterior probability computed for edge n°2600/3066
#> 2026-09-15 09:07:07.446231 - Posterior probability computed for edge n°2700/3066
#> 2026-09-15 09:07:08.976871 - Posterior probability computed for edge n°2800/3066
#> 2026-09-15 09:07:10.723595 - Posterior probability computed for edge n°2900/3066
#> 2026-09-15 09:07:12.230511 - Posterior probability computed for edge n°3000/3066
#> 2026-09-15 09:07:13.180778 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-15 09:07:15.947098 - Posterior probability computed for edge n°100/3066
#> 2026-09-15 09:07:16.695103 - Posterior probability computed for edge n°200/3066
#> 2026-09-15 09:07:17.603046 - Posterior probability computed for edge n°300/3066
#> 2026-09-15 09:07:19.402958 - Posterior probability computed for edge n°400/3066
#> 2026-09-15 09:07:20.036991 - Posterior probability computed for edge n°500/3066
#> 2026-09-15 09:07:20.293094 - Posterior probability computed for edge n°600/3066
#> 2026-09-15 09:07:20.619602 - Posterior probability computed for edge n°700/3066
#> 2026-09-15 09:07:20.891514 - Posterior probability computed for edge n°800/3066
#> 2026-09-15 09:07:21.197126 - Posterior probability computed for edge n°900/3066
#> 2026-09-15 09:07:21.532545 - Posterior probability computed for edge n°1000/3066
#> 2026-09-15 09:07:21.85067 - Posterior probability computed for edge n°1100/3066
#> 2026-09-15 09:07:22.417606 - Posterior probability computed for edge n°1200/3066
#> 2026-09-15 09:07:22.994102 - Posterior probability computed for edge n°1300/3066
#> 2026-09-15 09:07:23.872234 - Posterior probability computed for edge n°1400/3066
#> 2026-09-15 09:07:24.461196 - Posterior probability computed for edge n°1500/3066
#> 2026-09-15 09:07:25.468376 - Posterior probability computed for edge n°1600/3066
#> 2026-09-15 09:07:26.432855 - Posterior probability computed for edge n°1700/3066
#> 2026-09-15 09:07:27.129383 - Posterior probability computed for edge n°1800/3066
#> 2026-09-15 09:07:27.674557 - Posterior probability computed for edge n°1900/3066
#> 2026-09-15 09:07:28.141974 - Posterior probability computed for edge n°2000/3066
#> 2026-09-15 09:07:29.376208 - Posterior probability computed for edge n°2100/3066
#> 2026-09-15 09:07:30.042505 - Posterior probability computed for edge n°2200/3066
#> 2026-09-15 09:07:30.859213 - Posterior probability computed for edge n°2300/3066
#> 2026-09-15 09:07:31.424116 - Posterior probability computed for edge n°2400/3066
#> 2026-09-15 09:07:31.956002 - Posterior probability computed for edge n°2500/3066
#> 2026-09-15 09:07:32.606338 - Posterior probability computed for edge n°2600/3066
#> 2026-09-15 09:07:33.547673 - Posterior probability computed for edge n°2700/3066
#> 2026-09-15 09:07:35.024947 - Posterior probability computed for edge n°2800/3066
#> 2026-09-15 09:07:36.744907 - Posterior probability computed for edge n°2900/3066
#> 2026-09-15 09:07:38.227781 - Posterior probability computed for edge n°3000/3066
#> 2026-09-15 09:07:39.155128 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-15 09:07:41.959558 - Posterior probability computed for edge n°100/3066
#> 2026-09-15 09:07:42.698254 - Posterior probability computed for edge n°200/3066
#> 2026-09-15 09:07:43.597602 - Posterior probability computed for edge n°300/3066
#> 2026-09-15 09:07:45.244212 - Posterior probability computed for edge n°400/3066
#> 2026-09-15 09:07:45.924982 - Posterior probability computed for edge n°500/3066
#> 2026-09-15 09:07:46.17795 - Posterior probability computed for edge n°600/3066
#> 2026-09-15 09:07:46.496419 - Posterior probability computed for edge n°700/3066
#> 2026-09-15 09:07:46.76422 - Posterior probability computed for edge n°800/3066
#> 2026-09-15 09:07:47.085821 - Posterior probability computed for edge n°900/3066
#> 2026-09-15 09:07:47.403288 - Posterior probability computed for edge n°1000/3066
#> 2026-09-15 09:07:47.752697 - Posterior probability computed for edge n°1100/3066
#> 2026-09-15 09:07:48.330912 - Posterior probability computed for edge n°1200/3066
#> 2026-09-15 09:07:48.887678 - Posterior probability computed for edge n°1300/3066
#> 2026-09-15 09:07:49.76564 - Posterior probability computed for edge n°1400/3066
#> 2026-09-15 09:07:50.32741 - Posterior probability computed for edge n°1500/3066
#> 2026-09-15 09:07:51.341187 - Posterior probability computed for edge n°1600/3066
#> 2026-09-15 09:07:52.300519 - Posterior probability computed for edge n°1700/3066
#> 2026-09-15 09:07:52.988729 - Posterior probability computed for edge n°1800/3066
#> 2026-09-15 09:07:53.533127 - Posterior probability computed for edge n°1900/3066
#> 2026-09-15 09:07:53.995611 - Posterior probability computed for edge n°2000/3066
#> 2026-09-15 09:07:54.595692 - Posterior probability computed for edge n°2100/3066
#> 2026-09-15 09:07:55.330328 - Posterior probability computed for edge n°2200/3066
#> 2026-09-15 09:07:56.210664 - Posterior probability computed for edge n°2300/3066
#> 2026-09-15 09:07:56.831914 - Posterior probability computed for edge n°2400/3066
#> 2026-09-15 09:07:57.400274 - Posterior probability computed for edge n°2500/3066
#> 2026-09-15 09:07:58.112994 - Posterior probability computed for edge n°2600/3066
#> 2026-09-15 09:07:59.139106 - Posterior probability computed for edge n°2700/3066
#> 2026-09-15 09:08:01.358224 - Posterior probability computed for edge n°2800/3066
#> 2026-09-15 09:08:03.046046 - Posterior probability computed for edge n°2900/3066
#> 2026-09-15 09:08:04.502928 - Posterior probability computed for edge n°3000/3066
#> 2026-09-15 09:08:05.413893 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-09-15 09:08:08.313819 - Posterior probability computed for edge n°100/3066
#> 2026-09-15 09:08:09.044113 - Posterior probability computed for edge n°200/3066
#> 2026-09-15 09:08:09.931037 - Posterior probability computed for edge n°300/3066
#> 2026-09-15 09:08:11.607664 - Posterior probability computed for edge n°400/3066
#> 2026-09-15 09:08:12.218269 - Posterior probability computed for edge n°500/3066
#> 2026-09-15 09:08:12.467638 - Posterior probability computed for edge n°600/3066
#> 2026-09-15 09:08:12.794327 - Posterior probability computed for edge n°700/3066
#> 2026-09-15 09:08:13.060933 - Posterior probability computed for edge n°800/3066
#> 2026-09-15 09:08:13.361603 - Posterior probability computed for edge n°900/3066
#> 2026-09-15 09:08:13.676015 - Posterior probability computed for edge n°1000/3066
#> 2026-09-15 09:08:13.99384 - Posterior probability computed for edge n°1100/3066
#> 2026-09-15 09:08:14.55086 - Posterior probability computed for edge n°1200/3066
#> 2026-09-15 09:08:15.096921 - Posterior probability computed for edge n°1300/3066
#> 2026-09-15 09:08:15.970544 - Posterior probability computed for edge n°1400/3066
#> 2026-09-15 09:08:16.538916 - Posterior probability computed for edge n°1500/3066
#> 2026-09-15 09:08:17.525134 - Posterior probability computed for edge n°1600/3066
#> 2026-09-15 09:08:18.466546 - Posterior probability computed for edge n°1700/3066
#> 2026-09-15 09:08:19.150895 - Posterior probability computed for edge n°1800/3066
#> 2026-09-15 09:08:19.688034 - Posterior probability computed for edge n°1900/3066
#> 2026-09-15 09:08:20.143313 - Posterior probability computed for edge n°2000/3066
#> 2026-09-15 09:08:20.750242 - Posterior probability computed for edge n°2100/3066
#> 2026-09-15 09:08:21.462389 - Posterior probability computed for edge n°2200/3066
#> 2026-09-15 09:08:22.345352 - Posterior probability computed for edge n°2300/3066
#> 2026-09-15 09:08:22.941931 - Posterior probability computed for edge n°2400/3066
#> 2026-09-15 09:08:23.51843 - Posterior probability computed for edge n°2500/3066
#> 2026-09-15 09:08:24.220282 - Posterior probability computed for edge n°2600/3066
#> 2026-09-15 09:08:25.214309 - Posterior probability computed for edge n°2700/3066
#> 2026-09-15 09:08:26.79927 - Posterior probability computed for edge n°2800/3066
#> 2026-09-15 09:08:28.623318 - Posterior probability computed for edge n°2900/3066
#> 2026-09-15 09:08:30.199646 - Posterior probability computed for edge n°3000/3066
#> 2026-09-15 09:08:31.828287 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-15 09:08:34.469939 - Posterior probability computed for edge n°100/3066
#> 2026-09-15 09:08:35.204102 - Posterior probability computed for edge n°200/3066
#> 2026-09-15 09:08:36.084327 - Posterior probability computed for edge n°300/3066
#> 2026-09-15 09:08:37.848785 - Posterior probability computed for edge n°400/3066
#> 2026-09-15 09:08:38.470218 - Posterior probability computed for edge n°500/3066
#> 2026-09-15 09:08:38.720356 - Posterior probability computed for edge n°600/3066
#> 2026-09-15 09:08:39.036147 - Posterior probability computed for edge n°700/3066
#> 2026-09-15 09:08:39.301569 - Posterior probability computed for edge n°800/3066
#> 2026-09-15 09:08:39.601764 - Posterior probability computed for edge n°900/3066
#> 2026-09-15 09:08:39.930397 - Posterior probability computed for edge n°1000/3066
#> 2026-09-15 09:08:40.234291 - Posterior probability computed for edge n°1100/3066
#> 2026-09-15 09:08:40.786747 - Posterior probability computed for edge n°1200/3066
#> 2026-09-15 09:08:41.351621 - Posterior probability computed for edge n°1300/3066
#> 2026-09-15 09:08:42.226811 - Posterior probability computed for edge n°1400/3066
#> 2026-09-15 09:08:42.780621 - Posterior probability computed for edge n°1500/3066
#> 2026-09-15 09:08:43.761182 - Posterior probability computed for edge n°1600/3066
#> 2026-09-15 09:08:44.716444 - Posterior probability computed for edge n°1700/3066
#> 2026-09-15 09:08:45.398735 - Posterior probability computed for edge n°1800/3066
#> 2026-09-15 09:08:45.937493 - Posterior probability computed for edge n°1900/3066
#> 2026-09-15 09:08:46.392456 - Posterior probability computed for edge n°2000/3066
#> 2026-09-15 09:08:46.980838 - Posterior probability computed for edge n°2100/3066
#> 2026-09-15 09:08:47.698277 - Posterior probability computed for edge n°2200/3066
#> 2026-09-15 09:08:48.562063 - Posterior probability computed for edge n°2300/3066
#> 2026-09-15 09:08:49.172862 - Posterior probability computed for edge n°2400/3066
#> 2026-09-15 09:08:49.746935 - Posterior probability computed for edge n°2500/3066
#> 2026-09-15 09:08:50.449571 - Posterior probability computed for edge n°2600/3066
#> 2026-09-15 09:08:51.448632 - Posterior probability computed for edge n°2700/3066
#> 2026-09-15 09:08:53.027234 - Posterior probability computed for edge n°2800/3066
#> 2026-09-15 09:08:54.840137 - Posterior probability computed for edge n°2900/3066
#> 2026-09-15 09:08:56.416397 - Posterior probability computed for edge n°3000/3066
#> 2026-09-15 09:08:57.384721 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-15 09:09:00.363224 - Posterior probability computed for edge n°100/3066
#> 2026-09-15 09:09:01.919941 - Posterior probability computed for edge n°200/3066
#> 2026-09-15 09:09:02.710713 - Posterior probability computed for edge n°300/3066
#> 2026-09-15 09:09:04.177863 - Posterior probability computed for edge n°400/3066
#> 2026-09-15 09:09:04.736167 - Posterior probability computed for edge n°500/3066
#> 2026-09-15 09:09:04.966639 - Posterior probability computed for edge n°600/3066
#> 2026-09-15 09:09:05.254941 - Posterior probability computed for edge n°700/3066
#> 2026-09-15 09:09:05.497591 - Posterior probability computed for edge n°800/3066
#> 2026-09-15 09:09:05.773921 - Posterior probability computed for edge n°900/3066
#> 2026-09-15 09:09:06.06123 - Posterior probability computed for edge n°1000/3066
#> 2026-09-15 09:09:06.339499 - Posterior probability computed for edge n°1100/3066
#> 2026-09-15 09:09:06.847608 - Posterior probability computed for edge n°1200/3066
#> 2026-09-15 09:09:07.34519 - Posterior probability computed for edge n°1300/3066
#> 2026-09-15 09:09:08.125514 - Posterior probability computed for edge n°1400/3066
#> 2026-09-15 09:09:08.633499 - Posterior probability computed for edge n°1500/3066
#> 2026-09-15 09:09:09.533306 - Posterior probability computed for edge n°1600/3066
#> 2026-09-15 09:09:10.375234 - Posterior probability computed for edge n°1700/3066
#> 2026-09-15 09:09:11.004824 - Posterior probability computed for edge n°1800/3066
#> 2026-09-15 09:09:11.47726 - Posterior probability computed for edge n°1900/3066
#> 2026-09-15 09:09:11.892667 - Posterior probability computed for edge n°2000/3066
#> 2026-09-15 09:09:12.442609 - Posterior probability computed for edge n°2100/3066
#> 2026-09-15 09:09:13.08384 - Posterior probability computed for edge n°2200/3066
#> 2026-09-15 09:09:13.867863 - Posterior probability computed for edge n°2300/3066
#> 2026-09-15 09:09:14.409269 - Posterior probability computed for edge n°2400/3066
#> 2026-09-15 09:09:14.91981 - Posterior probability computed for edge n°2500/3066
#> 2026-09-15 09:09:15.5424 - Posterior probability computed for edge n°2600/3066
#> 2026-09-15 09:09:16.434563 - Posterior probability computed for edge n°2700/3066
#> 2026-09-15 09:09:17.847136 - Posterior probability computed for edge n°2800/3066
#> 2026-09-15 09:09:19.460085 - Posterior probability computed for edge n°2900/3066
#> 2026-09-15 09:09:20.870445 - Posterior probability computed for edge n°3000/3066
#> 2026-09-15 09:09:21.728714 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
