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
#> 2026-10-21 08:00:13.234784 - Fit 1 evolutionary model(s): ER.
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
#> 2026-10-21 08:00:15.696525 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-10-21 08:00:15.697396 - Run simulations for stochastic mapping.
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
#> 2026-10-21 08:01:28.155542 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-10-21 08:01:30.157222 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-10-21 08:01:32.909104 - Posterior probability computed for edge n°100/3066
#> 2026-10-21 08:01:33.479525 - Posterior probability computed for edge n°200/3066
#> 2026-10-21 08:01:34.108414 - Posterior probability computed for edge n°300/3066
#> 2026-10-21 08:01:35.308139 - Posterior probability computed for edge n°400/3066
#> 2026-10-21 08:01:35.746982 - Posterior probability computed for edge n°500/3066
#> 2026-10-21 08:01:35.925414 - Posterior probability computed for edge n°600/3066
#> 2026-10-21 08:01:36.163507 - Posterior probability computed for edge n°700/3066
#> 2026-10-21 08:01:36.360889 - Posterior probability computed for edge n°800/3066
#> 2026-10-21 08:01:36.580949 - Posterior probability computed for edge n°900/3066
#> 2026-10-21 08:01:36.804751 - Posterior probability computed for edge n°1000/3066
#> 2026-10-21 08:01:37.042223 - Posterior probability computed for edge n°1100/3066
#> 2026-10-21 08:01:37.434773 - Posterior probability computed for edge n°1200/3066
#> 2026-10-21 08:01:37.81753 - Posterior probability computed for edge n°1300/3066
#> 2026-10-21 08:01:38.420779 - Posterior probability computed for edge n°1400/3066
#> 2026-10-21 08:01:38.823348 - Posterior probability computed for edge n°1500/3066
#> 2026-10-21 08:01:39.521704 - Posterior probability computed for edge n°1600/3066
#> 2026-10-21 08:01:40.201763 - Posterior probability computed for edge n°1700/3066
#> 2026-10-21 08:01:40.700397 - Posterior probability computed for edge n°1800/3066
#> 2026-10-21 08:01:41.132274 - Posterior probability computed for edge n°1900/3066
#> 2026-10-21 08:01:41.459573 - Posterior probability computed for edge n°2000/3066
#> 2026-10-21 08:01:41.935773 - Posterior probability computed for edge n°2100/3066
#> 2026-10-21 08:01:42.484131 - Posterior probability computed for edge n°2200/3066
#> 2026-10-21 08:01:43.116631 - Posterior probability computed for edge n°2300/3066
#> 2026-10-21 08:01:43.573096 - Posterior probability computed for edge n°2400/3066
#> 2026-10-21 08:01:43.976305 - Posterior probability computed for edge n°2500/3066
#> 2026-10-21 08:01:44.48136 - Posterior probability computed for edge n°2600/3066
#> 2026-10-21 08:01:45.199374 - Posterior probability computed for edge n°2700/3066
#> 2026-10-21 08:01:46.345719 - Posterior probability computed for edge n°2800/3066
#> 2026-10-21 08:01:47.671908 - Posterior probability computed for edge n°2900/3066
#> 2026-10-21 08:01:48.807967 - Posterior probability computed for edge n°3000/3066
#> 2026-10-21 08:01:49.512623 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-21 08:01:51.653832 - Posterior probability computed for edge n°100/3066
#> 2026-10-21 08:01:52.224839 - Posterior probability computed for edge n°200/3066
#> 2026-10-21 08:01:53.494061 - Posterior probability computed for edge n°300/3066
#> 2026-10-21 08:01:54.57765 - Posterior probability computed for edge n°400/3066
#> 2026-10-21 08:01:55.00336 - Posterior probability computed for edge n°500/3066
#> 2026-10-21 08:01:55.177159 - Posterior probability computed for edge n°600/3066
#> 2026-10-21 08:01:55.394748 - Posterior probability computed for edge n°700/3066
#> 2026-10-21 08:01:55.574656 - Posterior probability computed for edge n°800/3066
#> 2026-10-21 08:01:55.781211 - Posterior probability computed for edge n°900/3066
#> 2026-10-21 08:01:55.995678 - Posterior probability computed for edge n°1000/3066
#> 2026-10-21 08:01:56.204875 - Posterior probability computed for edge n°1100/3066
#> 2026-10-21 08:01:56.606069 - Posterior probability computed for edge n°1200/3066
#> 2026-10-21 08:01:56.987677 - Posterior probability computed for edge n°1300/3066
#> 2026-10-21 08:01:57.595355 - Posterior probability computed for edge n°1400/3066
#> 2026-10-21 08:01:57.996027 - Posterior probability computed for edge n°1500/3066
#> 2026-10-21 08:01:58.686271 - Posterior probability computed for edge n°1600/3066
#> 2026-10-21 08:01:59.351026 - Posterior probability computed for edge n°1700/3066
#> 2026-10-21 08:01:59.824524 - Posterior probability computed for edge n°1800/3066
#> 2026-10-21 08:02:00.196052 - Posterior probability computed for edge n°1900/3066
#> 2026-10-21 08:02:00.514005 - Posterior probability computed for edge n°2000/3066
#> 2026-10-21 08:02:00.939118 - Posterior probability computed for edge n°2100/3066
#> 2026-10-21 08:02:01.432999 - Posterior probability computed for edge n°2200/3066
#> 2026-10-21 08:02:02.106737 - Posterior probability computed for edge n°2300/3066
#> 2026-10-21 08:02:02.534136 - Posterior probability computed for edge n°2400/3066
#> 2026-10-21 08:02:02.983408 - Posterior probability computed for edge n°2500/3066
#> 2026-10-21 08:02:03.483577 - Posterior probability computed for edge n°2600/3066
#> 2026-10-21 08:02:04.181229 - Posterior probability computed for edge n°2700/3066
#> 2026-10-21 08:02:05.295395 - Posterior probability computed for edge n°2800/3066
#> 2026-10-21 08:02:06.584347 - Posterior probability computed for edge n°2900/3066
#> 2026-10-21 08:02:07.692613 - Posterior probability computed for edge n°3000/3066
#> 2026-10-21 08:02:08.368009 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-21 08:02:10.627691 - Posterior probability computed for edge n°100/3066
#> 2026-10-21 08:02:11.204258 - Posterior probability computed for edge n°200/3066
#> 2026-10-21 08:02:11.92005 - Posterior probability computed for edge n°300/3066
#> 2026-10-21 08:02:13.767901 - Posterior probability computed for edge n°400/3066
#> 2026-10-21 08:02:14.180003 - Posterior probability computed for edge n°500/3066
#> 2026-10-21 08:02:14.354984 - Posterior probability computed for edge n°600/3066
#> 2026-10-21 08:02:14.565974 - Posterior probability computed for edge n°700/3066
#> 2026-10-21 08:02:14.748303 - Posterior probability computed for edge n°800/3066
#> 2026-10-21 08:02:14.959201 - Posterior probability computed for edge n°900/3066
#> 2026-10-21 08:02:15.178488 - Posterior probability computed for edge n°1000/3066
#> 2026-10-21 08:02:15.394361 - Posterior probability computed for edge n°1100/3066
#> 2026-10-21 08:02:15.78256 - Posterior probability computed for edge n°1200/3066
#> 2026-10-21 08:02:16.171832 - Posterior probability computed for edge n°1300/3066
#> 2026-10-21 08:02:16.771777 - Posterior probability computed for edge n°1400/3066
#> 2026-10-21 08:02:17.171041 - Posterior probability computed for edge n°1500/3066
#> 2026-10-21 08:02:17.85403 - Posterior probability computed for edge n°1600/3066
#> 2026-10-21 08:02:18.517579 - Posterior probability computed for edge n°1700/3066
#> 2026-10-21 08:02:19.054165 - Posterior probability computed for edge n°1800/3066
#> 2026-10-21 08:02:19.440083 - Posterior probability computed for edge n°1900/3066
#> 2026-10-21 08:02:19.760682 - Posterior probability computed for edge n°2000/3066
#> 2026-10-21 08:02:20.17943 - Posterior probability computed for edge n°2100/3066
#> 2026-10-21 08:02:20.67057 - Posterior probability computed for edge n°2200/3066
#> 2026-10-21 08:02:21.275472 - Posterior probability computed for edge n°2300/3066
#> 2026-10-21 08:02:21.686347 - Posterior probability computed for edge n°2400/3066
#> 2026-10-21 08:02:22.101166 - Posterior probability computed for edge n°2500/3066
#> 2026-10-21 08:02:22.596933 - Posterior probability computed for edge n°2600/3066
#> 2026-10-21 08:02:23.373951 - Posterior probability computed for edge n°2700/3066
#> 2026-10-21 08:02:24.504485 - Posterior probability computed for edge n°2800/3066
#> 2026-10-21 08:02:25.781085 - Posterior probability computed for edge n°2900/3066
#> 2026-10-21 08:02:26.93507 - Posterior probability computed for edge n°3000/3066
#> 2026-10-21 08:02:27.623478 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-10-21 08:02:29.83881 - Posterior probability computed for edge n°100/3066
#> 2026-10-21 08:02:30.557566 - Posterior probability computed for edge n°200/3066
#> 2026-10-21 08:02:31.276637 - Posterior probability computed for edge n°300/3066
#> 2026-10-21 08:02:32.596403 - Posterior probability computed for edge n°400/3066
#> 2026-10-21 08:02:33.810581 - Posterior probability computed for edge n°500/3066
#> 2026-10-21 08:02:33.993546 - Posterior probability computed for edge n°600/3066
#> 2026-10-21 08:02:34.217536 - Posterior probability computed for edge n°700/3066
#> 2026-10-21 08:02:34.40977 - Posterior probability computed for edge n°800/3066
#> 2026-10-21 08:02:34.619787 - Posterior probability computed for edge n°900/3066
#> 2026-10-21 08:02:34.833026 - Posterior probability computed for edge n°1000/3066
#> 2026-10-21 08:02:35.05522 - Posterior probability computed for edge n°1100/3066
#> 2026-10-21 08:02:35.454126 - Posterior probability computed for edge n°1200/3066
#> 2026-10-21 08:02:35.890049 - Posterior probability computed for edge n°1300/3066
#> 2026-10-21 08:02:36.574616 - Posterior probability computed for edge n°1400/3066
#> 2026-10-21 08:02:36.992968 - Posterior probability computed for edge n°1500/3066
#> 2026-10-21 08:02:37.705619 - Posterior probability computed for edge n°1600/3066
#> 2026-10-21 08:02:38.395056 - Posterior probability computed for edge n°1700/3066
#> 2026-10-21 08:02:38.885891 - Posterior probability computed for edge n°1800/3066
#> 2026-10-21 08:02:39.27655 - Posterior probability computed for edge n°1900/3066
#> 2026-10-21 08:02:39.608793 - Posterior probability computed for edge n°2000/3066
#> 2026-10-21 08:02:40.033547 - Posterior probability computed for edge n°2100/3066
#> 2026-10-21 08:02:40.53097 - Posterior probability computed for edge n°2200/3066
#> 2026-10-21 08:02:41.161487 - Posterior probability computed for edge n°2300/3066
#> 2026-10-21 08:02:41.590393 - Posterior probability computed for edge n°2400/3066
#> 2026-10-21 08:02:41.996874 - Posterior probability computed for edge n°2500/3066
#> 2026-10-21 08:02:42.513234 - Posterior probability computed for edge n°2600/3066
#> 2026-10-21 08:02:43.259747 - Posterior probability computed for edge n°2700/3066
#> 2026-10-21 08:02:44.49784 - Posterior probability computed for edge n°2800/3066
#> 2026-10-21 08:02:45.849227 - Posterior probability computed for edge n°2900/3066
#> 2026-10-21 08:02:47.001473 - Posterior probability computed for edge n°3000/3066
#> 2026-10-21 08:02:47.752464 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-21 08:02:50.100831 - Posterior probability computed for edge n°100/3066
#> 2026-10-21 08:02:50.723948 - Posterior probability computed for edge n°200/3066
#> 2026-10-21 08:02:51.640419 - Posterior probability computed for edge n°300/3066
#> 2026-10-21 08:02:53.026798 - Posterior probability computed for edge n°400/3066
#> 2026-10-21 08:02:53.547256 - Posterior probability computed for edge n°500/3066
#> 2026-10-21 08:02:53.773295 - Posterior probability computed for edge n°600/3066
#> 2026-10-21 08:02:54.042932 - Posterior probability computed for edge n°700/3066
#> 2026-10-21 08:02:54.923241 - Posterior probability computed for edge n°800/3066
#> 2026-10-21 08:02:55.138416 - Posterior probability computed for edge n°900/3066
#> 2026-10-21 08:02:55.364078 - Posterior probability computed for edge n°1000/3066
#> 2026-10-21 08:02:55.576763 - Posterior probability computed for edge n°1100/3066
#> 2026-10-21 08:02:55.968902 - Posterior probability computed for edge n°1200/3066
#> 2026-10-21 08:02:56.357194 - Posterior probability computed for edge n°1300/3066
#> 2026-10-21 08:02:56.965675 - Posterior probability computed for edge n°1400/3066
#> 2026-10-21 08:02:57.354635 - Posterior probability computed for edge n°1500/3066
#> 2026-10-21 08:02:58.057428 - Posterior probability computed for edge n°1600/3066
#> 2026-10-21 08:02:58.74088 - Posterior probability computed for edge n°1700/3066
#> 2026-10-21 08:02:59.231464 - Posterior probability computed for edge n°1800/3066
#> 2026-10-21 08:02:59.619614 - Posterior probability computed for edge n°1900/3066
#> 2026-10-21 08:02:59.982363 - Posterior probability computed for edge n°2000/3066
#> 2026-10-21 08:03:00.412088 - Posterior probability computed for edge n°2100/3066
#> 2026-10-21 08:03:00.935312 - Posterior probability computed for edge n°2200/3066
#> 2026-10-21 08:03:01.561891 - Posterior probability computed for edge n°2300/3066
#> 2026-10-21 08:03:02.001997 - Posterior probability computed for edge n°2400/3066
#> 2026-10-21 08:03:02.411134 - Posterior probability computed for edge n°2500/3066
#> 2026-10-21 08:03:02.934426 - Posterior probability computed for edge n°2600/3066
#> 2026-10-21 08:03:03.648797 - Posterior probability computed for edge n°2700/3066
#> 2026-10-21 08:03:04.758925 - Posterior probability computed for edge n°2800/3066
#> 2026-10-21 08:03:06.078664 - Posterior probability computed for edge n°2900/3066
#> 2026-10-21 08:03:07.276202 - Posterior probability computed for edge n°3000/3066
#> 2026-10-21 08:03:07.972721 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-21 08:03:10.262868 - Posterior probability computed for edge n°100/3066
#> 2026-10-21 08:03:10.888655 - Posterior probability computed for edge n°200/3066
#> 2026-10-21 08:03:11.827061 - Posterior probability computed for edge n°300/3066
#> 2026-10-21 08:03:13.204556 - Posterior probability computed for edge n°400/3066
#> 2026-10-21 08:03:13.723077 - Posterior probability computed for edge n°500/3066
#> 2026-10-21 08:03:13.93188 - Posterior probability computed for edge n°600/3066
#> 2026-10-21 08:03:14.197204 - Posterior probability computed for edge n°700/3066
#> 2026-10-21 08:03:14.420471 - Posterior probability computed for edge n°800/3066
#> 2026-10-21 08:03:15.369086 - Posterior probability computed for edge n°900/3066
#> 2026-10-21 08:03:15.60271 - Posterior probability computed for edge n°1000/3066
#> 2026-10-21 08:03:15.819104 - Posterior probability computed for edge n°1100/3066
#> 2026-10-21 08:03:16.206718 - Posterior probability computed for edge n°1200/3066
#> 2026-10-21 08:03:16.592193 - Posterior probability computed for edge n°1300/3066
#> 2026-10-21 08:03:17.190558 - Posterior probability computed for edge n°1400/3066
#> 2026-10-21 08:03:17.610142 - Posterior probability computed for edge n°1500/3066
#> 2026-10-21 08:03:18.313267 - Posterior probability computed for edge n°1600/3066
#> 2026-10-21 08:03:18.976007 - Posterior probability computed for edge n°1700/3066
#> 2026-10-21 08:03:19.45895 - Posterior probability computed for edge n°1800/3066
#> 2026-10-21 08:03:19.829885 - Posterior probability computed for edge n°1900/3066
#> 2026-10-21 08:03:20.141863 - Posterior probability computed for edge n°2000/3066
#> 2026-10-21 08:03:20.564842 - Posterior probability computed for edge n°2100/3066
#> 2026-10-21 08:03:21.067934 - Posterior probability computed for edge n°2200/3066
#> 2026-10-21 08:03:21.727011 - Posterior probability computed for edge n°2300/3066
#> 2026-10-21 08:03:22.139678 - Posterior probability computed for edge n°2400/3066
#> 2026-10-21 08:03:22.540497 - Posterior probability computed for edge n°2500/3066
#> 2026-10-21 08:03:23.04957 - Posterior probability computed for edge n°2600/3066
#> 2026-10-21 08:03:23.748573 - Posterior probability computed for edge n°2700/3066
#> 2026-10-21 08:03:24.84608 - Posterior probability computed for edge n°2800/3066
#> 2026-10-21 08:03:26.182768 - Posterior probability computed for edge n°2900/3066
#> 2026-10-21 08:03:27.272835 - Posterior probability computed for edge n°3000/3066
#> 2026-10-21 08:03:27.929428 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
