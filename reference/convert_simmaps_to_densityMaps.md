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
#> 2026-10-09 07:38:27.378039 - Fit 1 evolutionary model(s): ER.
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
#> 2026-10-09 07:38:30.866791 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-10-09 07:38:30.868202 - Run simulations for stochastic mapping.
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
#> 2026-10-09 07:40:51.554421 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-10-09 07:40:55.073487 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-10-09 07:40:58.931669 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 07:40:59.991895 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 07:41:01.27514 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 07:41:03.645948 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 07:41:04.555998 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 07:41:04.924563 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 07:41:05.903609 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 07:41:06.253637 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 07:41:06.645582 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 07:41:07.044786 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 07:41:07.432032 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 07:41:08.13534 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 07:41:08.823547 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 07:41:09.896638 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 07:41:10.594921 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 07:41:11.826453 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 07:41:12.992556 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 07:41:13.848098 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 07:41:14.524414 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 07:41:15.082403 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 07:41:15.842589 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 07:41:16.711269 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 07:41:17.793544 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 07:41:18.538215 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 07:41:19.243428 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 07:41:20.103134 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 07:41:21.342032 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 07:41:23.269579 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 07:41:25.500267 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 07:41:27.419743 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 07:41:28.628697 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-09 07:41:31.926075 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 07:41:32.87168 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 07:41:33.972664 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 07:41:36.126675 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 07:41:36.922649 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 07:41:37.245846 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 07:41:37.650449 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 07:41:37.9951 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 07:41:38.402006 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 07:41:38.807102 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 07:41:39.199102 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 07:41:39.927513 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 07:41:40.631225 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 07:41:41.736185 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 07:41:42.44897 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 07:41:43.70169 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 07:41:44.897614 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 07:41:45.764672 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 07:41:46.449537 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 07:41:47.030283 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 07:41:47.806603 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 07:41:49.208709 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 07:41:50.247992 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 07:41:51.000027 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 07:41:51.677898 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 07:41:52.502165 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 07:41:53.678288 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 07:41:55.547033 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 07:41:57.722176 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 07:41:59.606189 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 07:42:00.763299 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-09 07:42:04.099733 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 07:42:05.045631 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 07:42:06.147634 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 07:42:08.212394 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 07:42:09.012823 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 07:42:09.345787 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 07:42:09.745098 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 07:42:10.085333 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 07:42:10.4887 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 07:42:10.890696 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 07:42:11.285236 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 07:42:12.010733 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 07:42:12.711996 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 07:42:13.813515 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 07:42:14.519481 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 07:42:15.759605 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 07:42:16.945002 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 07:42:17.809756 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 07:42:18.511051 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 07:42:19.089589 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 07:42:19.849937 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 07:42:20.745629 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 07:42:21.86581 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 07:42:22.636518 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 07:42:23.354319 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 07:42:24.246647 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 07:42:25.522867 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 07:42:27.522852 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 07:42:30.271751 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 07:42:32.144008 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 07:42:33.318391 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-10-09 07:42:36.822587 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 07:42:37.76142 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 07:42:38.887944 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 07:42:40.933546 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 07:42:41.723174 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 07:42:42.048433 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 07:42:42.45121 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 07:42:42.794854 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 07:42:43.198905 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 07:42:43.601331 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 07:42:43.989667 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 07:42:44.712532 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 07:42:45.405518 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 07:42:46.499798 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 07:42:47.208594 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 07:42:48.466066 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 07:42:49.651177 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 07:42:50.513726 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 07:42:51.202214 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 07:42:51.774683 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 07:42:52.522663 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 07:42:53.405107 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 07:42:54.506336 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 07:42:55.275144 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 07:42:55.983932 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 07:42:56.866923 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 07:42:58.119591 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 07:43:00.088488 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 07:43:02.363625 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 07:43:04.334412 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 07:43:05.555209 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-09 07:43:09.250338 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 07:43:10.125801 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 07:43:11.173604 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 07:43:13.103353 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 07:43:13.844294 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 07:43:14.151368 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 07:43:14.631426 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 07:43:14.959933 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 07:43:15.330014 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 07:43:15.713922 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 07:43:16.086949 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 07:43:16.764432 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 07:43:17.428508 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 07:43:18.465176 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 07:43:19.142805 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 07:43:20.331685 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 07:43:21.449621 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 07:43:22.260199 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 07:43:22.911205 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 07:43:23.462836 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 07:43:24.18152 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 07:43:25.020303 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 07:43:26.059812 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 07:43:26.786691 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 07:43:27.467193 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 07:43:28.315214 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 07:43:29.507902 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 07:43:31.387556 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 07:43:33.549128 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 07:43:35.403031 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 07:43:36.568687 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-09 07:43:39.748813 - Posterior probability computed for edge n°100/3066
#> 2026-10-09 07:43:40.661374 - Posterior probability computed for edge n°200/3066
#> 2026-10-09 07:43:41.778209 - Posterior probability computed for edge n°300/3066
#> 2026-10-09 07:43:43.948603 - Posterior probability computed for edge n°400/3066
#> 2026-10-09 07:43:44.746212 - Posterior probability computed for edge n°500/3066
#> 2026-10-09 07:43:45.06672 - Posterior probability computed for edge n°600/3066
#> 2026-10-09 07:43:45.466001 - Posterior probability computed for edge n°700/3066
#> 2026-10-09 07:43:45.812747 - Posterior probability computed for edge n°800/3066
#> 2026-10-09 07:43:46.218204 - Posterior probability computed for edge n°900/3066
#> 2026-10-09 07:43:46.62331 - Posterior probability computed for edge n°1000/3066
#> 2026-10-09 07:43:47.009783 - Posterior probability computed for edge n°1100/3066
#> 2026-10-09 07:43:47.738274 - Posterior probability computed for edge n°1200/3066
#> 2026-10-09 07:43:48.433497 - Posterior probability computed for edge n°1300/3066
#> 2026-10-09 07:43:50.066752 - Posterior probability computed for edge n°1400/3066
#> 2026-10-09 07:43:50.74313 - Posterior probability computed for edge n°1500/3066
#> 2026-10-09 07:43:51.933559 - Posterior probability computed for edge n°1600/3066
#> 2026-10-09 07:43:53.055771 - Posterior probability computed for edge n°1700/3066
#> 2026-10-09 07:43:53.871975 - Posterior probability computed for edge n°1800/3066
#> 2026-10-09 07:43:54.525894 - Posterior probability computed for edge n°1900/3066
#> 2026-10-09 07:43:55.083363 - Posterior probability computed for edge n°2000/3066
#> 2026-10-09 07:43:55.80888 - Posterior probability computed for edge n°2100/3066
#> 2026-10-09 07:43:56.673101 - Posterior probability computed for edge n°2200/3066
#> 2026-10-09 07:43:57.728547 - Posterior probability computed for edge n°2300/3066
#> 2026-10-09 07:43:58.473988 - Posterior probability computed for edge n°2400/3066
#> 2026-10-09 07:43:59.162182 - Posterior probability computed for edge n°2500/3066
#> 2026-10-09 07:44:00.004272 - Posterior probability computed for edge n°2600/3066
#> 2026-10-09 07:44:01.211158 - Posterior probability computed for edge n°2700/3066
#> 2026-10-09 07:44:03.094702 - Posterior probability computed for edge n°2800/3066
#> 2026-10-09 07:44:05.282659 - Posterior probability computed for edge n°2900/3066
#> 2026-10-09 07:44:07.17509 - Posterior probability computed for edge n°3000/3066
#> 2026-10-09 07:44:08.347569 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
