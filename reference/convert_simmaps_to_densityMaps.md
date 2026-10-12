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
#> 2026-10-12 00:40:23.710518 - Fit 1 evolutionary model(s): ER.
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
#> 2026-10-12 00:40:26.998207 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-10-12 00:40:26.999994 - Run simulations for stochastic mapping.
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
#> 2026-10-12 00:42:43.272902 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-10-12 00:42:46.499769 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-10-12 00:42:50.217809 - Posterior probability computed for edge n°100/3066
#> 2026-10-12 00:42:51.109417 - Posterior probability computed for edge n°200/3066
#> 2026-10-12 00:42:52.17416 - Posterior probability computed for edge n°300/3066
#> 2026-10-12 00:42:54.129125 - Posterior probability computed for edge n°400/3066
#> 2026-10-12 00:42:54.900077 - Posterior probability computed for edge n°500/3066
#> 2026-10-12 00:42:55.210349 - Posterior probability computed for edge n°600/3066
#> 2026-10-12 00:42:55.604389 - Posterior probability computed for edge n°700/3066
#> 2026-10-12 00:42:55.93395 - Posterior probability computed for edge n°800/3066
#> 2026-10-12 00:42:56.307889 - Posterior probability computed for edge n°900/3066
#> 2026-10-12 00:42:56.700901 - Posterior probability computed for edge n°1000/3066
#> 2026-10-12 00:42:57.083455 - Posterior probability computed for edge n°1100/3066
#> 2026-10-12 00:42:57.76434 - Posterior probability computed for edge n°1200/3066
#> 2026-10-12 00:42:58.438644 - Posterior probability computed for edge n°1300/3066
#> 2026-10-12 00:42:59.488326 - Posterior probability computed for edge n°1400/3066
#> 2026-10-12 00:43:00.176783 - Posterior probability computed for edge n°1500/3066
#> 2026-10-12 00:43:01.373176 - Posterior probability computed for edge n°1600/3066
#> 2026-10-12 00:43:02.5088 - Posterior probability computed for edge n°1700/3066
#> 2026-10-12 00:43:03.33507 - Posterior probability computed for edge n°1800/3066
#> 2026-10-12 00:43:03.990773 - Posterior probability computed for edge n°1900/3066
#> 2026-10-12 00:43:04.545045 - Posterior probability computed for edge n°2000/3066
#> 2026-10-12 00:43:05.282811 - Posterior probability computed for edge n°2100/3066
#> 2026-10-12 00:43:06.139738 - Posterior probability computed for edge n°2200/3066
#> 2026-10-12 00:43:07.214322 - Posterior probability computed for edge n°2300/3066
#> 2026-10-12 00:43:07.942103 - Posterior probability computed for edge n°2400/3066
#> 2026-10-12 00:43:08.635194 - Posterior probability computed for edge n°2500/3066
#> 2026-10-12 00:43:09.484769 - Posterior probability computed for edge n°2600/3066
#> 2026-10-12 00:43:10.685087 - Posterior probability computed for edge n°2700/3066
#> 2026-10-12 00:43:12.587839 - Posterior probability computed for edge n°2800/3066
#> 2026-10-12 00:43:14.775119 - Posterior probability computed for edge n°2900/3066
#> 2026-10-12 00:43:16.666273 - Posterior probability computed for edge n°3000/3066
#> 2026-10-12 00:43:17.828233 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-12 00:43:21.045118 - Posterior probability computed for edge n°100/3066
#> 2026-10-12 00:43:21.956347 - Posterior probability computed for edge n°200/3066
#> 2026-10-12 00:43:23.604891 - Posterior probability computed for edge n°300/3066
#> 2026-10-12 00:43:25.479074 - Posterior probability computed for edge n°400/3066
#> 2026-10-12 00:43:26.201735 - Posterior probability computed for edge n°500/3066
#> 2026-10-12 00:43:26.497608 - Posterior probability computed for edge n°600/3066
#> 2026-10-12 00:43:26.875435 - Posterior probability computed for edge n°700/3066
#> 2026-10-12 00:43:27.192028 - Posterior probability computed for edge n°800/3066
#> 2026-10-12 00:43:27.5523 - Posterior probability computed for edge n°900/3066
#> 2026-10-12 00:43:27.928946 - Posterior probability computed for edge n°1000/3066
#> 2026-10-12 00:43:28.29052 - Posterior probability computed for edge n°1100/3066
#> 2026-10-12 00:43:28.961592 - Posterior probability computed for edge n°1200/3066
#> 2026-10-12 00:43:29.609768 - Posterior probability computed for edge n°1300/3066
#> 2026-10-12 00:43:30.635398 - Posterior probability computed for edge n°1400/3066
#> 2026-10-12 00:43:31.299818 - Posterior probability computed for edge n°1500/3066
#> 2026-10-12 00:43:32.461155 - Posterior probability computed for edge n°1600/3066
#> 2026-10-12 00:43:33.568692 - Posterior probability computed for edge n°1700/3066
#> 2026-10-12 00:43:34.36966 - Posterior probability computed for edge n°1800/3066
#> 2026-10-12 00:43:35.007744 - Posterior probability computed for edge n°1900/3066
#> 2026-10-12 00:43:35.547246 - Posterior probability computed for edge n°2000/3066
#> 2026-10-12 00:43:36.256225 - Posterior probability computed for edge n°2100/3066
#> 2026-10-12 00:43:37.086003 - Posterior probability computed for edge n°2200/3066
#> 2026-10-12 00:43:38.111648 - Posterior probability computed for edge n°2300/3066
#> 2026-10-12 00:43:38.825947 - Posterior probability computed for edge n°2400/3066
#> 2026-10-12 00:43:39.493729 - Posterior probability computed for edge n°2500/3066
#> 2026-10-12 00:43:40.321413 - Posterior probability computed for edge n°2600/3066
#> 2026-10-12 00:43:41.495546 - Posterior probability computed for edge n°2700/3066
#> 2026-10-12 00:43:43.335396 - Posterior probability computed for edge n°2800/3066
#> 2026-10-12 00:43:45.456296 - Posterior probability computed for edge n°2900/3066
#> 2026-10-12 00:43:47.292116 - Posterior probability computed for edge n°3000/3066
#> 2026-10-12 00:43:48.431767 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-12 00:43:51.713372 - Posterior probability computed for edge n°100/3066
#> 2026-10-12 00:43:52.626447 - Posterior probability computed for edge n°200/3066
#> 2026-10-12 00:43:53.7208 - Posterior probability computed for edge n°300/3066
#> 2026-10-12 00:43:56.266173 - Posterior probability computed for edge n°400/3066
#> 2026-10-12 00:43:56.987231 - Posterior probability computed for edge n°500/3066
#> 2026-10-12 00:43:57.282543 - Posterior probability computed for edge n°600/3066
#> 2026-10-12 00:43:57.654402 - Posterior probability computed for edge n°700/3066
#> 2026-10-12 00:43:57.967994 - Posterior probability computed for edge n°800/3066
#> 2026-10-12 00:43:58.324963 - Posterior probability computed for edge n°900/3066
#> 2026-10-12 00:43:58.702387 - Posterior probability computed for edge n°1000/3066
#> 2026-10-12 00:43:59.063392 - Posterior probability computed for edge n°1100/3066
#> 2026-10-12 00:43:59.715161 - Posterior probability computed for edge n°1200/3066
#> 2026-10-12 00:44:00.360185 - Posterior probability computed for edge n°1300/3066
#> 2026-10-12 00:44:01.3667 - Posterior probability computed for edge n°1400/3066
#> 2026-10-12 00:44:02.028097 - Posterior probability computed for edge n°1500/3066
#> 2026-10-12 00:44:03.182835 - Posterior probability computed for edge n°1600/3066
#> 2026-10-12 00:44:04.288789 - Posterior probability computed for edge n°1700/3066
#> 2026-10-12 00:44:05.086019 - Posterior probability computed for edge n°1800/3066
#> 2026-10-12 00:44:05.72489 - Posterior probability computed for edge n°1900/3066
#> 2026-10-12 00:44:06.264585 - Posterior probability computed for edge n°2000/3066
#> 2026-10-12 00:44:06.970284 - Posterior probability computed for edge n°2100/3066
#> 2026-10-12 00:44:07.801618 - Posterior probability computed for edge n°2200/3066
#> 2026-10-12 00:44:08.82874 - Posterior probability computed for edge n°2300/3066
#> 2026-10-12 00:44:09.543076 - Posterior probability computed for edge n°2400/3066
#> 2026-10-12 00:44:10.21215 - Posterior probability computed for edge n°2500/3066
#> 2026-10-12 00:44:11.041509 - Posterior probability computed for edge n°2600/3066
#> 2026-10-12 00:44:12.202896 - Posterior probability computed for edge n°2700/3066
#> 2026-10-12 00:44:14.053422 - Posterior probability computed for edge n°2800/3066
#> 2026-10-12 00:44:16.17244 - Posterior probability computed for edge n°2900/3066
#> 2026-10-12 00:44:18.007333 - Posterior probability computed for edge n°3000/3066
#> 2026-10-12 00:44:19.13134 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-10-12 00:44:22.421332 - Posterior probability computed for edge n°100/3066
#> 2026-10-12 00:44:23.340159 - Posterior probability computed for edge n°200/3066
#> 2026-10-12 00:44:24.569028 - Posterior probability computed for edge n°300/3066
#> 2026-10-12 00:44:26.596233 - Posterior probability computed for edge n°400/3066
#> 2026-10-12 00:44:27.374054 - Posterior probability computed for edge n°500/3066
#> 2026-10-12 00:44:27.692693 - Posterior probability computed for edge n°600/3066
#> 2026-10-12 00:44:28.632644 - Posterior probability computed for edge n°700/3066
#> 2026-10-12 00:44:28.948877 - Posterior probability computed for edge n°800/3066
#> 2026-10-12 00:44:29.30815 - Posterior probability computed for edge n°900/3066
#> 2026-10-12 00:44:29.688415 - Posterior probability computed for edge n°1000/3066
#> 2026-10-12 00:44:30.051565 - Posterior probability computed for edge n°1100/3066
#> 2026-10-12 00:44:30.708862 - Posterior probability computed for edge n°1200/3066
#> 2026-10-12 00:44:31.357651 - Posterior probability computed for edge n°1300/3066
#> 2026-10-12 00:44:32.369848 - Posterior probability computed for edge n°1400/3066
#> 2026-10-12 00:44:33.026784 - Posterior probability computed for edge n°1500/3066
#> 2026-10-12 00:44:34.171726 - Posterior probability computed for edge n°1600/3066
#> 2026-10-12 00:44:35.284651 - Posterior probability computed for edge n°1700/3066
#> 2026-10-12 00:44:36.086656 - Posterior probability computed for edge n°1800/3066
#> 2026-10-12 00:44:36.723662 - Posterior probability computed for edge n°1900/3066
#> 2026-10-12 00:44:37.26088 - Posterior probability computed for edge n°2000/3066
#> 2026-10-12 00:44:37.970661 - Posterior probability computed for edge n°2100/3066
#> 2026-10-12 00:44:38.802205 - Posterior probability computed for edge n°2200/3066
#> 2026-10-12 00:44:39.828686 - Posterior probability computed for edge n°2300/3066
#> 2026-10-12 00:44:40.540673 - Posterior probability computed for edge n°2400/3066
#> 2026-10-12 00:44:41.212787 - Posterior probability computed for edge n°2500/3066
#> 2026-10-12 00:44:42.037277 - Posterior probability computed for edge n°2600/3066
#> 2026-10-12 00:44:43.199686 - Posterior probability computed for edge n°2700/3066
#> 2026-10-12 00:44:45.039899 - Posterior probability computed for edge n°2800/3066
#> 2026-10-12 00:44:47.160692 - Posterior probability computed for edge n°2900/3066
#> 2026-10-12 00:44:48.999568 - Posterior probability computed for edge n°3000/3066
#> 2026-10-12 00:44:50.137645 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-12 00:44:53.32747 - Posterior probability computed for edge n°100/3066
#> 2026-10-12 00:44:54.242942 - Posterior probability computed for edge n°200/3066
#> 2026-10-12 00:44:55.351982 - Posterior probability computed for edge n°300/3066
#> 2026-10-12 00:44:57.469531 - Posterior probability computed for edge n°400/3066
#> 2026-10-12 00:44:58.249753 - Posterior probability computed for edge n°500/3066
#> 2026-10-12 00:44:58.563828 - Posterior probability computed for edge n°600/3066
#> 2026-10-12 00:44:58.970419 - Posterior probability computed for edge n°700/3066
#> 2026-10-12 00:44:59.311348 - Posterior probability computed for edge n°800/3066
#> 2026-10-12 00:44:59.69559 - Posterior probability computed for edge n°900/3066
#> 2026-10-12 00:45:00.102156 - Posterior probability computed for edge n°1000/3066
#> 2026-10-12 00:45:01.019552 - Posterior probability computed for edge n°1100/3066
#> 2026-10-12 00:45:01.673885 - Posterior probability computed for edge n°1200/3066
#> 2026-10-12 00:45:02.317814 - Posterior probability computed for edge n°1300/3066
#> 2026-10-12 00:45:03.327207 - Posterior probability computed for edge n°1400/3066
#> 2026-10-12 00:45:03.98936 - Posterior probability computed for edge n°1500/3066
#> 2026-10-12 00:45:05.134108 - Posterior probability computed for edge n°1600/3066
#> 2026-10-12 00:45:06.228703 - Posterior probability computed for edge n°1700/3066
#> 2026-10-12 00:45:07.037842 - Posterior probability computed for edge n°1800/3066
#> 2026-10-12 00:45:07.680609 - Posterior probability computed for edge n°1900/3066
#> 2026-10-12 00:45:08.218894 - Posterior probability computed for edge n°2000/3066
#> 2026-10-12 00:45:08.924617 - Posterior probability computed for edge n°2100/3066
#> 2026-10-12 00:45:09.753902 - Posterior probability computed for edge n°2200/3066
#> 2026-10-12 00:45:10.77983 - Posterior probability computed for edge n°2300/3066
#> 2026-10-12 00:45:11.491796 - Posterior probability computed for edge n°2400/3066
#> 2026-10-12 00:45:12.167563 - Posterior probability computed for edge n°2500/3066
#> 2026-10-12 00:45:12.98951 - Posterior probability computed for edge n°2600/3066
#> 2026-10-12 00:45:14.151909 - Posterior probability computed for edge n°2700/3066
#> 2026-10-12 00:45:15.995296 - Posterior probability computed for edge n°2800/3066
#> 2026-10-12 00:45:18.121488 - Posterior probability computed for edge n°2900/3066
#> 2026-10-12 00:45:19.960453 - Posterior probability computed for edge n°3000/3066
#> 2026-10-12 00:45:21.095299 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-12 00:45:24.378714 - Posterior probability computed for edge n°100/3066
#> 2026-10-12 00:45:25.301112 - Posterior probability computed for edge n°200/3066
#> 2026-10-12 00:45:26.407097 - Posterior probability computed for edge n°300/3066
#> 2026-10-12 00:45:28.452636 - Posterior probability computed for edge n°400/3066
#> 2026-10-12 00:45:29.228476 - Posterior probability computed for edge n°500/3066
#> 2026-10-12 00:45:29.547157 - Posterior probability computed for edge n°600/3066
#> 2026-10-12 00:45:29.948659 - Posterior probability computed for edge n°700/3066
#> 2026-10-12 00:45:30.291759 - Posterior probability computed for edge n°800/3066
#> 2026-10-12 00:45:30.679922 - Posterior probability computed for edge n°900/3066
#> 2026-10-12 00:45:31.081053 - Posterior probability computed for edge n°1000/3066
#> 2026-10-12 00:45:31.477334 - Posterior probability computed for edge n°1100/3066
#> 2026-10-12 00:45:32.701119 - Posterior probability computed for edge n°1200/3066
#> 2026-10-12 00:45:33.344831 - Posterior probability computed for edge n°1300/3066
#> 2026-10-12 00:45:34.349941 - Posterior probability computed for edge n°1400/3066
#> 2026-10-12 00:45:35.009784 - Posterior probability computed for edge n°1500/3066
#> 2026-10-12 00:45:36.152906 - Posterior probability computed for edge n°1600/3066
#> 2026-10-12 00:45:37.244297 - Posterior probability computed for edge n°1700/3066
#> 2026-10-12 00:45:38.048213 - Posterior probability computed for edge n°1800/3066
#> 2026-10-12 00:45:38.685437 - Posterior probability computed for edge n°1900/3066
#> 2026-10-12 00:45:39.221181 - Posterior probability computed for edge n°2000/3066
#> 2026-10-12 00:45:39.924369 - Posterior probability computed for edge n°2100/3066
#> 2026-10-12 00:45:40.75483 - Posterior probability computed for edge n°2200/3066
#> 2026-10-12 00:45:41.780333 - Posterior probability computed for edge n°2300/3066
#> 2026-10-12 00:45:42.48344 - Posterior probability computed for edge n°2400/3066
#> 2026-10-12 00:45:43.150412 - Posterior probability computed for edge n°2500/3066
#> 2026-10-12 00:45:43.969641 - Posterior probability computed for edge n°2600/3066
#> 2026-10-12 00:45:45.123827 - Posterior probability computed for edge n°2700/3066
#> 2026-10-12 00:45:46.956578 - Posterior probability computed for edge n°2800/3066
#> 2026-10-12 00:45:49.069443 - Posterior probability computed for edge n°2900/3066
#> 2026-10-12 00:45:50.894597 - Posterior probability computed for edge n°3000/3066
#> 2026-10-12 00:45:52.009096 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
