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
#> 2026-09-28 01:40:28.330905 - Fit 1 evolutionary model(s): ER.
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
#> 2026-09-28 01:40:33.79655 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-09-28 01:40:33.798517 - Run simulations for stochastic mapping.
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
#> 2026-09-28 01:42:58.448015 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-09-28 01:43:02.214246 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-09-28 01:43:06.180632 - Posterior probability computed for edge n°100/3066
#> 2026-09-28 01:43:07.313563 - Posterior probability computed for edge n°200/3066
#> 2026-09-28 01:43:08.636192 - Posterior probability computed for edge n°300/3066
#> 2026-09-28 01:43:11.076743 - Posterior probability computed for edge n°400/3066
#> 2026-09-28 01:43:12.507233 - Posterior probability computed for edge n°500/3066
#> 2026-09-28 01:43:12.864481 - Posterior probability computed for edge n°600/3066
#> 2026-09-28 01:43:13.276652 - Posterior probability computed for edge n°700/3066
#> 2026-09-28 01:43:13.621381 - Posterior probability computed for edge n°800/3066
#> 2026-09-28 01:43:14.016594 - Posterior probability computed for edge n°900/3066
#> 2026-09-28 01:43:14.423379 - Posterior probability computed for edge n°1000/3066
#> 2026-09-28 01:43:14.822193 - Posterior probability computed for edge n°1100/3066
#> 2026-09-28 01:43:15.56529 - Posterior probability computed for edge n°1200/3066
#> 2026-09-28 01:43:16.277666 - Posterior probability computed for edge n°1300/3066
#> 2026-09-28 01:43:17.380056 - Posterior probability computed for edge n°1400/3066
#> 2026-09-28 01:43:18.09632 - Posterior probability computed for edge n°1500/3066
#> 2026-09-28 01:43:19.336086 - Posterior probability computed for edge n°1600/3066
#> 2026-09-28 01:43:20.500756 - Posterior probability computed for edge n°1700/3066
#> 2026-09-28 01:43:21.362501 - Posterior probability computed for edge n°1800/3066
#> 2026-09-28 01:43:22.026097 - Posterior probability computed for edge n°1900/3066
#> 2026-09-28 01:43:22.602274 - Posterior probability computed for edge n°2000/3066
#> 2026-09-28 01:43:23.354496 - Posterior probability computed for edge n°2100/3066
#> 2026-09-28 01:43:24.248366 - Posterior probability computed for edge n°2200/3066
#> 2026-09-28 01:43:25.341809 - Posterior probability computed for edge n°2300/3066
#> 2026-09-28 01:43:26.09528 - Posterior probability computed for edge n°2400/3066
#> 2026-09-28 01:43:26.809629 - Posterior probability computed for edge n°2500/3066
#> 2026-09-28 01:43:27.680887 - Posterior probability computed for edge n°2600/3066
#> 2026-09-28 01:43:28.928619 - Posterior probability computed for edge n°2700/3066
#> 2026-09-28 01:43:30.891268 - Posterior probability computed for edge n°2800/3066
#> 2026-09-28 01:43:33.142416 - Posterior probability computed for edge n°2900/3066
#> 2026-09-28 01:43:35.090992 - Posterior probability computed for edge n°3000/3066
#> 2026-09-28 01:43:36.283999 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-28 01:43:39.635772 - Posterior probability computed for edge n°100/3066
#> 2026-09-28 01:43:40.596251 - Posterior probability computed for edge n°200/3066
#> 2026-09-28 01:43:41.721395 - Posterior probability computed for edge n°300/3066
#> 2026-09-28 01:43:43.933572 - Posterior probability computed for edge n°400/3066
#> 2026-09-28 01:43:44.721925 - Posterior probability computed for edge n°500/3066
#> 2026-09-28 01:43:45.052556 - Posterior probability computed for edge n°600/3066
#> 2026-09-28 01:43:45.479064 - Posterior probability computed for edge n°700/3066
#> 2026-09-28 01:43:45.82737 - Posterior probability computed for edge n°800/3066
#> 2026-09-28 01:43:46.224133 - Posterior probability computed for edge n°900/3066
#> 2026-09-28 01:43:46.63806 - Posterior probability computed for edge n°1000/3066
#> 2026-09-28 01:43:47.033875 - Posterior probability computed for edge n°1100/3066
#> 2026-09-28 01:43:47.774749 - Posterior probability computed for edge n°1200/3066
#> 2026-09-28 01:43:48.478605 - Posterior probability computed for edge n°1300/3066
#> 2026-09-28 01:43:49.611864 - Posterior probability computed for edge n°1400/3066
#> 2026-09-28 01:43:50.331156 - Posterior probability computed for edge n°1500/3066
#> 2026-09-28 01:43:51.619987 - Posterior probability computed for edge n°1600/3066
#> 2026-09-28 01:43:52.828955 - Posterior probability computed for edge n°1700/3066
#> 2026-09-28 01:43:53.707122 - Posterior probability computed for edge n°1800/3066
#> 2026-09-28 01:43:54.401175 - Posterior probability computed for edge n°1900/3066
#> 2026-09-28 01:43:55.594211 - Posterior probability computed for edge n°2000/3066
#> 2026-09-28 01:43:56.333791 - Posterior probability computed for edge n°2100/3066
#> 2026-09-28 01:43:57.190616 - Posterior probability computed for edge n°2200/3066
#> 2026-09-28 01:43:58.259118 - Posterior probability computed for edge n°2300/3066
#> 2026-09-28 01:43:59.014391 - Posterior probability computed for edge n°2400/3066
#> 2026-09-28 01:43:59.706775 - Posterior probability computed for edge n°2500/3066
#> 2026-09-28 01:44:00.552411 - Posterior probability computed for edge n°2600/3066
#> 2026-09-28 01:44:01.758706 - Posterior probability computed for edge n°2700/3066
#> 2026-09-28 01:44:03.6852 - Posterior probability computed for edge n°2800/3066
#> 2026-09-28 01:44:05.937259 - Posterior probability computed for edge n°2900/3066
#> 2026-09-28 01:44:07.878844 - Posterior probability computed for edge n°3000/3066
#> 2026-09-28 01:44:09.071975 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-28 01:44:12.513708 - Posterior probability computed for edge n°100/3066
#> 2026-09-28 01:44:13.462003 - Posterior probability computed for edge n°200/3066
#> 2026-09-28 01:44:14.619982 - Posterior probability computed for edge n°300/3066
#> 2026-09-28 01:44:16.77134 - Posterior probability computed for edge n°400/3066
#> 2026-09-28 01:44:17.582145 - Posterior probability computed for edge n°500/3066
#> 2026-09-28 01:44:17.909274 - Posterior probability computed for edge n°600/3066
#> 2026-09-28 01:44:18.323384 - Posterior probability computed for edge n°700/3066
#> 2026-09-28 01:44:18.672501 - Posterior probability computed for edge n°800/3066
#> 2026-09-28 01:44:19.066262 - Posterior probability computed for edge n°900/3066
#> 2026-09-28 01:44:19.498611 - Posterior probability computed for edge n°1000/3066
#> 2026-09-28 01:44:19.900604 - Posterior probability computed for edge n°1100/3066
#> 2026-09-28 01:44:20.626495 - Posterior probability computed for edge n°1200/3066
#> 2026-09-28 01:44:21.366598 - Posterior probability computed for edge n°1300/3066
#> 2026-09-28 01:44:22.489086 - Posterior probability computed for edge n°1400/3066
#> 2026-09-28 01:44:23.246728 - Posterior probability computed for edge n°1500/3066
#> 2026-09-28 01:44:24.520467 - Posterior probability computed for edge n°1600/3066
#> 2026-09-28 01:44:25.74046 - Posterior probability computed for edge n°1700/3066
#> 2026-09-28 01:44:26.637864 - Posterior probability computed for edge n°1800/3066
#> 2026-09-28 01:44:27.358994 - Posterior probability computed for edge n°1900/3066
#> 2026-09-28 01:44:27.948501 - Posterior probability computed for edge n°2000/3066
#> 2026-09-28 01:44:28.725261 - Posterior probability computed for edge n°2100/3066
#> 2026-09-28 01:44:29.63577 - Posterior probability computed for edge n°2200/3066
#> 2026-09-28 01:44:30.773366 - Posterior probability computed for edge n°2300/3066
#> 2026-09-28 01:44:31.562227 - Posterior probability computed for edge n°2400/3066
#> 2026-09-28 01:44:32.291981 - Posterior probability computed for edge n°2500/3066
#> 2026-09-28 01:44:33.198155 - Posterior probability computed for edge n°2600/3066
#> 2026-09-28 01:44:34.492429 - Posterior probability computed for edge n°2700/3066
#> 2026-09-28 01:44:37.086702 - Posterior probability computed for edge n°2800/3066
#> 2026-09-28 01:44:39.302542 - Posterior probability computed for edge n°2900/3066
#> 2026-09-28 01:44:41.230862 - Posterior probability computed for edge n°3000/3066
#> 2026-09-28 01:44:42.400728 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-09-28 01:44:45.957977 - Posterior probability computed for edge n°100/3066
#> 2026-09-28 01:44:46.913798 - Posterior probability computed for edge n°200/3066
#> 2026-09-28 01:44:48.058867 - Posterior probability computed for edge n°300/3066
#> 2026-09-28 01:44:50.204713 - Posterior probability computed for edge n°400/3066
#> 2026-09-28 01:44:51.011895 - Posterior probability computed for edge n°500/3066
#> 2026-09-28 01:44:51.337569 - Posterior probability computed for edge n°600/3066
#> 2026-09-28 01:44:51.749041 - Posterior probability computed for edge n°700/3066
#> 2026-09-28 01:44:52.095516 - Posterior probability computed for edge n°800/3066
#> 2026-09-28 01:44:52.504903 - Posterior probability computed for edge n°900/3066
#> 2026-09-28 01:44:52.915961 - Posterior probability computed for edge n°1000/3066
#> 2026-09-28 01:44:53.314625 - Posterior probability computed for edge n°1100/3066
#> 2026-09-28 01:44:54.052 - Posterior probability computed for edge n°1200/3066
#> 2026-09-28 01:44:54.758754 - Posterior probability computed for edge n°1300/3066
#> 2026-09-28 01:44:55.894019 - Posterior probability computed for edge n°1400/3066
#> 2026-09-28 01:44:56.636693 - Posterior probability computed for edge n°1500/3066
#> 2026-09-28 01:44:57.924075 - Posterior probability computed for edge n°1600/3066
#> 2026-09-28 01:44:59.136144 - Posterior probability computed for edge n°1700/3066
#> 2026-09-28 01:45:00.022975 - Posterior probability computed for edge n°1800/3066
#> 2026-09-28 01:45:00.722747 - Posterior probability computed for edge n°1900/3066
#> 2026-09-28 01:45:01.311242 - Posterior probability computed for edge n°2000/3066
#> 2026-09-28 01:45:02.10367 - Posterior probability computed for edge n°2100/3066
#> 2026-09-28 01:45:03.026973 - Posterior probability computed for edge n°2200/3066
#> 2026-09-28 01:45:04.171164 - Posterior probability computed for edge n°2300/3066
#> 2026-09-28 01:45:04.953456 - Posterior probability computed for edge n°2400/3066
#> 2026-09-28 01:45:05.706474 - Posterior probability computed for edge n°2500/3066
#> 2026-09-28 01:45:06.620529 - Posterior probability computed for edge n°2600/3066
#> 2026-09-28 01:45:07.903867 - Posterior probability computed for edge n°2700/3066
#> 2026-09-28 01:45:09.944959 - Posterior probability computed for edge n°2800/3066
#> 2026-09-28 01:45:12.283444 - Posterior probability computed for edge n°2900/3066
#> 2026-09-28 01:45:14.314914 - Posterior probability computed for edge n°3000/3066
#> 2026-09-28 01:45:16.167731 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-28 01:45:19.486275 - Posterior probability computed for edge n°100/3066
#> 2026-09-28 01:45:20.461933 - Posterior probability computed for edge n°200/3066
#> 2026-09-28 01:45:21.603163 - Posterior probability computed for edge n°300/3066
#> 2026-09-28 01:45:23.873548 - Posterior probability computed for edge n°400/3066
#> 2026-09-28 01:45:24.675845 - Posterior probability computed for edge n°500/3066
#> 2026-09-28 01:45:25.005512 - Posterior probability computed for edge n°600/3066
#> 2026-09-28 01:45:25.441259 - Posterior probability computed for edge n°700/3066
#> 2026-09-28 01:45:25.792592 - Posterior probability computed for edge n°800/3066
#> 2026-09-28 01:45:26.187613 - Posterior probability computed for edge n°900/3066
#> 2026-09-28 01:45:26.606289 - Posterior probability computed for edge n°1000/3066
#> 2026-09-28 01:45:27.030674 - Posterior probability computed for edge n°1100/3066
#> 2026-09-28 01:45:27.759607 - Posterior probability computed for edge n°1200/3066
#> 2026-09-28 01:45:28.492425 - Posterior probability computed for edge n°1300/3066
#> 2026-09-28 01:45:29.623584 - Posterior probability computed for edge n°1400/3066
#> 2026-09-28 01:45:30.371119 - Posterior probability computed for edge n°1500/3066
#> 2026-09-28 01:45:31.651575 - Posterior probability computed for edge n°1600/3066
#> 2026-09-28 01:45:32.874618 - Posterior probability computed for edge n°1700/3066
#> 2026-09-28 01:45:33.76569 - Posterior probability computed for edge n°1800/3066
#> 2026-09-28 01:45:34.486043 - Posterior probability computed for edge n°1900/3066
#> 2026-09-28 01:45:35.075671 - Posterior probability computed for edge n°2000/3066
#> 2026-09-28 01:45:35.866988 - Posterior probability computed for edge n°2100/3066
#> 2026-09-28 01:45:36.777654 - Posterior probability computed for edge n°2200/3066
#> 2026-09-28 01:45:37.928434 - Posterior probability computed for edge n°2300/3066
#> 2026-09-28 01:45:38.721917 - Posterior probability computed for edge n°2400/3066
#> 2026-09-28 01:45:39.456641 - Posterior probability computed for edge n°2500/3066
#> 2026-09-28 01:45:40.364513 - Posterior probability computed for edge n°2600/3066
#> 2026-09-28 01:45:41.669083 - Posterior probability computed for edge n°2700/3066
#> 2026-09-28 01:45:43.702009 - Posterior probability computed for edge n°2800/3066
#> 2026-09-28 01:45:46.04171 - Posterior probability computed for edge n°2900/3066
#> 2026-09-28 01:45:48.070987 - Posterior probability computed for edge n°3000/3066
#> 2026-09-28 01:45:49.317195 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-28 01:45:52.850163 - Posterior probability computed for edge n°100/3066
#> 2026-09-28 01:45:54.506571 - Posterior probability computed for edge n°200/3066
#> 2026-09-28 01:45:55.575321 - Posterior probability computed for edge n°300/3066
#> 2026-09-28 01:45:57.499608 - Posterior probability computed for edge n°400/3066
#> 2026-09-28 01:45:58.23593 - Posterior probability computed for edge n°500/3066
#> 2026-09-28 01:45:58.542131 - Posterior probability computed for edge n°600/3066
#> 2026-09-28 01:45:58.927786 - Posterior probability computed for edge n°700/3066
#> 2026-09-28 01:45:59.262664 - Posterior probability computed for edge n°800/3066
#> 2026-09-28 01:45:59.632729 - Posterior probability computed for edge n°900/3066
#> 2026-09-28 01:46:00.016554 - Posterior probability computed for edge n°1000/3066
#> 2026-09-28 01:46:00.389774 - Posterior probability computed for edge n°1100/3066
#> 2026-09-28 01:46:01.067408 - Posterior probability computed for edge n°1200/3066
#> 2026-09-28 01:46:01.73156 - Posterior probability computed for edge n°1300/3066
#> 2026-09-28 01:46:02.762953 - Posterior probability computed for edge n°1400/3066
#> 2026-09-28 01:46:03.43755 - Posterior probability computed for edge n°1500/3066
#> 2026-09-28 01:46:04.605816 - Posterior probability computed for edge n°1600/3066
#> 2026-09-28 01:46:05.740077 - Posterior probability computed for edge n°1700/3066
#> 2026-09-28 01:46:06.551006 - Posterior probability computed for edge n°1800/3066
#> 2026-09-28 01:46:07.20026 - Posterior probability computed for edge n°1900/3066
#> 2026-09-28 01:46:07.75062 - Posterior probability computed for edge n°2000/3066
#> 2026-09-28 01:46:08.476572 - Posterior probability computed for edge n°2100/3066
#> 2026-09-28 01:46:09.305784 - Posterior probability computed for edge n°2200/3066
#> 2026-09-28 01:46:10.344781 - Posterior probability computed for edge n°2300/3066
#> 2026-09-28 01:46:11.076006 - Posterior probability computed for edge n°2400/3066
#> 2026-09-28 01:46:11.752783 - Posterior probability computed for edge n°2500/3066
#> 2026-09-28 01:46:12.583944 - Posterior probability computed for edge n°2600/3066
#> 2026-09-28 01:46:13.771114 - Posterior probability computed for edge n°2700/3066
#> 2026-09-28 01:46:15.634364 - Posterior probability computed for edge n°2800/3066
#> 2026-09-28 01:46:17.79391 - Posterior probability computed for edge n°2900/3066
#> 2026-09-28 01:46:19.669865 - Posterior probability computed for edge n°3000/3066
#> 2026-09-28 01:46:20.801189 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
