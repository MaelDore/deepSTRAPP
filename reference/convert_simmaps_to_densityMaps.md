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
#> 2026-10-28 08:35:45.648207 - Fit 1 evolutionary model(s): ER.
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
#> 2026-10-28 08:35:49.25209 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-10-28 08:35:49.253884 - Run simulations for stochastic mapping.
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
#> 2026-10-28 08:38:05.391557 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-10-28 08:38:08.70325 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-10-28 08:38:12.429233 - Posterior probability computed for edge n°100/3066
#> 2026-10-28 08:38:13.324422 - Posterior probability computed for edge n°200/3066
#> 2026-10-28 08:38:14.404368 - Posterior probability computed for edge n°300/3066
#> 2026-10-28 08:38:16.377315 - Posterior probability computed for edge n°400/3066
#> 2026-10-28 08:38:17.166809 - Posterior probability computed for edge n°500/3066
#> 2026-10-28 08:38:17.479951 - Posterior probability computed for edge n°600/3066
#> 2026-10-28 08:38:17.872099 - Posterior probability computed for edge n°700/3066
#> 2026-10-28 08:38:18.202317 - Posterior probability computed for edge n°800/3066
#> 2026-10-28 08:38:18.577408 - Posterior probability computed for edge n°900/3066
#> 2026-10-28 08:38:18.978457 - Posterior probability computed for edge n°1000/3066
#> 2026-10-28 08:38:19.357042 - Posterior probability computed for edge n°1100/3066
#> 2026-10-28 08:38:20.041963 - Posterior probability computed for edge n°1200/3066
#> 2026-10-28 08:38:20.725504 - Posterior probability computed for edge n°1300/3066
#> 2026-10-28 08:38:21.775517 - Posterior probability computed for edge n°1400/3066
#> 2026-10-28 08:38:22.471596 - Posterior probability computed for edge n°1500/3066
#> 2026-10-28 08:38:23.675092 - Posterior probability computed for edge n°1600/3066
#> 2026-10-28 08:38:24.822622 - Posterior probability computed for edge n°1700/3066
#> 2026-10-28 08:38:25.659326 - Posterior probability computed for edge n°1800/3066
#> 2026-10-28 08:38:26.32107 - Posterior probability computed for edge n°1900/3066
#> 2026-10-28 08:38:26.880263 - Posterior probability computed for edge n°2000/3066
#> 2026-10-28 08:38:27.622871 - Posterior probability computed for edge n°2100/3066
#> 2026-10-28 08:38:28.488072 - Posterior probability computed for edge n°2200/3066
#> 2026-10-28 08:38:29.563137 - Posterior probability computed for edge n°2300/3066
#> 2026-10-28 08:38:30.308347 - Posterior probability computed for edge n°2400/3066
#> 2026-10-28 08:38:30.999673 - Posterior probability computed for edge n°2500/3066
#> 2026-10-28 08:38:31.864165 - Posterior probability computed for edge n°2600/3066
#> 2026-10-28 08:38:33.086872 - Posterior probability computed for edge n°2700/3066
#> 2026-10-28 08:38:35.008251 - Posterior probability computed for edge n°2800/3066
#> 2026-10-28 08:38:37.216846 - Posterior probability computed for edge n°2900/3066
#> 2026-10-28 08:38:39.12844 - Posterior probability computed for edge n°3000/3066
#> 2026-10-28 08:38:40.297255 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-28 08:38:43.541917 - Posterior probability computed for edge n°100/3066
#> 2026-10-28 08:38:44.466912 - Posterior probability computed for edge n°200/3066
#> 2026-10-28 08:38:46.153787 - Posterior probability computed for edge n°300/3066
#> 2026-10-28 08:38:48.064994 - Posterior probability computed for edge n°400/3066
#> 2026-10-28 08:38:48.802045 - Posterior probability computed for edge n°500/3066
#> 2026-10-28 08:38:49.099429 - Posterior probability computed for edge n°600/3066
#> 2026-10-28 08:38:49.475862 - Posterior probability computed for edge n°700/3066
#> 2026-10-28 08:38:49.793268 - Posterior probability computed for edge n°800/3066
#> 2026-10-28 08:38:50.155307 - Posterior probability computed for edge n°900/3066
#> 2026-10-28 08:38:50.537829 - Posterior probability computed for edge n°1000/3066
#> 2026-10-28 08:38:50.902065 - Posterior probability computed for edge n°1100/3066
#> 2026-10-28 08:38:51.562866 - Posterior probability computed for edge n°1200/3066
#> 2026-10-28 08:38:52.225448 - Posterior probability computed for edge n°1300/3066
#> 2026-10-28 08:38:53.261814 - Posterior probability computed for edge n°1400/3066
#> 2026-10-28 08:38:53.939139 - Posterior probability computed for edge n°1500/3066
#> 2026-10-28 08:38:55.113518 - Posterior probability computed for edge n°1600/3066
#> 2026-10-28 08:38:56.230901 - Posterior probability computed for edge n°1700/3066
#> 2026-10-28 08:38:57.045552 - Posterior probability computed for edge n°1800/3066
#> 2026-10-28 08:38:57.692586 - Posterior probability computed for edge n°1900/3066
#> 2026-10-28 08:38:58.239589 - Posterior probability computed for edge n°2000/3066
#> 2026-10-28 08:38:58.963436 - Posterior probability computed for edge n°2100/3066
#> 2026-10-28 08:38:59.800284 - Posterior probability computed for edge n°2200/3066
#> 2026-10-28 08:39:00.842524 - Posterior probability computed for edge n°2300/3066
#> 2026-10-28 08:39:01.579216 - Posterior probability computed for edge n°2400/3066
#> 2026-10-28 08:39:02.264231 - Posterior probability computed for edge n°2500/3066
#> 2026-10-28 08:39:03.108972 - Posterior probability computed for edge n°2600/3066
#> 2026-10-28 08:39:04.306081 - Posterior probability computed for edge n°2700/3066
#> 2026-10-28 08:39:06.169269 - Posterior probability computed for edge n°2800/3066
#> 2026-10-28 08:39:08.327246 - Posterior probability computed for edge n°2900/3066
#> 2026-10-28 08:39:10.192421 - Posterior probability computed for edge n°3000/3066
#> 2026-10-28 08:39:11.343071 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-28 08:39:14.552239 - Posterior probability computed for edge n°100/3066
#> 2026-10-28 08:39:15.621304 - Posterior probability computed for edge n°200/3066
#> 2026-10-28 08:39:16.744828 - Posterior probability computed for edge n°300/3066
#> 2026-10-28 08:39:18.796674 - Posterior probability computed for edge n°400/3066
#> 2026-10-28 08:39:20.12809 - Posterior probability computed for edge n°500/3066
#> 2026-10-28 08:39:20.424718 - Posterior probability computed for edge n°600/3066
#> 2026-10-28 08:39:20.800625 - Posterior probability computed for edge n°700/3066
#> 2026-10-28 08:39:21.118908 - Posterior probability computed for edge n°800/3066
#> 2026-10-28 08:39:21.478967 - Posterior probability computed for edge n°900/3066
#> 2026-10-28 08:39:21.85872 - Posterior probability computed for edge n°1000/3066
#> 2026-10-28 08:39:22.222866 - Posterior probability computed for edge n°1100/3066
#> 2026-10-28 08:39:22.89559 - Posterior probability computed for edge n°1200/3066
#> 2026-10-28 08:39:23.555763 - Posterior probability computed for edge n°1300/3066
#> 2026-10-28 08:39:24.575544 - Posterior probability computed for edge n°1400/3066
#> 2026-10-28 08:39:25.242292 - Posterior probability computed for edge n°1500/3066
#> 2026-10-28 08:39:26.415101 - Posterior probability computed for edge n°1600/3066
#> 2026-10-28 08:39:27.531512 - Posterior probability computed for edge n°1700/3066
#> 2026-10-28 08:39:28.343265 - Posterior probability computed for edge n°1800/3066
#> 2026-10-28 08:39:28.988693 - Posterior probability computed for edge n°1900/3066
#> 2026-10-28 08:39:29.535629 - Posterior probability computed for edge n°2000/3066
#> 2026-10-28 08:39:30.254762 - Posterior probability computed for edge n°2100/3066
#> 2026-10-28 08:39:31.091593 - Posterior probability computed for edge n°2200/3066
#> 2026-10-28 08:39:32.136379 - Posterior probability computed for edge n°2300/3066
#> 2026-10-28 08:39:32.859828 - Posterior probability computed for edge n°2400/3066
#> 2026-10-28 08:39:33.536556 - Posterior probability computed for edge n°2500/3066
#> 2026-10-28 08:39:34.366629 - Posterior probability computed for edge n°2600/3066
#> 2026-10-28 08:39:35.546917 - Posterior probability computed for edge n°2700/3066
#> 2026-10-28 08:39:37.405334 - Posterior probability computed for edge n°2800/3066
#> 2026-10-28 08:39:39.547628 - Posterior probability computed for edge n°2900/3066
#> 2026-10-28 08:39:41.397394 - Posterior probability computed for edge n°3000/3066
#> 2026-10-28 08:39:42.530864 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-10-28 08:39:45.869493 - Posterior probability computed for edge n°100/3066
#> 2026-10-28 08:39:46.791 - Posterior probability computed for edge n°200/3066
#> 2026-10-28 08:39:47.896689 - Posterior probability computed for edge n°300/3066
#> 2026-10-28 08:39:50.003016 - Posterior probability computed for edge n°400/3066
#> 2026-10-28 08:39:50.774283 - Posterior probability computed for edge n°500/3066
#> 2026-10-28 08:39:51.082931 - Posterior probability computed for edge n°600/3066
#> 2026-10-28 08:39:51.488431 - Posterior probability computed for edge n°700/3066
#> 2026-10-28 08:39:51.826268 - Posterior probability computed for edge n°800/3066
#> 2026-10-28 08:39:52.208586 - Posterior probability computed for edge n°900/3066
#> 2026-10-28 08:39:52.610769 - Posterior probability computed for edge n°1000/3066
#> 2026-10-28 08:39:53.526023 - Posterior probability computed for edge n°1100/3066
#> 2026-10-28 08:39:54.182003 - Posterior probability computed for edge n°1200/3066
#> 2026-10-28 08:39:54.821682 - Posterior probability computed for edge n°1300/3066
#> 2026-10-28 08:39:55.830121 - Posterior probability computed for edge n°1400/3066
#> 2026-10-28 08:39:56.48521 - Posterior probability computed for edge n°1500/3066
#> 2026-10-28 08:39:57.638572 - Posterior probability computed for edge n°1600/3066
#> 2026-10-28 08:39:58.728054 - Posterior probability computed for edge n°1700/3066
#> 2026-10-28 08:39:59.531617 - Posterior probability computed for edge n°1800/3066
#> 2026-10-28 08:40:00.168208 - Posterior probability computed for edge n°1900/3066
#> 2026-10-28 08:40:00.699769 - Posterior probability computed for edge n°2000/3066
#> 2026-10-28 08:40:01.416081 - Posterior probability computed for edge n°2100/3066
#> 2026-10-28 08:40:02.255508 - Posterior probability computed for edge n°2200/3066
#> 2026-10-28 08:40:03.281925 - Posterior probability computed for edge n°2300/3066
#> 2026-10-28 08:40:03.992489 - Posterior probability computed for edge n°2400/3066
#> 2026-10-28 08:40:04.657991 - Posterior probability computed for edge n°2500/3066
#> 2026-10-28 08:40:05.484048 - Posterior probability computed for edge n°2600/3066
#> 2026-10-28 08:40:06.656976 - Posterior probability computed for edge n°2700/3066
#> 2026-10-28 08:40:08.489111 - Posterior probability computed for edge n°2800/3066
#> 2026-10-28 08:40:10.599813 - Posterior probability computed for edge n°2900/3066
#> 2026-10-28 08:40:12.423677 - Posterior probability computed for edge n°3000/3066
#> 2026-10-28 08:40:13.543545 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-10-28 08:40:16.800534 - Posterior probability computed for edge n°100/3066
#> 2026-10-28 08:40:17.717475 - Posterior probability computed for edge n°200/3066
#> 2026-10-28 08:40:18.819997 - Posterior probability computed for edge n°300/3066
#> 2026-10-28 08:40:21.38196 - Posterior probability computed for edge n°400/3066
#> 2026-10-28 08:40:22.103055 - Posterior probability computed for edge n°500/3066
#> 2026-10-28 08:40:22.396309 - Posterior probability computed for edge n°600/3066
#> 2026-10-28 08:40:22.768559 - Posterior probability computed for edge n°700/3066
#> 2026-10-28 08:40:23.087433 - Posterior probability computed for edge n°800/3066
#> 2026-10-28 08:40:23.444026 - Posterior probability computed for edge n°900/3066
#> 2026-10-28 08:40:23.81437 - Posterior probability computed for edge n°1000/3066
#> 2026-10-28 08:40:24.173602 - Posterior probability computed for edge n°1100/3066
#> 2026-10-28 08:40:24.830014 - Posterior probability computed for edge n°1200/3066
#> 2026-10-28 08:40:25.469572 - Posterior probability computed for edge n°1300/3066
#> 2026-10-28 08:40:26.479878 - Posterior probability computed for edge n°1400/3066
#> 2026-10-28 08:40:27.142674 - Posterior probability computed for edge n°1500/3066
#> 2026-10-28 08:40:28.303346 - Posterior probability computed for edge n°1600/3066
#> 2026-10-28 08:40:29.403709 - Posterior probability computed for edge n°1700/3066
#> 2026-10-28 08:40:30.203231 - Posterior probability computed for edge n°1800/3066
#> 2026-10-28 08:40:30.838535 - Posterior probability computed for edge n°1900/3066
#> 2026-10-28 08:40:31.371686 - Posterior probability computed for edge n°2000/3066
#> 2026-10-28 08:40:32.078303 - Posterior probability computed for edge n°2100/3066
#> 2026-10-28 08:40:32.910471 - Posterior probability computed for edge n°2200/3066
#> 2026-10-28 08:40:33.934642 - Posterior probability computed for edge n°2300/3066
#> 2026-10-28 08:40:34.637875 - Posterior probability computed for edge n°2400/3066
#> 2026-10-28 08:40:35.308169 - Posterior probability computed for edge n°2500/3066
#> 2026-10-28 08:40:36.128996 - Posterior probability computed for edge n°2600/3066
#> 2026-10-28 08:40:37.288602 - Posterior probability computed for edge n°2700/3066
#> 2026-10-28 08:40:39.130783 - Posterior probability computed for edge n°2800/3066
#> 2026-10-28 08:40:41.243159 - Posterior probability computed for edge n°2900/3066
#> 2026-10-28 08:40:43.084082 - Posterior probability computed for edge n°3000/3066
#> 2026-10-28 08:40:44.220518 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-10-28 08:40:47.333805 - Posterior probability computed for edge n°100/3066
#> 2026-10-28 08:40:48.237578 - Posterior probability computed for edge n°200/3066
#> 2026-10-28 08:40:49.454491 - Posterior probability computed for edge n°300/3066
#> 2026-10-28 08:40:51.46461 - Posterior probability computed for edge n°400/3066
#> 2026-10-28 08:40:52.225224 - Posterior probability computed for edge n°500/3066
#> 2026-10-28 08:40:52.544987 - Posterior probability computed for edge n°600/3066
#> 2026-10-28 08:40:52.943036 - Posterior probability computed for edge n°700/3066
#> 2026-10-28 08:40:53.818252 - Posterior probability computed for edge n°800/3066
#> 2026-10-28 08:40:54.170851 - Posterior probability computed for edge n°900/3066
#> 2026-10-28 08:40:54.542951 - Posterior probability computed for edge n°1000/3066
#> 2026-10-28 08:40:54.899647 - Posterior probability computed for edge n°1100/3066
#> 2026-10-28 08:40:55.545612 - Posterior probability computed for edge n°1200/3066
#> 2026-10-28 08:40:56.184908 - Posterior probability computed for edge n°1300/3066
#> 2026-10-28 08:40:57.180222 - Posterior probability computed for edge n°1400/3066
#> 2026-10-28 08:40:57.830719 - Posterior probability computed for edge n°1500/3066
#> 2026-10-28 08:40:58.964425 - Posterior probability computed for edge n°1600/3066
#> 2026-10-28 08:41:00.056456 - Posterior probability computed for edge n°1700/3066
#> 2026-10-28 08:41:00.84942 - Posterior probability computed for edge n°1800/3066
#> 2026-10-28 08:41:01.490835 - Posterior probability computed for edge n°1900/3066
#> 2026-10-28 08:41:02.021203 - Posterior probability computed for edge n°2000/3066
#> 2026-10-28 08:41:02.716419 - Posterior probability computed for edge n°2100/3066
#> 2026-10-28 08:41:03.53505 - Posterior probability computed for edge n°2200/3066
#> 2026-10-28 08:41:04.552063 - Posterior probability computed for edge n°2300/3066
#> 2026-10-28 08:41:05.259919 - Posterior probability computed for edge n°2400/3066
#> 2026-10-28 08:41:05.917438 - Posterior probability computed for edge n°2500/3066
#> 2026-10-28 08:41:06.72893 - Posterior probability computed for edge n°2600/3066
#> 2026-10-28 08:41:07.876385 - Posterior probability computed for edge n°2700/3066
#> 2026-10-28 08:41:09.688311 - Posterior probability computed for edge n°2800/3066
#> 2026-10-28 08:41:11.783467 - Posterior probability computed for edge n°2900/3066
#> 2026-10-28 08:41:13.594406 - Posterior probability computed for edge n°3000/3066
#> 2026-10-28 08:41:14.701655 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
