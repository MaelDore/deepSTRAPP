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
#> 2026-09-09 05:40:18.425118 - Fit 1 evolutionary model(s): ER.
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
#> 2026-09-09 05:40:23.330756 - Compare model fits.
#> 
#>    model      logL k      AIC     AICc delta_AICc Akaike_weights rank
#> ER    ER -45.66509 1 93.33017 93.33278          0            100    1
#> 2026-09-09 05:40:23.331821 - Run simulations for stochastic mapping.
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
#> 2026-09-09 05:42:07.903723 - Extract ACE as posterior sampling from stochastic mapping.
#> 2026-09-09 05:42:10.704748 - Create densityMaps by summarizing simulations of evolutionary history (simmaps).
#> 
#> 2026-09-09 05:42:14.050789 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 05:42:14.972496 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 05:42:16.066481 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 05:42:18.071621 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 05:42:18.815122 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 05:42:19.678666 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 05:42:20.017638 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 05:42:20.283618 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 05:42:20.582539 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 05:42:20.895201 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 05:42:21.19938 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 05:42:21.756043 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 05:42:22.300607 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 05:42:23.153013 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 05:42:23.701548 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 05:42:24.674267 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 05:42:25.59773 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 05:42:26.287382 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 05:42:26.817214 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 05:42:27.262726 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 05:42:27.856626 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 05:42:28.559657 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 05:42:29.42643 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 05:42:30.023387 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 05:42:30.582794 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 05:42:31.268765 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 05:42:32.272791 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 05:42:33.819188 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 05:42:35.609349 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 05:42:37.153426 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 05:42:38.112243 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-09 05:42:40.905881 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 05:42:41.668922 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 05:42:42.580366 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 05:42:44.402906 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 05:42:45.026481 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 05:42:45.285061 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 05:42:45.626197 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 05:42:45.904511 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 05:42:46.215595 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 05:42:46.537424 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 05:42:46.866959 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 05:42:47.446462 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 05:42:48.006982 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 05:42:48.912981 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 05:42:49.485392 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 05:42:50.523632 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 05:42:51.49666 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 05:42:52.204089 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 05:42:52.758446 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 05:42:53.232298 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 05:42:54.501019 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 05:42:55.185721 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 05:42:56.015049 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 05:42:56.585621 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 05:42:57.123822 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 05:42:57.779577 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 05:42:58.739748 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 05:43:00.231222 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 05:43:01.967084 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 05:43:03.445898 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 05:43:04.3837 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-09 05:43:07.237971 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 05:43:07.991667 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 05:43:08.903048 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 05:43:10.581298 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 05:43:11.264691 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 05:43:11.526641 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 05:43:11.849726 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 05:43:12.122633 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 05:43:12.429395 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 05:43:12.773388 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 05:43:13.087277 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 05:43:13.657911 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 05:43:14.242139 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 05:43:15.137549 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 05:43:15.714879 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 05:43:16.740449 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 05:43:17.71989 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 05:43:18.429144 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 05:43:19.001555 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 05:43:19.467991 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 05:43:20.085416 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 05:43:20.811777 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 05:43:21.724119 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 05:43:22.350062 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 05:43:22.929567 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 05:43:23.655792 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 05:43:24.701937 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 05:43:26.981255 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 05:43:28.761328 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 05:43:30.235724 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 05:43:31.161098 - Posterior probabilities computed for State = terricolous - n°3/3

# Note that densityMaps are already produced by the [deepSTRAPP::prepare_trait_data()] function,
# but for the sake of example, we can convert the simmaps stored in the output into densityMaps,
# using [deepSTRAPP::convert_simmaps_to_densityMaps()].

## Convert simmaps to densityMaps as input format for deepSTRAPP
Ponerinae_densityMaps <- convert_simmaps_to_densityMaps(
   simmaps = Ponerinae_cat_3lvl_data_old_calib$simmaps)
#> 2026-09-09 05:43:34.087371 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 05:43:34.843148 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 05:43:35.749848 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 05:43:37.418806 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 05:43:38.082186 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 05:43:38.336417 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 05:43:38.659703 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 05:43:38.931446 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 05:43:39.258889 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 05:43:39.579535 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 05:43:39.894997 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 05:43:40.479095 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 05:43:41.037549 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 05:43:41.92833 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 05:43:42.499922 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 05:43:43.521653 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 05:43:44.484237 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 05:43:45.184817 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 05:43:45.740027 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 05:43:46.208605 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 05:43:46.895879 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 05:43:47.629697 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 05:43:48.527913 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 05:43:49.143745 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 05:43:49.741477 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 05:43:50.462072 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 05:43:51.479721 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 05:43:53.105807 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 05:43:54.97004 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 05:43:56.594472 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 05:43:58.297424 - Posterior probabilities computed for State = arboreal - n°1/3
#> 2026-09-09 05:44:01.016327 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 05:44:01.781187 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 05:44:02.687356 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 05:44:04.494105 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 05:44:05.110235 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 05:44:05.363417 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 05:44:05.701896 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 05:44:05.971902 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 05:44:06.275334 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 05:44:06.609591 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 05:44:06.923356 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 05:44:07.488712 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 05:44:08.062515 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 05:44:08.942375 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 05:44:09.53601 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 05:44:10.556015 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 05:44:11.512028 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 05:44:12.207613 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 05:44:12.772887 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 05:44:13.238313 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 05:44:13.845084 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 05:44:14.567376 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 05:44:15.466944 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 05:44:16.090613 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 05:44:16.65758 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 05:44:17.373106 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 05:44:18.416137 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 05:44:20.031718 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 05:44:21.882319 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 05:44:23.482071 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 05:44:24.468685 - Posterior probabilities computed for State = subterranean - n°2/3
#> 2026-09-09 05:44:27.475895 - Posterior probability computed for edge n°100/3066
#> 2026-09-09 05:44:29.064807 - Posterior probability computed for edge n°200/3066
#> 2026-09-09 05:44:29.869419 - Posterior probability computed for edge n°300/3066
#> 2026-09-09 05:44:31.363105 - Posterior probability computed for edge n°400/3066
#> 2026-09-09 05:44:31.931719 - Posterior probability computed for edge n°500/3066
#> 2026-09-09 05:44:32.164659 - Posterior probability computed for edge n°600/3066
#> 2026-09-09 05:44:32.456722 - Posterior probability computed for edge n°700/3066
#> 2026-09-09 05:44:32.704952 - Posterior probability computed for edge n°800/3066
#> 2026-09-09 05:44:32.985218 - Posterior probability computed for edge n°900/3066
#> 2026-09-09 05:44:33.277101 - Posterior probability computed for edge n°1000/3066
#> 2026-09-09 05:44:33.563093 - Posterior probability computed for edge n°1100/3066
#> 2026-09-09 05:44:34.079648 - Posterior probability computed for edge n°1200/3066
#> 2026-09-09 05:44:34.587703 - Posterior probability computed for edge n°1300/3066
#> 2026-09-09 05:44:35.382141 - Posterior probability computed for edge n°1400/3066
#> 2026-09-09 05:44:35.899922 - Posterior probability computed for edge n°1500/3066
#> 2026-09-09 05:44:36.818886 - Posterior probability computed for edge n°1600/3066
#> 2026-09-09 05:44:37.678375 - Posterior probability computed for edge n°1700/3066
#> 2026-09-09 05:44:38.299938 - Posterior probability computed for edge n°1800/3066
#> 2026-09-09 05:44:38.801125 - Posterior probability computed for edge n°1900/3066
#> 2026-09-09 05:44:39.227808 - Posterior probability computed for edge n°2000/3066
#> 2026-09-09 05:44:39.79069 - Posterior probability computed for edge n°2100/3066
#> 2026-09-09 05:44:40.44994 - Posterior probability computed for edge n°2200/3066
#> 2026-09-09 05:44:41.254313 - Posterior probability computed for edge n°2300/3066
#> 2026-09-09 05:44:41.809057 - Posterior probability computed for edge n°2400/3066
#> 2026-09-09 05:44:42.332052 - Posterior probability computed for edge n°2500/3066
#> 2026-09-09 05:44:42.971168 - Posterior probability computed for edge n°2600/3066
#> 2026-09-09 05:44:43.891397 - Posterior probability computed for edge n°2700/3066
#> 2026-09-09 05:44:45.349706 - Posterior probability computed for edge n°2800/3066
#> 2026-09-09 05:44:47.01131 - Posterior probability computed for edge n°2900/3066
#> 2026-09-09 05:44:48.465826 - Posterior probability computed for edge n°3000/3066
#> 2026-09-09 05:44:49.349368 - Posterior probabilities computed for State = terricolous - n°3/3

# Plot densityMaps one by one
plot(Ponerinae_densityMaps[[1]]) # densityMap for state n°1 ("arboreal")

plot(Ponerinae_densityMaps[[2]]) # densityMap for state n°1 ("subterranean")

plot(Ponerinae_densityMaps[[3]]) # densityMap for state n°1 ("terricolous")


# Plot overlay of all densityMaps
plot_densityMaps_overlay(densityMaps = Ponerinae_densityMaps)


```
