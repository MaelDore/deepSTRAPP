# Extract trait data from a trait_data object

Extract trait data from a trait_data_list object as produced within the
deepSTRAPP workflow for a specific time in the past (i.e. the
`focal_time`) and produce a melted dataset recording trait values across
(stochastic maps x) edges for the associated `focal_time`.

## Usage

``` r
extract_trait_data_melted_df_for_focal_time(trait_data_list)
```

## Arguments

- trait_data_list:

  Object summarizing trait data extracted for a given `focal_time`,
  including a `$trait_data` element storing trait data recorded across
  (stochastic maps x) edges for the associated `$focal_time`. Use
  [`extract_all_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_all_trait_values_for_focal_time.md)
  to obtain trait data across multiple stochastic maps. Use
  [`extract_most_likely_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_most_likely_trait_values_for_focal_time.md)
  to obtain only the most likely trait data values.

## Value

Returns a data.frame with five columns.

- `$focal_time` Integer. The time, in terms of time distance from the
  present, at which the trait data were extracted.

- `$Map_ID` Character string. ID of the stochastic map from which the
  trait data are extracted. If using 'rate_only' strategy to account for
  uncertainty in ancestral estimate, this is fixed to "Map_ML". If using
  'paired' or 'full' strategies, this records either true stochastic
  maps "Map_X", or dummy maps "Dummy_map_X", depending on whether you
  provided respectively stochastic maps (`contMaps`/`simmaps`) or simply
  `densityMaps` as inputs.

- `$tip_ID` Character string. Tip labels of the branches cut-off at
  `focal_time`.

  - If `keep_tip_labels = TRUE` was used in
    [`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md),
    cut-off branches with a single descendant tip retain their initial
    `tip.label`.

  - If `keep_tip_labels = FALSE` was used in
    [`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md),
    all cut-off branches are labeled using their tipward node ID.

- `$trait_data_type` Character string. Records the type of traits:
  "continuous", "categorical", or "biogeographic".

- `$trait_value` Numerical or Character string.

  - For "continuous" traits: numerical values.

  - For "categorical" traits: states.

  - For "biogeographic" data: ranges.

## Author

Maël Doré

## Examples

``` r
# ----- Example 1: Extract ML estimates data ----- #

## Load categorical trait data mapped on a phylogeny
data(eel_cat_3lvl_data, package = "deepSTRAPP")

# Explore data
str(eel_cat_3lvl_data, 1)
#> List of 6
#>  $ densityMaps       :List of 3
#>  $ trait_data_type   : chr "categorical"
#>  $ simmaps           :Class "multiPhylo"
#> List of 100
#>  $ ace               : num [1:60, 1:3] 0.26 0.3 0.39 0.43 0.43 0.44 0.43 0 0.47 0.49 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ best_model_fit    :List of 4
#>   ..- attr(*, "class")= chr [1:2] "gfit" "list"
#>  $ model_selection_df:'data.frame':  5 obs. of  6 variables:
eel_cat_3lvl_data$densityMaps # Three density maps: one per state
#> $Density_map_bite
#> Object of class "densityMap" containing:
#> 
#> (1) A phylogenetic tree with 61 tips and 60 internal nodes.
#> 
#> (2) The mapped posterior density of a discrete binary character
#>     with states (Not bite, bite).
#> 
#> 
#> $Density_map_kiss
#> Object of class "densityMap" containing:
#> 
#> (1) A phylogenetic tree with 61 tips and 60 internal nodes.
#> 
#> (2) The mapped posterior density of a discrete binary character
#>     with states (Not kiss, kiss).
#> 
#> 
#> $Density_map_suction
#> Object of class "densityMap" containing:
#> 
#> (1) A phylogenetic tree with 61 tips and 60 internal nodes.
#> 
#> (2) The mapped posterior density of a discrete binary character
#>     with states (Not suction, suction).
#> 
#> 

# Set focal time to 10 Mya
focal_time <- 10

## Extract ML estimates of ancestral states for the given focal_time

# Extract from the densityMaps
eel_cat_3lvl_data_10My <- extract_most_likely_trait_values_for_focal_time(
   densityMaps = eel_cat_3lvl_data$densityMaps,
   trait_data_type = "categorical",
   focal_time = focal_time)
#> WARNING: No ancestral character estimates (ace) for internal nodes have been provided. Using most likely states extracted from the densityMaps instead.
#> WARNING: No tip data have been provided. Using states extracted from the densityMaps instead.

## Format ML states data as a melted df

melted_df <- extract_trait_data_melted_df_for_focal_time(trait_data_list = eel_cat_3lvl_data_10My)
head(melted_df)
#>   focal_time Map_ID               tip_ID trait_data_type trait_value
#> 1         10 Map_ML    Moringua_edwardsi     categorical        bite
#> 2         10 Map_ML Kaupichthys_nuchalis     categorical        bite
#> 3         10 Map_ML                   69     categorical     suction
#> 4         10 Map_ML Venefica_proboscidea     categorical        bite
#> 5         10 Map_ML                   74     categorical     suction
#> 6         10 Map_ML     Anguilla_bicolor     categorical     suction


# ----- Example 2: Extract data across multiple stochastic maps ----- #

## Load biogeographic range data mapped on a phylogeny
data(eel_biogeo_data, package = "deepSTRAPP")

# Explore data
str(eel_biogeo_data, 1)
#> List of 9
#>  $ densityMaps           :List of 2
#>  $ densityMaps_all_ranges:List of 3
#>  $ trait_data_type       : chr "biogeographic"
#>  $ ace                   : num [1:60, 1:2] 0.361 0.524 0.623 0.334 0.47 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ ace_all_ranges        : num [1:60, 1:3] 0.174 0.373 0.53 0.244 0.411 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ BSM_output            :List of 2
#>  $ simmaps               :Class "multiPhylo"
#> List of 100
#>  $ best_model_fit        :List of 13
#>  $ model_selection_df    :'data.frame':  6 obs. of  8 variables:
length(eel_biogeo_data$simmaps) # 100 stochastic maps: one per simulated biogeographic history
#> [1] 100

# Set focal time to 10 Mya
focal_time <- 10

## Extract range data across stochastic maps for the given focal_time

# Extract from the simmaps
eel_biogeo_data_10My <- extract_all_trait_values_for_focal_time(
   simmaps = eel_biogeo_data$simmaps,
   trait_data_type = "biogeographic",
   focal_time = focal_time)
#> WARNING: No tip data have been provided. Using ranges extracted from the simmaps instead.

## Format range data as a melted df

melted_df <- extract_trait_data_melted_df_for_focal_time(trait_data_list = eel_biogeo_data_10My)
head(melted_df)
#>   focal_time Map_ID               tip_ID trait_data_type trait_value
#> 1         10  Map_1    Moringua_edwardsi   biogeographic          AB
#> 2         10  Map_1 Kaupichthys_nuchalis   biogeographic          AB
#> 3         10  Map_1                   69   biogeographic           B
#> 4         10  Map_1 Venefica_proboscidea   biogeographic           A
#> 5         10  Map_1                   74   biogeographic          AB
#> 6         10  Map_1     Anguilla_bicolor   biogeographic           B
```
