# Data summarizing the evolution of geographic ranges in eels

A list containing (fake) geographic ranges data of eels mapped on the
phylogeny, modeled with R package `BioGeoBEARS`. This object was
obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).
Initial data based on feeding habits was altered to be transformed into
ranges "A" and "B", and then adding arbitrarily multi-area "AB" ranges.
This is NOT real biogeographic data. Please refer to the initial article
for real data.

Original data source: Collar, D. C., P. C. Wainwright, M. E. Alfaro, L.
J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur
skull evolution in eels. Nature Communications, 5, 5505.
[doi:10.1038/ncomms6505](https://doi.org/10.1038/ncomms6505)

## Usage

``` r
data(eel_biogeo_data)
```

## Format

A list with 9 elements.

## Details

A list of 9 elements containing information on the evolution of
geographic ranges in eels. This object was obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).

- `$densityMaps` List of objects of class `"densityMap"` that contains a
  phylogenetic tree and associated mapping of probability to harbor a
  given range along branches. The list contains only a `"densityMap"`
  per unique area because `split_multi_area_ranges` was set to TRUE.

- `$densityMaps_all_ranges` List of objects of class `"densityMap"` that
  contains a phylogenetic tree and associated mapping of probability to
  harbor a given range along branches. The list contains one
  `"densityMap"` per range found along branches during the simulated
  biogeographic histories.

- `$trait_data_type` Character string. Record the type of trait data.
  Here: "biogeographic".

- `$ace` Numeric matrix. Record the posterior probabilities of ancestral
  ranges estimated at internal nodes. Only unique areas are considered
  among the ranges. Multi-area ranges have been split among unique
  ranges. Rows are internal nodes. Columns are ranges. Values are
  posterior probabilities of each range per node.

- `$ace_all_ranges` Numeric matrix. Record the posterior probabilities
  of ancestral ranges estimated at internal nodes. All ranges observed
  along branches during the simulated biogeographic histories are
  present. Rows are internal nodes. Columns are ranges. Values are
  posterior probabilities of each range per node.

- `$BSM_output` List of two lists that contain summary information on
  cladogenetic (`$RES_caldo_events_tables`) and anagenetic
  (`$RES_ana_events_tables`) events recorded across the 1000 simulations
  of biogeographic histories performed during Biogeographic Stochastic
  Mapping (BSM). Each element of the list is a data.frame recording
  events occurring during one simulation.

- `$simmaps` List of 1000 objects of class `"simmap"`. Each simmap
  object is a phylogeny with one simulated biogeographic history (i.e.,
  transitions in geographic ranges) mapped along branches.

- `$best_model_fit` List that provides the output of the best fitting
  model (Here: DEC+J model).

- `$model_selection_df` Data.frame that summarizes model comparisons
  used to select the best fitting model.

## References

Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S.
Mehta (2014) Biting disrupts integration to spur skull evolution in
eels. Nature Communications, 5, 5505.
[doi:10.1038/ncomms6505](https://doi.org/10.1038/ncomms6505)

## See also

[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md)
