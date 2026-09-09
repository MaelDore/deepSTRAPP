# Data summarizing the evolution of geographic ranges in Ponerinae ants using an old ill-calibrated phylogeny for illustrative purposes

A list containing geographic ranges data of Ponerinae mapped on the old
ill-calibrated phylogeny, modeled with R package `BioGeoBEARS`. Ranges
are labeled between "Old World" (O) and New World (N). This object was
obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).
The phylogeny used is also NOT a properly time-calibrated phylogeny. It
uses an ill-designed old calibration for illustrative purposes.

Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P.,
Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
Evolutionary history of ponerine ants highlights how the timing of
dispersal events shapes modern biodiversity. Nature Communications, 16,
8297.
[doi:10.1038/s41467-025-63709-3](https://doi.org/10.1038/s41467-025-63709-3)

## Usage

``` r
data(Ponerinae_biogeo_data_old_calib)
```

## Format

A list with 6 elements.

## Details

A list of six elements containing information on the biogeographic
history of Ponerinae ants. This object was obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).

- `$densityMaps` List of objects of class `"densityMap"` that contains a
  phylogenetic tree and associated mapping of probability to harbor a
  given range along branches. The list contains only a `"densityMap"`
  per unique area because `split_multi_area_ranges` was set to TRUE. In
  this case, unique areas are "N" (= "New World") and "O" (= "Old
  World")

- `$densityMaps_all_ranges` List of objects of class `"densityMap"` that
  contains a phylogenetic tree and associated mapping of probability to
  harbor a given range along branches. The list contains one
  `"densityMap"` per range found along branches during the simulated
  biogeographic histories. Here those ranges are "N" (= "New World"),
  "O" (= "Old World"), and "NO" for multi-area ranges encompassing both
  regions.

- `$trait_data_type` Character string. Record the type of trait data.
  Here: "biogeographic".

- `$ace` Numeric matrix. Record the posterior probabilities of ancestral
  ranges estimated at internal nodes. Only unique areas (i.e., "N" and
  "O") are considered among the ranges. Multi-area ranges (i.e., "NO")
  have been split among unique ranges. Rows are internal nodes. Columns
  are ranges. Values are posterior probabilities of each range per node.

- `$ace_all_ranges` Numeric matrix. Record the posterior probabilities
  of ancestral ranges estimated at internal nodes. All ranges observed
  along branches during the simulated biogeographic histories are
  present (i.e., "N", "O", and "NO"). Rows are internal nodes. Columns
  are ranges. Values are posterior probabilities of each range per node.

- `$model_selection_df` Data.frame that summarizes model comparisons
  used to select the best fitting model.

## References

Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher,
B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025). Evolutionary
history of ponerine ants highlights how the timing of dispersal events
shapes modern biodiversity. Nature Communications, 16, 8297.
[doi:10.1038/s41467-025-63709-3](https://doi.org/10.1038/s41467-025-63709-3)

## See also

[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md)
