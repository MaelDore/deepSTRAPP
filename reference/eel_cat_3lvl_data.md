# Data summarizing the evolution of feeding habits in eels using a 3-level factor as categorical trait

A list containing feeding habits data of eels mapped on the phylogeny,
modeled with
[geiger::fitDiscrete](https://rdrr.io/pkg/geiger/man/fitDiscrete.html).
This object was obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).
Initial data was altered arbitrarily to create three categories, adding
a "kiss" feeding habit to the initial "bite" and "suction" data. This is
NOT real biological data. Please refer to the initial article for real
data.

Original data source: Collar, D. C., P. C. Wainwright, M. E. Alfaro, L.
J. Revell, and R. S. Mehta (2014) Biting disrupts integration to spur
skull evolution in eels. Nature Communications, 5, 5505.
[doi:10.1038/ncomms6505](https://doi.org/10.1038/ncomms6505)

## Usage

``` r
data(eel_cat_3lvl_data)
```

## Format

A list with 6 elements.

## Details

A list of six elements containing information on the evolution of
feeding habits in eels. This object was obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).

- `$densityMaps` List of objects of class `"densityMap"` that contains a
  phylogenetic tree and associated mapping of probability to harbor a
  given state/range along branches. The list contains one `"densityMap"`
  per state/range found in the `tip_data`.

- `$trait_data_type` Character string. Record the type of trait data.
  Here: "categorical".

- `$simmaps` List of 100 stochastic mapping simulations for trait
  evolution. Each element is a `"simmap"` object (see
  [phytools::make.simmap](https://rdrr.io/pkg/phytools/man/make.simmap.html))
  representing a possible evolutionary history that fits states observed
  on tips, inferred ACE at internal nodes, and transition rates as
  estimated from the best fit model.

- `$ace` Numeric matrix. Record the posterior probabilities of ancestral
  states/ranges (characters) estimates (ACE) at internal nodes. Rows are
  internal nodes. Columns are states/ranges. Values are posterior
  probabilities of each state per node.

- `$best_model_fit` List that provides the output of the best fitting
  model (Here: ER model).

- `$model_selection_df` Data.frame that summarizes model comparisons
  used to select the best fitting model.

## References

Collar, D. C., P. C. Wainwright, M. E. Alfaro, L. J. Revell, and R. S.
Mehta (2014) Biting disrupts integration to spur skull evolution in
eels. Nature Communications, 5, 5505.
[doi:10.1038/ncomms6505](https://doi.org/10.1038/ncomms6505)

## See also

[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md)
