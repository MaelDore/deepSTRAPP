# Data summarizing the evolution of fake size data in Ponerinae ants using a 2-level factor as categorical trait

A list containing fake size data of Ponerinae ants mapped on the
phylogeny, modeled with
[geiger::fitDiscrete](https://rdrr.io/pkg/geiger/man/fitDiscrete.html).
This object was obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).
This is NOT real biological/ecological data. They were designed for
illustrative purposes only. The phylogeny used is also NOT a properly
time-calibrated phylogeny. It uses an ill-designed old calibration for
illustrative purposes.

## Usage

``` r
data(Ponerinae_cat_2lvl_data_old_calib)
```

## Format

A list with 5 elements.

## Details

A list of five objects containing information on the evolution of fake
size data in Ponerinae ants. This object was obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).

- `$densityMaps` List of two objects of class `"densityMap"` that
  contains a phylogenetic tree and associated mapping of probability to
  harbor a given state along branches. The list contains one
  `"densityMap"` per state found in the `tip_data` (i.e., "large" and
  "small").

- `$trait_data_type` Character string. Record the type of trait data.
  Here: "categorical".

- `$ace` Numeric matrix. Record the posterior probabilities of ancestral
  states (characters) estimates (ACE) at internal nodes. Rows are
  internal nodes. Columns are states (i.e., "large" and "small"). Values
  are posterior probabilities of each state per node.

- `$best_model_fit` List that provides the output of the best fitting
  model (Here: ARD model).

- `$model_selection_df` Data.frame that summarizes model comparisons
  used to select the best fitting model.

## See also

[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md)
