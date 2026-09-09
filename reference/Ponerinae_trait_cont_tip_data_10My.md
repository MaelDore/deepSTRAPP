# Data summarizing the evolution of a fake continuous trait in Ponerinae ants extracted for 10 Mya

A list containing estimated values for a fake continuous trait mapped on
the Ponerinae ant phylogeny, modeled with
[geiger::fitContinuous](https://rdrr.io/pkg/geiger/man/fitContinuous.html).
This object was obtained with
[`extract_most_likely_trait_values_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/extract_most_likely_trait_values_for_focal_time.md)
applied on trait evolution data obtained with
[`prepare_trait_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_trait_data.md).
The phylogeny used is NOT a properly time-calibrated phylogeny. It uses
an ill-designed old calibration for illustrative purposes.

## Usage

``` r
data(Ponerinae_trait_cont_tip_data_10My)
```

## Format

A list with 4 elements.

## Details

A list of four elements containing information on the evolution of a
continuous trait in Ponerinae ants extracted for 10 Mya.

- `$trait_data` Named numeric vector. Names are the taxa or internal
  tipward node ID associated with the values. Values are the continuous
  trait data estimated along branches for 10 Mya.

- `$focal_time` Numeric. Time in the past at which the trait data were
  extracted.

- `$trait_data_type` Character string. Record the type of trait data.
  Here: "continuous".

- `$contMap` A phylogenetic tree and associated mapping of estimated
  trait values. It was updated such that the tips correspond to lineages
  found 10 Mya (i.e., at the focal time in the past).
