# Extract diversification data from a BAMM_object

Extracts regime IDs and tip rates from a `BAMM_object` that has been
updated to provide diversification data for a specific time in the past
(i.e. the `focal_time`). Use
[`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md)
to obtain a `BAMM_object` updated for a given `focal_time`.

## Usage

``` r
extract_diversification_data_melted_df_for_focal_time(
  BAMM_object,
  verbose = TRUE
)
```

## Arguments

- BAMM_object:

  Object of class `"bammdata"`, typically generated with
  [`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md),
  that contains a phylogenetic tree and associated diversification rate
  mapping across selected posterior samples.

- verbose:

  Logical. Should progression be displayed? A message will be printed
  for every batch of 100 BAMM posterior samples extracted. Default is
  `TRUE`.

## Value

Returns a data.frame with six columns.

- `$focal_time` Integer. The time, in terms of time distance from the
  present, at which the trait data were extracted. Should be equal for
  all rows as a unique BAMM_object updated for a unique `focal_time` is
  being extracted.

- `$BAMM_sample_ID` Character string. ID of the posterior samples from
  which the diversification data are extracted named as 'BAMM_X' where X
  is the numerical ID.

- `$tip_ID` Character string. Tip labels of the branches cut-off at
  `focal_time`.

  - If `keep_tip_labels = TRUE` was used in
    [`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md),
    cut-off branches with a single descendant tip retain their initial
    `tip.label`.

  - If `keep_tip_labels = FALSE` was used in
    [`update_rates_and_regimes_for_focal_time()`](https://maeldore.github.io/deepSTRAPP/reference/update_rates_and_regimes_for_focal_time.md),
    all cut-off branches are labeled using their tipward node ID.

- `$regime_ID` Integer. The regime ID on tips at `focal_time`. IDs are
  integer. The root process equals "1", then they are incremented from
  oldest to youngest. Regime IDs are independent across posterior
  samples.

- `$rate_type` Character string. Type of rates: "lambda" for speciation
  rates, "mu" for extinction rates, and "net_diversification" for net
  diversification rates (lambda - mu).

- `$rates` Numeric. Rates in \\number of events / branch / evolutionary
  time\\.

## Author

Maël Doré

## Examples

``` r

if (deepSTRAPP::is_dev_version())
{
 # Load the BAMM_object summarizing 1000 posterior samples of BAMM updated for a focal_time of 10 My
 data(Ponerinae_BAMM_object_10My)
 ## This dataset is only available in development versions installed from GitHub.
 # It is not available in CRAN versions.
 # Use remotes::install_github(repo = "MaelDore/deepSTRAPP") to get the latest development version.

  # (May take several minutes to run)
 # Extract diversification data
 diversification_data_df <- extract_diversification_data_melted_df_for_focal_time(
     BAMM_object = Ponerinae_BAMM_object_10My,
     verbose = TRUE)
 # Print output
 head(diversification_data_df) 
}
```
