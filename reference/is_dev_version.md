# Check whether deepSTRAPP is a development version

Detect whether the current deepSTRAPP installation corresponds to a
development version (e.g., version sourced from GitHub), as opposed to a
CRAN release.

The development versions include deepSTRAPP outputs as additional
internal datasets to help produce examples and vignettes outputs. These
additional datasets are removed from the CRAN releases because their
size is not compatible with CRAN policies.

This function is only used to check if an example must be run or a
vignette chunk evaluated to produce output from data.

## Usage

``` r
is_dev_version(pkg = "deepSTRAPP")
```

## Arguments

- pkg:

  Character string. Name of the R package for which the version type is
  inspected. Default is "deepSTRAPP".

## Value

Logical. TRUE if running a development version.

## Author

Maël Doré

## Examples

``` r
# Check the current deepSTRAPP installation is a development version
is_dev_version(pkg = "deepSTRAPP")
#> [1] FALSE
```
