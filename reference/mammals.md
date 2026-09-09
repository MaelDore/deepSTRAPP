# Phylogeny and body mass data for extant and extinct mammal families/genera from Slater, 2013

A list containing two elements:

- `$mammal.mass` A data.frame with mean and standard error of mammal
  body masses in ln(mass in g).

- `$mammal.phy` A time-calibrated phylogeny of extinct and extant mammal
  families/genera.

Source: Slater, G. J. (2013). Phylogenetic evidence for a shift in the
mode of mammalian body size evolution at the Cretaceous-Palaeogene
boundary. Methods in Ecology and Evolution, 4(8), 734-744.
[doi:10.1111/2041-210X.12084](https://doi.org/10.1111/2041-210X.12084)

## Usage

``` r
data(mammals)
```

## Format

A list with 2 elements.

- `$mammal.mass` A data.frame of 213 observations and two columns.
  Values are numerical.

- `$mammal.phy` A time-calibrated phylogeny of 211 tips.

## Details

Initial dataset from Slater, G. J. (2013). The R object was imported
from R package `motmot` to include in deepSTRAPP.
