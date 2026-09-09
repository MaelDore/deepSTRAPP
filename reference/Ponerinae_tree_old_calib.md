# Dataset providing the extensive time-calibrated phylogeny of extant ponerine ants using an old calibration for illustrative purposes

A `phylo` object describing the time-calibrated phylogeny of the 1534
extant ponerine ants (Ponerinae subfamily). This is NOT a properly
time-calibrated phylogeny. It uses an ill-designed old calibration for
illustrative purposes. For a proper time-calibrated phylogeny of
ponerine ants, see the `Ponerinae_tree` object in deepSTRAPP.

Source: Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P.,
Fisher, B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025).
Evolutionary history of ponerine ants highlights how the timing of
dispersal events shapes modern biodiversity. Nature Communications, 16,
8297.
[doi:10.1038/s41467-025-63709-3](https://doi.org/10.1038/s41467-025-63709-3)

## Usage

``` r
data(Ponerinae_tree_old_calib)
```

## Format

A `phylo` object with 4 elements.

## Details

A time-calibrated phylogeny as a `phylo` object with 4 elements.

- `$edge` Integer matrix. Defines the tree topology by providing
  rootward and tipward node ID of each edge.

- `$edge.length` Numeric vector. Length of edges/branches.

- `$Nnode` Integer. Number of internal nodes.

- `$tip.label` Character vector. Labels of all tips.

## References

Doré, M., Borowiec, M. L., Branstetter, M. G., Camacho, G. P., Fisher,
B. L., Longino, J. T., Ward, P. S., Blaimer, B. B. (2025). Evolutionary
history of ponerine ants highlights how the timing of dispersal events
shapes modern biodiversity. Nature Communications, 16, 8297.
[doi:10.1038/s41467-025-63709-3](https://doi.org/10.1038/s41467-025-63709-3)
