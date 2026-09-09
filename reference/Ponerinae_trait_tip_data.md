# Dataset providing fake trait data for extant ponerine ants for illustrative purposes

A data.frame of fake trait data covering the 1534 extant ponerine ant
taxa (Ponerinae subfamily). This is NOT real biological/ecological data.
They were designed for illustrative purposes only.

## Usage

``` r
data(Ponerinae_trait_tip_data)
```

## Format

A data.frame with 1534 rows and 4 columns.

## Details

A data.frame of fake trait data covering the 1534 extant ponerine ant
taxa (Ponerinae subfamily).

- `$Taxa` Character string. Names of Ponerinae ant taxa.

- `$fake_cont_tip_data` Numeric. Fake continuous trait data.

- `$fake_cat_2lvl_tip_data` Character vector. Fake categorical size data
  with two levels: "large" and "small".

- `$fake_cat_3lvl_tip_data` Character vector. Fake categorical habitat
  data with three levels: "arboreal", "subterranean", and "terricolous".
