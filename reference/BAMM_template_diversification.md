# Template file for BAMM diversification analyses

Template file for BAMM diversification analyses provided as character
strings.

Source: bamm-2.5.0 References: http://bamm-project.org/;
https://github.com/macroevolution/bamm

## Usage

``` r
data(BAMM_template_diversification)
```

## Format

A vector of 260 character strings.

## Details

A vector of 260 character strings that can be displayed with
`print(BAMM_template_diversification)`. This is the template used to
generate the 'config_file.txt' controlling settings used for a BAMM run.
It provides detailed explanations of the `additional_BAMM_settings` that
can be used in
[`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md)
to customize the BAMM run. It is called internally by
[`prepare_diversification_data()`](https://maeldore.github.io/deepSTRAPP/reference/prepare_diversification_data.md)
to produce the custom 'config_file' used in the subsequent BAMM run.

## References

Authors: Daniel Rabosky (BAMM) & Pascal Title (`{BAMMtools}`). Modified
by Maël Doré for deepSTRAPP.

## See also

BAMM software website: <http://bamm-project.org/>
