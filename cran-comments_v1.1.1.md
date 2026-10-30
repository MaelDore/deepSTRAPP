# NA

## R CMD check results

0 errors \| 0 warning \| 2 notes

- This is a minor update of an existing R package from v1.1.0 to v1.1.1
  fixing a minor bug.

- Notes are similar to the previous release (v1.1.0)

- Note 1:

> Suggests or Enhances not in mainstream repositories: BioGeoBEARS,
> contsimmap Availability using Additional_repositories specification:
> BioGeoBEARS yes <https://maeldore.github.io/drat> contsimmap yes
> <https://maeldore.github.io/drat>

> Packages suggested but not available for checking: ‘BioGeoBEARS’,
> ‘contsimmap’

Core features of deepSTRAPP relies on the BioGeoBEARS and contsimmap
packages, which are not hosted on CRAN but are well-established and
actively maintained R packages widely used in macroevolutionary
research. These packages are indicated as Suggests and conditions are
implemented to check for its presence when running functions that rely
on them, ensuring the package is functional even if they are not
installed.

These suggested dependencies were already present in deepSTRAPP v1.1.0.

- Note 2:

> checking installed package size … NOTE installed size is 8.7Mb
> sub-directories of 1Mb or more: data 1.4Mb doc 5.7Mb

This package implements macroevolutionary modeling on large phylogenies,
which inherently generates sizable data objects. To provide meaningful
and reproducible examples, the package includes a few representative
datasets that reflect typical outputs of its workflow. These datasets,
which have been reduced to the bare minimum, account for the increased
overall package size and meet the 10Mb size limit for CRAN. Pre-rendered
visual outputs for vignettes, which replace even more massive datasets,
are also included and contribute to the overall size. As the package is
designed for downstream analyses of time-calibrated phylogenies and is
not intended as a dependency for other packages, its large size should
not pose practical issues.

This note was already present in deepSTRAPP v1.1.0.

- Note 3:

> Package unavailable to check Rd xrefs: ‘contsimmap’

The R package contsimmap is not on CRAN (yet); it is an optional
Suggests dependency, available from the repository declared in
Additional_repositories (<https://maeldore.github.io/drat>).

This suggested dependency was already present in deepSTRAPP v1.1.0.

## Replies to comments from the CRAN team

No comments yet.
