## R CMD check results

0 errors | 0 warning | 2 notes

* This is a new release.

* Note 1: 

> Suggests or Enhances not in mainstream repositories:
    BioGeoBEARS
Availability using Additional_repositories specification:
  BioGeoBEARS   yes   https://maeldore.github.io/drat

  A core feature of deepSTRAPP relies on the BioGeoBEARS package, which is not hosted on CRAN 
  but is a well-established and actively maintained R package widely used in macroevolutionary research.
  This package is indicated as Suggests and conditions are implemented to check for its presence when running functions that rely on BioGeoBEARS,
  ensuring the package is functional even if BioGeoBEARS is not installed.
  
> Found the following (possibly) invalid URLs:
  URL: https://support.posit.co/hc/en-us/articles/200486498-Package-Development-Prerequisites
  
  This URL is valid, but the server blocks automated requests.
  Please ignore this warning.
  
* Note 2:

> checking installed package size ... NOTE
    installed size is  8.3Mb
    sub-directories of 1Mb or more:
      data   1.4Mb
      doc    4.5Mb
      
This package implements macroevolutionary modeling on large phylogenies, which inherently generates sizable data objects.
To provide meaningful and reproducible examples, the package includes a few representative datasets that reflect typical outputs of its workflow.
These datasets, which have been reduced to the bare minimum, account for the increased overall package size and meet the 10Mb size limit for CRAN.
Pre-rendered visual outputs for vignettes, which replace even more massive datasets, are also included and contribute to the overall size.
As the package is designed for downstream analyses of time-calibrated phylogenies and is not intended as a dependency for other packages, its large size should not pose practical issues.


