## R CMD check results

0 errors | 0 warning | 3 notes

* This is a new release.

* Note 1: 

Suggests or Enhances not in mainstream repositories:
    BioGeoBEARS
Availability using Additional_repositories specification:
  BioGeoBEARS   yes   https://maeldore.github.io/drat

  A core feature of deepSTRAPP relies on the BioGeoBEARS package, which is not hosted on CRAN 
  but is a well-established and actively maintained R package widely used in macroevolutionary research.
  This package is indicated as Suggests and conditions are implemented to check for its presence when running functions that rely on BioGeoBEARS,
  ensuring the package is functional even if BioGeoBEARS is not installed.
  
Found the following (possibly) invalid URLs:
  URL: https://support.posit.co/hc/en-us/articles/200486498-Package-Development-Prerequisites
  
  This URL is valid, but the server blocks automated requests.
  Please ignore this warning.
  
Found the following URLs which should use \doi (with the DOI name only):

  All listed DOIs are valid. The \doi{} tag was not used because it did not generate valid links in the rendered documentation within RStudio.
  Please ignore this warning.
  
* Note 2:

checking installed package size ... NOTE
  installed size is 21.4Mb
    sub-directories of 1Mb or more:
      data  11.1Mb
      doc    7.8Mb
      help   1.1Mb
      
This package implements macroevolutionary modeling on large phylogenies, which inherently generates sizable data objects.
To provide meaningful and reproducible examples, the package includes a few representative datasets that reflect typical outputs of its workflow.
These datasets, which have been reduced to the bare minimum, account for the increased overall package size. 
Pre-rendered visual outputs for vignettes, which replace even more massive datasets, are also included and contribute to the overall size.
As the package is designed for downstream analyses of time-calibrated phylogenies and is not intended as a dependency for other packages, its large size should not pose practical issues.
    
* Note 3:

checking examples ... [178s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s

As with the package size, the examples involve large datasets and objects.
Most time-consuming steps have been shortened using pre-computed data, which account for the package size,
but some functions (e.g., plotting) inherently require longer runtimes that cannot be further reduced.


