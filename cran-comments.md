## R CMD check results

0 errors | 1 warning | 3 notes

* This is a new release.

* Warning: 

Strong dependencies not in mainstream repositories:  BioGeoBEARS
  
  A core feature of deepSTRAPP relies on the BioGeoBEARS package, which is not hosted on CRAN 
  but is a well-established and actively maintained R package widely used in macroevolutionary research.
  This dependency is essential for the package’s main functionality and should not present practical issues.
  
Found the following (possibly) invalid URLs:
  URL: https://support.posit.co/hc/en-us/articles/200486498-Package-Development-Prerequisites
  
  This URL is valid, but the server blocks automated requests.
  Please ignore this warning.
  
Found the following (possibly) invalid DOIs:
    DOI: TBA
      From: inst/CITATION
      Message: Invalid DOI
      
  This is the DOI of the future research paper. This will be updated as soon as this research is published.

Found the following URLs which should use \doi (with the DOI name only):

  All listed DOIs are valid. The \doi{} tag was not used because it did not generate valid links in the rendered documentation within RStudio.
  Please ignore this warning.
  
* Note 1:

checking installed package size ... NOTE
  installed size is 41.0Mb
      sub-directories of 1Mb or more:
        data  33.6Mb
        doc    4.8Mb
        help   1.2Mb
      
This package implements macroevolutionary modeling on large phylogenies, which inherently generates sizable data objects.
To provide meaningful and reproducible examples, the package includes a few representative datasets that reflect typical outputs of its workflow.
These datasets account for the increased overall package size.
As the package is designed for downstream analyses of time-calibrated phylogenies and is not intended as a dependency for other packages, its large size should not pose practical issues.
    
* Note 2:

checking examples ... [487s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s

As with the package size, the examples involve large datasets and objects.
Most time-consuming steps have been shortened using pre-computed data, which account for the package size,
but some functions (e.g., plotting) inherently require longer runtimes that cannot be further reduced.

* Note 3:

checking for future file timestamps ... NOTE
  unable to verify current time
