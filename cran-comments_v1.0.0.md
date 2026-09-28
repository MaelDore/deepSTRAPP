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
  
> Possibly misspelled words in DESCRIPTION:
    STRAPP (9:9)
    phylogenies (7:37, 10:11)
    
These words are not misspelled.
  
> Found the following (possibly) invalid URLs:
  URL: https://support.posit.co/hc/en-us/articles/200486498-Package-Development-Prerequisites
  
  This URL is valid, but the server blocks automated requests.
  Please ignore this warning.
  
* Note 2:

> checking installed package size ... NOTE
    installed size is  7.3Mb
    sub-directories of 1Mb or more:
      data   1.4Mb
      doc    4.4Mb
      
This package implements macroevolutionary modeling on large phylogenies, which inherently generates sizable data objects.
To provide meaningful and reproducible examples, the package includes a few representative datasets that reflect typical outputs of its workflow.
These datasets, which have been reduced to the bare minimum, account for the increased overall package size and meet the 10Mb size limit for CRAN.
Pre-rendered visual outputs for vignettes, which replace even more massive datasets, are also included and contribute to the overall size.
As the package is designed for downstream analyses of time-calibrated phylogenies and is not intended as a dependency for other packages, its large size should not pose practical issues.

## Replies to comments from the CRAN team

* Please add () behind and remove the single quotes around all function
names in the description texts (DESCRIPTION file). e.g: -->
BAMMtools::traitDependentBAMM()

Reply: Change for BAMMtools::traitDependentBAMM()

* If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'https:' and angle brackets for
auto-linking. (If you want to add a title as well please put it in
quotes: "Title")
For more details:
<https://contributor.r-project.org/cran-cookbook/description_issues.html#references>

Reply: The paper accompanying this R package will be submitted as soon as the CRAN submission is completed. 
I will update the DESCRIPTION to include the reference once the paper is published, as follows: 
"For more details see Doré & Blaimer (2026) <doi:10.XXX/XXXXX>". This information will also be updated in the inst/CITATION file.

* We see:
Warning: Unexecutable code in man/select_best_model_from_BioGeoBEARS.Rd:
   (May take:
I think you forgot to comment out a line there.

Reply: I commented the line out in all cases following the donttest tag.

* You have examples for unexported functions. Please either omit these
examples or export these functions.
"Using foo:::f instead of foo::f allows access to unexported objects.
This is generally not recommended, as the semantics of unexported
objects may be changed by the package author in routine maintenance."
Used ::: in documentation:
      man/compute_STRAPP_test_for_focal_time.Rd:
         if (deepSTRAPP:::is_dev_version()) {
Please omit one colon.

Reply: is_dev_version() is now an exported function. All references to the function have been updated to two colons such as deepSTRAPP::is_dev_version()

* \dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user. Does not seem necessary.
Please replace \dontrun with \donttest.
Please unwrap the examples if they are executable in < 5 sec, or replace
dontrun{} with \donttest{}.
For more details:
<https://contributor.r-project.org/cran-cookbook/general_issues.html#structuring-of-examples>
-> prepare_diversification_data.Rd;
update_rates_and_regimes_for_focal_time.Rd

Reply: The tag \dontrun is used in cases when the user (or the CRAN machine) may not have all needed components on their local machine such as the BAMM software
(ex: prepare_diversification_data). Other cases (i.e., time-consuming examples) are now dealt with the \donttest tag. 
In addition, examples using datasets that are only available on development version (due to size requirement on CRAN) are wrapped within the is_dev_version() function
and adequately commented to inform the user.

* Please ensure that your functions do not write by default or in your
examples/vignettes/tests in the user's home filespace (including the
package directory and getwd()). This is not allowed by CRAN policies.
Please omit any default path in writing functions. In your
examples/vignettes/tests you can write to tempdir().
For more details:
<https://contributor.r-project.org/cran-cookbook/code_issues.html#writing-files-and-directories-to-the-home-filespace>
-> R/prepare_diversification_data.R; R/prepare_trait_data.R

Reply: Destination folders have been set to NULL in writing functions (i.e.,  BAMM_output_directory_path in prepare_diversification_data and BioGeoBEARS_directory_path in prepare_trait_data).
Examples now use tempdir() as destination folder and remove files after usage.

* Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited.
If you're not familiar with the function, please check ?on.exit. This
function makes it possible to restore options before exiting a function
even if the function breaks. Therefore it needs to be called immediately
after the option change within a function.
For more details:
<https://contributor.r-project.org/cran-cookbook/code_issues.html#change-of-options-graphical-parameters-and-working-directory>
-> R/plot_BAMM_rates.R; R/plot_densityMaps_overlay.R;
R/plot_histogram_STRAPP_test_for_focal_time.R;
R/plot_histograms_STRAPP_tests_over_time.R;
R/plot_rates_vs_trait_data_for_focal_time.R;
R/plot_traits_vs_rates_on_phylogeny_for_focal_time.R;
R/prepare_diversification_data.R; R/prepare_trait_data.R

Reply: The two suggested lines (oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))) were added in preamble of each function involving plotting to ensure that initial user's options are preserved. 

* Please always make sure to reset to user's options(), working directory
or par() after you changed it in examples and vignettes and demos.
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)
-> inst/doc/cut_phylogenies.R
For more details:
<https://contributor.r-project.org/cran-cookbook/code_issues.html#change-of-options-graphical-parameters-and-working-directory>

Reply: The two suggested lines (old_par <- par(no.readonly = TRUE); par(old_par)) were added to wrap any plotting instance in vignettes.

* Please always add all authors, contributors and copyright holders in the
Authors@R field with the appropriate roles.
 From CRAN policies you agreed to:
"The ownership of copyright and intellectual property rights of all
components of the package must be clear and unambiguous (including from
the authors specification in the DESCRIPTION file). Where code is copied
(or derived) from the work of others (including from R itself), care
must be taken that any copyright/license statements are preserved and
authorship is not misrepresented.
Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles for the
authors of such code. Alternatively, the ‘Author’ field should list
these authors as contributors. Where copyrights are held by an entity
other than the package authors, this should preferably be indicated via
‘cph’ roles in the ‘Authors@R’ field, or using a ‘Copyright’ field (if
necessary referring to an inst/COPYRIGHTS file)."
e.g.: Dan Rabosky, Mike Grundler in "treetraverse.c"
Please explain in the submission comments what you did about this issue.
For more details:
<https://contributor.r-project.org/cran-cookbook/description_issues.html#using-authorsr>

Reply: A designated inst/COPYRIGHTS file has been added to clearly describe all copyright and licensing information for third-party code included in the deepSTRAPP package.
Reference to this file is added to DESCRIPTION. All copyright holders of third-party code are mentioned in Authors@R under the appropriate roles.

