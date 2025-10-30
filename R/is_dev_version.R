#' @title Check if package is in development version or not
#'
#' @description Detect if the deepSTRAPP package is in development or released version.
#'   The development versions include deepSTRAPP outputs as datasets to help
#'   produces examples and vignettes outputs.
#'
#'   These datasets are removed from the CRAN release because their size is not compatible
#'   with CRAN policies.
#'
#'   This function is used to check if an example must be ran or a vignette chunk evaluated
#'   to produce output from data.
#'
#' @return Logical. TRUE if development version or local check.
#'
#' @keywords internal
#' @noRd
#'

is_dev_version <- function (pkg = "deepSTRAPP")
{
  # # Check if ran on CRAN
  # not_cran <- identical(Sys.getenv("NOT_CRAN"), "true") # || interactive()

  # Version number check
  version <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "")
  dev_version <- grepl("\\.9000", version)

  # not_cran || dev_version

  return(dev_version)
}
