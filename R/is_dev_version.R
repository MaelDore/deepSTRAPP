
#' @title Check whether deepSTRAPP is a development version
#'
#' @description Detect whether the current deepSTRAPP installation corresponds to
#'   a development version (e.g., version sourced from GitHub), as opposed to a CRAN release.
#'
#'   The development versions include deepSTRAPP outputs as additional internal datasets to help
#'   produces examples and vignettes outputs. These additional datasets are removed from the CRAN releases
#'   because their size is not compatible with CRAN policies.
#'
#'   This function is only used to check if an example must be ran or a vignette chunk evaluated
#'   to produce output from data.
#'
#' @param pkg Character string. Name of the R package for which the version type is inspected.
#'   Default is "deepSTRAPP".
#'
#' @return Logical. TRUE if running a development version.
#'
#' @export
#' @importFrom utils packageVersion
#'
#' @author Maël Doré
#'
#' @examples
#' # Check the current deepSTRAPP installation is a development version
#' is_dev_version(pkg = "deepSTRAPP")
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
