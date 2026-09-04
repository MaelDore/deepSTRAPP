###########################################################################################
##   Does a deepSTRAPP plotting function leave any graphical parameter changed?           ##
##   Author: Maël Doré                                                                   ##
##                                                                                       ##
##   Two sources of leak are tested:                                                      ##
##    * parameters the function, or the functions it calls, set on their own.              ##
##    * parameters the USER passes through '...', which are forwarded to par().             ##
##   The second is the one an earlier version of this script missed entirely.               ##
###########################################################################################

devtools::load_all("."); library(BAMMtools)
data(whale_BAMM_object, package = "deepSTRAPP")

## Parameters that describe the plot just drawn. They are EXPECTED to change, and restoring
## them is what rewinds the panel counter and discards the coordinate system.
## /!\ This list MUST match 'plot_state_par_names' in plot_BAMM_rates(): the function decides
## what not to restore, this script decides what not to report. If the two drift apart, the
## script either hides a real leak, or reports one that is intentional.
## 'pin' is included because it is derived from the plot region, and moves with it.
## 'din', 'cra' and 'cxy' are read-only, so they never appear in par(no.readonly = TRUE) anyway.
## 'mai' is deliberately NOT excluded: it is 'mar' expressed in inches, so a leaked margin shows
## up twice, and silencing one of the two would only make the report harder to read, not safer.
PLOT_STATE <- c("usr", "plt", "fig", "mfg", "new", "xaxp", "yaxp", "pin")

compare_par <- function (label, call_expr, facets = FALSE)
{
  pdf(tempfile(fileext = ".pdf"), width = 10, height = 6)
  if (facets) { par(mfrow = c(1, 2)) }
  plot.new()
  before <- par(no.readonly = TRUE)
  eval(call_expr)
  after <- par(no.readonly = TRUE)
  dev.off()

  changed <- names(before)[!vapply(names(before),
                                   function (p) isTRUE(all.equal(before[[p]], after[[p]])),
                                   logical(1))]
  leaked <- setdiff(changed, PLOT_STATE)

  cat(sprintf("\n  %-46s %s\n", label, ifelse(facets, "mfrow(1,2)", "single")))
  cat(sprintf("    left changed: %s\n", ifelse(length(leaked) == 0, "none", paste(leaked, collapse = ", "))))
  for (p in leaked)
  {
    cat(sprintf("      %-6s %s  ->  %s\n", p,
                paste(format(before[[p]]), collapse = " "), paste(format(after[[p]]), collapse = " ")))
  }
  invisible(leaked)
}

cat("=== plot_BAMM_rates: no extra graphical arguments ===\n")
leaked <- list(
  compare_par("legend = TRUE",  quote(plot_BAMM_rates(whale_BAMM_object, legend = TRUE))),
  compare_par("legend = FALSE", quote(plot_BAMM_rates(whale_BAMM_object, legend = FALSE))),
  compare_par("legend = TRUE",  quote(plot_BAMM_rates(whale_BAMM_object, legend = TRUE)), facets = TRUE))

cat("\n=== plot_BAMM_rates: graphical parameters passed through '...' ===\n")
cat("    (these reach par() through BAMMtools::plot.bammdata(), and are NOT restored by it\n")
cat("     once par.reset defaults to FALSE)\n")
leaked <- c(leaked, list(
  compare_par("mar = c(2,2,2,2)",
              quote(plot_BAMM_rates(whale_BAMM_object, legend = TRUE, mar = c(2, 2, 2, 2)))),
  compare_par("cex.axis / col.lab / lend",
              quote(plot_BAMM_rates(whale_BAMM_object, legend = TRUE,
                                    cex.axis = 2, col.lab = "red", lend = "square"))),
  compare_par("oma / xpd / bty",
              quote(plot_BAMM_rates(whale_BAMM_object, legend = FALSE,
                                    oma = c(2, 2, 2, 2), xpd = NA, bty = "n"))),
  compare_par("las / tcl, faceted",
              quote(plot_BAMM_rates(whale_BAMM_object, legend = TRUE, las = 2, tcl = -0.8)),
              facets = TRUE)))

all_leaked <- unique(unlist(leaked))
cat("\n----------------------------------------------------------------\n")
if (length(all_leaked) == 0)
{
  cat("No graphical setting is left changed, including those passed through '...'.\n")
} else {
  cat("Still leaking: ", paste(all_leaked, collapse = ", "), "\n")
  cat("The restore is built from names(add_args_for_par), so a leak here means a parameter\n")
  cat("reached par() by another route, or was filtered out of 'add_args_for_par'.\n")
}
