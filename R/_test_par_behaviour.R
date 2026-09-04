###########################################################################################
##   Regression tests for graphical-parameter handling in deepSTRAPP plotting functions  ##
##   Author: Maël Doré                                                                   ##
##                                                                                       ##
##   Eight behaviours must hold for every plotting function:                             ##
##                                                                                       ##
##     1 leak         no graphical setting is left changed (CRAN policy)                 ##
##     2 facets       the function's own panels land in distinct device regions          ##
##     3 user layout  a user's par(mfrow = ...) set before the call survives it          ##
##     4 series       consecutive calls do not overlay one another                       ##
##     5 annotate     abline()/title() after the call hit the panel that was drawn       ##
##     6 dots         graphical args passed through '...' do not leak                    ##
##     7 pdf only     display_plot = FALSE leaves the active device untouched            ##
##     8 on error     an error mid-plot still restores the user's parameters             ##
##                                                                                       ##
##   Run as-is to test the package. Set SELFTEST <- TRUE to run the harness against      ##
##   three built-in fixtures instead (one correct, two broken in the two ways this       ##
##   package has actually been broken) and confirm the tests discriminate. Do that       ##
##   after any edit to this file: a test suite that cannot fail is worse than none.      ##
###########################################################################################

SELFTEST <- FALSE

###########################################################################################
##   Recorders                                                                           ##
##                                                                                        ##
##   trace() splices its expression into the traced function's body, so it is evaluated  ##
##   in that function's frame, whose enclosure is namespace:graphics. The only reachable ##
##   environment from there is globalenv(), so the logs live there.                      ##
##                                                                                        ##
##   plot.new()    -> the device region each panel was opened in                          ##
##   plot.window() -> the user coordinate system each panel ended up with. This is the    ##
##                    one that matters for annotation: a full par(no.readonly = TRUE)     ##
##                    restore silently resets usr to 0 1 0 1, which is why abline() and   ##
##                    title() used to land in the wrong place.                            ##
###########################################################################################

.panel_fig <- list()
.panel_usr <- list()
.panel_stop_at <- NA_integer_

start_recorder <- function (stop_at = NA_integer_)
{
  assign(".panel_fig", list(), envir = globalenv())
  assign(".panel_usr", list(), envir = globalenv())
  assign(".panel_stop_at", stop_at, envir = globalenv())

  invisible(suppressMessages(trace(what = "plot.new",
        exit = quote({
          assign(x = ".panel_fig",
                 value = c(get(x = ".panel_fig", envir = globalenv()), list(graphics::par("fig"))),
                 envir = globalenv())
          if (!is.na(get(x = ".panel_stop_at", envir = globalenv())) &&
              length(get(x = ".panel_fig", envir = globalenv())) >=
                get(x = ".panel_stop_at", envir = globalenv()))
          {
            stop("injected failure mid-plot")
          }
        }),
        print = FALSE, where = asNamespace("graphics"))))

  invisible(suppressMessages(trace(what = "plot.window",
        exit = quote({
          assign(x = ".panel_usr",
                 value = c(get(x = ".panel_usr", envir = globalenv()), list(graphics::par("usr"))),
                 envir = globalenv())
        }),
        print = FALSE, where = asNamespace("graphics"))))

  invisible(NULL)
}

stop_recorder <- function ()
{
  try(expr = invisible(suppressMessages(untrace(what = "plot.new",
                                                where = asNamespace("graphics")))), silent = TRUE)
  try(expr = invisible(suppressMessages(untrace(what = "plot.window",
                                                where = asNamespace("graphics")))), silent = TRUE)
  list(fig = get(x = ".panel_fig", envir = globalenv()),
       usr = get(x = ".panel_usr", envir = globalenv()))
}

n_panels_so_far <- function () { length(get(x = ".panel_fig", envir = globalenv())) }

close_devices <- function () { while (grDevices::dev.cur() > 1) { grDevices::dev.off() } }

format_region <- function (z) { sprintf("x[%.2f-%.2f] y[%.2f-%.2f]", z[1], z[2], z[3], z[4]) }

count_pdf_pages <- function (path)
{
  txt <- readLines(con = path, warn = FALSE)
  length(grep(pattern = "/Type\\s*/Page[^s]", x = txt))
}

## Graphical settings that describe where the last plot went, not a user preference.
## These are never restored: doing so rewinds the panel counter and discards the
## coordinate system, which is what breaks faceting and post-hoc annotation.
PLOT_STATE <- c("usr", "plt", "fig", "mfg", "new", "xaxp", "yaxp", "pin")

par_diff <- function (before, after)
{
  changed <- names(before)[!vapply(X = names(before),
                                   FUN = function (p) { isTRUE(all.equal(before[[p]], after[[p]])) },
                                   FUN.VALUE = logical(1))]
  setdiff(changed, PLOT_STATE)
}

###########################################################################################
##   Result collection                                                                   ##
###########################################################################################

RESULTS <- data.frame(target = character(0), test = character(0),
                      result = character(0), detail = character(0),
                      stringsAsFactors = FALSE)

record <- function (target, test, pass, detail = "")
{
  RESULTS[nrow(RESULTS) + 1L, ] <<- list(target, test,
                                         ifelse(is.na(pass), "SKIP", ifelse(pass, "PASS", "FAIL")),
                                         detail)
  invisible(NULL)
}

###########################################################################################
##   Tests                                                                               ##
###########################################################################################

## 1 / 6   no graphical setting left changed
test_leak <- function (target, extra = list(), test_name = "1 leak")
{
  on.exit({ stop_recorder() ; close_devices() }, add = TRUE)
  start_recorder()
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 12, height = 6)
  graphics::plot.new()
  before <- graphics::par(no.readonly = TRUE)
  err <- NULL
  tryCatch(expr = target$call(extra), error = function (e) { err <<- conditionMessage(e) })
  after <- graphics::par(no.readonly = TRUE)
  leaked <- par_diff(before, after)
  stop_recorder() ; close_devices()

  if (!is.null(err)) { record(target$id, test_name, NA, paste("call failed:", err)) ; return(invisible(NULL)) }
  record(target$id, test_name, length(leaked) == 0,
         ifelse(length(leaked) == 0, "nothing left changed",
                paste("leaked:", paste(leaked, collapse = ", "))))
}

## 2 / 3   the function's own facets, and the user's layout
test_facets_and_layout <- function (target)
{
  on.exit({ stop_recorder() ; close_devices() }, add = TRUE)
  start_recorder()
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 12, height = 6)
  user_mfrow <- c(2, 2)
  graphics::par(mfrow = user_mfrow)
  graphics::plot.new()
  n0 <- n_panels_so_far()
  err <- NULL
  tryCatch(expr = target$call(), error = function (e) { err <<- conditionMessage(e) })
  mfrow_after <- graphics::par("mfrow")
  log <- stop_recorder() ; close_devices()

  if (!is.null(err)) {
    record(target$id, "2 facets", NA, paste("call failed:", err))
    record(target$id, "3 user layout", NA, paste("call failed:", err))
    return(invisible(NULL))
  }
  panels <- log$fig[seq_along(log$fig) > n0]
  regions <- unique(vapply(X = panels, FUN = format_region, FUN.VALUE = character(1)))

  if (is.null(target$panels) || is.na(target$panels))
  {
    ## Nothing to assert against yet. Report what the function actually does so the
    ## number can be pinned in TARGETS, rather than inventing an expectation.
    record(target$id, "2 facets", NA,
           sprintf("observed %d panel(s) in %d distinct region(s) - pin panels = %d",
                   length(panels), length(regions), length(regions)))
  } else {
    record(target$id, "2 facets", length(regions) == target$panels,
           sprintf("%d panel(s) in %d distinct region(s), expected %d",
                   length(panels), length(regions), target$panels))
  }
  record(target$id, "3 user layout", isTRUE(all.equal(mfrow_after, user_mfrow)),
         sprintf("user mfrow 2x2 -> %s", paste(mfrow_after, collapse = "x")))
}

## 4   consecutive calls must not overlay
##
##     Region cycling alone cannot see this for a single-panel function (every call uses
##     the whole device either way), so the page count carries the test: one call should
##     consume exactly one page. A leftover par(new = TRUE) merges two calls onto one.
##     Nothing is annotated between the calls on purpose - title() clears `new` and would
##     mask the very bug this test is for.
test_series <- function (target, n_calls = NULL)
{
  if (is.null(n_calls)) { n_calls <- if (is.null(target$n_calls)) 3L else target$n_calls }
  on.exit({ stop_recorder() ; close_devices() }, add = TRUE)

  ## Calibration: how many pages does ONE call legitimately consume? Functions that loop
  ## over time slices produce many, so this cannot be assumed to be 1. A single call is
  ## never affected by the overlay bug, which only shows up between calls.
  calib_file <- tempfile(fileext = ".pdf")
  start_recorder()
  grDevices::pdf(calib_file, width = 12, height = 6, compress = FALSE)
  err <- NULL
  tryCatch(expr = target$call(), error = function (e) { err <<- conditionMessage(e) })
  calib_log <- stop_recorder() ; close_devices()
  if (!is.null(err)) { record(target$id, "4 series", NA, paste("call failed:", err)) ; return(invisible(NULL)) }
  pages_per_call <- count_pdf_pages(calib_file)

  ## How many distinct regions one call uses. Taken from the calibration run when the
  ## target does not declare it, so a new function can be tested without guessing.
  calib_regions <- unique(vapply(X = calib_log$fig, FUN = format_region, FUN.VALUE = character(1)))
  panels_used <- if (is.null(target$panels) || is.na(target$panels)) {
    max(1L, length(calib_regions))
  } else { target$panels }
  if (pages_per_call == 0) { record(target$id, "4 series", NA, "one call produced no page") ; return(invisible(NULL)) }

  device_file <- tempfile(fileext = ".pdf")
  start_recorder()
  grDevices::pdf(device_file, width = 12, height = 6, compress = FALSE)
  tryCatch(expr = { for (i in seq_len(n_calls)) { target$call() } },
           error = function (e) { err <<- conditionMessage(e) })
  new_after <- graphics::par("new")
  log <- stop_recorder() ; close_devices()

  if (!is.null(err)) { record(target$id, "4 series", NA, paste("call failed:", err)) ; return(invisible(NULL)) }

  n_pages <- count_pdf_pages(device_file)
  expected_pages <- n_calls * pages_per_call
  regions <- vapply(X = log$fig, FUN = format_region, FUN.VALUE = character(1))
  expected <- rep(regions[seq_len(min(panels_used, length(regions)))], length.out = length(regions))
  bad <- which(regions != expected)

  ok <- (n_pages == expected_pages) && (length(bad) == 0) && !isTRUE(new_after)
  record(target$id, "4 series", ok,
         sprintf("%d calls -> %d page(s) (expect %d = %d x %d/call), %d panel(s) repainted, par('new') left %s",
                 n_calls, n_pages, expected_pages, n_calls, pages_per_call, length(bad), new_after))
}

## 5   annotation after the call must reach the panel that was drawn
##
##     The decisive quantity is usr. A full par(no.readonly = TRUE) restore resets it to
##     0 1 0 1, so abline(v = 20) silently lands off-plot and title() attaches to the
##     wrong facet. fig is not a usable check: for a two-panel function that restores
##     mfrow it legitimately widens back to the whole device.
test_annotate <- function (target)
{
  on.exit({ stop_recorder() ; close_devices() }, add = TRUE)
  start_recorder()
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 12, height = 6)
  graphics::par(mfrow = c(1, 2))
  err <- NULL
  tryCatch(expr = target$call(), error = function (e) { err <<- conditionMessage(e) })

  drawn_usr <- get(x = ".panel_usr", envir = globalenv())
  n_before <- n_panels_so_far()
  usr_after <- graphics::par("usr")

  ## Bail out BEFORE annotating. If the call failed there is no plot on the device and
  ## abline() would raise "plot.new has not been called yet", aborting the whole run
  ## instead of recording one skipped test.
  if (!is.null(err))
  {
    stop_recorder() ; close_devices()
    record(target$id, "5 annotate", NA, paste("call failed:", err))
    return(invisible(NULL))
  }
  if (length(drawn_usr) == 0)
  {
    stop_recorder() ; close_devices()
    record(target$id, "5 annotate", NA, "the call drew nothing (no plot.window() seen)")
    return(invisible(NULL))
  }

  ## An error raised by the annotation itself is a genuine failure of this test:
  ## it means the call left the device in a state where nothing can be added.
  annotate_err <- NULL
  tryCatch(expr = {
    graphics::abline(v = usr_after[1] + 0.5 * diff(usr_after[1:2]), col = "red", lwd = 3)
    graphics::title(main = "annotation")
    graphics::text(x = usr_after[1] + 0.5 * diff(usr_after[1:2]),
                   y = usr_after[3] + 0.5 * diff(usr_after[3:4]), labels = "x")
  }, error = function (e) { annotate_err <<- conditionMessage(e) })

  n_after <- n_panels_so_far()
  stop_recorder() ; close_devices()

  if (!is.null(annotate_err))
  {
    record(target$id, "5 annotate", FALSE, paste("annotation raised:", annotate_err))
    return(invisible(NULL))
  }

  last_usr <- drawn_usr[[length(drawn_usr)]]
  same_usr <- isTRUE(all.equal(last_usr, usr_after))
  no_new_panel <- (n_after == n_before)

  record(target$id, "5 annotate", same_usr && no_new_panel,
         sprintf("usr %s, %d new panel(s) opened by abline/title/text",
                 ifelse(same_usr, "is the drawn panel's",
                        sprintf("RESET to %s (drawn panel had %s)",
                                paste(round(usr_after, 2), collapse = " "),
                                paste(round(last_usr, 2), collapse = " "))),
                 n_after - n_before))
}

## 7   display_plot = FALSE must leave the active device alone
test_pdf_only <- function (target)
{
  if (!isTRUE(target$has_pdf_arg))
  {
    record(target$id, "7 pdf only", NA, "no PDF_file_path argument")
    return(invisible(NULL))
  }
  on.exit({ close_devices() }, add = TRUE)
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 12, height = 6)
  graphics::plot(1:3)
  graphics::par(new = TRUE)                        # a flag the user set on purpose
  before <- graphics::par(no.readonly = TRUE)
  err <- NULL
  ## The flag that gates drawing is not always called display_plot: plot_histograms_*
  ## uses display_plots, prepare_trait_data uses plot_map.
  off_args <- list()
  off_args[[if (is.null(target$display_arg)) "display_plot" else target$display_arg]] <- FALSE
  off_args$PDF_file_path <- tempfile(fileext = ".pdf")

  tryCatch(expr = target$call(off_args),
           error = function (e) { err <<- conditionMessage(e) })
  after <- graphics::par(no.readonly = TRUE)
  new_after <- graphics::par("new")
  close_devices()

  if (!is.null(err)) { record(target$id, "7 pdf only", NA, paste("call failed:", err)) ; return(invisible(NULL)) }
  leaked <- par_diff(before, after)
  record(target$id, "7 pdf only", isTRUE(new_after) && length(leaked) == 0,
         sprintf("user par(new = TRUE) left %s%s", new_after,
                 ifelse(length(leaked) == 0, "", paste("; leaked:", paste(leaked, collapse = ", ")))))
}

## 8   an error mid-plot must still restore the user's parameters
test_error_recovery <- function (target)
{
  on.exit({ stop_recorder() ; close_devices() }, add = TRUE)
  start_recorder(stop_at = 2L)                     # abort inside the target's first panel
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 12, height = 6)
  user_mfrow <- c(2, 2)
  graphics::par(mfrow = user_mfrow)
  graphics::plot.new()                             # panel 1, counted by the recorder
  before <- graphics::par(no.readonly = TRUE)
  raised <- FALSE
  tryCatch(expr = target$call(), error = function (e) { raised <<- TRUE })
  after <- graphics::par(no.readonly = TRUE)
  mfrow_after <- graphics::par("mfrow")
  stop_recorder() ; close_devices()

  if (!raised)
  {
    record(target$id, "8 on error", NA, "injector did not fire")
    return(invisible(NULL))
  }
  leaked <- par_diff(before, after)
  record(target$id, "8 on error", length(leaked) == 0 && isTRUE(all.equal(mfrow_after, user_mfrow)),
         sprintf("mfrow %s%s", paste(mfrow_after, collapse = "x"),
                 ifelse(length(leaked) == 0, "", paste("; leaked:", paste(leaked, collapse = ", ")))))
}

## No single test may abort the run: a crash is recorded and the suite carries on.
guarded <- function (target, label, expr)
{
  tryCatch(expr = expr(),
           error = function (e) {
             stop_recorder() ; close_devices()
             record(target$id, label, FALSE, paste("harness aborted:", conditionMessage(e)))
           })
  invisible(NULL)
}

run_all_tests <- function (target)
{
  cat(sprintf("\n--- %s ---\n", target$id))
  guarded(target, "1 leak",        function () { test_leak(target) })
  guarded(target, "2 facets",      function () { test_facets_and_layout(target) })
  guarded(target, "4 series",      function () { test_series(target) })
  guarded(target, "5 annotate",    function () { test_annotate(target) })
  guarded(target, "6 dots",        function () { test_leak(target,
                                     extra = list(las = 2, tcl = -0.8, xpd = NA),
                                     test_name = "6 dots") })
  guarded(target, "7 pdf only",    function () { test_pdf_only(target) })
  guarded(target, "8 on error",    function () { test_error_recovery(target) })
  invisible(NULL)
}

###########################################################################################
##   Fixtures for the harness self-test                                                  ##
##                                                                                        ##
##   Identical drawing, three parameter policies:                                         ##
##     good      diff restore + guarded par(new = FALSE)                    -> all pass    ##
##     snapshot  the dd418ac policy, par(no.readonly = TRUE) restored       -> annotation  ##
##     no_new    diff restore without the `new` handling                    -> series      ##
###########################################################################################

make_fixture <- function (policy, n_panels = 2L)
{
  ## Both are promises captured by the returned closure. Without forcing them here they
  ## are evaluated at call time, by which point the loop below has moved on, and every
  ## fixture silently becomes the last one.
  force(policy) ; force(n_panels)

  function (display_plot = TRUE, PDF_file_path = NULL, ...)
  {
    entry_par <- graphics::par(no.readonly = TRUE)

    if (policy == "snapshot")
    {
      on.exit(graphics::par(entry_par), add = TRUE)
    } else {
      on.exit({
        exit_par <- graphics::par(no.readonly = TRUE)
        changed <- par_diff(entry_par, exit_par)
        if (length(changed) > 0) { graphics::par(entry_par[changed]) }
        if (policy == "good" &&
            (isTRUE(display_plot) || (isTRUE(exit_par$new) && !isTRUE(entry_par$new))))
        {
          graphics::par(new = FALSE)
        }
      }, add = TRUE)
    }

    draw <- function ()
    {
      initial_par <- graphics::par("oma", "mfrow")
      if (n_panels == 2L) { graphics::par(mfrow = c(1, 2), oma = c(0, 0, 3, 0)) }
      graphics::par(mar = c(6, 6, 4, 2))
      graphics::plot(1:3, main = "A", ...)
      if (n_panels == 2L) { graphics::plot(4:9, main = "B", ...) }
      graphics::par(new = TRUE)        # stand-in for a dependency leaving the flag up
      if (n_panels == 2L) { graphics::par(initial_par) }
    }

    if (!is.null(PDF_file_path))
    {
      grDevices::pdf(PDF_file_path, width = 12, height = 6)
      draw()
      grDevices::dev.off()
    }
    if (isTRUE(display_plot)) { draw() }
    invisible(NULL)
  }
}

if (isTRUE(SELFTEST))
{
  TARGETS <- list()
  for (n_panels in c(1L, 2L))
  {
    for (policy in c("good", "snapshot", "no_new"))
    {
      local({
        fx <- make_fixture(policy, n_panels)
        TARGETS[[length(TARGETS) + 1L]] <<- list(
          id = sprintf("%s / %d panel", policy, n_panels),
          panels = n_panels, has_pdf_arg = TRUE,
          call = function (extra = list()) { do.call(what = fx, args = extra) })
      })
    }
  }

} else {

###########################################################################################
##   Targets                                                                             ##
##                                                                                        ##
##   Add the remaining plotting functions here as you convert them. 'panels' is how many ##
##   distinct device regions one call is expected to use.                                ##
###########################################################################################

  devtools::load_all(".") ; library(BAMMtools) ; library(ape)

  deepSTRAPP_outputs <- readRDS(system.file("extdata",
    "Ponerinae_deepSTRAPP_cont_old_calib_0_40.rds", package = "deepSTRAPP"))

  ## Objects are found by class wherever they sit: several elements are LISTS of objects
  ## (note the plural names), so a direct $ pull hands the function a list and every test
  ## skips with an unhelpful "call failed".
  find_by_class <- function (x, cls, depth = 4L)
  {
    if (inherits(x, cls)) { return(x) }
    if (depth <= 0L || !is.list(x)) { return(NULL) }
    for (el in x) { f <- find_by_class(el, cls, depth - 1L) ; if (!is.null(f)) { return(f) } }
    NULL
  }
  ## a list whose elements are all of class `cls` (what plot_densityMaps_overlay wants)
  find_list_of_class <- function (x, cls, depth = 4L)
  {
    if (is.list(x) && !inherits(x, cls) && length(x) > 0L &&
        all(vapply(X = x, FUN = function (el) { inherits(el, cls) }, FUN.VALUE = logical(1))))
    {
      return(x)
    }
    if (depth <= 0L || !is.list(x)) { return(NULL) }
    for (el in x) { f <- find_list_of_class(el, cls, depth - 1L) ; if (!is.null(f)) { return(f) } }
    NULL
  }

  BAMM_object_test  <- find_by_class(deepSTRAPP_outputs, "bammdata")
  contMap_test      <- find_by_class(deepSTRAPP_outputs, "contMap")
  densityMaps_test  <- find_list_of_class(deepSTRAPP_outputs, "densityMap")

  cat(sprintf("bammdata found: %s | contMap found: %s | densityMaps found: %s\n",
              !is.null(BAMM_object_test), !is.null(contMap_test), !is.null(densityMaps_test)))
  if (is.null(BAMM_object_test))
  {
    stop(paste0("No 'bammdata' object inside deepSTRAPP_outputs. Elements: ",
                paste(names(deepSTRAPP_outputs), collapse = ", ")))
  }

###########################################################################################
##   Targets                                                                             ##
##                                                                                        ##
##   panels  = distinct device regions one call is expected to use. Leave NA on a first   ##
##             run: test 2 then REPORTS the observed number instead of asserting, and     ##
##             test 4 derives it from a calibration call. Pin it once you know it, so     ##
##             the test can actually fail afterwards.                                     ##
##   display_arg = the formal that gates whether anything is drawn. Defaults to           ##
##             "display_plot"; set it where the name differs.                             ##
##   n_calls = how many consecutive calls test 4 makes (default 3). Lower it for slow     ##
##             functions that loop over every time slice.                                 ##
###########################################################################################

  TARGETS <- list(

    list(id = "plot_STRAPP_pvalues_over_time", panels = NA, has_pdf_arg = TRUE,
         call = function (extra = list()) {
           do.call(what = plot_STRAPP_pvalues_over_time,
                   args = c(list(deepSTRAPP_outputs = deepSTRAPP_outputs), extra)) }),

    list(id = "plot_histogram_STRAPP_test_for_focal_time", panels = NA, has_pdf_arg = TRUE,
         call = function (extra = list()) {
           do.call(what = plot_histogram_STRAPP_test_for_focal_time,
                   args = c(list(deepSTRAPP_outputs = deepSTRAPP_outputs, focal_time = 10), extra)) }),

    ## display_plots, not display_plot
    list(id = "plot_histograms_STRAPP_tests_over_time", panels = NA, has_pdf_arg = TRUE,
         display_arg = "display_plots", n_calls = 2L,
         call = function (extra = list()) {
           do.call(what = plot_histograms_STRAPP_tests_over_time,
                   args = c(list(deepSTRAPP_outputs = deepSTRAPP_outputs), extra)) }),

    list(id = "plot_rates_through_time", panels = NA, has_pdf_arg = TRUE,
         call = function (extra = list()) {
           do.call(what = plot_rates_through_time,
                   args = c(list(deepSTRAPP_outputs = deepSTRAPP_outputs), extra)) }),

    list(id = "plot_rates_vs_trait_data_for_focal_time", panels = NA, has_pdf_arg = TRUE,
         call = function (extra = list()) {
           do.call(what = plot_rates_vs_trait_data_for_focal_time,
                   args = c(list(deepSTRAPP_outputs = deepSTRAPP_outputs, focal_time = 10), extra)) }),

    list(id = "plot_rates_vs_trait_data_over_time", panels = NA, has_pdf_arg = TRUE,
         n_calls = 2L,
         call = function (extra = list()) {
           do.call(what = plot_rates_vs_trait_data_over_time,
                   args = c(list(deepSTRAPP_outputs = deepSTRAPP_outputs), extra)) }),

    ## These two take a map object rather than the deepSTRAPP output. If the finders above
    ## reported FALSE, point them at a contMap / list of densityMaps you have to hand.
    list(id = "plot_contMap", panels = NA, has_pdf_arg = TRUE,
         call = function (extra = list()) {
           do.call(what = plot_contMap, args = c(list(contMap = contMap_test), extra)) }),

    list(id = "plot_densityMaps_overlay", panels = NA, has_pdf_arg = TRUE,
         call = function (extra = list()) {
           do.call(what = plot_densityMaps_overlay,
                   args = c(list(densityMaps = densityMaps_test), extra)) })
  )

###########################################################################################
##   Not included: prepare_diversification_data and prepare_trait_data                   ##
##                                                                                        ##
##   Both were patched, but neither can be driven from a deepSTRAPP output: the first     ##
##   needs a BAMM install directory and runs an MCMC, the second needs raw tip data and   ##
##   fits evolutionary models. Running them here would take hours and test the fitting,   ##
##   not the graphics. Their on.exit block is the same one the eight targets above        ##
##   exercise, with the draw gate set to plot_evaluations / plot_map respectively.        ##
##                                                                                        ##
##   To test them anyway, append a target with arguments you already have cached, e.g.    ##
##                                                                                        ##
##     list(id = "prepare_trait_data", panels = NA, has_pdf_arg = TRUE,                   ##
##          display_arg = "plot_map",                                                     ##
##          call = function (extra = list()) {                                            ##
##            do.call(what = prepare_trait_data,                                          ##
##                    args = c(list(tip_data = my_tip_data,                               ##
##                                  trait_data_type = "continuous",                       ##
##                                  phylo = my_phylo), extra)) })                         ##
###########################################################################################

}

###########################################################################################
##   Run                                                                                 ##
###########################################################################################

cat("===========================================================\n")
cat(ifelse(isTRUE(SELFTEST), " HARNESS SELF-TEST\n",
                             " deepSTRAPP graphical-parameter regression tests\n"))
cat("===========================================================\n")

for (target in TARGETS) { run_all_tests(target) }

cat("\n===========================================================\n")
for (i in seq_len(nrow(RESULTS)))
{
  cat(sprintf("  %-5s %-14s %-46s %s\n", RESULTS$result[i], RESULTS$test[i],
              RESULTS$target[i], RESULTS$detail[i]))
}
cat("===========================================================\n")

if (isTRUE(SELFTEST))
{
  fails <- function (pattern) { sum(RESULTS$result == "FAIL" & grepl(pattern, RESULTS$target)) }
  cat(sprintf("\n  good     failures: %d   (must be 0)\n", fails("^good")))
  cat(sprintf("  snapshot failures: %d   (must be > 0 - the annotation policy)\n", fails("^snapshot")))
  cat(sprintf("  no_new   failures: %d   (must be > 0 - the series policy)\n", fails("^no_new")))
  cat(sprintf("\n  HARNESS %s\n",
              ifelse(fails("^good") == 0 && fails("^snapshot") > 0 && fails("^no_new") > 0,
                     "TRUSTWORTHY", "NOT TRUSTWORTHY - it cannot tell good from broken")))
} else {
  cat(sprintf("\n  %d failure(s), %d skipped\n",
              sum(RESULTS$result == "FAIL"), sum(RESULTS$result == "SKIP")))
}
