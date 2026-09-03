###########################################################################################
##                                                                                       ##
##   Static audit of the @examples blocks of a package                                   ##
##                                                                                       ##
##   Author: Maël Doré                                                                   ##
##                                                                                       ##
##   Inspects example code WITHOUT running it, to anticipate the runtime errors that      ##
##   R CMD check would otherwise surface one at a time. Requires check_examples_parsing.R ##
##   to be sourced first, for prepare_example_code() and extract_examples_blocks().       ##
##                                                                                       ##
##   Usage, from the root of the package directory:                                       ##
##     source("dev/check_examples_parsing.R")                                             ##
##     source("dev/audit_examples.R")                                                     ##
##     audit_examples()                                                                   ##
##                                                                                       ##
##   Checks performed:                                                                    ##
##     ARGUMENT  a named argument that is not a formal of the function being called       ##
##     DATASET   data() referring to a dataset absent from data/                          ##
##     VARIABLE  an object used but never created within the block                        ##
##     PATH      an input file that is not shipped with the package                       ##
##                                                                                        ##
###########################################################################################

audit_examples <- function (path = ".", include_dontrun = TRUE)
{
  r_files <- list.files(file.path(path, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
  datasets <- sub("\\.rda$", "", list.files(file.path(path, "data"), pattern = "\\.rda$"))
  extdata <- list.files(file.path(path, "inst", "extdata"))

  ## Formals of every function defined in the package
  formals_of <- list()
  for (focal_file in r_files)
  {
    exprs <- tryCatch(parse(focal_file), error = function (e) NULL)
    for (e in exprs)
    {
      if (is.call(e) && (as.character(e[[1]])[1] %in% c("<-", "=")) &&
          is.call(e[[3]]) && (as.character(e[[3]][[1]])[1] == "function"))
      {
        formals_of[[as.character(e[[2]])]] <- names(as.list(e[[3]][[2]]))
      }
    }
  }

  findings <- data.frame(file = character(0), line = integer(0), severity = character(0),
                         kind = character(0), detail = character(0), stringsAsFactors = FALSE)
  add <- function (file, line, severity, kind, detail)
  {
    findings <<- rbind(findings, data.frame(file = basename(file), line = line, severity = severity,
                                            kind = kind, detail = detail, stringsAsFactors = FALSE))
  }

  for (focal_file in r_files)
  {
    for (block in extract_examples_blocks(focal_file))
    {
      raw <- block$code
      first <- block$first_code_line
      code <- prepare_example_code(raw)
      in_dontrun <- grepl("\\\\dontrun", paste(raw, collapse = "\n"))
      if (in_dontrun & !include_dontrun) { next }
      tag <- ifelse(in_dontrun, " [dontrun]", "")

      exprs <- tryCatch(parse(text = code), error = function (e) NULL)
      if (is.null(exprs)) { next }   # handled by check_examples_parsing()

      ## ---- datasets --------------------------------------------------------------
      ## 'data(' must not be preceded by a word character, to avoid matching calls such
      ## as prepare_trait_data(tip_data = ...)
      for (k in seq_along(raw))
      {
        hits <- regmatches(raw[k], gregexpr('(^|[^[:alnum:]._])data\\([^)]*\\)', raw[k]))[[1]]
        for (h in hits)
        {
          nm <- sub('^.*data\\(\\s*["\']?([A-Za-z0-9_.]+).*$', "\\1", h)
          pkg <- if (grepl('package\\s*=\\s*["\']([^"\']+)', h))
                   sub('.*package\\s*=\\s*["\']([^"\']+).*', "\\1", h) else NA_character_
          if ((is.na(pkg) | identical(pkg, "deepSTRAPP")) & !(nm %in% datasets))
          {
            add(focal_file, first + k - 1L, ifelse(is.na(pkg), "warning", "ERROR"), "DATASET",
                sprintf("data(%s)%s : not in data/%s", nm, tag,
                        ifelse(is.na(pkg), " (no package= given; fine if another package supplies it)", "")))
          }
        }
      }

      ## ---- named arguments, with R's partial matching rules -----------------------
      walk <- function (e)
      {
        if (!is.call(e)) { return(invisible(NULL)) }
        fn <- e[[1]]
        fname <- if (is.name(fn)) as.character(fn)
                 else if (is.call(fn) && (as.character(fn[[1]])[1] %in% c("::", ":::"))) as.character(fn[[3]])
                 else NA_character_

        if (!is.na(fname) && !is.null(formals_of[[fname]]))
        {
          fmls <- formals_of[[fname]]
          given <- names(as.list(e))[-1]
          given <- given[nzchar(given)]
          dots_position <- match("...", fmls, nomatch = length(fmls) + 1L)

          for (g in setdiff(given, fmls))
          {
            ## R matches a named argument exactly first, then by unique partial match,
            ## but only against formals that come before '...'
            candidates <- fmls[seq_len(dots_position - 1L)]
            candidates <- setdiff(candidates, intersect(given, fmls))
            partial <- candidates[startsWith(candidates, g)]

            if (length(partial) == 1L)
            {
              add(focal_file, first, "warning", "ARGUMENT",
                  sprintf("%s(%s = )%s : matches '%s' only by partial matching", fname, g, tag, partial))
            } else if ("..." %in% fmls) {
              add(focal_file, first, "note", "ARGUMENT",
                  sprintf("%s(%s = )%s : not a formal, absorbed by '...'", fname, g, tag))
            } else {
              add(focal_file, first, "ERROR", "ARGUMENT",
                  sprintf("%s(%s = )%s : not a formal and no '...'. Formals: %s",
                          fname, g, tag, paste(fmls, collapse = ", ")))
            }
          }
        }
        ## Recurse into the arguments. Empty arguments, as in x[, 1], are missing symbols:
        ## touching one throws, so the extraction itself is wrapped, not just the recursion.
        arg_list <- as.list(e)[-1]
        for (idx in seq_along(arg_list))
        {
          try(walk(arg_list[[idx]]), silent = TRUE)
        }
      }
      for (e in exprs) { walk(e) }

      ## ---- objects used but never created ------------------------------------------
      ## Only variables are checked. Functions are skipped, as an example may legitimately
      ## call a function from any package it attaches with library().
      loaded <- unlist(regmatches(raw, gregexpr('(^|[^[:alnum:]._])data\\(\\s*["\']?[A-Za-z0-9_.]+',
                                                raw, perl = TRUE)))
      loaded <- sub('^.*data\\(\\s*["\']?', "", loaded)
      known <- unique(c(loaded, datasets, "T", "F"))

      wrapped <- tryCatch(eval(parse(text = paste0("function() {\n", paste(code, collapse = "\n"), "\n}"))),
                          error = function (e) NULL)
      if (!is.null(wrapped))
      {
        msgs <- character(0)
        suppressWarnings(try(codetools::checkUsage(wrapped, name = "example",
                                                   report = function (x) { msgs <<- c(msgs, x) }),
                             silent = TRUE))
        for (m in msgs)
        {
          if (grepl("no visible binding for global variable", m))
          {
            v <- sub(".*global variable [\u2018'\"]([^\u2019'\"]+)[\u2019'\"].*", "\\1", m)
            if (!(v %in% known) && !exists(v, envir = baseenv()))
            {
              add(focal_file, first, ifelse(in_dontrun, "warning", "ERROR"), "VARIABLE",
                  sprintf("'%s'%s is used but never created in this block", v, tag))
            }
          }
        }
      }

      ## ---- input files ---------------------------------------------------------------
      ## Only files that are read. Output paths (.pdf, and anything passed to a *file_path
      ## argument) are skipped, and system.file() is resolved against inst/.
      for (k in seq_along(raw))
      {
        ## Only lines that actually read a file. This excludes output paths, and dataset
        ## names such as data("eel.tree", package = "phytools") that merely look like filenames.
        if (!grepl("read\\.|readRDS|\\bload\\(|\\bsource\\(|\\bscan\\(|\\bfile\\(", raw[k])) { next }
        if (grepl("file_path\\s*=|PDF_file_path|\\.pdf", raw[k])) { next }
        quoted <- regmatches(raw[k], gregexpr('"[^"]+\\.(txt|tree|csv|nex|rds)"', raw[k]))[[1]]
        for (q in quoted)
        {
          rel <- gsub('^"|"$', "", q)
          base <- basename(rel)
          ok <- file.exists(file.path(path, rel)) | (base %in% extdata)
          if (!ok)
          {
            add(focal_file, first + k - 1L, ifelse(in_dontrun, "note", "ERROR"), "PATH",
                sprintf("'%s'%s is not shipped in the package sources", rel, tag))
          }
        }
      }
    }
  }

  ## ---- report -------------------------------------------------------------------------
  findings <- findings[!duplicated(paste(findings$file, findings$kind, findings$detail)), ]
  findings$severity <- factor(findings$severity, levels = c("ERROR", "warning", "note"))
  findings <- findings[order(findings$severity, findings$file, findings$line), ]

  for (lv in levels(findings$severity))
  {
    sub_findings <- findings[findings$severity == lv, ]
    if (nrow(sub_findings) == 0) { next }
    cat(sprintf("\n======== %s (%d) ========\n", lv, nrow(sub_findings)))
    for (i in seq_len(nrow(sub_findings)))
    {
      cat(sprintf("  %-46s L%-6s %-9s %s\n", sub_findings$file[i], sub_findings$line[i],
                  sub_findings$kind[i], sub_findings$detail[i]))
    }
  }
  cat(sprintf("\n%d finding(s).\n\n", nrow(findings)))

  return(invisible(findings))
}
