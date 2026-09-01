###########################################################################################
##                                                                                       ##
##   Check that all roxygen2 @examples blocks parse as valid R code                      ##
##                                                                                       ##
##   Author: Maël Doré                                                                   ##
##                                                                                       ##
##   'R CMD check' parses the examples of every .Rd file at once, in the step             ##
##   "checking for unstated dependencies in examples". A syntax error there is reported    ##
##   against a line number in that concatenated file, which is impossible to trace back    ##
##   to a source file. This script runs the same parsing, but one @examples block at a     ##
##   time, and reports the offending line in the R/*.R file where it must be fixed.        ##
##                                                                                        ##
##   A typical cause is a line of prose inside @examples that lost its '#', so that        ##
##   roxygen2 emits it as R code.                                                          ##
##                                                                                        ##
##   Usage, from the root of the package directory:                                        ##
##     source("dev/check_examples_parsing.R")                                              ##
##     check_examples_parsing()                # scan R/*.R, before running document()      ##
##     check_examples_parsing(source = "man")  # scan man/*.Rd, exactly as R CMD check does ##
##                                                                                        ##
###########################################################################################


## Helper function to turn an @examples block into parseable R code ####

#' @title Convert the content of an @examples block into parseable R code
#'
#' @description Replace the Rd macros that are allowed inside `@examples`, but that are not
#'   valid R syntax, so that the block can be handed to [base::parse()].
#'
#'   `\dontrun{`, `\donttest{`, `\dontshow{` and `\testonly{` are removed together with their
#'   matching closing brace, splicing their content back inline. The matching brace is found by
#'   counting every brace, including those that sit inside an R comment, which is how the Rd parser
#'   itself matches them. `\%` becomes `%`, and `\dots` becomes `...`.
#'
#'   Substitutions are made in place, so that the returned vector holds exactly as many lines
#'   as the input, and a parse error can be traced back to its line in the source file.
#'
#' @param code_lines Vector of character strings. The content of an `@examples` block,
#'   with the leading `#'` already removed.
#'
#' @return A vector of character strings of the same length as `code_lines`.
#'
#' @author Maël Doré
#'
#' @noRd
#'

prepare_example_code <- function (code_lines)
{
  ## Work on the whole block as a single string, so that a macro and its closing brace can be
  ## matched across lines. Characters are replaced by spaces rather than deleted, so that the
  ## number of lines, and the position of everything else, is preserved.
  block_text <- paste(code_lines, collapse = "\n")
  block_chars <- strsplit(block_text, "", fixed = TRUE)[[1]]

  ## Rd wrapper macros. Their braces are matched by the Rd parser, which counts every brace,
  ## including those that sit inside an R comment. The macro name, its opening brace, and its
  ## matching closing brace are therefore blanked out, splicing the content back inline.
  repeat
  {
    macro_match <- regexpr("\\\\(dontrun|donttest|dontshow|testonly)[[:space:]]*\\{",
                           paste(block_chars, collapse = ""))
    if (macro_match == -1) { break }

    macro_start <- as.integer(macro_match)
    macro_stop <- macro_start + attr(macro_match, "match.length") - 1L   # position of the '{'

    ## Find the brace matching the opening one, counting every brace as the Rd parser does
    depth <- 0L
    closing_position <- NA_integer_
    for (k in macro_stop:length(block_chars))
    {
      if (block_chars[k] == "{") { depth <- depth + 1L }
      if (block_chars[k] == "}")
      {
        depth <- depth - 1L
        if (depth == 0L) { closing_position <- k ; break }
      }
    }

    ## Blank out the macro name and its opening brace
    block_chars[macro_start:macro_stop] <- " "
    ## Blank out the matching closing brace, when the block holds one
    if (!is.na(closing_position)) { block_chars[closing_position] <- " " }
  }

  ## Other Rd escapes that can legitimately appear inside examples
  code_lines <- strsplit(paste(block_chars, collapse = ""), "\n", fixed = TRUE)[[1]]
  code_lines <- gsub("\\\\%", "%", code_lines)
  code_lines <- gsub("\\\\dots", "...", code_lines)

  return(code_lines)
}


## Helper function to extract the @examples blocks of an R source file ####

#' @title Extract the @examples blocks of an R source file
#'
#' @description Locate every `@examples` block in a roxygen2-documented R file, and return its
#'   content alongside the line at which it starts, so that parse errors can be reported
#'   against the source file.
#'
#' @param file_path Character string. Path to the R source file.
#'
#' @return A list of lists, one per `@examples` block, each holding:
#'   * `$start_line` Integer. Line of the `@examples` tag in the source file.
#'   * `$code` Vector of character strings. Content of the block, with `#'` removed.
#'   * `$first_code_line` Integer. Line of the source file matching `$code[1]`.
#'
#' @author Maël Doré
#'
#' @noRd
#'

extract_examples_blocks <- function (file_path)
{
  source_lines <- readLines(file_path, warn = FALSE)
  nb_lines <- length(source_lines)
  is_roxygen <- grepl("^[[:space:]]*#'", source_lines)

  blocks <- list()
  i <- 1L
  while (i <= nb_lines)
  {
    ## '@examples' and '@examplesIf' both open a block of example code
    if (is_roxygen[i] & grepl("^[[:space:]]*#'[[:space:]]*@examples", source_lines[i]))
    {
      block_start <- i
      j <- i + 1L
      ## The block runs until the next roxygen tag, or the end of the roxygen comment
      while ((j <= nb_lines) & is_roxygen[j] & !grepl("^[[:space:]]*#'[[:space:]]*@[a-zA-Z]", source_lines[j]))
      {
        j <- j + 1L
      }

      if (j > (block_start + 1L))
      {
        blocks[[length(blocks) + 1L]] <- list(
          start_line = block_start,
          first_code_line = block_start + 1L,
          code = sub("^[[:space:]]*#'[[:space:]]?", "", source_lines[(block_start + 1L):(j - 1L)]))
      }
      i <- j
    } else {
      i <- i + 1L
    }
  }

  return(blocks)
}


## Main function to check that all @examples blocks parse ####

#' @title Check that all @examples blocks of a package parse as valid R code
#'
#' @description Parse every `@examples` block of a package, one at a time, and report those
#'   that hold a syntax error, pointing at the file and line where the error must be fixed.
#'
#'   This reproduces the parsing done by `R CMD check` in its step
#'   "checking for unstated dependencies in examples", but it runs in a fraction of a second,
#'   and it reports a usable location.
#'
#' @param path Character string. Path to the root of the package directory. Default = `"."`.
#' @param source Character string. Where to read the examples from.
#'   * `"R"` (default) reads the roxygen2 `@examples` blocks in `R/*.R`. This works before
#'     [devtools::document()] has been run, and reports the line in the R source file.
#'   * `"man"` reads the examples extracted from `man/*.Rd` with [tools::Rd2ex()]. This is
#'     exactly what `R CMD check` parses, but it requires `man/` to be up to date, and it
#'     reports the line in the `.Rd` file.
#' @param verbose Logical. Whether to list every block that was checked. Default = `FALSE`.
#'
#' @return Invisibly, a data.frame with one row per failing block, holding `$file`, `$line`,
#'   `$message` and `$offending_line`. Returns an empty data.frame when every block parses.
#'
#' @author Maël Doré
#'
#' @noRd
#'

check_examples_parsing <- function (path = ".", source = c("R", "man"), verbose = FALSE)
{
  source <- match.arg(source)
  failures <- data.frame(file = character(0), line = integer(0),
                         message = character(0), offending_line = character(0),
                         stringsAsFactors = FALSE)
  nb_blocks <- 0L

  if (source == "R")
  {
    files <- list.files(file.path(path, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
    if (length(files) == 0)
    {
      stop(paste0("No R source file found in '", file.path(path, "R"), "'.\n",
                  "Please run this function from the root of the package directory."))
    }

    for (focal_file in files)
    {
      for (block in extract_examples_blocks(focal_file))
      {
        nb_blocks <- nb_blocks + 1L
        code <- prepare_example_code(block$code)
        parse_error <- tryCatch({ parse(text = code) ; NULL },
                                error = function (e) { e })

        if (!is.null(parse_error))
        {
          ## Recover the line of the error within the block, from the "<text>:LINE:COL:" prefix
          message_text <- conditionMessage(parse_error)
          line_in_block <- suppressWarnings(as.integer(sub("^<text>:([0-9]+):.*$", "\\1",
                                                           strsplit(message_text, "\n")[[1]][1])))
          if (is.na(line_in_block)) { line_in_block <- 1L }
          source_line <- block$first_code_line + line_in_block - 1L

          failures <- rbind(failures, data.frame(
            file = focal_file, line = source_line,
            message = strsplit(message_text, "\n")[[1]][1],
            offending_line = trimws(block$code[min(line_in_block, length(block$code))]),
            stringsAsFactors = FALSE))
        }

        if (verbose)
        {
          cat(sprintf("  %-46s @examples line %5d  %s\n", basename(focal_file), block$start_line,
                      ifelse(is.null(parse_error), "OK", "PARSE ERROR")))
        }
      }
    }

  } else {

    files <- list.files(file.path(path, "man"), pattern = "\\.Rd$", full.names = TRUE)
    if (length(files) == 0)
    {
      stop(paste0("No .Rd file found in '", file.path(path, "man"), "'.\n",
                  "Please run devtools::document() first, or use source = 'R'."))
    }

    for (focal_file in files)
    {
      example_file <- tempfile(fileext = ".R")
      extracted <- tryCatch({ tools::Rd2ex(focal_file, out = example_file) ; TRUE },
                            error = function (e) { FALSE })
      if (!extracted | !file.exists(example_file)) { next }   # no \examples section

      nb_blocks <- nb_blocks + 1L
      parse_error <- tryCatch({ parse(file = example_file) ; NULL }, error = function (e) { e })

      if (!is.null(parse_error))
      {
        message_text <- strsplit(conditionMessage(parse_error), "\n")[[1]][1]
        failures <- rbind(failures, data.frame(
          file = focal_file, line = NA_integer_,
          message = message_text, offending_line = NA_character_,
          stringsAsFactors = FALSE))
      }

      if (verbose)
      {
        cat(sprintf("  %-46s  %s\n", basename(focal_file),
                    ifelse(is.null(parse_error), "OK", "PARSE ERROR")))
      }
      unlink(example_file)
    }
  }

  ## Report
  cat(sprintf("\nChecked %d @examples block(s) from '%s/'.\n", nb_blocks, source))
  if (nrow(failures) == 0)
  {
    cat("All example blocks parse correctly.\n\n")
  } else {
    cat(sprintf("%d block(s) failed to parse:\n\n", nrow(failures)))
    for (k in seq_len(nrow(failures)))
    {
      cat(sprintf("  %s:%s\n", failures$file[k],
                  ifelse(is.na(failures$line[k]), "", failures$line[k])))
      cat(sprintf("    %s\n", failures$message[k]))
      if (!is.na(failures$offending_line[k]))
      {
        cat(sprintf("    offending line: %s\n", failures$offending_line[k]))
      }
      cat("\n")
    }
    cat("A frequent cause is a line of prose inside @examples that lost its leading '#'.\n\n")
  }

  return(invisible(failures))
}
