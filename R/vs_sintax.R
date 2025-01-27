#' Taxonomic classification
#'
#' @param fasta_input A FASTA file path or a FASTA object with reads to
#' classify. See Details.
#' @param db A FASTA/FASTQ file path or object containing the reference database
#' in FASTQ/FASTA format. The sequences need to be annotated with taxonomy.
#' @param tabbedout The name of the tab-separated output file. See Details.
#' @param cutoff The minimum level of bootstrap support for the taxonomic ranks
#' that will be included in column 4 of the output file.
#' @param randseed The seed for the random number generator used in the Sintax
#' algorithm. Defaults to \code{NULL}.
#' @param strand \code{"plus"} (default) or \code{"both"}. When comparing
#' sequences only check the plus strand or both strands.
#' @param threads The number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param vsearch_options A character string of additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See Details.
#'
#' @details The sequences in the input file are classified according to the
#' Sintax algorithm, using \code{VSEARCH}.
#'
#' \code{fasta_input} can either be a FASTA file or object. FASTA objects are
#' tibbles that contain the columns \code{Header} and \code{Sequence}.
#'
#' \code{tabbedout} is a tab-separated text format file. Column 1 contains the
#' query label. Column 2 contains the predicted taxonomy in the same format as
#' for the reference data, with bootstrap support indicated in parentheses after
#' each rank. Column 3 contains the strand. If \code{cutoff} option is specified,
#' the predicted taxonomy will be repeated in column 4 while omitting the
#' bootstrap values and including only the ranks with support at or above the
#' threshold.
#'
#' \code{vsearch_options} can be used to pass additional arguments to
#' \code{VSEARCH}, that are not implemented in \code{Rsearch}. See the
#' \code{VSEARCH} manual for additional arguments, and how to use them.
#'
#' @returns A tibble or \code{NULL}.
#'
#' If \code{tabbedout} is unspecified, a tibble containing the results is
#' returned. If \code{tabbedout} is specified, results are written to file, and
#' nothing is returned.
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_sintax sintax classify
#'
#' @export
#'
vs_sintax <- function(fasta_input,
                      db,
                      tabbedout = NULL,
                      cutoff = NULL,
                      randseed = NULL,
                      threads = 1,
                      strand = "plus",
                      vsearch_options = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate strand
  if (!strand %in% c("plus", "both")) {
    stop("Invalid value for 'strand'. Choose from 'plus' or 'both'.")
  }

  # Create empty vector for collecting temporary files
  temp_files <- character()

  # Set up cleanup of temporary files
  on.exit({
    if (length(temp_files) > 0 && is.character(temp_files)) {
      existing_files <- temp_files[file.exists(temp_files)]
      if (length(existing_files) > 0) {
        file.remove(existing_files)
      }
    }
  }, add = TRUE)

  # Check if FASTA input is file or tibble
  if (!is.character(fasta_input)){
    temp_file <- tempfile(pattern = "input", fileext = ".fa")
    temp_files <- c(temp_files, temp_file)
    microseq::writeFasta(fasta_input, temp_file)
    fasta_file <- temp_file

    # Capture original name for statistics table later
    fasta_input_name <- as.character(substitute(fasta_input))
  } else {
    fasta_file <- fasta_input

    # Capture original name for statistics table later
    fasta_input_name <- basename(fasta_input)
  }

  # Check is input file exists at given path
  if (!file.exists(fasta_file)) stop("Cannot find input file: ", fasta_file)

  # Handle input data base sequences
  if (!is.character(db)){
    if ("Quality" %in% colnames(db)){

      # Validate tibble
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(db))) {
        stop("FASTQ data base must contain columns: Header, Sequence, Quality")
      }

      temp_file_db <- tempfile(pattern = "db_input", fileext = ".fq")
      temp_files <- c(temp_files, temp_file_db)
      microseq::writeFastq(db, temp_file_db)

      db_file <- temp_file_db

    } else {

      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(db))) {
        stop("FASTA data base must contain columns: Header and Sequence")
      }

      temp_file_db <- tempfile(pattern = "db_input", fileext = ".fa")
      temp_files <- c(temp_files, temp_file_db)
      microseq::writeFasta(db, temp_file_db)

      db_file <- temp_file_db

    }
  } else {
    if (!file.exists(db)) stop("Cannot find input file: ", db)

    db_file <- db
  }

  # Normalize file paths
  fasta_file <- normalizePath(fasta_file)
  db_file <- normalizePath(db_file)

  # Determine tabbedout file
  if (is.null(tabbedout)) {
    outfile <- tempfile(pattern = "output", fileext = ".txt")
    temp_files <- c(temp_files, outfile)
  } else {
    outfile <- tabbedout
  }

  # Build argument string for command line
  args <- c("--sintax", shQuote(fasta_file),
            "--db", shQuote(db_file),
            "--threads", threads,
            "--strand", strand,
            "--tabbedout", outfile)

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Add log file if specified
  if (!is.null(log_file)){
    args <- c(args, "--log", log_file)
  }

  # Add cutoff if specified
  if (!is.null(cutoff)){
    args <- c(args, "--sintax_cutoff", cutoff)
  }

  # Add random seed for sintax if specified
  if (!is.null(randseed)){
    args <- c(args, "--randseed", randseed)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  if (is.null(tabbedout)) {

    outfile_colnames <- c("query", "predicted_taxonomy", "strand")

    if(!is.null(cutoff)){
      outfile_colnames <- c(outfile_colnames, "predicted_taxonomy_cutoff")
    }

    # Read output into R tibble
    tabbedout.tbl <- utils::read.delim(outfile,
                                       header = FALSE,
                                       col.names = outfile_colnames)
  }

  # Return results
  if (is.null(tabbedout)) { # Return tibble
    return(tabbedout.tbl)
  } else {
    return(invisible(NULL)) # No return when output file is written
  }
}
