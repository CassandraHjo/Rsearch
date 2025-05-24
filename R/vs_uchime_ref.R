#' Detect chimeras by comparing sequences to a reference database
#'
#' @description \code{vs_uchime_ref} detects chimeras present in the FASTA
#' sequences in using \code{VSEARCH}'s \code{uchime_ref} algorithm.
#'
#' @param fasta_input (Required). A FASTA file path or a FASTA object with reads.
#' See \emph{Details}.
#' @param database (Required). A FASTA file path or FASTA tibble object
#' containing the reference sequences. These sequences are assumed to be
#' chimera-free.
#' @param nonchimeras (Optional). Name of the FASTA output file for the
#' non-chimeric sequences. If \code{NULL} (default), no output is written to
#' file.
#' @param chimeras (Optional). Name of the FASTA output file for the chimeric
#' sequences. If \code{NULL} (default), no output is written to file.
#' @param sizein (Optional). If \code{TRUE} (default), abundance annotations
#' present in sequence headers are taken into account.
#' @param sizeout (Optional). If \code{TRUE} (default), abundance annotations
#' are added to FASTA headers.
#' @param relabel (Optional). Relabel sequences using the given prefix and a
#' ticker to construct new headers. Defaults to \code{NULL}.
#' @param relabel_sha1 (Optional). If \code{TRUE} (default), relabel sequences
#' using the SHA1 message digest algorithm. Defaults to \code{FALSE}.
#' @param fasta_width (Optional). Number of characters per line in the output
#' FASTA file. Defaults to \code{0}, which eliminates wrapping.
#' @param sample (Optional). Add the given sample identifier string to sequence
#' headers. For instance, if the given string is "ABC", the text ";sample=ABC"
#' will be added to the header. If \code{NULL} (default), no identifier is added.
#' @param log_file (Optional). Name of the log file to capture messages from
#' \code{VSEARCH}. If \code{NULL} (default), no log file is created.
#' @param threads (Optional). Number of computational threads to be used by
#' \code{VSEARCH}. Defaults to \code{1}.
#' @param vsearch_options (Optional). Additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See \emph{Details}.
#' @param tmpdir (Optional). Path to the directory where temporary files should
#' be written when tables are used as input or output. Defaults to
#' \code{NULL}, which resolves to the session-specific temporary directory
#' (\code{tempdir()}).
#'
#' @details
#' Chimeras in the input FASTA sequences are detected using \code{VSEARCH}´s
#' \code{uchime_ref}.
#'
#' \code{fasta_input} can either be a FASTA file or a FASTA object. FASTA objects
#' are tibbles that contain the columns \code{Header} and \code{Sequence}, see
#' \code{\link{readFasta}}.
#'
#' \code{database} must be a FASTA file or a FASTA object with high-quality
#' non-chimeric sequences.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{nonchimeras} and \code{chimeras} are specified, the resulting
#' sequences after chimera detection written directly to the specified files in
#' FASTA format, and no tibbles are returned.
#'
#' If \code{nonchimeras} and \code{chimeras} are \code{NULL}, A FASTA object
#' containing non-chimeric sequences with an attribute \code{"chimeras"}
#' containing a tibble of chimeric sequences is returned. If no chimeras are
#' found, the \code{"chimeras"} attribute is an empty data frame.
#'
#' Additionally, the returned tibble (when applicable) has an attribute
#' \code{"statistics"} containing a tibble with chimera detection statistics.
#'
#' The statistics tibble has the following columns:
#' \itemize{
#'   \item \code{num_nucleotides}: Total number of nucleotides used as input
#'   for chimera detection.
#'   \item \code{num_sequences}: Total number of sequences used as input for
#'   chimera detection.
#'   \item \code{min_length_input_seq}: Length of the shortest sequence used
#'   as input for chimera detection.
#'   \item \code{max_length_input_seq}: Length of the longest sequence used as
#'   input for chimera detection.
#'   \item \code{avg_length_input_seq}: Average length of the sequences used as
#'   input for chimera detection.
#'   \item \code{num_non_chimeras}: Number of non-chimeric sequences.
#'   \item \code{num_chimeras}: Number of chimeric sequences.
#'   \item \code{input}: Name of the input file/object for the chimera
#'   detection.
#' }
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' query_file <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                         "small.fasta")
#' db <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                 "sintax_db.fasta")
#'
#' # Detect chimeras with default parameters and return FASTA files
#' vs_uchime_ref(fasta_input = query_file,
#'               database = db,
#'               nonchimeras = "nonchimeras.fa",
#'               chimeras = "chimeras.fa")
#'
#' # Detect chimeras with default parameters and return a FASTA tibble
#' nonchimeras.tbl <- vs_uchime_ref(fasta_input = query_file,
#'                                  database = db,
#'                                  nonchimeras = NULL,
#'                                  chimeras = NULL)
#'
#' # Get chimeras tibble
#' chimeras.tbl <- attr(nonchimeras.tbl, "chimeras")
#'
#' # Get statistics tibble
#' statistics.tbl <- attr(nonchimeras.tbl, "statistics")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_uchime_ref uchime_ref
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_uchime_ref <- function(fasta_input,
                          database,
                          nonchimeras = NULL,
                          chimeras = NULL,
                          sizein = TRUE,
                          sizeout = TRUE,
                          relabel = NULL,
                          relabel_sha1 = FALSE,
                          fasta_width = 0,
                          sample = NULL,
                          log_file = NULL,
                          threads = 1,
                          vsearch_options = NULL,
                          tmpdir = NULL) {

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Set temporary directory if not provided
  if (is.null(tmpdir)) tmpdir <- tempdir()

  # Check if both output files are specified, or both unspecified
  if (is.null(nonchimeras) != is.null(chimeras)) {
    stop("nonchimeras and chimeras must either both be specified or both unspecified.")
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
    temp_file <- tempfile(pattern = "input",
                          tmpdir = tmpdir,
                          fileext = ".fa")
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

  # Check if database is file or tibble
  if (!is.character(database)){
    temp_file_db <- tempfile(pattern = "db",
                             tmpdir = tmpdir,
                             fileext = ".fa")
    temp_files <- c(temp_files, temp_file_db)
    microseq::writeFasta(database, temp_file_db)
    db_file <- temp_file_db

  } else {
    db_file <- database
  }

  # Determine nonchimeras file
  if (is.null(nonchimeras)) {
    nonchimeras_file <- tempfile(pattern = "nonchimeras",
                                 tmpdir = tmpdir,
                                 fileext = ".fa")
    temp_files <- c(temp_files, nonchimeras_file)
  } else {
    nonchimeras_file <- nonchimeras
  }

  # Determine chimeras file
  if (is.null(chimeras)) {
    chimeras_file <- tempfile(pattern = "chimeras",
                              tmpdir = tmpdir,
                              fileext = ".fa")
    temp_files <- c(temp_files, chimeras_file)
  } else {
    chimeras_file <- chimeras
  }

  # Normalize file paths
  fasta_file <- normalizePath(fasta_file)
  db_file <- normalizePath(db_file)

  # Build argument string for command line
  args <- c("--uchime_ref", shQuote(fasta_file),
            "--db", shQuote(db_file),
            "--fasta_width", fasta_width,
            "--nonchimeras", shQuote(nonchimeras_file),
            "--chimeras", shQuote(chimeras_file),
            "--threads", threads
  )

  if (sizein) {
    args <- c(args, "--sizein", "")
  }

  if (sizeout) {
    args <- c(args, "--sizeout", "")
  }

  # Add relabeling arguments if specified
  if (!is.null(relabel)){
    args <- c(args, "--relabel", relabel)
  }

  if (relabel_sha1){
    args <- c(args, "--relabel_sha1", "")
  }

  # Add sample identifier if specified
  if (!is.null(sample)) {
    args <- c(args, "--sample", sample)
  }

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Add log file if specified
  if (!is.null(log_file)){
    args <- c(args, "--log", log_file)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Check for VSEARCH failure
  check_vsearch_status(vsearch_output, args)

  # Handle output
  if (is.null(nonchimeras) && is.null(chimeras)) {

    # Read output into FASTA object (tbl)
    nonchimeras.tbl <- microseq::readFasta(nonchimeras_file)

    # Create empty data frame
    chimeras.tbl <- data.frame()

    # Check if chimeras file contains something
    if (file.info(chimeras_file)$size > 0){

      chimeras.tbl <- microseq::readFasta(chimeras_file)
    }

    # Add additional table as attribute to the primary table
    attr(nonchimeras.tbl, "chimeras") <- chimeras.tbl

    # Add statistics
    statistics <- calculate_uchime_statistics(fasta_file,
                                              fasta_input_name,
                                              nonchimeras.tbl,
                                              attr(nonchimeras.tbl, "chimeras"))

    attr(nonchimeras.tbl, "statistics") <- statistics

  }

  # Return results
  if (is.null(nonchimeras) && is.null(chimeras)) { # Return tibble
    return(nonchimeras.tbl)
  } else {
    return(invisible(NULL)) # No return when output file is written
  }
}
