#' Detect chimeras without external references (i.e. de novo)
#'
#' @description \code{vs_uchime_denovo} detects chimeras present in the FASTA
#' sequences in using \code{VSEARCH}'s \code{uchime_denovo} algorithm.
#' Automatically sorts sequences by decreasing abundance to enhance chimera
#' detection accuracy.
#'
#' @param fasta_input (Required). A FASTA file path or a FASTA object with reads.
#' If a tibble is provided, any columns in addition to \code{Header} and
#' \code{Sequence} will be preserved in the output. See \emph{Details}.
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
#' @param vsearch_options (Optional). Additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See \emph{Details}.
#' @param tmpdir (Optional). Path to the directory where temporary files should
#' be written when tables are used as input or output. Defaults to
#' \code{NULL}, which resolves to the session-specific temporary directory
#' (\code{tempdir()}).
#'
#' @details
#' Chimeras in the input FASTA sequences are detected using \code{VSEARCH}´s
#' \code{uchime_denovo}. In de novo mode, input FASTA file/object must present
#' abundance annotations (i.e. a pattern [;]size=integer[;] in the header).
#' Input order matters for chimera detection, so it is recommended to sort
#' sequences by decreasing abundance.
#'
#' \code{fasta_input} can either be a FASTA file or a FASTA object. FASTA objects
#' are tibbles that contain the columns \code{Header} and \code{Sequence}, see
#' \code{\link[microseq]{readFasta}}.
#'
#' When providing a tibble as \code{fasta_input}, you can include additional
#' columns with metadata (e.g., OTU IDs, sample origins). The function will
#' preserve these columns by joining them back to the results based on the
#' DNA sequence. This allows you to keep your metadata associated with your
#' sequences throughout the chimera detection process.
#'
#' If \code{nonchimeras} and \code{chimeras} are specified, resulting
#' non-chimeric and chimeric sequences are written to these files in FASTA
#' format.
#'
#' If \code{nonchimeras} and \code{chimeras} are \code{NULL}, results are
#' returned as a FASTA-objects.
#'
#' \code{nonchimeras} and \code{chimeras} must either both be specified or both
#' be \code{NULL}.
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
#' If \code{nonchimeras} and \code{chimeras} are \code{NULL}, a FASTA object
#' containing non-chimeric sequences is returned. This output tibble will
#' include any additional columns that were present in the \code{fasta_input}
#' tibble. An attribute named \code{"chimeras"} will contain a tibble of the
#' chimeric sequences, also with the additional columns preserved.
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
#' fasta_input <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                          "small_R1.fq")
#' nonchimeras <- "nonchimeras.fa"
#' chimeras <- "chimeras.fa"
#'
#' # Detect chimeras with default parameters and return FASTA files
#' vs_uchime_denovo(fasta_input = fasta_input,
#'                  nonchimeras = nonchimeras,
#'                  chimeras = chimeras)
#'
#' # Detect chimeras with default parameters and return a FASTA tibble
#' nonchimeras.tbl <- vs_uchime_denovo(fasta_input = fasta_input,
#'                                     nonchimeras = NULL,
#'                                     chimeras = NULL)
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
#' @aliases vs_uchime_denovo uchime_denovo chimera
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_uchime_denovo <- function(fasta_input,
                             nonchimeras = NULL,
                             chimeras = NULL,
                             sizein = TRUE,
                             sizeout = TRUE,
                             relabel = NULL,
                             relabel_sha1 = FALSE,
                             fasta_width = 0,
                             sample = NULL,
                             log_file = NULL,
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
    fasta_input_vsearch <- dplyr::select(fasta_input, Header, Sequence)
    microseq::writeFasta(fasta_input_vsearch, temp_file)
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

  # Normalize file path
  fasta_file <- normalizePath(fasta_file)

  # Build argument string for command line
  args <- c("--uchime_denovo", shQuote(fasta_file),
            "--fasta_width", fasta_width,
            "--nonchimeras", shQuote(nonchimeras_file),
            "--chimeras", shQuote(chimeras_file)
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

  # Run VSEARCH
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Check for VSEARCH failure
  check_vsearch_status(vsearch_output, args)

  if (is.null(nonchimeras) && is.null(chimeras)) {

    # Read output into FASTA object (tbl)
    nonchimeras.tbl <- microseq::readFasta(nonchimeras_file)|>
      dplyr::mutate(Sequence = toupper(Sequence))

    # Join with input table if possible
    if (!is.character(fasta_input)){
        fasta_input_join <- fasta_input |>
          dplyr::select(-Header)

        nonchimeras.tbl <- dplyr::left_join(nonchimeras.tbl,
                                            fasta_input_join,
                                            by = "Sequence")
    }

    # Create empty table
    chimeras.tbl <- data.frame()

    # Check if chimeras file contains something
    if (file.info(chimeras_file)$size > 0){

      chimeras.tbl <- microseq::readFasta(chimeras_file) |>
        dplyr::mutate(Sequence = toupper(Sequence))

      # Join with input table if possible
      if (!is.character(fasta_input)){
          chimeras.tbl <- dplyr::left_join(chimeras.tbl,
                                           fasta_input_join,
                                           by = "Sequence")
        }
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

#' Calculate chimera detection statistics
#'
#' @description Calculates important chimera detection statistics after running
#' \code{vs_uchime_denovo()} or \code{vs_uchime_ref()}, including the number of chimeric and non-chimeric
#' sequences.
#'
#' @param fasta_file The FASTA file containing the input sequences to the
#' chimera detection.
#' @param fasta_input_name The name of the file/object with the input sequences
#' that was used in the chimera detection.
#' @param nonchimeras.tbl The output tibble from chimera detection with the
#' non-cimeric sequences. Contains the columns: Header and Sequence.
#' @param chimeras.tbl The output tibble from chimera detection with the
#' chimeric sequences. Contains the columns: Header and Sequence. If the table
#' is \code{NULL}, it means that no chimeras were found.
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item \code{num_nucleotides}: The total number of nucleotides used as input
#'   for chimera detection.
#'   \item \code{num_sequences}: The total number of sequences used as input for
#'   chimera detection.
#'   \item \code{min_length_input_seq}: The length of the shortest sequence used
#'   as input for chimera detection.
#'   \item \code{max_length_input_seq}: The length of the longest sequence used
#'   as input for chimera detection.
#'   \item \code{avg_length_input_seq}: The average length of the sequences used
#'   as input for chimera detection.
#'   \item \code{num_non_chimeras}: The number of non-chimeric sequences.
#'   \item \code{num_chimeras}: The number of chimeric sequences.
#'   \item \code{input}: The name of the input file/object for the chimera
#'   detection.
#' }
#'
#' @return A tibble with chimera detection statistics.
#'
#' @noRd
#'
calculate_uchime_statistics <- function(fasta_file,
                                        fasta_input_name,
                                        nonchimeras.tbl,
                                        chimeras.tbl) {

  # Make tibble from input sequences to the clustering
  input.df <- microseq::readFasta(fasta_file)

  # Calculate statistics
  num_nucleotides <- sum(nchar(input.df$Sequence))
  num_sequences <- nrow(input.df)
  min_length_input_seq <- min(nchar(input.df$Sequence))
  max_length_input_seq <- max(nchar(input.df$Sequence))
  avg_length_input_seq <- mean(nchar(input.df$Sequence))
  num_non_chimeras <- nrow(nonchimeras.tbl)

  num_chimeras <- nrow(chimeras.tbl)

  # Create table
  result_table <- data.frame(
    num_nucleotides = num_nucleotides,
    num_sequences = num_sequences,
    min_length_input_seq = min_length_input_seq,
    max_length_input_seq = max_length_input_seq,
    avg_length_input_seq = avg_length_input_seq,
    num_non_chimeras = num_non_chimeras,
    num_chimeras = num_chimeras,
    input = fasta_input_name
  )

  return(result_table)
}
