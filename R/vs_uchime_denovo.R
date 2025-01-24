#' Detect chimeras
#'
#' @description Detects chimeras present in the FASTA sequences in the given
#' file or object.
#'
#' @param fasta_input A FASTA file path or a FASTA object with reads. See
#' Details.
#' @param nonchimeras Name of the FASTA output file for the non-chimeric
#' sequences. If \code{NULL} (default) no output will be written to file.
#' @param chimeras Name of the FASTA output file for the chimeric sequences.
#' If \code{NULL} (default) no output will be written to file.
#' @param sizein Decides if abundance annotations present in sequence headers
#' should be taken into account. Defaults to \code{TRUE}.
#' @param sizeout Decides if abundance annotations should be added to
#' FASTA headers. Defaults to \code{TRUE}.
#' @param relabel Relabel sequences using the given prefix and a ticker to
#' construct new headers. Defaults to \code{NULL}.
#' @param relabel_sha1 Relabel sequences using the SHA1 message digest
#' algorithm. Defaults to \code{FALSE}.
#' @param fasta_width The number of characters in the width of sequences in the
#' output FASTA file. Defaults to \code{0}. See Details.
#' @param log_file Name of the log file to capture messages from \code{VSEARCH}.
#' If \code{NULL}, no log file is created. Defaults to \code{NULL}.
#' @param vsearch_options A character string of additional arguments to pass to
#' \code{VSEARCH}. Defaults to \code{NULL}. See Details.
#'
#' @details Detects chimeras present in the FASTA-formated input, without
#' external references (i.e. de novo). Automatically sort the sequences in the
#' input by decreasing abundance beforehand.
#'
#' Chimeras in the input FASTA sequences are detected using \code{VSEARCH}´s
#' \code{uchime_denovo}. In de novo mode, input FASTA file/object must present
#' abundance annotations (i.e. a pattern [;]size=integer[;] in the header).
#' Input order matters for chimera detection, so it is recommended to sort
#' sequences by decreasing abundance.
#'
#' \code{fasta_input} can either be a FASTA file or object. FASTA objects are
#' tibbles that contain the columns \code{Header} and \code{Sequence}.
#'
#' If \code{nonchimeras} and \code{chimeras} are specified, the resulting
#' sequences after chimera detection are output to these files in FASTA format.
#' If unspecified (\code{NULL}) the results are returned as a FASTA-objects.
#' \code{nonchimeras} and \code{chimeras} must either both be specified, or both
#' unspecified.
#'
#' FASTA files produced by \code{VSEARCH} are wrapped
#' (sequences are written on lines of integer nucleotides).\code{fasta_width} is
#' by default set to zero to eliminate the wrapping.
#'
#' \code{vsearch_options} can be used to pass additional arguments to \code{VSEARCH},
#' that are not implemented in \code{Rsearch}. See the \code{VSEARCH} manual for
#' additional arguments, and how to use them.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{nonchimeras} and \code{chimeras} are specified, the resulting
#' sequences after chimera detection are output to these files in FASTA format,
#' and nothing is returned.
#' If unspecified (\code{NULL}) the results are returned as a FASTA-objects.
#' The tibble containing chimeric sequences is an attribute, called
#' \code{"chimera"}, to the primary table with non-chimeric sequences. If no
#' chimeric sequences are found, the attribute is empty (\code{NULL}).
#'
#' When a FASTA object is returned, the statistics from the chimera detection,
#' \code{statistics}, is an attribute of the non-chimeras tibble.
#' The statistics tibble has the following columns:
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
                             log_file = NULL,
                             vsearch_options = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

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

  # Determine nonchimeras file
  if (is.null(nonchimeras)) {
    nonchimeras_file <- tempfile(pattern = "nonchimeras", fileext = ".fa")
    temp_files <- c(temp_files, nonchimeras_file)
  } else {
    nonchimeras_file <- normalizePath(nonchimeras)
  }

  # Determine chimeras file
  if (is.null(chimeras)) {
    chimeras_file <- tempfile(pattern = "chimeras", fileext = ".fa")
    temp_files <- c(temp_files, chimeras_file)
  } else {
    chimeras_file <- normalizePath(chimeras)
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

  if (is.null(nonchimeras) && is.null(chimeras)) {

    # Read output into FASTA object (tbl)
    nonchimeras.tbl <- microseq::readFasta(nonchimeras_file)

    # Check if chimeras file contains something
    if (file.info(chimeras_file)$size > 0){

      chimeras.tbl <- microseq::readFasta(chimeras_file)

      # Add additional table as attribute to the primary table
      attr(nonchimeras.tbl, "chimeras") <- chimeras.tbl

    } else {
      # Add empty attribute to the primary table because no chimeras were found
      attr(nonchimeras.tbl, "chimeras") <- NULL
    }

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
#' \code{vs_uchime_denovo()}, including the number of chimeric and non-chimeric
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

  if (!is.null(chimeras.tbl)) {
    num_chimeras <- nrow(chimeras.tbl)
  } else {
    num_chimeras <- 0
  }

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
