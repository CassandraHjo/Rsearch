#' Trim and/or filter sequences in FASTA/FASTQ format
#'
#' @description Trim and/or filters biological sequences (DNA) in the FASTA/FASTQ format.
#'
#' @param fastx_input A FASTA/FASTQ file path or object containing (forward) reads. See Details.
#' @param reverse An optional FASTA/FASTQ file path or object, if the input consists of paired sequences, containing reverse reads. If provided, it will be processed alongside \code{fastx_input}. Defaults to \code{NULL}. See Details.
#' @param file_format Format of input files \code{fastx_input} (and \code{reverse}), and desired output format for file/tibble: \code{"fasta"} or \code{"fastq"} (default).
#' @param maxee_rate Threshold for average expected error. Numeric value ranging form \code{0.0} to \code{1.0}. Defaults to \code{NULL}. See Details.
#' @param minlen The minimum number of bases a sequence must have to be retained. Defaults to \code{0}. See Details.
#' @param maxlen The maximum number of bases for a given sequence. Sequences with more bases than the specified number are discarded. Defaults to \code{NULL}.
#' @param maxns The maximum number of N's for a given sequence. Sequences with more N's than the specified number are discarded. Defaults to \code{NULL}.
#' @param maxsize The maximum abundance for a given sequence. Sequences with abundance higher than the specified value are discarded. Defaults to \code{NULL}.
#' @param minsize The minimum abundance for a given sequence. Sequences with abundance lower than the specified value are discarded. Defaults to \code{NULL}.
#' @param trunclen The sequences are truncated to the specified length. Shorter sequences are discarded. Defaults to \code{NULL}.
#' @param truncqual The sequences are truncated starting from the first base with the specified base quality score value or lower. Defaults to \code{NULL}.
#' @param truncee The sequences are truncated so that their total expected error is not higher than the specified value. Defaults to \code{NULL}.
#' @param fastaout Name of the FASTA output file for the sequences given in \code{fastx_input}. If \code{NULL} no FASTA sequences will be written to file. Defaults to \code{NULL}. See Details.
#' @param fastqout Name of the FASTQ output file for the sequences given in \code{fastx_input}. If \code{NULL} no FASTQ sequences will be written to file. Defaults to \code{NULL}. See Details.
#' @param fastaout_rev Name of the FASTA output file for the sequences given in \code{reverse}. If \code{NULL} no FASTA sequences will be written to file. Defaults to \code{NULL}. See Details.
#' @param fastqout_rev Name of the FASTQ output file for the sequences given in \code{reverse}. If \code{NULL} no FASTQ sequences will be written to file. Defaults to \code{NULL}. See Details.
#' @param fasta_width Number of characters per line in the output FASTA file. Only applies if the output file is in FASTA format. Defaults to \code{0}. See Details.
#' @param log_file Name of the log file to capture messages from \code{vsearch}. If \code{NULL}, no log file is created. Defaults to \code{NULL}.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details The function trims and/or filters sequences from the input FASTA/FASTQ file or object based on the specified options, using \code{vsearch}.
#' If a \code{reverse} input is provided, these reads are processed in the same way. The the format of the output is determined by \code{file_format}.
#'
#' \code{fastx_input} and \code{reverse} can either be FASTA/FASTQ files or objects. FASTA objects are tibbles that contain the columns \code{Header} and \code{Sequence}. FASTQ objects are tibbles that contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#' \code{reverse} is an optional argument to the function. If provided, it will be processed alongside \code{fastx_input}, meaning the same specified options will be used for both inputs.
#'
#' If \code{fastaout} and \code{fastaout_rev} or \code{fastqout} and \code{fastqout_rev} are specified, the remaining sequences after trimming and/or filtering are output to these files in either FASTA or FASTQ format.
#' If unspecified (\code{NULL}), results are returned as a tibble, and no output is written to file. \code{file_format} has to match the desired output files/objects.
#'
#' Sequences with an average expected error greater than the specified \code{maxee_rate} are discarded.
#' For a given sequence, the average expected error is the sum of error probabilities for all the positions in the sequence, divided by the length of the sequence.
#'
#' Any input sequence with fewer bases than the value set in \code{minlen} will be discarded. By default, \code{minlen} is set to 0, which means that no sequences are removed.
#' However, using the default value may allow empty sequences to remain in the results.
#'
#' FASTA files produced by \code{vsearch} are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#'
#' @return Tibble or \code{NULL}.
#'
#' If output files are not specified, a tibble containing the filtered reads from \code{fastx_input} in the format specified by \code{file_format} is returned. If output files are specified, results are written to file and nothing is returned.
#'
#' If \code{reverse} is specified, the resulting tibble containing the trimmed and/or filtered reverse reads in the format specified by \code{file_format} is an attribute, called \code{"filt_reverse"}, to the primary table.
#'
#' When a FASTA/FASTQ object is returned, the statistics from the filtering, \code{statistics}, is an attribute, called \code{"statistics"}, to the primary filtering tibble.
#' This tibble contains statistics, including number of kept and discarded sequences, and the names of the FASTA/FASTQ files or objects that were used as input.
#'
#' @examples
#' \dontrun{
#' # Read example FASTQ files
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fq")
#' reverse <- file.path(file.path(path.package("Rsearch"), "extdata"), "R2_sample1_small.fq")
#'
#' # Define other arguments
#' file_format <- "fastq"
#' maxee_rate <- 0.01
#' minlen <- 0
#'
#' # Execute filtering, with tibble as output
#' filt_seqs <- vs_fastx_trim_filt(fastx_input = fastx_input,
#'                                 reverse = reverse,
#'                                 file_format = file_format,
#'                                 maxee_rate = maxee_rate,
#'                                 minlen = minlen)
#'
#' # Extract tibbles with filtered sequences
#' R1_filt <- filt_seqs
#' R2_filt <- attr(filt_seqs, "filt_reverse")
#'
#' # Extract filtering statistics
#' statistics <- attr(filt_seqs, "statistics")
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastx_trim_filt <- function(fastx_input,
                               reverse = NULL,
                               file_format = "fastq",
                               maxee_rate = NULL,
                               minlen = 0,
                               maxlen = NULL,
                               maxns = NULL,
                               maxsize = NULL,
                               minsize = NULL,
                               trunclen = NULL,
                               truncqual = NULL,
                               truncee = NULL,
                               fastaout = NULL,
                               fastqout = NULL,
                               fastaout_rev = NULL,
                               fastqout_rev = NULL,
                               fasta_width = 0,
                               log_file = NULL,
                               threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate file_format
  if (!file_format %in% c("fasta", "fastq")) {
    stop("Invalid file_format. Choose from fasta or fastq.")
  }

  # If file_format is "fasta", fastqout and fastqout_rev can not be defined
  if (file_format == "fasta") {
    if (!is.null(fastqout) || !is.null(fastqout_rev)) {
      stop("When file_format is defined as 'fasta', 'fastqout' and 'fastqout_rev' cannot be used. Use 'fastaout' and 'fastaout_rev' instead.")
    }
  }

  # If file_format is "fastq", fastaout and fastaout_rev can not be defined
  if (file_format == "fastq") {
    if (!is.null(fastaout) || !is.null(fastaout_rev)) {
      stop("When file_format is defined as 'fastq', 'fastaout' and 'fastaout_rev' cannot be used. Use 'fastqout' and 'fastqout_rev' instead.")
    }
  }

  # If reverse is specified, ensure paired output parameters are both NULL or both character strings
  if (!is.null(reverse)) {
    if (file_format == "fasta") {
      # Check that both fastaout and fastaout_rev are NULL or both are character strings
      if ((is.null(fastaout) && !is.null(fastaout_rev)) ||
          (!is.null(fastaout) && is.null(fastaout_rev))) {
        stop("When 'reverse' is specified and file_format is 'fasta', both 'fastaout' and 'fastaout_rev' must be NULL or both specified as character strings.")
      }
    }

    if (file_format == "fastq") {
      # Check that both fastqout and fastqout_rev are NULL or both are character strings
      if ((is.null(fastqout) && !is.null(fastqout_rev)) ||
          (!is.null(fastqout) && is.null(fastqout_rev))) {
        stop("When 'reverse' is specified and file_format is 'fastq', both 'fastqout' and 'fastqout_rev' must be NULL or both specified as character strings.")
      }
    }
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

  # Handle input when file_format = "fasta"
  if (file_format == "fasta") {

    # Handle input for primary sequences: file or tibble
    if (!is.character(fastx_input)){
      # Ensure required columns exist
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }
      temp_file_primary <- tempfile(pattern = "primary_input", fileext = ".fa")
      microseq::writeFasta(fastx_input, temp_file_primary)
      temp_files <- c(temp_files, temp_file_primary)

      fastx_file <- temp_file_primary

      # Capture original name for statistics table later
      fastx_input_name <- as.character(substitute(fastx_input))

    } else {
      if (!file.exists(fastx_input)) stop("Cannot find input FASTA file: ", fastx_input)

      fastx_file <- fastx_input

      # Capture original name for statistics table later
      fastx_input_name <- basename(fastx_input)
    }

    # Handle input for reverse sequences: file or tibble
    if (!is.null(reverse)){
      if (!is.character(reverse)){
        # Ensure required columns exist
        required_cols_rev <- c("Header", "Sequence")
        if (!all(required_cols_rev %in% colnames(reverse))) {
          stop("Reverse FASTA object must contain columns: Header and Sequence")
        }
        temp_reverse_file <- tempfile(pattern = "reverse_temp_", fileext = ".fa")
        microseq::writeFasta(reverse, temp_reverse_file)
        temp_files <- c(temp_files, temp_reverse_file)

        reverse_file <- temp_reverse_file

        # Capture original name for statistics table later
        reverse_name <- as.character(substitute(reverse))

      } else {
        if (!file.exists(reverse)) stop("Cannot find reverse FASTA file: ", reverse)

        reverse_file <- reverse

        # Capture original name for statistics table later
        reverse_name <- basename(reverse)
      }
    }
  }

  # Handle input when file_format = "fastq"
  if (file_format == "fastq") {

    # Handle input for primary sequences: file or tibble
    if (!is.character(fastx_input)){
      # Ensure required columns exist
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }
      temp_file_primary <- tempfile(pattern = "primary_input", fileext = ".fq")
      microseq::writeFastq(fastx_input, temp_file_primary)
      temp_files <- c(temp_files, temp_file_primary)

      fastx_file <- temp_file_primary

      # Capture original name for statistics table later
      fastx_input_name <- as.character(substitute(fastx_input))

    } else {
      if (!file.exists(fastx_input)) stop("Cannot find input FASTQ file: ", fastx_input)

      fastx_file <- fastx_input

      # Capture original name for statistics table later
      fastx_input_name <- basename(fastx_input)
    }

    # Handle input for reverse sequences: file or tibble
    if (!is.null(reverse)){
      if (!is.character(reverse)){
        # Ensure required columns exist
        required_cols_rev <- c("Header", "Sequence", "Quality")
        if (!all(required_cols_rev %in% colnames(reverse))) {
          stop("Reverse FASTQ object must contain columns: Header, Sequence, Quality")
        }
        temp_reverse_file <- tempfile(pattern = "reverse_temp_", fileext = ".fq")
        microseq::writeFastq(reverse, temp_reverse_file)
        temp_files <- c(temp_files, temp_reverse_file)

        reverse_file <- temp_reverse_file

        # Capture original name for statistics table later
        reverse_name <- as.character(substitute(reverse))

      } else {
        if (!file.exists(reverse)) stop("Cannot find reverse FASTQ file: ", reverse)

        reverse_file <- reverse

        # Capture original name for statistics table later
        reverse_name <- basename(reverse)
      }
    }
  }

  # Handle output for primary sequences
  if (file_format == "fasta") {
    if (is.null(fastaout)) {
      outfile_fasta <- tempfile(pattern = "filtered_primary_", fileext = ".fa")
      temp_files <- c(temp_files, outfile_fasta)
    } else {
      outfile_fasta <- fastaout
    }
  }

  if (file_format == "fastq") {
    if (is.null(fastqout)) {
      outfile_fastq <- tempfile(pattern = "filtered_primary_", fileext = ".fq")
      temp_files <- c(temp_files, outfile_fastq)
    } else {
      outfile_fastq <- fastqout
    }
  }

  # Handle output for reverse sequences
  if (!is.null(reverse)) {
    if (file_format == "fasta") {
      if (is.null(fastaout_rev)) {
        outfile_fasta_rev <- tempfile(pattern = "filtered_reverse_", fileext = ".fa")
        temp_files <- c(temp_files, outfile_fasta_rev)
      } else {
        outfile_fasta_rev <- fastaout_rev
      }
    }

    if (file_format == "fastq") {
      if (is.null(fastqout_rev)) {
        outfile_fastq_rev <- tempfile(pattern = "filtered_reverse_", fileext = ".fq")
        temp_files <- c(temp_files, outfile_fastq_rev)
      } else {
        outfile_fastq_rev <- fastqout_rev
      }
    }
  }

  # Normalize file paths
  fastx_file <- normalizePath(fastx_file)
  if (!is.null(reverse)) {
    reverse_file <- normalizePath(reverse_file)
  }

  # Build argument string for command line
  args <- c("--fastx_filter", fastx_file,
            "--fastq_minlen", minlen,
            "--threads", threads)

  # Add reverse to arguments if provided
  if (!is.null(reverse)) {
    args <- c(args, "--reverse", reverse_file)
  }

  # Add trimming and filtering arguments if provided
  if (!is.null(maxee_rate)) {
    args <- c(args, "--fastq_maxee_rate", maxee_rate)
  }

  if (!is.null(maxlen)) {
    args <- c(args, "--fastq_maxlen", maxlen)
  }

  if (!is.null(maxns)) {
    args <- c(args, "--fastq_maxns", maxns)
  }

  if (!is.null(maxsize)) {
    args <- c(args, "--maxsize", maxsize)
  }

  if (!is.null(minsize)) {
    args <- c(args, "--minsize", minsize)
  }

  if (!is.null(trunclen)) {
    args <- c(args, "--fastq_trunclen", trunclen)
  }

  if (!is.null(truncqual)) {
    args <- c(args, "--fastq_truncqual", truncqual)
  }

  if (!is.null(truncee)) {
    args <- c(args, "--fastq_truncee", truncee)
  }

  # Add output files based on file_format
  if (file_format == "fastq") {
    args <- c(args, "--fastqout", outfile_fastq)
    if (!is.null(reverse)) {
      args <- c(args, "--fastqout_rev", outfile_fastq_rev)
    }
  } else if(file_format == "fasta") {
    args <- c(args, "--fastaout", outfile_fasta, "--fasta_width", fasta_width)
    if (!is.null(reverse)) {
      args <- c(args, "--fastaout_rev", outfile_fasta_rev)
    }
  }

  # Add log file if specified
  if (!is.null(log_file)) {
    args <- c(args, "--log", log_file)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Handle output if output files are NULL
  if ((file_format == "fasta" && is.null(fastaout)) ||
      (file_format == "fastq" && is.null(fastqout))) {

    # Extract statistics
    if (!is.null(reverse)){
      statistics <- parse_trim_filt_statistics(vsearch_output, fastx_input_name, reverse_name)
    } else {
      statistics <- parse_trim_filt_statistics(vsearch_output, fastx_input_name)
    }

    # Process primary sequences
    if (file_format == "fasta") {
      filt_seqs <- microseq::readFasta(outfile_fasta)
    } else if (file_format == "fastq") {
      filt_seqs <- microseq::readFastq(outfile_fastq)
    }

    # Process reverse sequences if provided
    if (!is.null(reverse)) {
      if (file_format == "fasta") {
        filt_reverse <- microseq::readFasta(outfile_fasta_rev)
      } else if (file_format == "fastq") {
        filt_reverse <- microseq::readFastq(outfile_fastq_rev)
      }
    }

    # Add additional tables as attributes to the primary table
    attr(filt_seqs, "statistics") <- statistics
    if (!is.null(reverse)) {
      attr(filt_seqs, "filt_reverse") <- filt_reverse
    }
  }

  # Return results
  if ((file_format == "fasta" && is.null(fastaout)) ||
      (file_format == "fastq" && is.null(fastqout))) {
    return(filt_seqs)
  } else {
    return(invisible(NULL))
  }
}


#' Parse filtering statistics from string to tibble
#'
#' @description This function transforms the output from \code{vsearch} when running \code{vs_fastx_trim_filt()} into a tibble. The most important statistics are included in the tibble such as kept, truncated, and discarded sequences.
#'
#' @param output A string of output from filtering reads based on quality with \code{vsearch}.
#' @param fastq The name of the file/object with R1 reads.
#' @param reverse The name of the file/object with R2 reads
#'
#' @return A tibble with filtering metrics, including number of kept, truncated and discarded sequences after filtering.
#'
#' @noRd
parse_trim_filt_statistics <- function(output, fastx, reverse = NULL) {

  # Find line with statistics
  stats_line <- stringr::str_subset(output, "sequences kept")

  # Extract number of kept sequences
  kept <- as.numeric(stringr::str_extract(stats_line, "^\\d+"))

  # Extract number of truncated sequences
  truncated <- as.numeric(stringr::str_extract(stats_line, "(?<=of which )\\d+(?= truncated)"))

  # Extract number of discarded sequences
  discarded <- as.numeric(stringr::str_extract(stats_line, "(?<=, )\\d+(?= sequences discarded)"))

  # Create table
  result_table <- data.frame(
    Kept_Sequences = kept,
    Truncated_Sequences = truncated,
    Discarded_Sequences = discarded,
    fastx_source = fastx
  )

  # Add reverse column if provided
  if (!is.null(reverse)){
    result_table$reverse_source <- reverse
  }

  return(result_table)
}
