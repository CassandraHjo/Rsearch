#' Quality filtering
#'
#' @description Quality filtering of FASTQ sequences.
#'
#' @param fastq_input A FASTQ file path containing (forward) reads or a FASTQ object (tibble), see Details.
#' @param reverse An optional FASTQ file path containing reverse reads or a FASTQ object (tibble), see Details. If provided, it will be processed alongside \code{fastq_input}.
#' @param output_format Desired output format of tibble: \code{"fasta"} or \code{"fastq"}. Determines the format for both forward and reverse outputs.
#' @param fastq_maxee_rate Threshold for average expected error. Numeric value ranging form \code{0.0} to \code{1.0}. See Details.
#' @param fastaout Name of the FASTA output file for primary sequences (forward reads). If \code{NULL} no FASTA sequences will be written to file. See Details.
#' @param fastqout Name of the FASTQ output file for primary sequences (forward reads). If \code{NULL} no FASTQ sequences will be written to file. See Details.
#' @param fastaout_rev Name of the FASTA output file for reverse reads. If \code{NULL} no FASTA sequences will be written to file. See Details.
#' @param fastqout_rev Name of the FASTQ output file for reverse reads. If \code{NULL} no FASTQ sequences will be written to file. See Details.
#' @param fasta_width Number of characters per line in the output FASTA file. Only applies if the output file is in FASTA format. See Details.
#' @param minlen The minimum number of bases a sequence must have to be retained. Default is 0. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}.
#' @param log_file Name of the log file to capture messages from \code{vsearch}. If \code{NULL}, no log file is created.
#'
#' @details The function filters sequences from the input FASTQ file or object based on the average expected error rate using \code{vsearch}.
#' If a \code{reverse} input is provided, it filters the reverse reads similarly. The output format for both primary and reverse sequences
#' is determined by the \code{output_format} parameter.
#'
#' \code{fastq_input} and \code{reverse} can either be FASTQ files or FASTQ objects. If provided as tibbles, they must contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#' \code{reverse} is an optional argument to the function. If provided, it will be processed alongside \code{fastq_input}, meaning the same \code{fastq_maxee_rate} will be used for both FASTQ objects.
#'
#' If \code{fastaout}, \code{fastqout}, \code{fastaout_rev}, or \code{fastqout_rev} are specified, the remaining sequences after quality filtering are output to these files in either FASTA or FASTQ format.
#' If unspecified (\code{NULL}) no output is written to file. \code{output_format} has to match the desired output files.
#'
#' Sequences with an average expected error greater than the specified \code{fastq_maxee_rate} are discarded.
#' For a given sequence, the average expected error is the sum of error probabilities for all the positions in the sequence, divided by the length of the sequence.
#'
#' FASTA files produced by \code{vsearch} are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' Any input sequence with fewer bases than the value set in \code{minlen} will be discarded. By default, \code{minlen} is set to 0, which means that no sequences are removed.
#' However, using the default value may allow empty sequences to remain in the results.
#'
#' @return If output files are not specified, a tibble containing the filtered reads from \code{fastq_input} in the format specified by \code{output_format} is returned. If output files are specified, nothing is returned.
#'
#' If \code{reverse} is specified, the resulting tibble (\code{filt_reverse}) containing the filtered reverse reads in the format specified by \code{output_format} is an attribute to the primary table (\code{filt_seqs}).
#' This table can be accessed by running \code{attributes(filt_seqs)$filt_reverse} or \code{attr(filt_seqs, "filt_reverse")}.
#'
#' When a FASTA/FASTQ object is returned, the statistics from the filtering, \code{statistics}, is an attribute of the filtering tibble (\code{filt_seqs}).
#' This tibble contains filtering statistics, including number of kept and discarded sequences, and the names of the FASTQ files or objects that were filtered.
#' The statistics can be accessed by running \code{attributes(filt_seqs)$statistics} or \code{attr(filt_seqs, "statistics")}.
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastq_filter <- function(fastq_input,
                            reverse = NULL,
                            output_format = "fasta",
                            fastq_maxee_rate = 0.01,
                            fastaout = NULL,
                            fastqout = NULL,
                            fastaout_rev = NULL,
                            fastqout_rev = NULL,
                            fasta_width = 0,
                            minlen = 0,
                            threads = 1,
                            log_file = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate output_format
  if (!output_format %in% c("fasta", "fastq")) {
    stop("Invalid output_format. Choose from fasta or fastq.")
  }

  # If output_format is "fasta", fastqout and fastqout_rev can not be defined
  if (output_format == "fasta") {
    if (!is.null(fastqout) || !is.null(fastqout_rev)) {
      stop("When output_format is defined as 'fasta', 'fastqout' and 'fastqout_rev' cannot be used. Use 'fastaout' and 'fastaout_rev' instead.")
    }
  }

  # If output_format is "fastq", fastaout and fastaout_rev can not be defined
  if (output_format == "fastq") {
    if (!is.null(fastaout) || !is.null(fastaout_rev)) {
      stop("When output_format is defined as 'fastq', 'fastaout' and 'fastaout_rev' cannot be used. Use 'fastqout' and 'fastqout_rev' instead.")
    }
  }

  # If reverse is specified, ensure paired output parameters are both NULL or both character strings
  if (!is.null(reverse)) {
    if (output_format == "fasta") {
      # Check that both fastaout and fastaout_rev are NULL or both are character strings
      if ((is.null(fastaout) && !is.null(fastaout_rev)) ||
          (!is.null(fastaout) && is.null(fastaout_rev))) {
        stop("When 'reverse' is specified and output_format is 'fasta', both 'fastaout' and 'fastaout_rev' must be NULL or both specified as character strings.")
      }
    }

    if (output_format == "fastq") {
      # Check that both fastqout and fastqout_rev are NULL or both are character strings
      if ((is.null(fastqout) && !is.null(fastqout_rev)) ||
          (!is.null(fastqout) && is.null(fastqout_rev))) {
        stop("When 'reverse' is specified and output_format is 'fastq', both 'fastqout' and 'fastqout_rev' must be NULL or both specified as character strings.")
      }
    }
  }

  # Create empty vector for collecting temporary files
  temp_files <- c()

  # Set up cleanup of temporary files
  on.exit({
    existing_files <- temp_files[file.exists(temp_files)]
    if (length(existing_files) > 0) {
      file.remove(existing_files)
    }
  }, add = TRUE)

  # Handle output for primary sequences
  if (output_format == "fasta") {
    if (is.null(fastaout)) {
      outfile_fasta <- tempfile(pattern = "filtered_primary_", fileext = ".fa")
      temp_files <- c(temp_files, outfile_fasta)
    } else {
      outfile_fasta <- fastaout
    }
  }

  if (output_format == "fastq") {
    if (is.null(fastqout)) {
      outfile_fastq <- tempfile(pattern = "filtered_primary_", fileext = ".fq")
      temp_files <- c(temp_files, outfile_fastq)
    } else {
      outfile_fastq <- fastqout
    }
  }

  # Handle output for reverse sequences
  if (!is.null(reverse)) {
    if (output_format == "fasta") {
      if (is.null(fastaout_rev)) {
        outfile_fasta_rev <- tempfile(pattern = "filtered_reverse_", fileext = ".fa")
        temp_files <- c(temp_files, outfile_fasta_rev)
      } else {
        outfile_fasta_rev <- fastaout_rev
      }
    }

    if (output_format == "fastq") {
      if (is.null(fastqout_rev)) {
        outfile_fastq_rev <- tempfile(pattern = "filtered_reverse_", fileext = ".fq")
        temp_files <- c(temp_files, outfile_fastq_rev)
      } else {
        outfile_fastq_rev <- fastqout_rev
      }
    }
  }

  # Handle input for primary sequences: file or tibble
  if (!is.character(fastq_input)){
    # Ensure required columns exist
    required_cols <- c("Header", "Sequence", "Quality")
    if (!all(required_cols %in% colnames(fastq_input))) {
      stop("FASTQ object must contain columns: Header, Sequence, Quality")
    }
    temp_fastq_file <- tempfile(pattern = "fastq_input_temp_", fileext = ".fq")
    microseq::writeFastq(fastq_input, temp_fastq_file)
    fastq_file <- temp_fastq_file
    temp_files <- c(temp_files, temp_fastq_file)

    # Capture original name for statistics table later
    fastq_input_name <- as.character(substitute(fastq_input))

  } else {
    fastq_file <- fastq_input

    # Capture original name for statistics table later
    fastq_input_name <- basename(fastq_input)
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
      reverse_file <- temp_reverse_file
      temp_files <- c(temp_files, temp_reverse_file)

      # Capture original name for statistics table later
      reverse_name <- as.character(substitute(reverse))

    } else {
      reverse_file <- reverse

      # Capture original name for statistics table later
      reverse_name <- basename(reverse)
    }
  }

  # Check if input files exists
  if (!file.exists(fastq_file)) stop("Cannot find input FASTQ file: ", fastq_file)
  if (!is.null(reverse) && !file.exists(reverse_file)) stop("Cannot find reverse FASTQ file: ", reverse_file)

  # Normalize file paths
  fastq_file <- normalizePath(fastq_file)
  if (!is.null(reverse)) {
    reverse_file <- normalizePath(reverse_file)
  }

  # Build argument string for command line
  args <- c("--fastq_filter", fastq_file,
            "--threads", threads,
            "--fastq_maxee_rate", fastq_maxee_rate,
            "--fastq_minlen", minlen)

  # Add reverse to arguments if provided
  if (!is.null(reverse)) {
    args <- c(args, "--reverse", reverse_file)
  }

  # Add output files based on output_format
  if (output_format == "fastq") {
    args <- c(args, "--fastqout", outfile_fastq)
    if (!is.null(reverse)) {
      args <- c(args, "--fastqout_rev", outfile_fastq_rev)
    }
  } else if(output_format == "fasta") {
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
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {

    # Extract statistics
    if (!is.null(reverse)){
      statistics <- parse_filter_statistics(vsearch_output, fastq_input_name, reverse_name)
    } else {
      statistics <- parse_filter_statistics(vsearch_output, fastq_input_name)
    }

    # Process primary sequences
    if (output_format == "fasta") {
      filt_seqs <- microseq::readFasta(outfile_fasta)
    } else if (output_format == "fastq") {
      filt_seqs <- microseq::readFastq(outfile_fastq)
    }

    # Process reverse sequences if provided
    if (!is.null(reverse)) {
      if (output_format == "fasta") {
        filt_reverse <- microseq::readFasta(outfile_fasta_rev)
      } else if (output_format == "fastq") {
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
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {
    return(filt_seqs)
  } else {
    return(invisible(NULL))
  }
}


#' Parse filtering statistics from string to tibble
#'
#' @description This function transforms the output from \code{vsearch} when running \code{vs_fastq_filter()} into a tibble. The most important statistics are included in the tibble such as kept, truncated, and discarded sequences.
#'
#' @param output A string of output from filtering reads based on quality with \code{vsearch}.
#' @param fastq The name of the file/object with R1 reads.
#' @param reverse The name of the file/object with R2 reads
#'
#' @return A tibble with filtering metrics, including number of kept, truncated and discarded sequences after filtering.
#'
#' @noRd
parse_filter_statistics <- function(output, fastq, reverse = NULL) {

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
    fastq_source = fastq
  )

  # Add reverse column if provided
  if (!is.null(reverse)){
    result_table$reverse_source <- reverse
  }

  return(result_table)
}
