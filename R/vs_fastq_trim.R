#' Trimming sequences in FASTQ format
#'
#' @description Removes given number of nucleotides from the ends of biological sequences (DNA) in the FASTQ format.
#'
#' @param fastq_input A FASTQ file path or FASTQ object containing (forward) reads. See Details.
#' @param reverse An optional FASTQ file path or FASTQ object, if the input consists of paired sequences, containing reverse reads. If provided, it will be processed alongside \code{fastq_input}. Defaults to \code{NULL}. See Details.
#' @param output_format Desired output format of file or tibble: \code{"fasta"} or \code{"fastq"}. Determines the format for both forward and reverse (if provided) outputs.Defaults to \code{"fasta"}.
#' @param stripright The number of bases stripped from the right end of the reads. Defaults to \code{0}.
#' @param stripleft The number of bases stripped from the left end of the reads. Defaults to \code{0}.
#' @param fastaout Name of the FASTA output file for the sequences given in \code{fastq_input}. If \code{NULL} no FASTA sequences will be written to file. Defaults to \code{NULL}. See Details.
#' @param fastqout fastqout Name of the FASTQ output file for the sequences given in \code{fastq_input}. If \code{NULL} no FASTQ sequences will be written to file. Defaults to \code{NULL}. See Details.
#' @param fastaout_rev Name of the FASTA output file for the sequences given in \code{reverse}. If \code{NULL} no FASTA sequences will be written to file. Defaults to \code{NULL}. See Details.
#' @param fastqout_rev Name of the FASTQ output file for the sequences given in \code{reverse}. If \code{NULL} no FASTQ sequences will be written to file. Defaults to \code{NULL}. See Details.
#' @param fasta_width Number of characters per line in the output FASTA file. Only applies if the output file is in FASTA format. Defaults to \code{0}. See Details.
#' @param minlen The minimum number of bases a sequence must have to be retained. Defaults to \code{0}. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details The reads in the input FASTQ-file (\code{fastq_input}) are trimmed based on the specified number of bases for each end of the read, using \code{vsearch}.
#' If a \code{reverse} input is provided, it trims the reverse reads similarly. The output format for both primary and reverse sequences is determined by \code{output_format}.
#'
#' \code{fastq_input} and \code{reverse} can either be FASTQ files or FASTQ objects. FASTQ objects are tibbles that contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#' \code{reverse} is an optional argument to the function. If provided, it will be processed alongside \code{fastq_input}, meaning the same \code{stripright} and \code{stripleft} will be used for both inputs.
#'
#' If \code{fastaout} and \code{fastaout_rev} or \code{fastqout} and \code{fastqout_rev} are specified, the remaining sequences after trimming are output to these files in either FASTA or FASTQ format.
#' If unspecified (\code{NULL}) no output is written to file. \code{output_format} has to match the desired output files/objects.
#'
#' FASTA files produced by \code{vsearch} are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' Any input sequence with fewer bases than the value set in \code{minlen} will be discarded. By default, \code{minlen} is set to 0, which means that no sequences are removed.
#' However, using the default value may allow empty sequences to remain in the results.
#'
#' @return If output files are not specified, a tibble containing the trimmed reads from \code{fastq_input} in the format specified by \code{output_format} is returned. If output files are specified, results are written to file and nothing is returned.
#'
#' If \code{reverse} is specified, the resulting tibble (\code{trimmed_reverse}) containing the filtered reverse reads in the format specified by \code{output_format} is an attribute to the primary table (\code{trimmed_seqs}).
#' This table can be accessed by running \code{attributes(trimmed_seqs)$trimmed_reverse} or \code{attr(trimmed_seqs, "trimmed_reverse")}.
#'
#' @examples
#' \dontrun{
#' # Read example FASTQ files
#' fastq_input <- file.path(file.path(path.package("Rsearch"), "inst/extdata"), "R1_sample1_small.fq")
#' reverse <- file.path(file.path(path.package("Rsearch"), "inst/extdata"), "R2_sample1_small.fq")
#'
#' # Define other arguments
#' output_format <- "fastq"
#' stripright <- 10
#' stripleft <- 10
#' minlen <- 0
#'
#' # Execute trimming, with tibble as output
#' trim_seqs <- vs_fastq_trim(fastq_input = fastq_input,
#'                            reverse = reverse,
#'                            output_format = output_format,
#'                            stripright = stripright,
#'                            stripleft = stripleft,
#'                            minlen = minlen)
#'
#' # Extract tibbles with trimmed sequences
#' R1_trim <- trim_seqs
#' R2_trim <- attr(trim_seqs, "trimmed_reverse")
#'
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastq_trim <- function(fastq_input,
                          reverse = NULL,
                          output_format = "fasta",
                          stripright = 0,
                          stripleft = 0,
                          fastaout = NULL,
                          fastqout = NULL,
                          fastaout_rev = NULL,
                          fastqout_rev = NULL,
                          minlen = 0,
                          fasta_width = 0,
                          threads = 1){


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
    if (length(temp_files) > 0 && is.character(temp_files)) {
      existing_files <- temp_files[file.exists(temp_files)]
      if (length(existing_files) > 0) {
        file.remove(existing_files)
      }
    }
  }, add = TRUE)

  # Handle output for primary sequences
  if (output_format == "fasta") {
    if (is.null(fastaout)) {
      outfile_fasta <- tempfile(pattern = "trimmed_primary_", fileext = ".fa")
      temp_files <- c(temp_files, outfile_fasta)
    } else {
      outfile_fasta <- fastaout
    }
  }

  if (output_format == "fastq") {
    if (is.null(fastqout)) {
      outfile_fastq <- tempfile(pattern = "trimmed_primary_", fileext = ".fq")
      temp_files <- c(temp_files, outfile_fastq)
    } else {
      outfile_fastq <- fastqout
    }
  }

  # Handle output for reverse sequences
  if (!is.null(reverse)) {
    if (output_format == "fasta") {
      if (is.null(fastaout_rev)) {
        outfile_fasta_rev <- tempfile(pattern = "trimmed_reverse_", fileext = ".fa")
        temp_files <- c(temp_files, outfile_fasta_rev)
      } else {
        outfile_fasta_rev <- fastaout_rev
      }
    }

    if (output_format == "fastq") {
      if (is.null(fastqout_rev)) {
        outfile_fastq_rev <- tempfile(pattern = "trimmed_reverse_", fileext = ".fq")
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
  } else {
    fastq_file <- fastq_input
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
    } else {
      reverse_file <- reverse
    }
  }

  # Check is input files exists
  if (!file.exists(fastq_file)) stop("Cannot find input FASTQ file: ", fastq_file)
  if (!is.null(reverse) && !file.exists(reverse_file)) stop("Cannot find reverse FASTQ file: ", reverse_file)


  # Normalize file paths
  fastq_file <- normalizePath(fastq_file)
  if (!is.null(reverse)) {
    reverse_file <- normalizePath(reverse_file)
  }

  # Build argument string for command line
  args <- c("--fastq_filter", fastq_file,
            "--fastq_stripright", stripright,
            "--fastq_stripleft", stripleft,
            "--fastq_minlen", minlen,
            "--threads", threads)

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

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  # Handle output if output files are NULL
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {

    # Process primary sequences
    if (output_format == "fasta") {
      trimmed_seqs <- microseq::readFasta(outfile_fasta)
    } else if (output_format == "fastq") {
      trimmed_seqs <- microseq::readFastq(outfile_fastq)
    }

    # Process reverse sequences if provided
    if (!is.null(reverse)) {
      if (output_format == "fasta") {
        trimmed_reverse <- microseq::readFasta(outfile_fasta_rev)
      } else if (output_format == "fastq") {
        trimmed_reverse <- microseq::readFastq(outfile_fastq_rev)
      }
    }

    # Add additional tables as attributes to the primary table
    if (!is.null(reverse)) {
      attr(trimmed_seqs, "trimmed_reverse") <- trimmed_reverse
    }
  }

  # Return results
  if ((output_format == "fasta" && is.null(fastaout)) ||
      (output_format == "fastq" && is.null(fastqout))) {
    return(trimmed_seqs)
  } else {
    return(invisible(NULL))
  }
}
