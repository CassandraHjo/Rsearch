#' Subsample sequences
#'
#' @description Subsample sequences in FASTA/FASTQ file or object by randomly extracting sequences based on number or precentage.
#'
#' @param fastx_input A FASTQ/FASTA file path or object. See Details.
#' @param fastx_output Name of the output file for dereplicated reads from \code{fastx_input}. File can be in either FASTA or FASTQ format, depending on \code{file_format}. If \code{NULL} (default) no sequences will be written to file. See Details.
#' @param file_format Format of input file \code{fastx_input}, and desired output format for file/tibble: \code{"fasta"} or \code{"fastq"} (default).
#' @param sample_pct The given percentage of the input sequences to be subsampled. Numeric value ranging from \code{0.0} to \code{100.0}. Defaults to \code{NULL}.
#' @param sample_size The given number of sequences to extract. Must be a positive integer if specified. Defaults to \code{NULL}.
#' @param randseed Random seed. Must be a positive integer. A given seed always produces the same output, which is useful for replicability. Defaults to \code{NULL}.
#' @param sizein Decides if abundance annotations present in sequence headers in \code{fastx_input} should be taken into account. Defaults to \code{TRUE}.
#' @param sizeout Decides if abundance annotations should be added to headers in output table or output file (\code{fastx_output}). Defaults to \code{TRUE}.
#' @param relabel_sha1 Relabel sequences using the SHA1 message digest algorithm. Defaults to \code{FALSE}.
#' @param relabel Relabel sequences using the given prefix and a ticker to construct new headers. Defaults to \code{NULL}.
#' @param fasta_width Number of characters per line in the output FASTA file. Only applies if the output file is in FASTA format. Defaults to \code{0}. See Details.
#' @param threads Number of computational threads to be used by \code{vsearch}. Defaults to \code{1}.
#'
#' @details The reads in the input file/object (\code{fastx_input}) are subsampled by randomly extracting a certain number or a certain percentage of the sequences in the input, using \code{vsearch}.
#' The extraction is performed as a random sampling with a uniform distribution among the input sequences and is performed without replacement.
#'
#' \code{fastx_input} can either be a FASTA/FASTQ file or object with reads. FASTA objects are tibbles that contain the columns \code{Header} and \code{Sequence}. FASTQ objects are tibbles that contain the columns \code{Header}, \code{Sequence}, and \code{Quality}.
#'
#' If \code{fastx_output} is specified, the sampled sequences are output to this file in format given by \code{file_format}.
#' If unspecified (\code{NULL}) the result is returned as a FASTA/FASTQ object, depending on \code{file_format}.
#'
#' FASTA files produced by\code{vsearch} are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' @return Tibble or \code{NULL}.
#'
#' If \code{fastx_output} is not specified, a tibble containing the dereplicated reads is returned. If \code{fastx_output} is specified nothing is returned.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fq")
#' fastx_output <- NULL
#' file_format <- "fastq"
#' sample_size <- 1
#'
#' # Subsample sequences, with tibble as output
#' subsample_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
#'                                    fastx_output = fastx_output,
#'                                    file_format = file_format,
#'                                    sample_size = sample_size)
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @export
#'
vs_fastx_subsample <- function(fastx_input,
                               fastx_output = NULL,
                               file_format = "fastq",
                               sample_pct = NULL,
                               sample_size = NULL,
                               randseed = NULL,
                               sizein = TRUE,
                               sizeout = TRUE,
                               relabel_sha1 = FALSE,
                               relabel = NULL,
                               fasta_width = 0,
                               threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate file_format
  if (!file_format %in% c("fasta", "fastq")) {
    stop("Invalid file_format. Choose from fasta or fastq.")
  }

  # Validate that only sample_pct or sample_size is specified

  # At least one of them must be specified
  if (is.null(sample_size) && is.null(sample_pct)) {
    stop("Either sample_size or sample_pct must be specified.")
  }

  # Only one option can be specified at a time
  if (!is.null(sample_size) && !is.null(sample_pct)) {
    stop("Only specify one of the following parameters, not both: sample_size, sample_pct ")
  }

  # Create empty vector for collecting temporary files
  temp_files <- character()

  # Set up cleanup of temporary files
  on.exit({
    existing_files <- temp_files[file.exists(temp_files)]
    if (length(existing_files) > 0) {
      file.remove(existing_files)
    }
  }, add = TRUE)

  # Handle file_format = "fasta"
  if (file_format == "fasta") {
    # Handle input: file or tibble
    if (!is.character(fastx_input)){
      # Validate tibble
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }
      temp_file <- tempfile(pattern = "input", fileext = ".fa")
      temp_files <- c(temp_files, temp_file)
      microseq::writeFasta(fastx_input, temp_file)
      input_file <- temp_file
    } else {
      input_file <- fastx_input
    }

    # Handle output
    if (is.null(fastx_output)) {
      output_file <- tempfile(pattern = "subsample", fileext = ".fa")
      temp_files <- c(temp_files, output_file)
    } else {
      output_file <- fastx_output
    }
  }

  # Handle file_format = "fastq"
  if (file_format == "fastq") {
    # Handle input: file or tibble
    if (!is.character(fastx_input)){
      # Validate tibble
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }
      temp_file <- tempfile(pattern = "input", fileext = ".fq")
      temp_files <- c(temp_files, temp_file)
      microseq::writeFastq(fastx_input, temp_file)
      input_file <- temp_file
    } else {
      input_file <- fastx_input
    }

    # Handle output
    if (is.null(fastx_output)) {
      output_file <- tempfile(pattern = "subsample", fileext = ".fq")
      temp_files <- c(temp_files, output_file)
    } else {
      output_file <- fastx_output
    }
  }

  # Check is input file exists at given path
  if (!file.exists(input_file)) stop("Cannot find input file: ", input_file)

  # Normalize file paths
  input_file <- normalizePath(input_file) |>
    shQuote()

  # Build argument string for command line
  args <- c("--fastx_subsample", input_file,
            "--threads", threads)

  if (file_format == "fasta") {
    args <- c(args,
              "--fasta_width", fasta_width,
              "--fastaout", output_file)
  }

  if (file_format == "fastq") {
    args <- c(args, "--fastqout", output_file)
  }

  if (!is.null(sample_size)) {
    args <- c(args, "--sample_size", sample_size)
  }

  if (!is.null(sample_pct)) {
    args <- c(args, "--sample_pct", sample_pct)
  }

  if (!is.null(randseed)) {
    args <- c(args, "--randseed", randseed)
  }

  if (sizein) {
    args <- c(args, "--sizein", "")
  }

  if (sizeout) {
    args <- c(args, "--sizeout", "")
  }

  if (relabel_sha1) {
    args <- c(args, "--relabel_sha1", "")
  }

  if (!is.null(relabel)) {
    args <- c(args, "--relabel", relabel)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  if (is.null(fastx_output)) {
    if (file_format == "fasta"){
      # Read output into FASTA object
      derep_tbl <- microseq::readFasta(output_file)
    }
    if (file_format == "fastq"){
      # Read output into FASTQ object
      derep_tbl <- microseq::readFastq(output_file)
    }
  }

  # Return results
  if (is.null(fastx_output)) { # Return tibble
    return(derep_tbl)
  } else {
    return(invisible(NULL)) # No return when output file is written
  }
}
