#' Subsample sequences
#'
#' @description \code{vs_fastx_subsample} subsamples sequences in FASTA/FASTQ
#' file or object by randomly extracting sequences based on number or percentage
#' using \code{VSEARCH}.
#'
#' @param fastx_input A FASTA/FASTQ file path or FASTA/FASTQ object. See
#' \emph{Details}.
#' @param output_format Desired output format of file or tibble: \code{"fasta"}
#' or \code{"fastq"} (default). If \code{fastx_input} is a FASTA file path or a
#' FASTA object, \code{output_format} cannot be \code{"fastq"}.
#' @param fastx_output Name of the output file for subsampled reads from
#' \code{fastx_input}. File can be in either FASTA or FASTQ format, depending on
#' \code{output_format}. If \code{NULL} (default), no sequences are written to
#' file. See \emph{Details}.
#' @param sample_pct Percentage of the input sequences to be subsampled.
#' Numeric value ranging from \code{0.0} to \code{100.0}. Defaults to \code{NULL}.
#' @param sample_size The given number of sequences to extract.
#' Must be a positive integer if specified. Defaults to \code{NULL}.
#' @param sizein If \code{TRUE} (default), abundance annotations present in
#' sequence headers are taken into account.
#' @param sizeout If \code{TRUE} (default), abundance annotations are added to
#' FASTA headers.
#' @param relabel Relabel sequences using the given prefix and a ticker to
#' construct new headers. Defaults to \code{NULL}.
#' @param relabel_sha1 If \code{TRUE} (default), relabel sequences using the
#' SHA1 message digest algorithm. Defaults to \code{FALSE}.
#' @param randseed Random seed. Must be a positive integer. A given seed always
#' produces the same output, which is useful for replicability. Defaults to
#' \code{NULL}.
#' @param fasta_width Number of characters per line in the output FASTA
#' file. Defaults to \code{0}, which eliminates wrapping.
#' @param threads Number of computational threads to be used by \code{VSEARCH}.
#' Defaults to \code{1}.
#' @param vsearch_options Additional arguments to pass to \code{VSEARCH}.
#' Defaults to \code{NULL}. See \emph{Details}.
#'
#' @details
#' Sequences in the input file/object (\code{fastx_input}) are subsampled by
#' randomly extracting a specified number or percentage of sequences. Extraction
#' is performed as random sampling with a uniform distribution among the input
#' sequences and without replacement.
#'
#' \code{fastx_input} can either be a FASTA/FASTQ file or a FASTA/FASTQ object.
#' FASTA objects are tibbles that contain the columns \code{Header} and
#' \code{Sequence}.FASTQ objects are tibbles that contain the columns
#' \code{Header}, \code{Sequence}, and \code{Quality}.
#'
#' Specify either \code{sample_size} or \code{sample_pct} to determine the
#' number or percentage of sequences to subsample. Only one of these parameters
#' can be specified at a time. If neither is specified, an error is thrown.
#'
#' If \code{fastx_output} is specified, the sampled sequences are output to this
#' file in format given by \code{output_format}.
#' If \code{fastx_output} is \code{NULL}, the sample sequences are returned as a
#' FASTA or FASTQ object, depending on \code{output_format}.
#'
#' \code{vsearch_options} allows users to pass additional command-line arguments
#' to \code{VSEARCH} that are not directly supported by this function. Refer to
#' the \code{VSEARCH} manual for more details.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{fastx_output} is specified, the subsampled sequences are written to
#' the specified output file, and no tibble is returned.
#'
#' If \code{fastx_output} \code{NULL}, a tibble containing the subsampled reads
#' in the format specified by \code{output_format} is returned.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' fastx_input <- file.path(file.path(path.package("Rsearch"), "extdata"), "R1_sample1_small.fq")
#' fastx_output <- NULL
#' output_format <- "fastq"
#' sample_size <- 10
#'
#' # Subsample sequences and return a FASTQ tibble
#' subsample_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
#'                                    fastx_output = fastx_output,
#'                                    output_format = output_format,
#'                                    sample_size = sample_size)
#'
#' # Subsample sequences and write subsampled sequences to a file
#' vs_fastx_subsample(fastx_input = fastx_input,
#'                    fastx_output = "subsample.fq",
#'                    output_format = output_format,
#'                    sample_size = sample_size)
#' }
#'
#' @references \url{https://github.com/torognes/vsearch}
#'
#' @aliases vs_fastx_subsample vs_fastq_subsample vs_fasta_subsample
#' vs_subsample fastx_subsample subsample
#' @export
#'
vs_fastx_subsample <- function(fastx_input,
                               output_format = "fastq",
                               fastx_output = NULL,
                               sample_pct = NULL,
                               sample_size = NULL,
                               sizein = TRUE,
                               sizeout = TRUE,
                               relabel = NULL,
                               relabel_sha1 = FALSE,
                               randseed = NULL,
                               fasta_width = 0,
                               threads = 1,
                               vsearch_options = NULL){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Validate output_format
  if (!output_format %in% c("fasta", "fastq")) {
    stop("Invalid output_format. Choose from fasta or fastq.")
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

  # Handle input
  if (!is.character(fastx_input)){
    if ("Quality" %in% colnames(fastx_input)){

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

      if (output_format == "fastq") {
        stop("Invalid output_format when input tibble is of type 'fasta'")
      }

      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(fastx_input))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_file <- tempfile(pattern = "input", fileext = ".fa")
      temp_files <- c(temp_files, temp_file)
      microseq::writeFasta(fastx_input, temp_file)
      input_file <- temp_file
    }
  } else {
    input_file <- fastx_input
  }

  # Handle output_format = "fasta"
  if (output_format == "fasta") {
    if (is.null(fastx_output)) {
      output_file <- tempfile(pattern = "subsample", fileext = ".fa")
      temp_files <- c(temp_files, output_file)
    } else {
      output_file <- fastx_output
    }
  }

  # Handle output_format = "fastq"
  if (output_format == "fastq") {
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
  input_file <- normalizePath(input_file)

  # Build argument string for command line
  args <- c("--fastx_subsample", shQuote(input_file),
            "--threads", threads)

  if (output_format == "fasta") {
    args <- c(args,
              "--fasta_width", fasta_width,
              "--fastaout", output_file)
  }

  if (output_format == "fastq") {
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

  # Add additional arguments if specified
  if (!is.null(vsearch_options)) {
    args <- c(args, vsearch_options)
  }

  # Run vsearch
  vsearch_output <- system2(command = vsearch_executable,
                            args = args,
                            stdout = TRUE,
                            stderr = TRUE)

  if (is.null(fastx_output)) {
    if (output_format == "fasta"){
      # Read output into FASTA object
      derep_tbl <- microseq::readFasta(output_file)
    }
    if (output_format == "fastq"){
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
