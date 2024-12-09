#' Synchronize FASTQ and FASTA files and objects
#'
#' @param file1 A FASTQ/FASTA file path or a FASTQ/FASTA object (tibble), see Details.
#' @param file2 A FASTQ/FASTA file path or a FASTQ/FASTA object (tibble), see Details.
#' @param file_format Format of input files \code{file1} and \code{file2}, and desired output format of tibbles: \code{"fasta"} or \code{"fastq"}. Determines the format for both outputs.
#' @param file1_out Name of the output file for synchronized reads from \code{file1}. File can be in either FASTA or FASTQ format, depending on \code{file_format}. If \code{NULL} no sequences will be written to file. See Details.
#' @param file2_out Name of the output file for synchronized reads from \code{file2}. File can be in either FASTA or FASTQ format, depending on \code{file_format}. If \code{NULL} no sequences will be written to file. See Details.
#'
#' @description The function synchronizes sequences in two FASTA/FASTQ files or objects, by retaining the common sequences.
#'
#' @details
#' \code{file1} and \code{file2} can either be FASTQ/FASTA files or FASTQ/FASTA objects. If provided as tibbles, they must contain the columns \code{Header}, \code{Sequence}, and \code{Quality} or the columns \code{Header} and \code{Sequence}, depending on the \code{file_format}.
#'
#' #' If \code{file1_out} or \code{file2_out} are specified, the remaining sequences after synchronizing are output to these files in either FASTA or FASTQ format.
#' If unspecified (\code{NULL}) no output is written to file.
#'
#' @return A list with two tibbles:
#' \describe{
#'   \item{sync_file1}{A tibble containing the synchronized reads from \code{file1}.}
#'   \item{sync_file2}{A tibble containing the synchronized reads from \code{file2}.}
#'   }
#'
#' @export
#'
fastx_synchronize <- function(file1,
                              file2,
                              file_format = "fastq",
                              file1_out = NULL,
                              file2_out = NULL) {

  # Validate file_format
  if (!file_format %in% c("fasta", "fastq")) {
    stop("Invalid file_format. Choose from fasta or fastq.")
  }

  # Create empty vector for collecting temporary files
  temp_files <- c()

  # Handle input file1: file or tibble
  if (!is.character(file1)){
    if (file_format == "fastq") {
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(file1))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }

      temp_fastq_file <- tempfile(pattern = "file1_temp_", fileext = ".fq")
      microseq::writeFastq(file1, temp_fastq_file)
      file1 <- temp_fastq_file
      temp_files <- c(temp_files, temp_fastq_file)
    }
    if (file_format == "fasta") {
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(file1))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_fasta_file <- tempfile(pattern = "file1_temp_", fileext = ".fa")
      microseq::writeFasta(file1, temp_fasta_file)
      file1 <- temp_fasta_file
      temp_files <- c(temp_files, temp_fasta_file)
    }
  }

  # Handle input file2: file or tibble
  if (!is.character(file2)){
    if (file_format == "fastq") {
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(file2))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }

      temp_fastq_file <- tempfile(pattern = "file2_temp_", fileext = ".fq")
      microseq::writeFastq(file2, temp_fastq_file)
      file2 <- temp_fastq_file
      temp_files <- c(temp_files, temp_fastq_file)
    }
    if (file_format == "fasta") {
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(file2))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }

      temp_fasta_file <- tempfile(pattern = "file2_temp_", fileext = ".fa")
      microseq::writeFasta(file1, temp_fasta_file)
      file2 <- temp_fasta_file
      temp_files <- c(temp_files, temp_fasta_file)
    }
  }

  # Check is input files exists
  if (!file.exists(file1)) stop("Cannot find input file: ", file1)
  if (!file.exists(file2)) stop("Cannot find input file: ", file2)


  # Normalize file paths
  file1 <- normalizePath(file1)
  file2 <- normalizePath(file2)

  # Create tibbles from files
  if (file_format == "fastq"){
    R1 <- microseq::readFastq(file1)
    R2 <- microseq::readFastq(file2)
  }

  if (file_format == "fasta"){
    R1 <- microseq::readFasta(file1)
    R2 <- microseq::readFasta(file2)
  }

  # Create tag column with sequence id
  R1 <- R1 %>%
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+"))

  R2 <- R2 %>%
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+"))

  # Find common tags
  common_tags <- intersect(R1$tag, R2$tag)

  # Keep only sequences from common tags
  sync_file1 <- R1 %>%
    dplyr::filter(tag %in% common_tags) %>%
    dplyr::arrange(tag) %>%
    dplyr::select(-tag)

  sync_file2 <- R2 %>%
    dplyr::filter(Header %in% common_tags) %>%
    dplyr::arrange(Header) %>%
    dplyr::select(-tag)

  # Write output files if specified
  if (file_format == "fastq" && !is.null(file1_out)) {
    microseq::writeFastq(sync_file1, file1_out)
  }

  if (file_format == "fastq" && !is.null(file2_out)) {
    microseq::writeFastq(sync_file2, file2_out)
  }

  if (file_format == "fasta" && !is.null(file1_out)) {
    microseq::writeFasta(sync_file1, file1_out)
  }

  if (file_format == "fasta" && !is.null(file2_out)) {
    microseq::writeFasta(sync_file2, file2_out)
  }

  # Cleanup temporary files
  if (length(temp_files) > 0) {
    on.exit(
      file.remove(temp_files),
      add = TRUE)
  }

  return(list(sync_file1 = sync_file1, sync_file2 = sync_file2))
}
