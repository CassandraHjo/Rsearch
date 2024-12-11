#' Synchronize FASTQ and FASTA files and objects
#'
#' @param file1 A FASTQ/FASTA file path or a FASTQ/FASTA object (tibble), see Details.
#' @param file2 A FASTQ/FASTA file path or a FASTQ/FASTA object (tibble), see Details.
#' @param file_format Format of input files \code{file1} and \code{file2}, and desired output format: \code{"fasta"} or \code{"fastq"}. Determines the format for both outputs.
#' @param file1_out Name of the output file for synchronized reads from \code{file1}. File can be in either FASTA or FASTQ format, depending on \code{file_format}. If \code{NULL} no sequences will be written to file. See Details.
#' @param file2_out Name of the output file for synchronized reads from \code{file2}. File can be in either FASTA or FASTQ format, depending on \code{file_format}. If \code{NULL} no sequences will be written to file. See Details.
#'
#' @description The function synchronizes sequences in two FASTA/FASTQ files or objects, by retaining the common sequences.
#'
#' @details
#' \code{file1} and \code{file2} can either be FASTQ/FASTA files or FASTQ/FASTA objects. If provided as tibbles, they must contain the columns \code{Header}, \code{Sequence}, and \code{Quality} or the columns \code{Header} and \code{Sequence}, depending on \code{file_format}.
#' In order for the synchronizing to work, it is necessary that the sequence IDs in the \code{Header}s are identical for each read pair in the two files.
#'
#' If \code{file1_out} or \code{file2_out} are specified, the remaining sequences after synchronizing are output to these files in either FASTA or FASTQ format depending on \code{file_format}.
#' If unspecified (\code{NULL}) no output is written to file, and the synchronized reads are returned as a FASTQ/FASTA object (tibble). \code{file1_out} or \code{file2_out} must either both be \code{NULL} or both \code{charachter}.
#'
#' @return If output files are not specified, a tibble containing the synchronized reads from \code{file1} is returned. The tibble containing the synchronized reads from \code{file2} is an attribute to the primary table (\code{sync_file1}).
#' This table can be accessed by running \code{attributes(sync_file1)$sync_file2} or \code{attr(sync_file1, "sync_file2")}.
#' If output files, \code{file1_out} or \code{file2_out}, are specified nothing is returned.
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

  # Validate output files
  if ((is.null(file1_out) && !is.null(file2_out)) ||
      (!is.null(file1_out) && is.null(file2_out))) {
    stop("Either both file1_out and file2_out must be NULL, or both must be specified.")
  }

  if (!is.null(file1_out) && !is.character(file1_out)) {
    stop("file1_out must be a character string specifying the output file path.")
  }

  if (!is.null(file2_out) && !is.character(file2_out)) {
    stop("file2_out must be a character string specifying the output file path.")
  }

  # Handle input file1: file or tibble
  if (!is.character(file1)){ # If tibble
    if (file_format == "fastq") {
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(file1))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }
    }
    if (file_format == "fasta") {
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(file1))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }
    }
  } else {
    # Check if file 1 exists
    if (!file.exists(file1)) stop("Cannot find input file: ", file1)
    # Normalize file path
    file1 <- normalizePath(file1)

    if (file_format == "fastq") {
      file1 <- microseq::readFastq(file1)
    }

    if (file_format == "fasta") {
      file1 <- microseq::readFasta(file1)
    }
  }

  # Handle input file2: file or tibble
  if (!is.character(file2)){ # If tibble
    if (file_format == "fastq") {
      required_cols <- c("Header", "Sequence", "Quality")
      if (!all(required_cols %in% colnames(file2))) {
        stop("FASTQ object must contain columns: Header, Sequence, Quality")
      }
    }
    if (file_format == "fasta") {
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(file2))) {
        stop("FASTA object must contain columns: Header and Sequence")
      }
    }
  } else {
    # Check if file2 exists
    if (!file.exists(file2)) {
      stop("Cannot find input file: ", file2)
    }
    # Normalize file paths
    file2 <- normalizePath(file2)

    if (file_format == "fastq") {
      file2 <- microseq::readFastq(file2)
    }

    if (file_format == "fasta") {
      file2 <- microseq::readFasta(file2)
    }
  }

  # Create tag column with sequence id
  file1 <- file1 %>%
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) %>%
    dplyr::mutate(tag = stringr::str_remove(tag, "/1$")) %>%
    dplyr::mutate(tag = stringr::str_remove(tag, "/2$"))

  file2 <- file2 %>%
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) %>%
    dplyr::mutate(tag = stringr::str_remove(tag, "/1$")) %>%
    dplyr::mutate(tag = stringr::str_remove(tag, "/2$"))

  # Find common tags
  common_tags <- intersect(file1$tag, file2$tag)

  # Keep only sequences from common tags
  sync_file1 <- file1 %>%
    dplyr::filter(tag %in% common_tags) %>%
    dplyr::arrange(tag) %>%
    dplyr::select(-tag)

  sync_file2 <- file2 %>%
    dplyr::filter(tag %in% common_tags) %>%
    dplyr::arrange(tag) %>%
    dplyr::select(-tag)

  # Write output files if specified
  if (file_format == "fastq" && !is.null(file1_out) && !is.null(file2_out)) {
    microseq::writeFastq(sync_file1, file1_out)
    microseq::writeFastq(sync_file2, file2_out)
  }

  if (file_format == "fasta" && !is.null(file1_out) && !is.null(file2_out)) {
    microseq::writeFasta(sync_file1, file1_out)
    microseq::writeFasta(sync_file2, file2_out)
  }

  # Return results
  if (is.null(file1_out) && is.null(file2_out)) { # Return tibble
    # Add 'sync_file2' as attribute
    attr(sync_file1, "sync_file2") <- sync_file2
    return(sync_file1)
  } else {
    return(invisible(NULL)) # No return when files are written
  }
}
