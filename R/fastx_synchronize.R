#' Synchronize FASTA and FASTQ files or objects
#'
#' @description \code{fastx_synchronize} synchronizes sequences between two
#' FASTA/FASTQ files or objects by retaining only the common sequences present
#' in both.
#'
#' @param file1 A FASTA/FASTQ file path or a FASTA/FASTQ tibble object. See
#' \emph{Details}.
#' @param file2 A FASTA/FASTQ file path or a FASTA/FASTQ tibble object. See
#' \emph{Details}.
#' @param file_format Format of the input (\code{file1} and \code{file2})
#' and the desired output format: \code{"fasta"} or \code{"fastq"} (default).
#' This determines the format for both outputs.
#' @param file1_out Name of the output file for synchronized reads from
#' \code{file1}. The file is in either FASTA or FASTQ format, depending on
#' \code{file_format}. If \code{NULL} (default), no sequences are written to a
#' file. See \emph{Details}.
#' @param file2_out Name of the output file for synchronized reads from
#' \code{file2}. The file is in either FASTA or FASTQ format, depending on
#' \code{file_format}. If \code{NULL} (default), no sequences are written to a
#' file. See \emph{Details}.
#'
#' @details
#' \code{file1} and \code{file2} can either be paths to FASTA/FASTQ files or
#' tibble objects containing the sequences.
#' FASTA objects are tibbles that contain the columns \code{Header} and
#' \code{Sequence}.
#' FASTQ objects are tibbles that contain the columns \code{Header},
#' \code{Sequence}, and \code{Quality}.
#'
#' Sequence IDs in the \code{Header} fields must be identical for each read pair
#' in both \code{file1} and \code{file2} for synchronization to work correctly.
#'
#' If \code{file1_out} and \code{file2_out} are specified, the synchronized
#' sequences are written to these files in the format specified by
#' \code{file_format}.
#'
#' If \code{file1_out} and \code{file2_out} are \code{NULL}, the function
#' returns a FASTA/FASTQ object containing synchronized reads from \code{file1}.
#' The synchronized reads from \code{file2} are included as an attribute named
#' \code{"reverse"} in the returned tibble.
#'
#' Both \code{file1_out} and \code{file2_out} must either be \code{NULL} or both
#' must be character strings specifying the file paths.
#'
#' @return A tibble or \code{NULL}.
#'
#' If both \code{file1_out} and \code{file2_out} are \code{NULL}, a tibble
#' containing the synchronized reads from \code{file1} is returned. The
#' synchronized reads from \code{file2} are accessible via the \code{"reverse"}
#' attribute of the returned tibble.
#'
#' If both \code{file1_out} and \code{file2_out} are specified, the synchronized
#' sequences are written to the specified output files, and no tibble is
#' returned.
#'
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' file1 <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                    "R1_sample1_small.fq")
#' file2 <- file.path(file.path(path.package("Rsearch"), "extdata"),
#'                    "R2_sample1_small.fq")
#' file_format <- "fastq"
#' file1_out <- NULL
#' file2_out <- NULL
#'
#' # Synchronize files and return as a tibble
#' sync_seqs <- fastx_synchronize(file1 = file1,
#'                                file2 = file2,
#'                                file_format = file_format,
#'                                file1_out = file1_out,
#'                                file2_out = file2_out)
#'
#' # Extract tibbles with synchronized sequences
#' R1_sync <- sync_seqs
#' R2_sync <- attr(sync_seqs, "reverse")
#'
#' # Synchronize files and write to output files
#' fastx_synchronize(file1 = file1,
#'                   file2 = file2,
#'                   file_format = file_format,
#'                   file1_out = "synchronized_R1.fastq",
#'                   file2_out = "synchronized_R2.fastq")
#' }
#'
#' @aliases fastx_synchronize fastq_synchronize fasta_synchronize
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
        stop("file1 FASTQ object must contain columns: Header, Sequence, Quality")
      }
    }
    if (file_format == "fasta") {
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(file1))) {
        stop("file1 FASTA object must contain columns: Header and Sequence")
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
        stop("file2 FASTQ object must contain columns: Header, Sequence, Quality")
      }
    }
    if (file_format == "fasta") {
      required_cols <- c("Header", "Sequence")
      if (!all(required_cols %in% colnames(file2))) {
        stop("file2 FASTA object must contain columns: Header and Sequence")
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
  file1 <- file1 |>
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/1$")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/2$"))

  file2 <- file2 |>
    dplyr::mutate(tag = stringr::str_extract(Header, "^\\S+")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/1$")) |>
    dplyr::mutate(tag = stringr::str_remove(tag, "/2$"))

  # Find common tags
  common_tags <- intersect(file1$tag, file2$tag)

  # Keep only sequences from common tags
  sync_file1 <- file1 |>
    dplyr::filter(tag %in% common_tags) |>
    dplyr::arrange(tag) |>
    dplyr::select(-tag)

  sync_file2 <- file2 |>
    dplyr::filter(tag %in% common_tags) |>
    dplyr::arrange(tag) |>
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
    # Add sync_file2 as "reverse" attribute
    attr(sync_file1, "reverse") <- sync_file2
    return(sync_file1)
  } else {
    return(invisible(NULL)) # No return when files are written
  }
}
