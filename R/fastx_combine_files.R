#' Combine FASTA/FASTQ files in a directory into a single file or object
#'
#' @description \code{fastx_combine_files} combines all FASTA or FASTQ files
#' within a specified directory into a single output file or a tibble object.
#'
#' @param files_dir A character string specifying the path to the directory
#' containing the files to be combined. Files must be uncompressed.
#' @param output_file A character string specifying the name of the output file.
#' If \code{NULL} (default), the combined data is returned as a FASTA/FASTQ
#' object depending on \code{file_format} instead of being written to a file.
#' @param file_ext File extension of the files to be combined. Defaults to
#' \code{".fq"}.
#' @param file_format Format of files to be combined and the desired output
#' format: either \code{"fasta"} or \code{"fastq"} (default). See \emph{Details}.
#'
#' @details
#' \code{files_dir} must contain uncompressed FASTA or FASTQ files matching the
#' specified \code{file_ext}.
#'
#' All files with the specified \code{file_ext} in \code{files_dir} are
#' concatenated into a single output file or tibble.
#'
#' A FASTA object is a tibble containing the columns \code{Header} and
#' \code{Sequence}.
#' A FASTQ object is a tibble containing the columns \code{Header},
#' \code{Sequence}, and \code{Quality}.
#'
#' If \code{output_file} is specified, the combined sequences are written to
#' this file in the format specified by \code{file_format}.
#'
#' If \code{output_file} is \code{NULL}, the combined sequences are returned as
#' a tibble in the format specified by \code{file_format}, and no file is
#' written.
#'
#' @return A tibble or \code{NULL}.
#'
#' If \code{output_file} is specified, the combined sequences are written to the
#' specified file.
#'
#' If \code{output_file} is \code{NULL}, the combined sequences are returned as
#' a tibble in the format specified by \code{file_format}.
#'
#' @examples
#' \dontrun{
#' # Define arguments
#' files_dir <- file.path(path.package("Rsearch"), "extdata")
#' output_file <- NULL
#' file_ext <- ".fq"
#' file_format <- "fastq"
#'
#' # Combine files and return tibble object
#' combined_files <- fastx_combine_files(files_dir = files_dir,
#'                                       output_file = output_file,
#'                                       file_ext = file_ext,
#'                                       file_format = file_format)
#'
#' # Combine files and write to an output file
#' fastx_combine_files(files_dir = files_dir,
#'                     output_file = "combined.fastq",
#'                     file_ext = file_ext,
#'                     file_format = file_format)
#' }
#'
#' @aliases fastx_combine_files fastq_combine_files fasta_combine_files
#' @export
#'
fastx_combine_files <- function(files_dir,
                                output_file = NULL,
                                file_ext = ".fq",
                                file_format = "fastq") {

  # Check if input directory exists
  if (!dir.exists(files_dir)) {
    stop("Directory does not exist: ", files_dir)
  }

  # Validate file_format
  if (!file_format %in% c("fasta", "fastq")) {
    stop("Invalid file_format. Choose from fasta or fastq.")
  }

  # Create empty vector for collecting temporary files
  temp_files <- c()

  # Set up cleanup of temporary files
  on.exit({
    if (length(temp_files) > 0) {
      file.remove(temp_files)
    }
  }, add = TRUE)

  # Find all files in the folder
  comb_files <- list.files(files_dir,
                           pattern = paste0("\\", file_ext, "$"),
                           full.names = TRUE)

  # Check if any files are found
  if (length(comb_files) == 0) {
    stop("No ", file_ext, " files found in the specified folder: ", files_dir)
  }

  # Handle output file if NULL
  if (is.null(output_file)) {
    out_file <- tempfile(pattern = "all_fasta", fileext = file_ext)
    temp_files <- c(temp_files, out_file)
  } else {
    out_file <- output_file
  }

  # Full path to the combined output file
  out_file <- file.path(out_file)

  # Combine content of all files
  file.append(out_file, comb_files)

  # Create tibble if output_file not specified
  if (is.null(output_file)) {
    # Create FASTA object
    if (file_format == "fasta") {
      all_seq_tbl <- microseq::readFasta(out_file)
    }
    # Create FASTQ object
    if (file_format == "fastq") {
      all_seq_tbl <- microseq::readFastq(out_file)
    }
  }

  # Return results
  if (is.null(output_file)) { # Return tibble
    return(all_seq_tbl)
  } else {
    return(invisible(NULL)) # No return when output file is written
  }
}
