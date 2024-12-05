#' Combine FASTA files
#'
#' @description Combines all FASTA-files in given directory into one FASTA file or object.
#'
#' @param fasta_dir Path to directory with FASTA files.
#' @param output_file Name of output FASTA-file or \code{NULL}. If not specified, the function will only return a FASTA object (a tibble with the columns \code{Header} and \code{Sequence}).
#'
#' @return A string with path to combined FASTA file. If \code{output_file} is unspecified, a FASTA object will be returned.
#'
#' @return A FASTA object with the columns \code{Header} and \code{Sequence}. If \code{output_file} is specified the path to the resulting FASTA file will be an attribute (\code{file_path})of the FASTA object.
#' This attribute can be accessed by running \code{attributes(all_fasta_tbl)$file_path} or \code{attr(all_fasta_tbl, "file_path")}.
#'
#' @export
combine_fasta_files <- function(fasta_dir, output_file = NULL) {

  # Check if input directory exists
  if (!dir.exists(fasta_dir)) {
    stop("Directory does not exist: ", fasta_dir)
  }

  # Create empty vector for collecting temporary files
  temp_files <- c()

  # Find all .fa files in the folder
  fa_files <- list.files(fasta_dir, pattern = "\\.fa$", full.names = TRUE)

  # Check if any .fa files are found
  if (length(fa_files) == 0) {
    stop("No .fa files found in the specified folder: ", fasta_dir)
  }

  # Read and combine the content of all .fa files
  fasta_content <- lapply(fa_files, readLines)
  combined_fasta <- unlist(fasta_content)

  if (!is.null(output_file)){

    # Full path to the combined output file
    all_fasta <- file.path(output_file)

    # Remove all.fasta if it already exists
    if (file.exists(all_fasta)) {
      file.remove(all_fasta)
    }
    # Write the combined content to all_fasta
    writeLines(combined_fasta, all_fasta)
    message("Sequences written to: ", all_fasta)

    # Create FASTA object
    all_fasta_tbl <- microseq::readFasta(all_fasta)

    # Create file path attribute to FASTA object
    attr(all_fasta_tbl, "file_path") <- all_fasta
  } else {
    temp_file_all_fasta <- tempfile(pattern = "all_fasta", fileext = ".fa")
    temp_files <- c(temp_files, temp_file_all_fasta)
    writeLines(combined_fasta, temp_file_all_fasta)

    all_fasta_tbl <- microseq::readFasta(temp_file_all_fasta)
  }

  # Cleanup temporary files
  if (length(temp_files) > 0) {
    on.exit(
      file.remove(temp_files),
      add = TRUE)
  }

  # Return path to combined file
  return(all_fasta_tbl)
}
