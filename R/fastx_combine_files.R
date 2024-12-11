#' Combine content files in one directory
#'
#' @description Combines all files of given type in given directory into one file/object.
#'
#' @param files_dir Path to directory with files to combine.
#' @param output_file Name of output file or \code{NULL}. If not specified, the function will only return a FASTA/FASTQ object depending on \code{file_format}.
#' @param file_ext The file extension for the files you want to combine. Must be written with a "." in front of the extension.
#' @param file_format Format of files you want to combine, and desired output format: \code{"fasta"} or \code{"fastq"}.
#'
#' @return If \code{output_file} is not specified, a tibble containing the combined reads in the format specified in \code{file_format} is returned. If \code{output_file} is specified nothing is returned.
#'
#' @export
fastx_combine_files <- function(files_dir, output_file = NULL, file_ext = ".fa", file_format = "fasta") {

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
  comb_files <- list.files(files_dir, pattern = paste0("\\", file_ext, "$"), full.names = TRUE)

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
