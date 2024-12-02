#' Combine fasta files into one file
#'
#' @param folder_path path to where the fasta files are stored
#' @param output_file name of output file
#'
#' @export
combine_fasta_files <- function(folder_path, output_file) {

  # Check if input directory exists
  if (!dir.exists(folder_path)) {
    stop("Directory does not exist: ", folder_path)
  }

  # Full path to the combined output file
  all_fasta <- file.path(folder_path, output_file)

  # Remove all.fasta if it already exists
  if (file.exists(all_fasta)) {
    file.remove(all_fasta)
  }

  # Find all .fa files in the folder
  fa_files <- list.files(folder_path, pattern = "\\.fa$", full.names = TRUE)

  # Check if any .fa files are found
  if (length(fa_files) == 0) {
    stop("No .fa files found in the specified folder: ", folder_path)
  }

  # Read and combine the content of all .fa files
  fasta_content <- lapply(fa_files, readLines)
  combined_fasta <- unlist(fasta_content)

  # Write the combined content to all.fasta
  writeLines(combined_fasta, all_fasta)

  # Return path to combined file
  return(all_fasta)
}
