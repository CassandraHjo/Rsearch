
# TODO: Må finne en måte å gjøre dette for mange samples
# TODO: Skrive dokumentasjon
# TODO: Skrive tester

#' Title
#'
#' @param fastq_file a FASTQ-file with forward reads (R1)
#' @param reverse a FASTQ-file with reverse reads (R2)
#' @param threads number of computational threads to use
#' @param fastqout name of the FASTQ-file with the output.
#' When defined as NULL, no file is written.
#'
#' @return A tibble with merged fastq sequences
#' @export
#'
fastq_mergepairs <- function(fastq_file,
                             reverse,
                             #log_file = NULL,
                             #maxseqlength = 50000,
                             #minseqlength = 32,
                             #sample = NULL,
                             threads = 1,
                             fastqout = NULL){

  if (!file.exists(fastq_file)) stop("Cannot find input file: ", fastq_file)
  if (!file.exists(reverse)) stop("Cannot find reverse file: ", reverse)

  fastq_file <- normalizePath(fastq_file)
  reverse <- normalizePath(reverse)

  if (is.null(fastqout) || fastqout == "-") {
    # Checks if vsearch should write output to file
    message("No filename for output file. No output file will be written.")
    merged_fastq_text <- system2(command = "vsearch",
                                 args = c("--fastq_mergepairs", fastq_file,
                                          "--reverse", reverse,
                                          "--threads", threads,
                                          "--fastqout", "-"),
                                 stdout = TRUE)

    # Transforming output string to tibble
    if (length(merged_fastq_text) %% 4 != 0) {
      stop("FASTQ-text is not valid or incomplete.")
    }

    headers <- merged_fastq_text[seq(1, length(merged_fastq_text), by = 4)]
    sequences <- merged_fastq_text[seq(2, length(merged_fastq_text), by = 4)]
    qualities <- merged_fastq_text[seq(4, length(merged_fastq_text), by = 4)]

    merged_fastq <- tibble::tibble(
      Header = stringr::str_remove(headers, "^@"),
      Sequence = sequences,
      Quality = qualities
    )
    return(merged_fastq)

  } else {
    message("Writing output to file:", fastqout)
    system2("vsearch",
            args = c("--fastq_mergepairs", fastq_file,
                     "--reverse", reverse,
                     "--threads", threads,
                     "--fastqout", fastqout))

    merged_fastq <- microseq::readFastq(fastqout)
    return(merged_fastq)
  }
}



