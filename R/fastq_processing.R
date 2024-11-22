#' Parse output text from read merging to tibble
#'
#' @param output_text string of output from running fastq_mergepairs in vsearch
#'
#' @return table with merging metrics
#' @export
#'
parse_merge_pairs_output <- function(output_text) {

  # Ekstraher linjer og verdier fra output_text
  pairs <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Pairs$"), "\\d+"))
  merged <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Merged \\("), "\\d+"))
  merged_percent <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Merged \\("), "\\d+\\.\\d+"))
  not_merged <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Not merged"), "\\d+"))
  not_merged_percent <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Not merged"), "\\d+\\.\\d+"))

  too_many_diff <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "too many differences"), "\\d+"))
  align_score_low <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "alignment score too low"), "\\d+"))

  mean_read_length <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Mean read length"), "\\d+\\.\\d+"))

  mean_frag_length <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Mean fragment length"), "\\d+\\.\\d+"))
  sd_frag_length <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Standard deviation of fragment length"), "\\d+\\.\\d+"))

  mean_expected_error_forward <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Mean expected error in forward sequences"), "\\d+\\.\\d+"))
  mean_expected_error_reverse <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Mean expected error in reverse sequences"), "\\d+\\.\\d+"))
  mean_expected_error_merged <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Mean expected error in merged sequences"), "\\d+\\.\\d+"))

  mean_observed_error_forward <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Mean observed errors in merged region of forward sequences"), "\\d+\\.\\d+"))
  mean_observed_error_reverse <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Mean observed errors in merged region of reverse sequences"), "\\d+\\.\\d+"))
  mean_observed_error_merged <- as.numeric(stringr::str_extract(stringr::str_subset(output_text, "Mean observed errors in merged region$"), "\\d+\\.\\d+"))

  metrics <- tibble::tibble(
    Metric = c(
      "Total Pairs",
      "Merged Pairs",
      "Merged Percentage",
      "Not Merged Pairs",
      "Not Merged Percentage",
      "Too Many Differences",
      "Alignment Score Too Low",
      "Mean Read Length",
      "Mean Fragment Length",
      "SD Fragment Length",
      "Mean Expected Error Forward",
      "Mean Expected Error Reverse",
      "Mean Expected Error Merged",
      "Mean Observed Error Forward",
      "Mean Observed Error Reverse",
      "Mean Observed Error Merged"
    ),
    Value = c(
      pairs,
      merged,
      merged_percent,
      not_merged,
      not_merged_percent,
      too_many_diff,
      align_score_low,
      mean_read_length,
      mean_frag_length,
      sd_frag_length,
      mean_expected_error_forward,
      mean_expected_error_reverse,
      mean_expected_error_merged,
      mean_observed_error_forward,
      mean_observed_error_reverse,
      mean_observed_error_merged
    )
  )

  return(metrics)
}


# TODO: Må finne en måte å gjøre dette for mange samples
# TODO: Skrive dokumentasjon
# TODO: Skrive tester

#' Merge FASTQ files
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
    output_text <- system2("vsearch",
                           args = c("--fastq_mergepairs", fastq_file,
                                    "--reverse", reverse,
                                    "--threads", threads,
                                    "--fastqout", fastqout),
                           stdout = TRUE,
                           stderr = TRUE)

    # output statistics in table
    metrics <- parse_merge_pairs_output(output_text)

    # merged reads in table
    merged_fastq <- microseq::readFastq(fastqout)

    return(list(metrics = metrics, merged_fastq = merged_fastq))
  }
}
