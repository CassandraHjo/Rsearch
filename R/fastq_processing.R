
# TODO: Må finne en måte å gjøre dette for mange samples

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
                                     #"--threads", threads,
                                     "--fastqout", "-"),
                            stdout = TRUE)

    # Transforming output string to tibble
    if (length(merged_fastq_text) %% 4 != 0) {
      stop("FASTQ-text is not valid or incomplete.")
    }

    headers <- merged_fastq_text[seq(1, length(merged_fastq_text), by = 4)]
    sequences <- fastq_text[seq(2, length(fastq_text), by = 4)]
    qualities <- fastq_text[seq(4, length(fastq_text), by = 4)]

    merged_fastq <- tibble(
      Header = str_remove(headers, "^@"),
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

    merged_fastq <- readFastq(fastqout) # fra microseq
    return(merged_fastq)
  }
}



