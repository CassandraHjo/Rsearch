


fastq_mergepairs <- function(fastq_file,
                             reverse,
                             #log_file = NULL,
                             #maxseqlength = 50000,
                             #minseqlength = 32,
                             #sample = NULL,
                             threads = 1,
                             fastqout = NULL){

  if (is.null(fastqout) || fastqout == "-") {
    # Checks if vsearch should write output to file
    print("No filename for output file. No output file will be written.")
    merged_fastq <- system2(command = "vsearch",
                            args = paste("vsearch",
                                         "--fastq_mergepairs", fastq_file,
                                         "--reverse", reverse,
                                         "--threads", threads,
                                         "--fastqout", "-"),
                            stdout = TRUE)
    # Må enten returneres som en stor tekst, eller bruke microseq-pakken for å oversette til tabell i R
    # return()

  } else {
    paste("Writing output to file:", fastqout)
    cmd <- paste("vsearch",
                 "--fastq_mergepairs", fastq_file,
                 "--reverse", reverse,
                 "--threads", threads,
                 "--fastqout", fastqout)
    system(cmd)
    merged_fastq <- readFastq(fastqout) # fra microseq
    return(merged_fastq)
  }
}


