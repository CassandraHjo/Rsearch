#' fastq_filter
#'
#' @param fastq_file A merged fastq file
#' @param fastq_maxee_rate numeric
#' @param fasta_width numeric
#' @param fastaout output file in FASTA format
#' @param threads number of threads
#'
#' @return a fasta file
#' @export
#'
vs_fastq_filter <- function(fastq_file,
                            fastq_maxee_rate,
                            fasta_width = 0,
                            fastaout = NULL,
                            threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

  # Her skal det skrives kode
}
