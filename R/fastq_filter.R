#' Filter sequences
#'
#' @description Trim and/or filter sequences in the given FASTQ file
#'
#' @param fastq_file a FASTQ file with reads
#' @param fastaout name of the FASTA-file with the output or NULL, see Details.
#' @param fastq_maxee_rate threshold for average expected error. Value ranging form 0.0 to 1.0. See Details.
#' @param fasta_width number of characters in the width of sequences in the output FASTA file. See Detalis.
#' @param threads number of computational threads to be used by vsearch.
#' @param log_file name of the log file with messages from running vsearch or NULL, see Details.
#'
#' @details The reads in the input FASTQ-file (\code{fastq_file}) are filtered based on ..., using vsearch.
#'
#' If \code{fastaout} is specified, the remaining sequences after trimming/filtering are output to this file in FASTA-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTA-object, i.e. a tibble with
#' columns \code{Header} and \code{Sequence}.
#'
#' Sequences with an average expected error greater than the specified \code{fastq_maxee_rate} are discarded.
#' For a given sequence, the average expected error is the sum of error probabilities for all the positions in the sequence,
#' divided by the length of the sequence.
#'
#' FASTA files produced by vsearch are wrapped (sequences are written on lines of integer nucleotides).
#' \code{fasta_width} is by default set to zero to eliminate the wrapping.
#'
#' #' If \code{log_file} is specified, the messages are output to this file.
#' If unspecified (\code{NULL}) no log file is written.
#'
#' @return a tibble with FASTA sequences
#' @export
#'
vs_fastq_filter <- function(fastq_file,
                            fastaout = NULL,
                            fastq_maxee_rate,
                            fasta_width = 0,
                            threads = 1){

  # Check if vsearch is available
  vsearch_executable <- options("Rsearch.vsearch_executable")[[1]]
  vsearch_available(vsearch_executable)

}

#' Merge read-pairs
#'
#' @description Merging read-pairs with overlapping regions.
#'
#' @param fastq_file a FASTQ-file with forward reads (R1).
#' @param reverse a FASTQ-file with reverse reads (R2).
#' @param fastqout name of the FASTQ-file with the output or NULL, see Details.
#' @param log_file name of the log file with messages from running vsearch or NULL, see Details.
#' @param threads number of computational threads to be used by vsearch.
#'
#' @details The read-pairs in the input FASTQ-files (\code{fastq_file} and \code{reverse})
#' are merged if they have sufficient overlap, using vsearch.
#'
#' If \code{fastqout} is specified, the merged reads are output to this file in FASTQ-format.
#' If unspecified (\code{NULL}) the result is returned as a FASTQ-object, i.e. a tibble with
#' columns \code{Header}, \code{Sequence} and \code{Quality}.
#'
#' #' If \code{log_file} is specified, the messages are output to this file.
#' If unspecified (\code{NULL}) no log file is written.
#'
#' @return A list with two tibbles, one with merged FASTQ sequences and one with merging statistics.
#'
#' @export
#'
