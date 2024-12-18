test_that("error when wrong output_format", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_fastq_filter(fastq_input = R1,
                               reverse = R2,
                               output_format = "fastx"),
               "Invalid output_format. Choose from fasta or fastq.")
})

test_that("error when output_format is 'fasta', and fastqout and fastqout_rev are defined", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  output_format <- "fasta"
  fastqout <- "some_file.fq"
  fastqout_rev <- "some_other_file.fq"

  expect_error(vs_fastq_filter(fastq_input = R1,
                               reverse = R2,
                               output_format = output_format,
                               fastqout = fastqout,
                               fastqout_rev = fastqout_rev),
               "When output_format is defined as 'fasta', 'fastqout' and 'fastqout_rev' cannot be used. Use 'fastaout' and 'fastaout_rev' instead.")
})

test_that("error when output_format is 'fastq', and fastaout and fastaout_rev are defined", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  output_format <- "fastq"
  fastaout <- "some_file.fa"
  fastaout_rev <- "some_other_file.fa"

  expect_error(vs_fastq_filter(fastq_input = R1,
                               reverse = R2,
                               output_format = output_format,
                               fastaout = fastaout,
                               fastaout_rev = fastaout_rev),
               "When output_format is defined as 'fastq', 'fastaout' and 'fastaout_rev' cannot be used. Use 'fastqout' and 'fastqout_rev' instead.")
})

test_that("error when reverse is specified, but output files are not both NULL or both character strings", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  output_format <- "fasta"
  fastaout <- "some_file.fa"
  fastaout_rev <- NULL

  expect_error(vs_fastq_filter(fastq_input = R1,
                               reverse = R2,
                               output_format = output_format,
                               fastaout = fastaout,
                               fastaout_rev = fastaout_rev),
               "When 'reverse' is specified and output_format is 'fasta', both 'fastaout' and 'fastaout_rev' must be NULL or both specified as character strings.")
})

test_that("error when reverse is specified, but output files are not both NULL or both character strings", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  output_format <- "fastq"
  fastaout <- "some_file.fa"
  fastaout_rev <- NULL

  expect_error(vs_fastq_filter(fastq_input = R1,
                               reverse = R2,
                               output_format = output_format,
                               fastqout = fastaout,
                               fastqout_rev = fastaout_rev),
               "When 'reverse' is specified and output_format is 'fastq', both 'fastqout' and 'fastqout_rev' must be NULL or both specified as character strings.")
})

test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_fastq_filter(fastq_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  expect_error(vs_fastq_filter(fastq_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when input file does not exist", {

  fastq_input <- "some_file.fq"
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")

  expect_error(vs_fastq_filter(fastq_input = fastq_input,
                             reverse = reverse),
               paste("Cannot find input FASTQ file:", fastq_input))
})

test_that("error when reverse file does not exist", {

  fastq_input <- test_path("testdata", "sample1", "R2_sample1.fq")
  reverse <- "some_file.fq"

  expect_error(vs_fastq_filter(fastq_input = fastq_input,
                             reverse = reverse),
               paste("Cannot find reverse FASTQ file:", reverse))
})

test_that("filter fastq sequences from two files, and return two fastq files", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastqout <- withr::local_tempfile()
  fastqout_rev <- withr::local_tempfile()
  output_format <- "fastq"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  log_file <- withr::local_tempfile()
  threads <- 1

  return_value <- vs_fastq_filter(fastq_input = fastq_input,
                                  reverse = reverse,
                                  fastqout = fastqout,
                                  fastqout_rev = fastqout_rev,
                                  output_format = output_format,
                                  fastq_maxee_rate = fastq_maxee_rate,
                                  minlen = minlen,
                                  log_file = log_file,
                                  threads = threads)

  expect_null(return_value)

  expect_true(file.exists(log_file))

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "R1_filtered_sample1.fq")))

  expect_equal(microseq::readFastq(fastqout_rev),
               microseq::readFastq(test_path("testdata", "output", "R2_filtered_sample1.fq")))
})

test_that("filter fastq sequences from two files, and return fastq tibble", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  output_format <- "fastq"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  threads <- 1

  filtered_sample1 <- vs_fastq_filter(fastq_input = fastq_input,
                                      reverse = reverse,
                                      output_format = output_format,
                                      fastq_maxee_rate = fastq_maxee_rate,
                                      minlen = minlen,
                                      threads = threads)

  expect_equal(filtered_sample1,
               readRDS(test_path("testdata", "output", "filtered_sample1_fastq_files.rds")))

})

test_that("filter fastq sequences from two files, and return two fasta files", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastaout <- withr::local_tempfile()
  fastaout_rev <- withr::local_tempfile()
  output_format <- "fasta"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  threads <- 1

  return_value <- vs_fastq_filter(fastq_input = fastq_input,
                                  reverse = reverse,
                                  fastaout = fastaout,
                                  fastaout_rev = fastaout_rev,
                                  output_format = output_format,
                                  fastq_maxee_rate = fastq_maxee_rate,
                                  minlen = minlen,
                                  threads = threads)

  expect_null(return_value)

  expect_equal(microseq::readFasta(fastaout),
               microseq::readFasta(test_path("testdata", "output", "R1_filtered_sample1.fa")))

  expect_equal(microseq::readFasta(fastaout_rev),
               microseq::readFasta(test_path("testdata", "output", "R2_filtered_sample1.fa")))
})

test_that("filter fastq sequences from two files, and return fasta tibble", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  output_format <- "fasta"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  threads <- 1

  filtered_sample1 <- vs_fastq_filter(fastq_input = fastq_input,
                                      reverse = reverse,
                                      output_format = output_format,
                                      fastq_maxee_rate = fastq_maxee_rate,
                                      minlen = minlen,
                                      threads = threads)

  expect_equal(filtered_sample1,
               readRDS(test_path("testdata", "output", "filtered_sample1_fasta_files.rds")))
})

test_that("filter fastq sequences from two tibbles, and return two fastq files", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))
  fastqout <- withr::local_tempfile()
  fastqout_rev <- withr::local_tempfile()
  output_format <- "fastq"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  threads <- 1

  return_value <- vs_fastq_filter(fastq_input = fastq_input,
                                  reverse = reverse,
                                  fastqout = fastqout,
                                  fastqout_rev = fastqout_rev,
                                  output_format = output_format,
                                  fastq_maxee_rate = fastq_maxee_rate,
                                  minlen = minlen,
                                  threads = threads)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "R1_filtered_sample1.fq")))

  expect_equal(microseq::readFastq(fastqout_rev),
               microseq::readFastq(test_path("testdata", "output", "R2_filtered_sample1.fq")))
})

test_that("filter fastq sequences from two tibbles, and return fastq tibble", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))
  output_format <- "fastq"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  threads <- 1

  filtered_sample1 <- vs_fastq_filter(fastq_input = fastq_input,
                                      reverse = reverse,
                                      output_format = output_format,
                                      fastq_maxee_rate = fastq_maxee_rate,
                                      minlen = minlen,
                                      threads = threads)


  expect_equal(filtered_sample1,
               readRDS(test_path("testdata", "output", "filtered_sample1_fastq_tibbles.rds")))
})

test_that("filter fastq sequences from two tibbles, and return two fasta files", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))
  fastaout <- withr::local_tempfile()
  fastaout_rev <- withr::local_tempfile()
  output_format <- "fasta"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  threads <- 1

  return_value <- vs_fastq_filter(fastq_input = fastq_input,
                                  reverse = reverse,
                                  fastaout = fastaout,
                                  fastaout_rev = fastaout_rev,
                                  output_format = output_format,
                                  fastq_maxee_rate = fastq_maxee_rate,
                                  minlen = minlen,
                                  threads = threads)

  expect_null(return_value)

  expect_equal(microseq::readFasta(fastaout),
               microseq::readFasta(test_path("testdata", "output", "R1_filtered_sample1.fa")))

  expect_equal(microseq::readFasta(fastaout_rev),
               microseq::readFasta(test_path("testdata", "output", "R2_filtered_sample1.fa")))
})

test_that("filter fastq sequences from two tibbles, and return fasta tibble", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))
  output_format <- "fasta"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  threads <- 1

  filtered_sample1 <- vs_fastq_filter(fastq_input = fastq_input,
                                      reverse = reverse,
                                      output_format = output_format,
                                      fastq_maxee_rate = fastq_maxee_rate,
                                      minlen = minlen,
                                      threads = threads)


  expect_equal(filtered_sample1,
               readRDS(test_path("testdata", "output", "filtered_sample1_fasta_tibbles.rds")))
})

test_that("filter fastq sequences from one tibble, and return fasta tibble", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  output_format <- "fasta"
  fastq_maxee_rate <- 0.01
  minlen <- 0
  threads <- 1

  filtered_sample1 <- vs_fastq_filter(fastq_input = fastq_input,
                                      output_format = output_format,
                                      fastq_maxee_rate = fastq_maxee_rate,
                                      minlen = minlen,
                                      threads = threads)


  expect_equal(filtered_sample1,
               readRDS(test_path("testdata", "output", "filtered_sample1_R1_fasta_tibbles.rds")))
})
