test_that("error when wrong output_format", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  output_format <- "fastx"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                               reverse = R2,
                               output_format = output_format),
               "Invalid output_format. Choose from fasta or fastq.")
})

test_that("error when output_format is 'fasta', and fastqout is defined", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  output_format <- "fasta"
  fastqout <- "some_file.fq"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                               reverse = R2,
                               output_format = output_format,
                               fastqout = fastqout),
               "When output_format is defined as 'fasta', 'fastqout' cannot be used. Use 'fastaout' instead.")
})

test_that("error when output_format is 'fastq', and fastaout is defined", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))
  output_format <- "fastq"
  fastaout <- "some_file.fa"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                               reverse = R2,
                               output_format = output_format,
                               fastaout = fastaout),
               "When output_format is defined as 'fastq', 'fastaout' cannot be used. Use 'fastqout' instead.")
})

test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  output_format <- "fastq"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                                   reverse = R2,
                                   output_format = output_format),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  output_format <- "fastq"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                                   reverse = R2,
                                   output_format = output_format),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when input file does not exist", {

  fastq_input <- "some_file.fq"
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  output_format <- "fastq"

  expect_error(vs_fastq_mergepairs(fastq_input = fastq_input,
                                   reverse = reverse,
                                   output_format = output_format),
               paste("Cannot find input FASTQ file:", fastq_input))
})

test_that("error when reverse file does not exist", {

  fastq_input <- test_path("testdata", "sample1", "R2_sample1.fq")
  reverse <- "some_file.fq"
  output_format <- "fastq"

  expect_error(vs_fastq_mergepairs(fastq_input = fastq_input,
                                   reverse = reverse,
                                   output_format = output_format),
               paste("Cannot find reverse FASTQ file:", reverse))
})

test_that("fastq_input and reverse can be merged when files, and results written to fastq file", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastqout <- withr::local_tempfile()
  output_format <- "fastq"
  minlen <- 0
  log_file <- NULL
  threads <- 1

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      output_format = output_format,
                                      minlen = minlen,
                                      log_file = log_file,
                                      threads = threads)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "merged_sample1.fq")))
})

test_that("fastq_input and reverse can be merged when files, and results written to fasta file", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastaout <- withr::local_tempfile()
  output_format <- "fasta"
  minlen <- 0
  log_file <- NULL
  threads <- 1

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastaout = fastaout,
                                      output_format = output_format,
                                      minlen = minlen,
                                      log_file = log_file,
                                      threads = threads)

  expect_null(return_value)

  expect_equal(microseq::readFasta(fastaout),
               microseq::readFasta(test_path("testdata", "output", "merged_sample1.fa")))
})

test_that("fastq_input and reverse can be merged when files, and results given as fastq tibble", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastqout <- NULL
  output_format <- "fastq"
  minlen <- 0
  log_file <- NULL
  threads <- 1

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      output_format = output_format,
                                      minlen = minlen,
                                      log_file = log_file,
                                      threads = threads)

  expect_equal(merged_sample1,
               readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds")))

})

test_that("fastq_input and reverse can be merged when files, and results given as fasta tibble", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastaout <- NULL
  output_format <- "fasta"
  minlen <- 0
  log_file <- NULL
  threads <- 1

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                        reverse = reverse,
                                        fastaout = fastaout,
                                        output_format = output_format,
                                        minlen = minlen,
                                        log_file = log_file,
                                        threads = threads)

  expect_equal(merged_sample1,
               readRDS(test_path("testdata", "output", "merged_sample1_fasta_files.rds")))

})

test_that("fastq_input and reverse can be merged when tibbles, and results written to fastq file", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))
  fastqout <- withr::local_tempfile()
  output_format <- "fastq"
  minlen <- 0
  log_file <- NULL
  threads <- 1

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      output_format = output_format,
                                      minlen = minlen,
                                      log_file = log_file,
                                      threads = threads)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "merged_sample1.fq")))
})

test_that("fastq_input and reverse can be merged when tibbles, and results given as tibble", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))
  fastqout <- NULL
  output_format <- "fastq"
  minlen <- 0
  log_file <- NULL
  threads <- 1

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                        reverse = reverse,
                                        fastqout = fastqout,
                                        output_format = output_format,
                                        minlen = minlen,
                                        log_file = log_file,
                                        threads = threads)

  expect_equal(merged_sample1,
               readRDS(test_path("testdata", "output", "merged_sample1_fastq_tibbles.rds")))
})

test_that("log file exists when specified", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastqout <- withr::local_tempfile()
  output_format <- "fastq"
  minlen <- 0
  log_file <- withr::local_tempfile()
  threads <- 1

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      output_format = output_format,
                                      minlen = minlen,
                                      log_file = log_file,
                                      threads = threads)
  expect_null(return_value)

  expect_true(file.exists(log_file))

})




