test_that("error when wrong output_format", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
  output_format <- "fastx"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                                   reverse = R2,
                                   output_format = output_format),
               "Invalid output_format. Choose from fasta or fastq.")
})

test_that("error when output_format is 'fasta', and fastqout is defined", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
  output_format <- "fasta"
  fastqout <- "some_file.fq"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                                   reverse = R2,
                                   output_format = output_format,
                                   fastqout = fastqout),
               "When output_format is 'fasta', 'fastqout' cannot be used. Use 'fastaout' instead.")
})

test_that("error when output_format is 'fastq', and fastaout is defined", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))
  output_format <- "fastq"
  fastaout <- "some_file.fa"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                                   reverse = R2,
                                   output_format = output_format,
                                   fastaout = fastaout),
               "When output_format is 'fastq', 'fastaout' cannot be used. Use 'fastqout' instead.")
})

test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds"))

  output_format <- "fastq"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                                   reverse = R2,
                                   output_format = output_format),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))

  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds")) |>
    dplyr::select(-Header)

  output_format <- "fastq"

  expect_error(vs_fastq_mergepairs(fastq_input = R1,
                                   reverse = R2,
                                   output_format = output_format),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when input file does not exist", {

  fastq_input <- "some_file.fq"
  reverse <- test_path("testdata", "R2.fastq")
  output_format <- "fastq"

  expect_error(vs_fastq_mergepairs(fastq_input = fastq_input,
                                   reverse = reverse,
                                   output_format = output_format),
               paste("Cannot find input FASTQ file:", fastq_input))
})

test_that("error when reverse file does not exist", {

  fastq_input <- test_path("testdata", "R2.fastq")
  reverse <- "some_file.fq"
  output_format <- "fastq"

  expect_error(vs_fastq_mergepairs(fastq_input = fastq_input,
                                   reverse = reverse,
                                   output_format = output_format),
               paste("Cannot find reverse FASTQ file:", reverse))
})

test_that("fastq_input and reverse can be merged when files, and results written to fastq file", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  fastqout <- withr::local_tempfile()
  output_format <- "fastq"

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      output_format = output_format)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "merged_fq_files.fastq")))
})

test_that("fastq_input and reverse can be merged when files, and results written to fasta file", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  fastaout <- withr::local_tempfile()
  output_format <- "fasta"

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastaout = fastaout,
                                      output_format = output_format)

  expect_null(return_value)

  expect_equal(microseq::readFasta(fastaout),
               microseq::readFasta(test_path("testdata", "output", "merged_fq_files.fasta")))
})

test_that("fastq_input and reverse can be merged when files, and results given as fastq tibble", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  output_format <- "fastq"

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                        reverse = reverse,
                                        output_format = output_format)

  expect_equal(merged_sample1,
               readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds")))
})

test_that("fastq_input and reverse can be merged when files, and results given as fasta tibble", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  fastaout <- NULL
  output_format <- "fasta"

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                        reverse = reverse,
                                        fastaout = fastaout,
                                        output_format = output_format)

  expect_equal(merged_sample1,
               readRDS(test_path("testdata", "output", "merged_fq_files_fa_tibble.rds")))
})

test_that("fastq_input and reverse can be merged when tibbles, and results written to fastq file", {

  fastq_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  reverse <- microseq::readFastq(test_path("testdata", "R2.fastq"))
  fastqout <- withr::local_tempfile()
  output_format <- "fastq"

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      output_format = output_format)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "merged_fq_files.fastq")))
})

test_that("fastq_input and reverse can be merged when tibbles, and results given as tibble", {

  fastq_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  reverse <- microseq::readFastq(test_path("testdata", "R2.fastq"))
  output_format <- "fastq"

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                        reverse = reverse,
                                        output_format = output_format)

  expect_equal(merged_sample1,
               readRDS(test_path("testdata", "output", "merged_fq_tibbles_fq_tibble.rds")))
})

test_that("log file exists when specified", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  fastqout <- withr::local_tempfile()
  output_format <- "fastq"
  log_file <- withr::local_tempfile()

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      output_format = output_format,
                                      log_file = log_file)
  expect_null(return_value)

  expect_true(file.exists(log_file))
})

test_that("fastq_input and reverse can be merged when files, and results given as fastq tibble with vsearch_options", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")
  output_format <- "fastq"
  sample <- "sample1"
  vsearch_options <- c("--relabel", "OTU")

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                        reverse = reverse,
                                        output_format = output_format,
                                        sample = sample,
                                        vsearch_options = vsearch_options)

  expect_equal(merged_sample1,
               readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble_vs.rds")))
})
