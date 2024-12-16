test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_fastq_mergepairs(fastq_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  expect_error(vs_fastq_mergepairs(fastq_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when input file does not exist", {

  fastq_input <- "some_file.fq"
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")

  expect_error(vs_fastq_mergepairs(fastq_input = fastq_input,
                                   reverse = reverse),
               paste("Cannot find input FASTQ file:", fastq_input))
})

test_that("error when reverse file does not exist", {

  fastq_input <- test_path("testdata", "sample1", "R2_sample1.fq")
  reverse <- "some_file.fq"

  expect_error(vs_fastq_mergepairs(fastq_input = fastq_input,
                                   reverse = reverse),
               paste("Cannot find reverse FASTQ file:", reverse))
})

test_that("fastq_input and reverse can be merged when files, and results written to file", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastqout <- withr::local_tempfile()
  minlen <- 0
  log_file <- NULL
  threads <- 1

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      minlen = minlen,
                                      log_file = log_file,
                                      threads = threads)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastqout),
               microseq::readFastq(test_path("testdata", "output", "merged_sample1.fq")))
})

test_that("fastq_input and reverse can be merged when files, and results given as tibble", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")
  fastqout <- NULL
  minlen <- 0
  log_file <- NULL
  threads <- 1

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      minlen = minlen,
                                      log_file = log_file,
                                      threads = threads)

  expect_equal(merged_sample1,
               readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds")))

})

test_that("fastq_input and reverse can be merged when tibbles, and results written to file", {

  fastq_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  reverse <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))
  fastqout <- withr::local_tempfile()
  minlen <- 0
  log_file <- NULL
  threads <- 1

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
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
  minlen <- 0
  log_file <- NULL
  threads <- 1

  merged_sample1 <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                        reverse = reverse,
                                        fastqout = fastqout,
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
  minlen <- 0
  log_file <- withr::local_tempfile()
  threads <- 1

  return_value <- vs_fastq_mergepairs(fastq_input = fastq_input,
                                      reverse = reverse,
                                      fastqout = fastqout,
                                      minlen = minlen,
                                      log_file = log_file,
                                      threads = threads)
  expect_null(return_value)

  expect_true(file.exists(log_file))

})




