test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_fastq_mergepairs(fastq_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  expect_error(vs_fastq_mergepairs(fastq_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when input file does not exist", {

  fastq_input <- "some_file.fq"
  reverse <- test_path("testdata", "R2_sample1.fq")

  expect_error(vs_fastq_mergepairs(fastq_input = fastq_input,
                                   reverse = reverse),
               paste("Cannot find input FASTQ file:", fastq_input))
})

test_that("error when reverse file does not exist", {

  fastq_input <- test_path("testdata", "R2_sample1.fq")
  reverse <- "some_file.fq"

  expect_error(vs_fastq_mergepairs(fastq_input = fastq_input,
                                   reverse = reverse),
               paste("Cannot find reverse FASTQ file:", reverse))
})
