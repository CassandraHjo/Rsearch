test_that("error when wrong output_format", {

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_fastq_filter(fastq_input = R1,
                               reverse = R2,
                               output_format = "fastx"),
               "Invalid output_format. Choose from fasta or fastq.")
})

test_that("error when output_format is 'fasta', and fastqout and fastqout_rev are defined", {

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds"))
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

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds"))
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

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds"))
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

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds"))
  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds"))
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

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds"))

  expect_error(vs_fastq_filter(fastq_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "R2_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  expect_error(vs_fastq_filter(fastq_input = R1, reverse = R2),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when input file does not exist", {

  fastq_input <- "some_file.fq"
  reverse <- test_path("testdata", "R2_sample1.fq")

  expect_error(vs_fastq_trim(fastq_input = fastq_input,
                             reverse = reverse),
               paste("Cannot find input FASTQ file:", fastq_input))
})

test_that("error when reverse file does not exist", {

  fastq_input <- test_path("testdata", "R2_sample1.fq")
  reverse <- "some_file.fq"

  expect_error(vs_fastq_trim(fastq_input = fastq_input,
                             reverse = reverse),
               paste("Cannot find reverse FASTQ file:", reverse))
})
