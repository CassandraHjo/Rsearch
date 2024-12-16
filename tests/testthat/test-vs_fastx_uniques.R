test_that("error when wrong file_format", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file_format <- "fastx"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input,
                             file_format = file_format),
               "Invalid file_format. Choose from fasta or fastq.")
})

test_that("error when wrong strand", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file_format <- "fastq"
  strand <- "something"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input,
                                file_format = file_format,
                                strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when fastx_input has incorrect columns if input is tibble and file_format = 'fastq'", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(-Header)

  file_format <- "fastq"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input, file_format = file_format),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when fastx_input has incorrect columns if input is tibble and file_format = 'fasta'", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) %>%
    dplyr::select(Header)

  file_format <- "fasta"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input, file_format = file_format),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when input file does not exist", {

  fastx_input <- "some_file.fq"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input),
               paste("Cannot find input file:", fastx_input))
})

