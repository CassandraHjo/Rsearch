test_that("error when wrong strand", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userout <- withr::local_tempfile()
  strand <- "wrong_input"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               db = db,
                               userout = userout,
                               strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when wrong columns in fastx_input fastq", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(Quality)
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               db = db,
                               userout = userout,
                               strand = strand),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when wrong columns in db fastq", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))|>
    dplyr::select(Quality)
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               db = db,
                               userout = userout,
                               strand = strand),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})




test_that("error when wrong columns in fastx_input fasta", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(Header)
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds")) |>
    dplyr::select(-Quality)
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               db = db,
                               userout = userout,
                               strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when wrong columns in db fasta", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Quality)
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))|>
    dplyr::select(Header)
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               db = db,
                               userout = userout,
                               strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when fastx_input does not exist", {

  fastx_input <- "some_file.fq"
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userout <- withr::local_tempfile()

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               db = db,
                               userout = userout),
               paste0("Cannot find input file: ", fastx_input))
})

test_that("error when db does not exist", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- "some_file.fq"
  userout <- withr::local_tempfile()

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               db = db,
                               userout = userout),
               paste0("Cannot find input file: ", db))
})

test_that("search with default values with fastq tibbles as input", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))[1:500, ]
  userout <- withr::local_tempfile()
  vsearch_options <- c("")

  return_value <- vs_search_exact(fastx_input = fastx_input,
                                  db = db,
                                  userout = userout,
                                  vsearch_options = vsearch_options)

  actual <- read.delim(userout,
                       sep = "\t",
                       header = FALSE)

  colnames(actual) <- c("query", "target", "id", "alnlen",
                        "mism", "opens", "qlo", "qhi",
                        "tlo", "thi", "evalue", "bits")

  expect_null(return_value)


  expect_equal(actual,
               readRDS(test_path("testdata", "output", "sample1_search_exact_default.rds")))
})

test_that("search with default values with fasta files as input", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- test_path("testdata", "sample1", "R1_sample1.fa")
  userout <- withr::local_tempfile()

  return_value <- vs_search_exact(fastx_input = fastx_input,
                                  db = db,
                                  userout = userout)

  actual <- read.delim(userout,
                       sep = "\t",
                       header = FALSE)

  colnames(actual) <- c("query", "target", "id", "alnlen",
                        "mism", "opens", "qlo", "qhi",
                        "tlo", "thi", "evalue", "bits")

  expect_null(return_value)


  expect_equal(actual,
               readRDS(test_path("testdata", "output", "sample1_search_exact_default_fasta_files.rds")))
})

test_that("search with default values with fasta file and tibble as input", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))[1:500, ] |>
    dplyr::select(-Quality)
  userout <- withr::local_tempfile()

  return_value <- vs_search_exact(fastx_input = fastx_input,
                                  db = db,
                                  userout = userout)

  actual <- read.delim(userout,
                       sep = "\t",
                       header = FALSE)

  colnames(actual) <- c("query", "target", "id", "alnlen",
                        "mism", "opens", "qlo", "qhi",
                        "tlo", "thi", "evalue", "bits")

  expect_null(return_value)


  expect_equal(actual,
               readRDS(test_path("testdata", "output", "sample1_search_exact_default.rds")))
})
