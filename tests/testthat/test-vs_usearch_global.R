test_that("error when wrong strand", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  blast6out <- withr::local_tempfile()
  strand <- "wrong_input"

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 db = db,
                                 blast6out = blast6out,
                                 strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when wrong columns in fastx_input fastq", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(Quality)
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  blast6out <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 db = db,
                                 blast6out = blast6out,
                                 strand = strand),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when wrong columns in db fastq", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))|>
    dplyr::select(Quality)
  blast6out <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 db = db,
                                 blast6out = blast6out,
                                 strand = strand),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})




test_that("error when wrong columns in fastx_input fasta", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(Header)
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds")) |>
    dplyr::select(-Quality)
  blast6out <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 db = db,
                                 blast6out = blast6out,
                                 strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when wrong columns in db fasta", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Quality)
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))|>
    dplyr::select(Header)
  blast6out <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 db = db,
                                 blast6out = blast6out,
                                 strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when fastx_input does not exist", {

  fastx_input <- "some_file.fq"
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  blast6out <- withr::local_tempfile()

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 db = db,
                                 blast6out = blast6out),
               paste0("Cannot find input file: ", fastx_input))
})

test_that("error when db does not exist", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- "some_file.fq"
  blast6out <- withr::local_tempfile()

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 db = db,
                                 blast6out = blast6out),
               paste0("Cannot find input file: ", db))
})

test_that("allignment with default values with fastq tibbles as input", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  blast6out <- withr::local_tempfile()
  vsearch_options <- c("")

  return_value <- vs_usearch_global(fastx_input = fastx_input,
                                    db = db,
                                    blast6out = blast6out,
                                    vsearch_options = vsearch_options)

  actual <- read.delim(blast6out,
                       sep = "\t",
                       header = FALSE)

  colnames(actual) <- c("query", "target", "id", "alnlen",
                        "mism", "opens", "qlo", "qhi",
                        "tlo", "thi", "evalue", "bits")

  expect_null(return_value)


  expect_equal(actual,
               readRDS(test_path("testdata", "output", "sample1_usearch_global_default.rds")))
})

test_that("allignment with default values with fasta files as input", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- test_path("testdata", "output", "merged_sample1.fa")
  blast6out <- withr::local_tempfile()

  return_value <- vs_usearch_global(fastx_input = fastx_input,
                                    db = db,
                                    blast6out = blast6out)

  actual <- read.delim(blast6out,
                       sep = "\t",
                       header = FALSE)

  colnames(actual) <- c("query", "target", "id", "alnlen",
                        "mism", "opens", "qlo", "qhi",
                        "tlo", "thi", "evalue", "bits")

  expect_null(return_value)


  expect_equal(actual,
               readRDS(test_path("testdata", "output", "sample1_usearch_global_default_fasta_files.rds")))
})

test_that("allignment with default values with fasta file and tibble as input", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds")) |>
    dplyr::select(-Quality)
  blast6out <- withr::local_tempfile()

  return_value <- vs_usearch_global(fastx_input = fastx_input,
                                    db = db,
                                    blast6out = blast6out)

  actual <- read.delim(blast6out,
                       sep = "\t",
                       header = FALSE)

  colnames(actual) <- c("query", "target", "id", "alnlen",
                        "mism", "opens", "qlo", "qhi",
                        "tlo", "thi", "evalue", "bits")

  expect_null(return_value)


  expect_equal(actual,
               readRDS(test_path("testdata", "output", "sample1_usearch_global_default_fasta_files.rds")))
})
