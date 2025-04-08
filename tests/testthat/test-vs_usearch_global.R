test_that("error when wrong strand", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userout <- withr::local_tempfile()
  strand <- "wrong_input"

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 database = db,
                                 userout = userout,
                                 strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when both outputs are specified", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userout <- withr::local_tempfile()
  otutabout <- withr::local_tempfile()

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 database = db,
                                 userout = userout,
                                 otutabout = otutabout),
               "Only one of 'userout' or 'otutabout' can be specified.")
})

test_that("error when wrong columns in fastx_input fastq", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(Quality)
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 database = db,
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

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 database = db,
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

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 database = db,
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

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 database = db,
                                 userout = userout,
                                 strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when fastx_input does not exist", {

  fastx_input <- "some_file.fq"
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userout <- withr::local_tempfile()

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 database = db,
                                 userout = userout),
               paste0("Cannot find input file: ", fastx_input))
})

test_that("error when db does not exist", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- "some_file.fq"
  userout <- withr::local_tempfile()

  expect_error(vs_usearch_global(fastx_input = fastx_input,
                                 database = db,
                                 userout = userout),
               paste0("Cannot find input file: ", db))
})

test_that("allignment with default values with fastq tibbles as input", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userout <- withr::local_tempfile()
  vsearch_options <- c("")

  return_value <- vs_usearch_global(fastx_input = fastx_input,
                                    database = db,
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
               readRDS(test_path("testdata", "output", "sample1_usearch_global_default.rds")))
})

test_that("allignment with default values with fasta files as input", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- test_path("testdata", "output", "merged_sample1.fa")
  userout <- withr::local_tempfile()

  return_value <- vs_usearch_global(fastx_input = fastx_input,
                                    database = db,
                                    userout = userout)

  actual <- read.delim(userout,
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
  userout <- withr::local_tempfile()

  return_value <- vs_usearch_global(fastx_input = fastx_input,
                                    database = db,
                                    userout = userout)

  actual <- read.delim(userout,
                       sep = "\t",
                       header = FALSE)

  colnames(actual) <- c("query", "target", "id", "alnlen",
                        "mism", "opens", "qlo", "qhi",
                        "tlo", "thi", "evalue", "bits")

  expect_null(return_value)


  expect_equal(actual,
               readRDS(test_path("testdata", "output", "sample1_usearch_global_default_fasta_files.rds")))
})

test_that("vs_usearch_global returns OTU table tibble when otutabout = TRUE", {

  fasta_input <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa")) |>
    dplyr::mutate(Header = paste0(Header, ";sample=sample1"))

  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds")) |>
    dplyr::select(-Quality)

  otu_tbl <- vs_usearch_global(fastx_input = fasta_input,
                               database = db,
                               otutabout = TRUE)

  expect_equal(otu_tbl,
               readRDS(test_path("testdata",
                                 "output",
                                 "sample1_usearch_global_otu.rds")))
})

test_that("vs_usearch_global writes OTU table to file when otutabout is path", {

  fasta_input <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa")) |>
    dplyr::mutate(Header = paste0(Header, ";sample=sample1"))

  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds")) |>
    dplyr::select(-Quality)

  otutable_out <- withr::local_tempfile(fileext = ".tsv")

  return_val <- vs_usearch_global(fastx_input = fasta_input,
                                  database = db,
                                  otutabout = otutable_out)

  expect_null(return_val)
  expect_true(file.exists(otutable_out))
  expect_equal(suppressMessages(readr::read_delim(otutable_out)),
               readRDS(test_path("testdata",
                                 "output",
                                 "sample1_usearch_global_otu.rds")))
})

test_that("vs_usearch_global returns tibble with default userfields when no output file is specified", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))

  result <- vs_usearch_global(fastx_input = fastx_input,
                              database = db)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("query", "target", "id", "alnlen", "mism", "opens",
                         "qlo", "qhi", "tlo", "thi", "evalue", "bits"))
  expect_equal(result,
               readRDS(test_path("testdata", "output", "sample1_usearch_global_userfields.rds")))
})

test_that("vs_usearch_global respects custom userfields", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))
  userfields <- "query+target+id+alnlen"

  result <- vs_usearch_global(fastx_input = fastx_input,
                              database = db,
                              userfields = userfields)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("query", "target", "id", "alnlen"))
  expect_equal(result,
               readRDS(test_path("testdata", "output", "sample1_usearch_global_userfields_custom.rds")))
})

test_that("vs_usearch_global runs when strand is 'both' (without maxaccepts)", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_sample1_fastq_files.rds"))

  result <- vs_usearch_global(fastx_input = fastx_input,
                              database = db,
                              strand = "both")

  expect_s3_class(result, "tbl_df")
  expect_equal(result,
               readRDS(test_path("testdata", "output", "sample1_usearch_global_strand_both.rds")))
})

