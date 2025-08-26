test_that("error when wrong strand", {

  fastx_input <- test_path("testdata", "R1.fasta")
  db <- test_path("testdata", "sintax_db.fasta")
  lcaout <- withr::local_tempfile()
  strand <- "wrong_input"

  expect_error(vs_alignment_classification(fastx_input = fastx_input,
                                 database = db,
                                 lcaout = lcaout,
                                 strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when lca_cutoff value is out of range", {

  fastx_input <- test_path("testdata", "R1.fasta")
  db <- test_path("testdata", "sintax_db.fasta")
  lcaout <- withr::local_tempfile()
  lca_cutoff <- 1.5

  expect_error(vs_alignment_classification(fastx_input = fastx_input,
                                           database = db,
                                           lcaout = lcaout,
                                           lca_cutoff = lca_cutoff),
               "Invalid value for 'lca_cutoff'. Must be between 0.5 and 1.0.")
})

test_that("error when wrong columns in fastx_input fastq", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(Quality)
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))
  lcaout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_alignment_classification(fastx_input = fastx_input,
                                 database = db,
                                 lcaout = lcaout,
                                 strand = strand),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when wrong columns in db fastq", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))|>
    dplyr::select(Quality)
  lcaout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_alignment_classification(fastx_input = fastx_input,
                                 database = db,
                                 lcaout = lcaout,
                                 strand = strand),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when wrong columns in fastx_input fasta", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(Header)
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds")) |>
    dplyr::select(-Quality)
  lcaout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_alignment_classification(fastx_input = fastx_input,
                                 database = db,
                                 lcaout = lcaout,
                                 strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when wrong columns in db fasta", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Quality)
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))|>
    dplyr::select(Header)
  lcaout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_alignment_classification(fastx_input = fastx_input,
                                 database = db,
                                 lcaout = lcaout,
                                 strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when fastx_input does not exist", {

  fastx_input <- "some_file.fq"
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))
  lcaout <- withr::local_tempfile()

  expect_error(vs_alignment_classification(fastx_input = fastx_input,
                                 database = db,
                                 lcaout = lcaout),
               paste0("Cannot find input file: ", fastx_input))
})

test_that("error when db does not exist", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  db <- "some_file.fq"
  lcaout <- withr::local_tempfile()

  expect_error(vs_alignment_classification(fastx_input = fastx_input,
                                 database = db,
                                 lcaout = lcaout),
               paste0("Cannot find input file: ", db))
})

test_that("allignment with default values with fasta files as input", {

  fastx_input <- test_path("testdata", "R1.fasta")
  db <- test_path("testdata", "sintax_db.fasta")
  lcaout <- withr::local_tempfile()
  vsearch_options <- c("")

  return_value <- vs_alignment_classification(fastx_input = fastx_input,
                                    database = db,
                                    lcaout = lcaout,
                                    vsearch_options = vsearch_options)

  actual <- read.delim(lcaout,
                       sep = "\t",
                       header = FALSE)

  expected <- read.delim(test_path("testdata", "output", "alignment_classification_default.tsv"),
                         sep = "\t",
                         header = FALSE)

  expect_null(return_value)
  expect_equal(actual, expected)
})

test_that("allignment with default values with fasta tables as input, and table as output", {

  fastx_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  db <- microseq::readFasta(test_path("testdata", "sintax_db.fasta"))
  lcaout <- NULL
  vsearch_options <- c("")

  return_value <- vs_alignment_classification(fastx_input = fastx_input,
                                              database = db,
                                              lcaout = lcaout,
                                              vsearch_options = vsearch_options)

  expected <- readRDS(test_path("testdata", "output", "alignment_classification_default.rds"))

  expect_equal(return_value, expected)
})

test_that("alignment returns tibble with top_hits_only = TRUE", {
  fastx_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  db <- microseq::readFasta(test_path("testdata", "sintax_db.fasta"))
  tmpd <- withr::local_tempdir()

  res <- vs_alignment_classification(
    fastx_input = fastx_input,
    database = db,
    lcaout = NULL,
    top_hits_only = TRUE,
    tmpdir = tmpd,
    vsearch_options = c("")
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("Header","Sequence","domain","phylum","class","order","family","genus","species") %in% names(res)))
  expect_equal(nrow(res), nrow(fastx_input))
})

test_that("alignment runs with strand = 'both'", {
  fastx_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  db <- microseq::readFasta(test_path("testdata", "sintax_db.fasta"))

  res <- vs_alignment_classification(
    fastx_input = fastx_input,
    database = db,
    lcaout = NULL,
    strand = "both",
    vsearch_options = c("")
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("Header","Sequence") %in% names(res)))
})

test_that("character FASTQ input path is handled and joined back to sequences", {
  fq_tbl <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  fq_file <- withr::local_tempfile(fileext = ".fq")
  microseq::writeFastq(fq_tbl, fq_file)

  db <- test_path("testdata", "sintax_db.fasta")

  res <- vs_alignment_classification(
    fastx_input = fq_file,   # karaktersti -> dekker readFastq + header-rydding
    database = db,
    lcaout = NULL,
    vsearch_options = c("")
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("Header","Sequence") %in% names(res)))
  expect_equal(nrow(res), nrow(fq_tbl))
})

test_that("non-NULL tmpdir is accepted and used", {
  fastx_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  db <- microseq::readFasta(test_path("testdata", "sintax_db.fasta"))
  tmpd <- withr::local_tempdir()

  res <- vs_alignment_classification(
    fastx_input = fastx_input,
    database = db,
    lcaout = NULL,
    tmpdir = tmpd,
    vsearch_options = c("")
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(fastx_input))
})

test_that("alignment accepts FASTQ tibble as database (covers writeFastq branch)", {
  fastx_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  db_fastq_tbl <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))
  tmpd <- withr::local_tempdir()

  res <- vs_alignment_classification(
    fastx_input = fastx_input,
    database = db_fastq_tbl,
    lcaout = NULL,
    tmpdir = tmpd,
    vsearch_options = c("")
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(fastx_input))
})

