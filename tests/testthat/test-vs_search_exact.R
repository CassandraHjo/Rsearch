test_that("error when wrong strand", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))
  userout <- withr::local_tempfile()
  strand <- "wrong_input"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               database = db,
                               userout = userout,
                               strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when both outputs are specified", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))
  userout <- withr::local_tempfile()
  otutabout <- withr::local_tempfile()

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               database = db,
                               userout = userout,
                               otutabout = otutabout),
               "Only one of 'userout' or 'otutabout' can be specified.")
})

test_that("error when wrong columns in fastx_input fastq", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(Quality)
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               database = db,
                               userout = userout,
                               strand = strand),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when wrong columns in db fastq", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))|>
    dplyr::select(Quality)
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               database = db,
                               userout = userout,
                               strand = strand),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when wrong columns in fastx_input fasta", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(Header)
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds")) |>
    dplyr::select(-Quality)
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               database = db,
                               userout = userout,
                               strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when wrong columns in db fasta", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Quality)
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))|>
    dplyr::select(Header)
  userout <- withr::local_tempfile()
  strand <- "plus"

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               database = db,
                               userout = userout,
                               strand = strand),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when fastx_input does not exist", {

  fastx_input <- "some_file.fq"
  db <- readRDS(test_path("testdata", "output", "merged_fq_files_fq_tibble.rds"))
  userout <- withr::local_tempfile()

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               database = db,
                               userout = userout),
               paste0("Cannot find input file: ", fastx_input))
})

test_that("error when db does not exist", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  db <- "some_file.fq"
  userout <- withr::local_tempfile()

  expect_error(vs_search_exact(fastx_input = fastx_input,
                               database = db,
                               userout = userout),
               paste0("Cannot find input file: ", db))
})

test_that("search with default values with fastq tibbles as input", {

  fastx_input <- readRDS(test_path("testdata", "R1_fastq_df.rds"))
  db <- readRDS(test_path("testdata", "R1_fastq_df.rds"))[1:500, ]
  userout <- withr::local_tempfile()
  vsearch_options <- c("")

  return_value <- vs_search_exact(fastx_input = fastx_input,
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
               readRDS(test_path("testdata", "output", "search_exact_default.rds")))
})

test_that("search with default values with fasta files as input", {

  fastx_input <- test_path("testdata", "R1.fasta")
  db <- test_path("testdata", "R1.fasta")
  userout <- withr::local_tempfile()

  return_value <- vs_search_exact(fastx_input = fastx_input,
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
               readRDS(test_path("testdata", "output", "search_exact_fa_files.rds")))
})

test_that("search with default values with fasta file and tibble as input", {

  fastx_input <- test_path("testdata", "R1.fasta")
  db <- readRDS(test_path("testdata", "R1_fastq_df.rds"))[1:500, ] |>
    dplyr::select(-Quality)
  userout <- withr::local_tempfile()

  return_value <- vs_search_exact(fastx_input = fastx_input,
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
               readRDS(test_path("testdata", "output", "search_exact_default.rds")))
})

test_that("vs_search_exact returns OTU table tibble when otutabout = TRUE", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta")) |>
    dplyr::mutate(Header = paste0(Header, ";sample=sample1"))

  db <- fasta_input[1:10, ]  # ensure exact matches exist

  otu_tbl <- vs_search_exact(fastx_input = fasta_input,
                             database = db,
                             otutabout = TRUE)

  expect_s3_class(otu_tbl, "tbl_df")
  expect_true("otu_id" %in% names(otu_tbl))
  expect_equal(otu_tbl, readRDS(test_path("testdata", "output", "search_exact_otu.rds")))
})

test_that("vs_search_exact writes OTU table to file when otutabout is path", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta")) |>
    dplyr::mutate(Header = paste0(Header, ";sample=sample1"))

  db <- fasta_input[1:10, ]

  otu_outfile <- withr::local_tempfile(fileext = ".tsv")

  return_val <- vs_search_exact(fastx_input = fasta_input,
                                database = db,
                                otutabout = otu_outfile)

  actual <- suppressMessages(readr::read_delim(otu_outfile))
  expected <- suppressMessages(readr::read_delim(test_path("testdata",
                                                           "output",
                                                           "search_exact_otu.txt")))

  expect_null(return_val)
  expect_true(file.exists(otu_outfile))
  expect_equal(actual, expected)
})

test_that("vs_search_exact returns alignment tibble with default userfields", {

  fasta_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  db <- fasta_input[1:10, ]
  out_tbl <- vs_search_exact(fastx_input = fasta_input,
                             database = db)

  expect_s3_class(out_tbl, "tbl_df")
  expect_named(out_tbl, c("query", "target", "id", "alnlen", "mism", "opens",
                          "qlo", "qhi", "tlo", "thi", "evalue", "bits"))
  expect_equal(out_tbl,
               readRDS(test_path("testdata", "output", "search_exact_userfields.rds")))
})
