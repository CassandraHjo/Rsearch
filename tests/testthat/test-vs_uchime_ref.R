test_that("error when input file does not exist", {

  expect_error(vs_uchime_ref(fasta_input = "missing_file.fa", db = "db.fa"),
               "Cannot find input file: missing_file.fa")
})

test_that("error when only one of nonchimeras/chimeras is provided", {

  fasta <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- test_path("testdata", "sample1", "R1_sample1.fa")

  expect_error(
    vs_uchime_ref(fasta_input = fasta, db = db, nonchimeras = "nonchim.fa"),
    "nonchimeras and chimeras must either both be specified or both unspecified."
  )
})

test_that("vs_uchime_ref returns nonchimeras tibble with attributes", {

  fasta <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- test_path("testdata", "sample1", "R1_sample1.fa") # Needs to be updated with a real e

  out <- vs_uchime_ref(fasta_input = fasta, db = db)

  expect_s3_class(out, "tbl_df")
  expect_true("statistics" %in% names(attributes(out)))
  expect_true("chimeras" %in% names(attributes(out)))
})

test_that("vs_uchime_ref writes files with relabel and log", {

  fasta <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- test_path("testdata", "sample1", "R1_sample1.fa")
  nonchimeras <- withr::local_tempfile()
  chimeras <- withr::local_tempfile()
  log_file <- withr::local_tempfile()

  expect_invisible(vs_uchime_ref(fasta_input = fasta,
                                 db = db,
                                 nonchimeras = nonchimeras,
                                 chimeras = chimeras,
                                 relabel = "relabeled",
                                 log_file = log_file))

  expect_true(file.exists(nonchimeras))
  expect_true(file.exists(chimeras))
  expect_true(file.exists(log_file))
})

test_that("vs_uchime_ref works with FASTA tibble as input", {

  fasta_tbl <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  db_tbl <- fasta_tbl


  out <- vs_uchime_ref(fasta_input = fasta_tbl,
                       db = db_tbl,
                       relabel_sha1 = TRUE,
                       vsearch_options = c("--threads", "1"))

  expect_s3_class(out, "tbl_df")
  expect_true("statistics" %in% names(attributes(out)))
  expect_true("chimeras" %in% names(attributes(out)))
})

test_that("vs_uchime_ref returns empty data frame for chimeras if none found", {

  fasta <- test_path("testdata", "sample1", "R1_sample1.fa")
  db <- test_path("testdata", "sample1", "R1_sample1.fa")

  out <- vs_uchime_ref(fasta_input = fasta,
                       db = db,
                       sample = "sample1")

  expect_s3_class(attr(out, "chimeras"), "data.frame")
  expect_equal(nrow(attr(out, "chimeras")), 0)
})

