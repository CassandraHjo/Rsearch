test_that("vs_uchime_denovo errors if only one output file is given", {

  fasta_input <- test_path("testdata", "sample1", "R2_sample1.fa")
  nonchimeras <- "some_file.fa"
  chimeras <- NULL

  expect_error(vs_uchime_denovo(fasta_input = fasta_input,
                                nonchimeras = nonchimeras,
                                chimeras = chimeras),
               "nonchimeras and chimeras must either both be specified or both unspecified.")
})

test_that("error when input file is missing", {

  expect_error(vs_uchime_denovo(fasta_input = "some_fasta_file.fa"),
               paste0("Cannot find input file: ", "some_fasta_file.fa"))
})

test_that("vs_uchime_denovo returns tibble with chimeras and stats", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")

  out <- vs_uchime_denovo(fasta_input = R1)

  expected_df <- readRDS(test_path(
    "testdata", "output", "sample1_uchime_denovo_default_fasta_file.rds"))

  expect_s3_class(out, "tbl_df")
  expect_true("chimeras" %in% names(attributes(out)))
  expect_true("statistics" %in% names(attributes(out)))
  expect_equal(out, expected_df)
})

test_that("vs_uchime_denovo works with tibble input", {

  R1 <- microseq::readFasta(test_path("testdata", "output", "derep_R1_sample1.fa"))

  out <- vs_uchime_denovo(fasta_input = R1)

  expected_df <- readRDS(test_path(
    "testdata", "output", "sample1_uchime_denovo_default_fasta_tibble.rds"))

  expect_s3_class(out, "tbl_df")
  expect_true("chimeras" %in% names(attributes(out)))
  expect_true("statistics" %in% names(attributes(out)))
  expect_equal(out, expected_df)
})

test_that("vs_uchime_denovo writes output files when paths are specified", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")

  tmp_nonchimeras <- withr::local_tempfile()
  tmp_chimeras <- withr::local_tempfile()

  expect_invisible(vs_uchime_denovo(fasta_input = R1,
                                    nonchimeras = tmp_nonchimeras,
                                    chimeras = tmp_chimeras,
                                    relabel_sha1 = TRUE))

  expect_true(file.exists(tmp_nonchimeras))
  expect_true(file.exists(tmp_chimeras))
})

test_that("vs_uchime_denovo handles optional arguments", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")
  log_file <- withr::local_tempfile()

  out <- vs_uchime_denovo(
    fasta_input = R1,
    relabel = "OTU",
    sample = "sample1",
    log_file = log_file,
    vsearch_options = c("")
  )

  expected_df <- readRDS(test_path(
    "testdata", "output", "sample1_uchime_denovo_optional_arguments.rds"))

  expect_s3_class(out, "tbl_df")
  expect_true("chimeras" %in% names(attributes(out)))
  expect_true("statistics" %in% names(attributes(out)))
  expect_equal(out, expected_df)
  expect_true(file.exists(log_file))
})
