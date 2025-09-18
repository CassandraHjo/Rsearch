test_that("vs_uchime_denovo errors if only one output file is given", {

  fasta_input <- test_path("testdata", "R2.fasta")
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

  R1 <- test_path("testdata", "output", "R1_derep.fasta")

  out <- vs_uchime_denovo(fasta_input = R1)

  expected_df <- readRDS(test_path(
    "testdata", "output", "uchime_denovo_fa_file.rds"))

  expect_s3_class(out, "tbl_df")
  expect_true("chimeras" %in% names(attributes(out)))
  expect_true("statistics" %in% names(attributes(out)))
  expect_equal(out, expected_df)
})

test_that("vs_uchime_denovo works with tibble input", {

  R1 <- microseq::readFasta(test_path("testdata", "output", "R1_derep.fasta"))

  out <- vs_uchime_denovo(fasta_input = R1)

  expected_df <- readRDS(test_path(
    "testdata", "output", "uchime_denovo_fa_tibble.rds"))

  expect_s3_class(out, "tbl_df")
  expect_true("chimeras" %in% names(attributes(out)))
  expect_true("statistics" %in% names(attributes(out)))
  expect_equal(out, expected_df)
})

test_that("vs_uchime_denovo writes output files when paths are specified", {

  R1 <- test_path("testdata", "output", "R1_derep.fasta")

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

  R1 <- test_path("testdata", "output", "R1_derep.fasta")
  log_file <- withr::local_tempfile()

  out <- vs_uchime_denovo(
    fasta_input = R1,
    relabel = "OTU",
    sample = "sample1",
    log_file = log_file,
    vsearch_options = c("")
  )

  expected_df <- readRDS(test_path(
    "testdata", "output", "uchime_denovo_optional_args.rds"))

  expect_s3_class(out, "tbl_df")
  expect_true("chimeras" %in% names(attributes(out)))
  expect_true("statistics" %in% names(attributes(out)))
  expect_equal(out, expected_df)
  expect_true(file.exists(log_file))
})

test_that("vs_uchime_denovo joins otu_id when present in input", {
  fasta_input <- microseq::readFasta(test_path("testdata", "output", "R1_derep.fasta"))
  fasta_input$otu_id <- paste0("OTU", seq_len(nrow(fasta_input)))

  out <- vs_uchime_denovo(fasta_input = fasta_input)

  expect_true("otu_id" %in% colnames(out))
})

test_that("vs_uchime_denovo returns empty chimera table when none found", {
  fasta_input <- microseq::readFasta(test_path("testdata", "output", "R1_derep.fasta"))

  out <- vs_uchime_denovo(fasta_input = fasta_input)

  chimeras_tbl <- attr(out, "chimeras")
  expect_s3_class(chimeras_tbl, "data.frame")
  expect_equal(nrow(chimeras_tbl), 0)
})

test_that("vs_uchime_denovo preserves all extra columns from input tibble", {
  fasta_input <- microseq::readFasta(test_path("testdata", "chimera_test.fasta"))
  fasta_input$otu_id <- c("ParentA_otu", "ParentB_otu", "Chimera_otu")
  fasta_input$sample_source <- c("Sample1", "Sample1", "Sample1")
  fasta_input$qc_passed <- c(TRUE, TRUE, TRUE)

  out <- vs_uchime_denovo(fasta_input = fasta_input)

  # 3. Verifisering: Hent ut begge resultat-tabellene
  nonchimeras_tbl <- out
  chimeras_tbl <- attr(out, "chimeras")

  # Definer de forventede kolonnene
  expected_cols <- c("Header", "Sequence", "otu_id", "sample_source", "qc_passed")

  # Sjekk at alle ekstra kolonner finnes i nonchimeras-tabellen
  expect_true(all(expected_cols %in% colnames(nonchimeras_tbl)))

  # Sjekk at det finnes minst én kimære for å validere neste sjekk
  expect_gt(nrow(chimeras_tbl), 0)

  # Sjekk at alle ekstra kolonner finnes i chimeras-tabellen
  expect_true(all(expected_cols %in% colnames(chimeras_tbl)))
})
