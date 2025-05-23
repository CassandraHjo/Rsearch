test_that("throws error if input FASTA file doesn't exist", {
  expect_error(
    vs_cluster_subseq("nonexistent_file.fasta"),
    "Cannot find input file:"
  )
})

test_that("throws error for invalid strand input", {
  dummy_fasta <- tibble::tibble(Header = "seq1;size=1", Sequence = "ATGC")
  expect_error(
    vs_cluster_subseq(dummy_fasta, strand = "invalid_strand"),
    "Invalid value for 'strand'"
  )
})

test_that("returns a tibble with expected columns", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta")) |>
    dplyr::mutate(Header = paste0(Header, ";size=", nchar(Sequence)))

  result <- vs_cluster_subseq(fasta_input)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("Header", "Sequence", "members", "size") %in% names(result)))
})

test_that("writes output to file if centroids is specified", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta")) |>
    dplyr::mutate(Header = paste0(Header, ";size=", nchar(Sequence)))

  outfile <- tempfile(fileext = ".fa")
  log_file <- tempfile(fileext = ".log")

  result <- vs_cluster_subseq(fasta_input,
                              centroids = outfile,
                              log_file = log_file,
                              vsearch_options = "")
  expect_null(result)
  expect_true(file.exists(outfile))
  expect_true(file.exists(log_file))
})

test_that("respects sizein = FALSE", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))

  result <- vs_cluster_subseq(fasta_input, sizein = FALSE)
  expect_s3_class(result, "tbl_df")
  expect_true("members" %in% names(result))
  expect_false("size" %in% names(result))
})
