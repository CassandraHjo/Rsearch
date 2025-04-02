test_that("error when wrong strand", {

  fasta_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  centroids <- withr::local_tempfile()
  strand <- "wrong_input"

  expect_error(vs_cluster_unoise(fasta_input = fasta_input,
                               centroids = centroids,
                               strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when input fasta_input does not exist", {

  fasta_input <- test_path("testdata", "some_file.fa")

  expect_error(vs_cluster_unoise(fasta_input = fasta_input),
               paste("Cannot find input file:", fasta_input))
})

test_that("error when both outputs are specified", {

  fasta_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  centroids <- withr::local_tempfile()
  otutabout <- withr::local_tempfile()

  expect_error(vs_cluster_unoise(fasta_input = fasta_input,
                               centroids = centroids,
                               otutabout = otutabout),
               "Only one of 'centroids' or 'otutabout' can be specified.")
})

test_that("vs_cluster_unoise works with relabel_sha1", {

  fasta_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  centroids <- withr::local_tempfile()

  expect_invisible(
    vs_cluster_unoise(fasta_input = fasta_input,
                      centroids = centroids,
                      relabel_sha1 = TRUE)
  )

  expect_true(file.exists(centroids))
})
