test_that("error when input fasta_input does not exist", {

  fasta_input <- test_path("testdata", "some_file.fa")

  expect_error(vs_cluster_size(fasta_input = fasta_input),
               paste("Cannot find input file:", fasta_input))
})

test_that("cluster sequences from fasta file, and return fasta file", {

  fasta_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  centroids <- withr::local_tempfile()

  return_value <- vs_cluster_size(fasta_input = fasta_input,
                                  centroids = centroids)

  expect_null(return_value)

  expect_equal(microseq::readFasta(centroids),
               microseq::readFasta(test_path("testdata", "output", "cluster_R1_sample1_file.fa")))
})

test_that("cluster sequences from fasta file, and return fasta tibble", {

  fasta_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  centroids <- NULL

  cluster_sample1_R1 <- vs_cluster_size(fasta_input = fasta_input,
                                        centroids = centroids)

  expect_equal(cluster_sample1_R1,
               readRDS(test_path("testdata", "output", "cluster_R1_sample1_file.rds")))
})

test_that("cluster sequences from fasta tibble, and return fasta file", {

  fasta_input <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  centroids <- withr::local_tempfile()

  return_value <- vs_cluster_size(fasta_input = fasta_input,
                                  centroids = centroids)

  expect_null(return_value)

  expect_equal(microseq::readFasta(centroids),
               microseq::readFasta(test_path("testdata", "output", "cluster_R1_sample1_file.fa")))
})

test_that("cluster sequences from fasta tibble, and return fasta tibble", {

  fasta_input <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  centroids <- NULL

  cluster_sample1_R1 <- vs_cluster_size(fasta_input = fasta_input,
                                        centroids = centroids)

  expect_equal(cluster_sample1_R1,
               readRDS(test_path("testdata", "output", "cluster_R1_sample1_tibble.rds")))
})
