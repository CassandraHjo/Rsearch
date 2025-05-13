test_that("error when wrong strand", {

  fasta_input <- test_path("testdata", "R1.fasta")
  centroids <- withr::local_tempfile()
  strand <- "wrong_input"

  expect_error(vs_cluster_size(fasta_input = fasta_input,
                               centroids = centroids,
                               strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when input fasta_input does not exist", {

  fasta_input <- test_path("testdata", "some_file.fa")

  expect_error(vs_cluster_size(fasta_input = fasta_input),
               paste("Cannot find input file:", fasta_input))
})

test_that("error when both outputs are specified", {

  fasta_input <- test_path("testdata", "R1.fasta")
  centroids <- withr::local_tempfile()
  otutabout <- withr::local_tempfile()

  expect_error(vs_cluster_size(fasta_input = fasta_input,
                               centroids = centroids,
                               otutabout = otutabout),
               "Only one of 'centroids' or 'otutabout' can be specified.")
})

test_that("cluster sequences from fasta file, and return fasta file", {

  fasta_input <- test_path("testdata", "R1.fasta")
  centroids <- withr::local_tempfile()

  return_value <- vs_cluster_size(fasta_input = fasta_input,
                                  centroids = centroids,
                                  relabel_sha1 = TRUE)

  expect_null(return_value)

  expect_equal(microseq::readFasta(centroids),
               microseq::readFasta(test_path("testdata", "output", "cluster.fasta")))
})

test_that("cluster sequences from fasta file, and return fasta file", {

  fasta_input <- test_path("testdata", "R1.fasta")
  centroids <- withr::local_tempfile()
  relabel <- "OTU"

  return_value <- vs_cluster_size(fasta_input = fasta_input,
                                  centroids = centroids,
                                  relabel = relabel)

  expect_null(return_value)

  expect_equal(microseq::readFasta(centroids),
               microseq::readFasta(test_path("testdata", "output", "cluster_relabel.fasta")))
})

test_that("cluster sequences from fasta file, and return fasta tibble", {

  fasta_input <- test_path("testdata", "R1.fasta")
  vsearch_options <- c("")

  cluster_sample1_R1 <- vs_cluster_size(fasta_input = fasta_input,
                                        relabel = "OTU",
                                        vsearch_options = vsearch_options)

  expect_equal(cluster_sample1_R1,
               readRDS(test_path("testdata", "output", "cluster_fa_file_fa.rds")))
})

test_that("cluster sequences from fasta tibble, and return fasta file", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  centroids <- withr::local_tempfile()
  log_file <- withr::local_tempfile()

  return_value <- vs_cluster_size(fasta_input = fasta_input,
                                  centroids = centroids,
                                  relabel_sha1 = TRUE,
                                  log_file = log_file)

  expect_null(return_value)

  expect_true(file.exists(log_file))

  expect_equal(microseq::readFasta(centroids),
               microseq::readFasta(test_path("testdata", "output", "cluster.fasta")))
})

test_that("cluster sequences from fasta tibble, and return fasta tibble", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  centroids <- NULL

  cluster_sample1_R1 <- vs_cluster_size(fasta_input = fasta_input,
                                        centroids = centroids,
                                        size_column = TRUE)

  expect_equal(cluster_sample1_R1,
               readRDS(test_path("testdata", "output", "cluster_fa_tibble_fa.rds")))
})

test_that("vs_cluster_size returns OTU table tibble when otutabout = TRUE", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta")) |>
    dplyr::mutate(Header = paste0(Header, ";sample=sample1"))

  otu_tbl <- vs_cluster_size(fasta_input = fasta_input,
                             otutabout = TRUE,
                             relabel = "OTU",
                             sample = "sample1")

  expect_equal(otu_tbl,
               readRDS(test_path("testdata",
                                 "output",
                                 "cluster_fa_tibble_otu.rds")))
})

test_that("vs_cluster_size writes OTU table to file when otutabout is path", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta")) |>
    dplyr::mutate(Header = paste0(Header, ";sample=sample1"))

  otutable_out <- withr::local_tempfile(fileext = ".tsv")

  return_val <- vs_cluster_size(fasta_input = fasta_input,
                                otutabout = otutable_out,
                                relabel = "OTU")

  expect_null(return_val)
  expect_true(file.exists(otutable_out))
  expect_equal(suppressMessages(readr::read_delim(otutable_out)),
               suppressMessages(readr::read_delim(test_path("testdata",
                                                            "output",
                                                            "cluster_fa_tibble_otu.txt"))))
})

test_that("cluster sequences from empty fasta tibble, and return fasta tibble", {

  fasta_input <- data.frame(
    Header = "seq1",
    Sequence = "ATCGGCTA"
  )

  expect_warning(vs_cluster_size(fasta_input = fasta_input),
                 "No centroid sequences were returned by VSEARCH. Check input quality or parameters.")
})
