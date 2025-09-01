test_that("error when input fasta_input does not exist", {

  fasta_input <- test_path("testdata", "some_file.fa")

  expect_error(vs_cluster_unoise(fasta_input = fasta_input),
               paste("Cannot find input file:", fasta_input))
})

test_that("cluster sequences from fasta file, and return file", {

  fasta_input <- test_path("testdata", "R1.fasta")
  otutabout <- withr::local_tempfile()
  relabel <- "OTU"

  return_value <- vs_cluster_unoise(fasta_input = fasta_input,
                                    otutabout = otutabout,
                                    minsize = 1,
                                    relabel = relabel)

  expect_null(return_value)

  expect_equal(read.delim(otutabout),
               read.delim(test_path("testdata", "output", "cluster_unoise_otutabout.tsv")))
})

test_that("cluster sequences from fasta file, and return tibble", {

  fasta_input <- test_path("testdata", "R1.fasta")
  vsearch_options <- c("")

  cluster_sample1_R1 <- vs_cluster_unoise(fasta_input = fasta_input,
                                          vsearch_options = vsearch_options,
                                          minsize = 1,
                                          relabel_sha1 = TRUE)

  expect_equal(cluster_sample1_R1,
               readRDS(test_path("testdata", "output", "cluster_unoise_otutabout.rds")))
})

test_that("cluster sequences from fasta tibble, and return file", {

  fasta_input <- microseq::readFasta(test_path("testdata", "R1.fasta"))
  otutabout <- withr::local_tempfile()
  log_file <- withr::local_tempfile()
  relabel <- "OTU"

  return_value <- vs_cluster_unoise(fasta_input = fasta_input,
                                    otutabout = otutabout,
                                    log_file = log_file,
                                    minsize = 1,
                                    relabel = relabel)

  expect_null(return_value)

  expect_true(file.exists(log_file))

  expect_equal(read.delim(otutabout),
               read.delim(test_path("testdata", "output", "cluster_unoise_otutabout.tsv")))
})

test_that("vs_cluster_unoise returns NULL and warning when no clusters found", {

  fasta_input <- microseq::readFasta(test_path("testdata", "output", "R1_derep.fasta"))

  expect_warning(
    result <- vs_cluster_unoise(fasta_input = fasta_input, minsize = 20),
    "No clusters found, try to lower minsize"
  )

  expect_null(result)
})
