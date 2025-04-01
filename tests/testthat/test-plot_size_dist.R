test_that("plot_size_dist errors if tibble is missing required columns", {

  R1 <- microseq::readFastq(test_path("testdata", "output", "derep_R1_sample1.fq")) |>
    dplyr::select(-Header)

  expect_error(plot_size_dist(fastx_input = R1),
               "FASTX object must contain columns: Header and Sequence")
})

test_that("plot_size_dist errors if input file path is invalid", {

  expect_error(plot_size_dist(fastx_input = "nonexistent.fa",
                              input_format = "fasta"),
               "Cannot find input file")
})

test_that("plot_size_dist errors if input_format is missing or invalid", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fq")

  expect_error(plot_size_dist(fastx_input = R1),
               "Input format must be specified")

  expect_error(plot_size_dist(fastx_input = R1, input_format = "txt"),
               "Input format must be specified")
})

test_that("plot_size_dist works with FASTA file input and default args", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")

  size_plot <- plot_size_dist(fastx_input = R1,
                              input_format = "fasta")

  expect_s3_class(size_plot, "ggplot")
})

test_that("plot_size_dist works with FASTQ file input", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fq")

  size_plot <- plot_size_dist(fastx_input = R1,
                              input_format = "fastq")
  expect_s3_class(size_plot, "ggplot")
})

test_that("plot_size_dist works with in-memory tibble", {

  R1 <- microseq::readFasta(test_path("testdata", "output", "derep_R1_sample1.fa"))

  size_plot <- plot_size_dist(fastx_input = R1)

  expect_s3_class(size_plot, "ggplot")
})

test_that("plot_size_dist works with cutoff parameter", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")

  size_plot <- plot_size_dist(fastx_input = R1,
                              input_format = "fasta",
                              cutoff = 100)

  expect_s3_class(size_plot, "ggplot")
})

test_that("plot_size_dist works with no log scale", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")

  size_plot <- plot_size_dist(fastx_input = R1,
                              input_format = "fasta",
                              log_scale_y = FALSE)
  expect_s3_class(size_plot, "ggplot")
})

test_that("plot_size_dist works with custom y-axis breaks", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")

  size_plot <- plot_size_dist(fastx_input = R1,
                              input_format = "fasta",
                              y_breaks = c(1, 10, 100, 1000))

  expect_s3_class(size_plot, "ggplot")
})

test_that("plot_size_dist works with empty plot title", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")

  size_plot <- plot_size_dist(fastx_input = R1,
                              input_format = "fasta",
                              plot_title = "")

  expect_s3_class(size_plot, "ggplot")
})

test_that("plot_size_dist applies log scale with custom y_breaks", {

  R1 <- test_path("testdata", "output", "derep_R1_sample1.fa")

  size_plot <- plot_size_dist(fastx_input = R1,
                              input_format = "fasta",
                              cutoff = 100,
                              log_scale_y = TRUE,
                              y_breaks = c(1, 10, 100, 1000))

  expect_s3_class(size_plot, "ggplot")
})

