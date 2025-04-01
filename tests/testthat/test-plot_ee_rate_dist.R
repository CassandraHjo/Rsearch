test_that("plot_ee_rate_dist errors when tibble is missing required columns", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  expect_error(plot_ee_rate_dist(fastq_input = R1),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("plot_ee_rate_dist works with FASTQ file path", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  ee_plot <- plot_ee_rate_dist(fastq_input = R1)

  expect_s3_class(ee_plot, "ggplot")
})

test_that("plot_ee_rate_dist works with in-memory tibble", {

  R1 <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))

  ee_plot <- plot_ee_rate_dist(fastq_input = R1)

  expect_s3_class(ee_plot, "ggplot")
})

test_that("plot_ee_rate_dist uses custom bin count", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")
  # This test only checks that it doesn't error
  ee_plot <- plot_ee_rate_dist(fastq_input = R1, n_bins = 10)

  expect_s3_class(ee_plot, "ggplot")
})

test_that("plot_ee_rate_dist works with no title", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  ee_plot <- plot_ee_rate_dist(fastq_input = R1, plot_title = "")

  expect_s3_class(ee_plot, "ggplot")
})
