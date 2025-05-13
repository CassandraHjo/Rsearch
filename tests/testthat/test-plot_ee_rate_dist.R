test_that("plot_ee_rate_dist errors when tibble is missing required columns", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Header)

  expect_error(plot_ee_rate_dist(fastq_input = R1),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("plot_ee_rate_dist works with FASTQ file path", {

  R1 <- test_path("testdata", "R1.fastq")

  ee_plot <- plot_ee_rate_dist(fastq_input = R1)

  expect_s3_class(ee_plot, "ggplot")
})

test_that("plot_ee_rate_dist works with in-memory tibble", {

  R1 <- microseq::readFastq(test_path("testdata", "R1.fastq"))

  ee_plot <- plot_ee_rate_dist(fastq_input = R1)

  expect_s3_class(ee_plot, "ggplot")
})

test_that("plot_ee_rate_dist uses custom bin count", {

  R1 <- test_path("testdata", "R1.fastq")
  # This test only checks that it doesn't error
  ee_plot <- plot_ee_rate_dist(fastq_input = R1, n_bins = 10)

  expect_s3_class(ee_plot, "ggplot")
})

test_that("plot_ee_rate_dist works with no title", {

  R1 <- test_path("testdata", "R1.fastq")

  ee_plot <- plot_ee_rate_dist(fastq_input = R1, plot_title = "")

  expect_s3_class(ee_plot, "ggplot")
})

test_that("plot_ee_rate_dist downsamples if more than 10 000 reads", {

  R1 <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  # Create a larger datasets by repeating the original data
  R1_large <- R1[rep(1:nrow(R1), length.out = 10001), ]

  ee_plot <- plot_ee_rate_dist(fastq_input = R1_large)

  expect_s3_class(ee_plot, "ggplot")
})
