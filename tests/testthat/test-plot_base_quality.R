test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  expect_error(plot_base_quality(fastq_input = R1),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))

  R2 <- readRDS(test_path("testdata", "sample1", "R2_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  expect_error(plot_base_quality(fastq_input = R1,
                                 reverse = R2),
               "Reverse FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error if invalid quantile range", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")
  R2 <- test_path("testdata", "sample1", "R2_sample1.fq")

  expect_error(plot_base_quality(fastq_input = R1,
                                 reverse = R2,
                                 quantile_lower = 0.25,
                                 quantile_upper = 1.1),
               "Invalid values for quantile range. Choose values between 0 and 1.")

  expect_error(plot_base_quality(fastq_input = R1,
                                 reverse = R2,
                                 quantile_lower = 0.25,
                                 quantile_upper = 0.10),
               "Invalid quantile range: 'quantile_lower' must be smaller than 'quantile_upper'.")
})

test_that("plot_base_quality handles forward FASTQ input as file", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  quality_plot <- plot_base_quality(fastq_input = R1)

  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality handles reverse input correctly", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")
  R2 <- test_path("testdata", "sample1", "R2_sample1.fq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    reverse = R2)

  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality disables title correctly", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    plot_title = "")
  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality only shows mean line", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    show_median = FALSE)

  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality works with in-memory tibble input", {

  R1 <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  R2 <- microseq::readFastq(test_path("testdata", "sample1", "R2_sample1.fq"))

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    reverse = R2)

  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality returns ggplot with neither median nor mean", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    show_median = FALSE,
                                    show_mean = FALSE)
  expect_s3_class(quality_plot, "ggplot")
})
