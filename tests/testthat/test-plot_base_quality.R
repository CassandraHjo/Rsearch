test_that("error when fastq_input has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds")) |>
    dplyr::select(-Header)

  expect_error(plot_base_quality(fastq_input = R1),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when reverse has incorrect columns if input is tibble", {

  R1 <- readRDS(test_path("testdata", "R1_fastq_df.rds"))

  R2 <- readRDS(test_path("testdata", "R2_fastq_df.rds")) |>
    dplyr::select(-Header)

  expect_error(plot_base_quality(fastq_input = R1,
                                 reverse = R2),
               "Reverse FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error if invalid quantile range", {

  R1 <- test_path("testdata", "R1.fastq")
  R2 <- test_path("testdata", "R2.fastq")

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

  R1 <- test_path("testdata", "R1.fastq")

  quality_plot <- plot_base_quality(fastq_input = R1)

  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality handles reverse input correctly", {

  R1 <- test_path("testdata", "R1.fastq")
  R2 <- test_path("testdata", "R2.fastq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    reverse = R2)

  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality disables title correctly", {

  R1 <- test_path("testdata", "R1.fastq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    plot_title = "")
  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality only shows mean line", {

  R1 <- test_path("testdata", "R1.fastq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    show_median = FALSE)

  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality works with in-memory tibble input", {

  R1 <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  R2 <- microseq::readFastq(test_path("testdata", "R2.fastq"))

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    reverse = R2)

  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality returns ggplot with neither median nor mean", {

  R1 <- test_path("testdata", "R1.fastq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    show_median = FALSE,
                                    show_mean = FALSE)
  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality returns ggplot with overlap box", {

  R1 <- test_path("testdata", "R1.fastq")
  R2 <- test_path("testdata", "R2.fastq")

  quality_plot <- plot_base_quality(fastq_input = R1,
                                    reverse = R2,
                                    show_overlap_box = TRUE)
  expect_s3_class(quality_plot, "ggplot")
})

test_that("plot_base_quality downsamples if more than 10 000 reads", {

  R1 <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  R2 <- microseq::readFastq(test_path("testdata", "R2.fastq"))

  # Create a larger datasets by repeating the original data
  R1_large <- R1[rep(1:nrow(R1), length.out = 10001), ]
  R2_large <- R2[rep(1:nrow(R2), length.out = 10001), ]

  # Generate unique headers for both datasets
  unique_ids <- sprintf("read_%05d", seq_len(nrow(R1_large)))
  R1_large$Header <- unique_ids
  R2_large$Header <- unique_ids

  quality_plot <- plot_base_quality(fastq_input = R1_large,
                                    reverse = R2_large,
                                    show_overlap_box = TRUE)
  expect_s3_class(quality_plot, "ggplot")
})
