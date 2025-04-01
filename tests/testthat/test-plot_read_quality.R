test_that("error when fastq_input tibble is missing required columns", {

  R1 <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  expect_error(plot_read_quality(R1),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("plot_read_quality works with FASTQ file path (default args)", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  read_plot <- plot_read_quality(fastq_input = R1)

  expect_s3_class(read_plot, "ggExtraPlot")
})

test_that("plot_read_quality works with in-memory FASTQ tibble", {

  R1 <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))

  read_plot <- plot_read_quality(fastq_input = R1)

  expect_s3_class(read_plot, "ggExtraPlot")
})

test_that("plot_read_quality switches to expected error rate if specified", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  read_plot <- plot_read_quality(fastq_input = R1,
                                 use_ee_rate = TRUE)

  expect_s3_class(read_plot, "ggExtraPlot")
})

test_that("plot_read_quality works with plot_title = FALSE", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  read_plot <- plot_read_quality(fastq_input = R1,
                                 plot_title = FALSE)

  expect_s3_class(read_plot, "ggExtraPlot")
})

test_that("plot_read_quality respects alpha transparency setting", {

  R1 <- test_path("testdata", "sample1", "R1_sample1.fq")

  read_plot <- plot_read_quality(fastq_input = R1,
                                 alpha = 0.1)

  expect_s3_class(read_plot, "ggExtraPlot")
})
