test_that("optimizing truncee_rate with default values and files as input", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")

  optimize.tbl <- vs_optimize_truncee_rate(fastq_input = fastq_input,
                                           reverse = reverse)

  expected_df <- readRDS(test_path("testdata", "output", "optimize_truncee_rate.rds"))

  expect_s3_class(attr(optimize.tbl, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(optimize.tbl, "plot") <- NULL
  attr(expected_df, "plot") <- NULL

  expect_equal(optimize.tbl, expected_df)

})

test_that("optimizing truncee_rate with tibbles as input", {

  fastq_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  reverse <- microseq::readFastq(test_path("testdata", "R2.fastq"))

  optimize.tbl <- vs_optimize_truncee_rate(fastq_input = fastq_input,
                                           reverse = reverse,
                                           plot_title = FALSE)

  expected_df <- readRDS(test_path("testdata", "output", "optimize_truncee_rate.rds"))

  expect_s3_class(attr(optimize.tbl, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(optimize.tbl, "plot") <- NULL
  attr(expected_df, "plot") <- NULL

  expect_equal(optimize.tbl, expected_df)

})

test_that("optimizing truncee_rate with pe_df tibble as input", {

  fastq_input <- readRDS(test_path("testdata", "pe_df.rds"))

  optimize.tbl <- vs_optimize_truncee_rate(fastq_input = fastq_input,
                                           min_size = 1,
                                           maxee_rate = 1.0)

  expected_df <- readRDS(test_path("testdata", "output", "optimize_truncee_rate_tibble.rds"))

  expect_s3_class(attr(optimize.tbl, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(optimize.tbl, "plot") <- NULL
  attr(expected_df, "plot") <- NULL

  expect_equal(optimize.tbl, expected_df)

  attr(fastq_input, "reverse") <- NULL
  expect_error(vs_optimize_truncee_rate(fastq_input = fastq_input),
               "fastq_input has class 'pe_df' but no 'reverse' attribute found.")
})

test_that("vs_optimize_truncee_rate handles vs_fastq_mergepairs failure gracefully", {
  fastq_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))

  tbl <- vs_optimize_truncee_rate(
    fastq_input = fastq_input,
    reverse = fastq_input,
    truncee_rate_range = c(0.002, 0.04),
    minovlen = 20,
    min_size = 1,
    maxee_rate = 0.01,
    threads = 1,
    plot_title = FALSE
  )

  expect_equal(tbl$merged_read_pairs[1], 0)
  expect_true(is.numeric(tbl$R1_length[1]))
  expect_true(is.numeric(tbl$R2_length[1]))
  expect_s3_class(attr(tbl, "plot"), "gg")
})
