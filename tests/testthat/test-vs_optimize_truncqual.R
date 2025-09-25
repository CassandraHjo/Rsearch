test_that("optimizing truncqual with default values and files as input", {

  fastq_input <- test_path("testdata", "R1.fastq")
  reverse <- test_path("testdata", "R2.fastq")

  optimize.tbl <- vs_optimize_truncqual(fastq_input = fastq_input,
                                        reverse = reverse)

  expected_df <- readRDS(test_path("testdata", "output", "optimize_truncqual.rds"))

  expect_s3_class(attr(optimize.tbl, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(optimize.tbl, "plot") <- NULL
  attr(expected_df, "plot") <- NULL
  attr(optimize.tbl, "optimal_truncqual") <- NULL
  attr(expected_df, "optimal_truncqual") <- NULL


  expect_equal(optimize.tbl, expected_df)
})

test_that("optimizing truncqual with tibbles as input", {

  fastq_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))
  reverse <- microseq::readFastq(test_path("testdata", "R2.fastq"))

  optimize.tbl <- vs_optimize_truncqual(fastq_input = fastq_input,
                                        reverse = reverse,
                                        plot_title = FALSE)

  expected_df <- readRDS(test_path("testdata", "output", "optimize_truncqual.rds"))

  expect_s3_class(attr(optimize.tbl, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(optimize.tbl, "plot") <- NULL
  attr(expected_df, "plot") <- NULL
  attr(optimize.tbl, "optimal_truncqual") <- NULL
  attr(expected_df, "optimal_truncqual") <- NULL

  expect_equal(optimize.tbl, expected_df)
})

test_that("optimizing truncqual with pe_df tibble as input", {

  fastq_input <- readRDS(test_path("testdata", "pe_df.rds"))

  optimize.tbl <- vs_optimize_truncqual(fastq_input = fastq_input,
                                        min_size = 1,
                                        maxee_rate = 1.0)

  expected_df <- readRDS(test_path("testdata", "output", "optimize_truncqual_tibble.rds"))

  expect_s3_class(attr(optimize.tbl, "plot"), "ggplot")

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(optimize.tbl, "plot") <- NULL
  attr(expected_df, "plot") <- NULL
  attr(optimize.tbl, "optimal_truncqual") <- NULL
  attr(expected_df, "optimal_truncqual") <- NULL

  expect_equal(optimize.tbl, expected_df)

  attr(fastq_input, "reverse") <- NULL
  expect_error(vs_optimize_truncqual(fastq_input = fastq_input),
               "fastq_input has class 'pe_df' but no 'reverse' attribute found.")
})

test_that("vs_optimize_truncqual handles vs_fastq_mergepairs failure gracefully", {
  fastq_input <- microseq::readFastq(test_path("testdata", "R1.fastq"))

  tbl <- vs_optimize_truncqual(
    fastq_input = fastq_input,
    reverse = fastq_input,
    truncqual_range = c(5, 10),
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

