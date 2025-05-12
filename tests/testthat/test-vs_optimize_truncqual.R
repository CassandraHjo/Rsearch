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


  expect_equal(optimize.tbl, expected_df)

})
