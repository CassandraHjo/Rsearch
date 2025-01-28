test_that("optimizing truncqual with default values", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")

  optimize.tbl <- vs_optimize_truncqual(fastq_input = fastq_input,
                                        reverse = reverse)

  expected_df <- readRDS(test_path("testdata", "output", "optimize_truncqual.rds"))

  # Remove 'plot' attribute before comparison due to errors with ggplot
  attr(optimize.tbl, "sum_size_plot") <- NULL
  attr(optimize.tbl, "read_lengths_plot") <- NULL
  attr(expected_df, "sum_size_plot") <- NULL
  attr(expected_df, "read_lengths_plot") <- NULL


  expect_equal(optimize.tbl, expected_df)

})
