test_that("optimizing truncqual with default values", {

  fastq_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  reverse <- test_path("testdata", "sample1", "R2_sample1.fq")

  optimize.tbl <- vs_optimize_truncqual(fastq_input = fastq_input,
                                        reverse = reverse)

  expect_equal(optimize.tbl,
               readRDS(test_path("testdata", "output", "optimize_truncqual.rds")))

})
