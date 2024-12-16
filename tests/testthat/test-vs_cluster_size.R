test_that("error when input fasta_input does not exist", {

  fasta_input <- test_path("testdata", "some_file.fa")

  expect_error(vs_cluster_size(fasta_input = fasta_input),
               paste("Cannot find input file:", fasta_input))
})
