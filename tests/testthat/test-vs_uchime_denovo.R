test_that("error when wrong specification of chimeras and nonchimeras", {

  fasta_input <- test_path("testdata", "sample1", "R2_sample1.fa")
  nonchimeras <- "some_file.fa"
  chimeras <- NULL

  expect_error(vs_uchime_denovo(fasta_input = fasta_input,
                                nonchimeras = nonchimeras,
                                chimeras = chimeras),
               "nonchimeras and chimeras must either both be specified or both unspecified.")
})

test_that("error when input file does not exist", {

  fasta_input <- "some_fasta_file.fa"
  nonchimeras <- NULL
  chimeras <- NULL

  expect_error(vs_uchime_denovo(fasta_input = fasta_input,
                                nonchimeras = nonchimeras,
                                chimeras = chimeras),
               paste0("Cannot find input file: ", fasta_input))
})
