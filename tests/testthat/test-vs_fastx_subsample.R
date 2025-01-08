test_that("error when wrong file_format", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file_format <- "fastx"
  sample_size <- 10

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  file_format = file_format,
                                  sample_size = sample_size),
               "Invalid file_format. Choose from fasta or fastq.")
})

test_that("error if neither sample_size or sample_pct is specified", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file_format <- "fastq"
  sample_size <- NULL
  sample_pct <- NULL

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  file_format = file_format,
                                  sample_size = sample_size,
                                  sample_pct = sample_pct),
               "Either sample_size or sample_pct must be specified.")
})

test_that("error if both sample_size and sample_pct are specified", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  file_format <- "fastq"
  sample_size <- 10
  sample_pct <- 10.0

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  file_format = file_format,
                                  sample_size = sample_size,
                                  sample_pct = sample_pct),
               "Only specify one of the following parameters, not both: sample_size, sample_pct ")
})

test_that("error when fastx_input has incorrect columns if input is tibble and file_format = 'fastq'", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  file_format <- "fastq"
  sample_size <- 10

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  file_format = file_format,
                                  sample_size = sample_size),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when fastx_input has incorrect columns if input is tibble and file_format = 'fasta'", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(Header)

  file_format <- "fasta"
  sample_size <- 10

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  file_format = file_format,
                                  sample_size = sample_size),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when input file does not exist", {

  fastx_input <- "some_file.fq"
  sample_size <- 10

  expect_error(vs_fastx_subsample(fastx_input = fastx_input,
                                  sample_size = sample_size),
               paste("Cannot find input file:", fastx_input))
})

# -------------------------------
test_that("subsample fastq file with size, and return fastq file", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- withr::local_tempfile()
  file_format <- "fastq"
  sample_size <- 100
  randseed <- 1

  return_value <- vs_fastx_subsample(fastx_input = fastx_input,
                                     fastx_output = fastx_output,
                                     file_format = file_format,
                                     sample_size = sample_size,
                                     randseed = randseed
  )

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "subsample_R1_sample1.fq")))

})

test_that("subsample fastq file with size, and return fastq tibble", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- NULL
  file_format <- "fastq"
  sample_size <- 100
  randseed <- 1

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             file_format = file_format,
                                             sample_size = sample_size,
                                             randseed = randseed
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_sample1_R1_fastq.rds")))

})

test_that("subsample fasta file with pct, and return fasta file", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  fastx_output <- withr::local_tempfile()
  file_format <- "fasta"
  sample_pct <- 10.0
  randseed <- 1

  return_value <- vs_fastx_subsample(fastx_input = fastx_input,
                                     fastx_output = fastx_output,
                                     file_format = file_format,
                                     sample_pct = sample_pct,
                                     randseed = randseed
  )

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "subsample_R1_sample1.fa")))

})

test_that("subsample fasta file with pct, and return fasta tibble", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  fastx_output <- NULL
  file_format <- "fasta"
  sample_pct <- 10.0
  randseed <- 1

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             file_format = file_format,
                                             sample_pct = sample_pct,
                                             randseed = randseed
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_sample1_R1_fasta.rds")))

})

test_that("subsample fastq file with size, and return fastq file with relabeling", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- withr::local_tempfile()
  file_format <- "fastq"
  sample_size <- 100
  randseed <- 1
  relabel <- "OTU"

  return_value <- vs_fastx_subsample(fastx_input = fastx_input,
                                     fastx_output = fastx_output,
                                     file_format = file_format,
                                     sample_size = sample_size,
                                     randseed = randseed,
                                     relabel = relabel
  )

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "subsample_R1_sample1_relabel.fq")))

})

test_that("subsample fastq file with size, and return fastq file with sha1 relabeling", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- withr::local_tempfile()
  file_format <- "fastq"
  sample_size <- 100
  randseed <- 1
  relabel_sha1 <- TRUE

  return_value <- vs_fastx_subsample(fastx_input = fastx_input,
                                     fastx_output = fastx_output,
                                     file_format = file_format,
                                     sample_size = sample_size,
                                     randseed = randseed,
                                     relabel_sha1 = relabel_sha1
  )

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "subsample_R1_sample1_relabel_sha1.fq")))

})

test_that("subsample fastq tibble with size, and return fastq tibble", {

  fastx_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  fastx_output <- NULL
  file_format <- "fastq"
  sample_size <- 100
  randseed <- 1

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             file_format = file_format,
                                             sample_size = sample_size,
                                             randseed = randseed
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_sample1_R1_fastq_tibble.rds")))

})

test_that("subsample fasta tibble with size, and return fasta tibble", {

  fastx_input <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  fastx_output <- NULL
  file_format <- "fasta"
  sample_size <- 100
  randseed <- 1

  subsample_sample1_R1 <- vs_fastx_subsample(fastx_input = fastx_input,
                                             fastx_output = fastx_output,
                                             file_format = file_format,
                                             sample_size = sample_size,
                                             randseed = randseed
  )

  expect_equal(subsample_sample1_R1,
               readRDS(test_path("testdata", "output", "subsample_sample1_R1_fasta_tibble.rds")))

})
