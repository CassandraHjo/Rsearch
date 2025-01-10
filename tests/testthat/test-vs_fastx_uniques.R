test_that("error when wrong input_format", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  input_format <- "fastx"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input,
                                input_format = input_format),
               "Invalid input_format. Choose from fasta or fastq.")
})

test_that("error when wrong output_format", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  output_format <- "fastx"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input,
                                output_format = output_format),
               "Invalid output_format. Choose from fasta or fastq.")
})

test_that("error when wrong strand", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds"))
  strand <- "something"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input,
                                strand = strand),
               "Invalid value for 'strand'. Choose from 'plus' or 'both'.")
})

test_that("error when fastx_input has incorrect columns if input is tibble and input_format = 'fastq'", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(-Header)

  input_format <- "fastq"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input, input_format = input_format),
               "FASTQ object must contain columns: Header, Sequence, Quality")
})

test_that("error when fastx_input has incorrect columns if input is tibble and input_format = 'fasta'", {

  fastx_input <- readRDS(test_path("testdata", "sample1", "R1_sample1_fastq_dataframe.rds")) |>
    dplyr::select(Header)

  input_format <- "fasta"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input, input_format = input_format),
               "FASTA object must contain columns: Header and Sequence")
})

test_that("error when input file does not exist", {

  fastx_input <- "some_file.fq"

  expect_error(vs_fastx_uniques(fastx_input = fastx_input),
               paste("Cannot find input file:", fastx_input))
})

test_that("dereplicate fastq file, and return fastq file", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- withr::local_tempfile()
  input_format <- "fastq"
  output_format <- "fastq"

  return_value <- vs_fastx_uniques(fastx_input = fastx_input,
                                   fastx_output = fastx_output,
                                   input_format = input_format,
                                   output_format = output_format)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "derep_R1_sample1.fq")))

})

test_that("dereplicate fastq file, and return fastq tibble", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- NULL
  input_format <- "fastq"
  output_format <- "fastq"

  derep_sample1_R1 <- vs_fastx_uniques(fastx_input = fastx_input,
                                       fastx_output = fastx_output,
                                       input_format = input_format,
                                       output_format = output_format)

  expect_equal(derep_sample1_R1,
               readRDS(test_path("testdata", "output", "derep_sample1_R1_fastq.rds")))

})

test_that("dereplicate fastq tibble, and return fastq file", {

  fastx_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  fastx_output <- withr::local_tempfile()
  input_format <- "fastq"
  output_format <- "fastq"

  return_value <- vs_fastx_uniques(fastx_input = fastx_input,
                                   fastx_output = fastx_output,
                                   input_format = input_format,
                                   output_format = output_format)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "derep_R1_sample1.fq")))
})

test_that("dereplicate fastq tibble, and return fastq tibble", {

  fastx_input <- microseq::readFastq(test_path("testdata", "sample1", "R1_sample1.fq"))
  fastx_output <- NULL
  input_format <- "fastq"
  output_format <- "fastq"

  derep_sample1_R1 <- vs_fastx_uniques(fastx_input = fastx_input,
                                       fastx_output = fastx_output,
                                       input_format = input_format,
                                       output_format = output_format)

  expect_equal(derep_sample1_R1,
               readRDS(test_path("testdata", "output", "derep_sample1_R1_fastq.rds")))
})

test_that("dereplicate fasta file, and return fasta file", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  fastx_output <- withr::local_tempfile()
  input_format <- "fasta"
  output_format <- "fasta"

  return_value <- vs_fastx_uniques(fastx_input = fastx_input,
                                   fastx_output = fastx_output,
                                   input_format = input_format,
                                   output_format = output_format)

  expect_null(return_value)

  expect_equal(microseq::readFasta(fastx_output),
               microseq::readFasta(test_path("testdata", "output", "derep_R1_sample1.fa")))
})

test_that("dereplicate fasta file, and return fasta tibble", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fa")
  fastx_output <- NULL
  input_format <- "fasta"
  output_format <- "fasta"

  derep_sample1_R1 <- vs_fastx_uniques(fastx_input = fastx_input,
                                       fastx_output = fastx_output,
                                       input_format = input_format,
                                       output_format = output_format)

  expect_equal(derep_sample1_R1,
               readRDS(test_path("testdata", "output", "derep_sample1_R1_fasta.rds")))
})

test_that("dereplicate fasta tibble, and return fasta file", {

  fastx_input <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  fastx_output <- withr::local_tempfile()
  input_format <- "fasta"
  output_format <- "fasta"

  return_value <- vs_fastx_uniques(fastx_input = fastx_input,
                                   fastx_output = fastx_output,
                                   input_format = input_format,
                                   output_format = output_format)

  expect_null(return_value)

  expect_equal(microseq::readFasta(fastx_output),
               microseq::readFasta(test_path("testdata", "output", "derep_R1_sample1.fa")))
})

test_that("dereplicate fasta tibble, and return fasta tibble", {

  fastx_input <- microseq::readFasta(test_path("testdata", "sample1", "R1_sample1.fa"))
  fastx_output <- NULL
  input_format <- "fasta"
  output_format <- "fasta"

  derep_sample1_R1 <- vs_fastx_uniques(fastx_input = fastx_input,
                                       fastx_output = fastx_output,
                                       input_format = input_format,
                                       output_format = output_format)

  expect_equal(derep_sample1_R1,
               readRDS(test_path("testdata", "output", "derep_sample1_R1_fasta.rds")))
})

test_that("dereplicate fastq file, and return fastq file with relabeling", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- withr::local_tempfile()
  input_format <- "fastq"
  output_format <- "fastq"
  relabel <- "OTU"

  return_value <- vs_fastx_uniques(fastx_input = fastx_input,
                                   fastx_output = fastx_output,
                                   input_format = input_format,
                                   output_format = output_format,
                                   relabel = relabel)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "derep_R1_sample1_relabel.fq")))
})

test_that("dereplicate fastq file, and return fastq file with relbeling sha1", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- withr::local_tempfile()
  input_format <- "fastq"
  output_format <- "fastq"
  relabel_sha1 <- TRUE

  return_value <- vs_fastx_uniques(fastx_input = fastx_input,
                                   fastx_output = fastx_output,
                                   input_format = input_format,
                                   output_format = output_format,
                                   relabel_sha1 = relabel_sha1)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "derep_R1_sample1_relabel_sha1.fq")))
})

test_that("dereplicate fastq file, and return fastq file with fastq_qout_max", {

  fastx_input <- test_path("testdata", "sample1", "R1_sample1.fq")
  fastx_output <- withr::local_tempfile()
  input_format <- "fastq"
  output_format <- "fastq"
  fastq_qout_max <- TRUE

  return_value <- vs_fastx_uniques(fastx_input = fastx_input,
                                   fastx_output = fastx_output,
                                   input_format = input_format,
                                   output_format = output_format,
                                   fastq_qout_max = fastq_qout_max)

  expect_null(return_value)

  expect_equal(microseq::readFastq(fastx_output),
               microseq::readFastq(test_path("testdata", "output", "derep_R1_sample1_fastq_qout_max.fq")))
})
