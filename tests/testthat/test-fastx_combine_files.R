test_that("error when wrong file_format", {

  files_dir <- test_path("testdata")
  file_format <- "fastx"

  expect_error(fastx_combine_files(files_dir = files_dir,
                                   output_file = NULL,
                                   file_ext = ".fa",
                                   file_format = file_format),
               "Invalid file_format. Choose from fasta or fastq.")
})

test_that("error when directory does not exist", {

  files_dir <- "some_directory"

  expect_error(fastx_combine_files(files_dir = files_dir),
               paste("Directory does not exist:", files_dir))
})

test_that("error when no files are found in specified folder", {

  files_dir <- test_path("testdata", "empty_dir")
  file_ext <- ".fa"

  expect_error(fastx_combine_files(files_dir = files_dir),
               paste("No", file_ext, "files found in the specified folder:", files_dir))
})

test_that("two fastq files can be combined, and written to fastq file", {

  files_dir <- test_path("testdata", "sample1")
  output_file <- withr::local_tempfile()
  file_ext <- ".fq"
  file_format <- "fastq"

  fastx_combine_files(files_dir = files_dir,
                      output_file = output_file,
                      file_ext = file_ext,
                      file_format = file_format)

  expect_equal(microseq::readFastq(output_file),
               microseq::readFastq(test_path("testdata", "output", "combine_sample1.fq")))

})

test_that("two fastq files can be combined, and returned as fastq tibble", {

  files_dir <- test_path("testdata", "sample1")
  output_file <- NULL
  file_ext <- ".fq"
  file_format <- "fastq"

  expect_equal(fastx_combine_files(files_dir = files_dir,
                                   output_file = output_file,
                                   file_ext = file_ext,
                                   file_format = file_format),
               microseq::readFastq(test_path("testdata", "output", "combine_sample1.fq")))

})

test_that("two fasta files can be combined, and written to fasta file", {

  files_dir <- test_path("testdata", "sample1")
  output_file <- withr::local_tempfile()
  file_ext <- ".fa"
  file_format <- "fasta"

  fastx_combine_files(files_dir = files_dir,
                      output_file = output_file,
                      file_ext = file_ext,
                      file_format = file_format)

  expect_equal(microseq::readFasta(output_file),
               microseq::readFasta(test_path("testdata", "output", "combine_sample1.fa")))

})

test_that("two fasta files can be combined, and returned as fasta tibble", {

  files_dir <- test_path("testdata", "sample1")
  output_file <- NULL
  file_ext <- ".fa"
  file_format <- "fasta"

  expect_equal(fastx_combine_files(files_dir = files_dir,
                                   output_file = output_file,
                                   file_ext = file_ext,
                                   file_format = file_format),
               microseq::readFasta(test_path("testdata", "output", "combine_sample1.fa")))

})
