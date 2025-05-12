test_that("error when wrong file_format", {

  files_dir <- test_path("testdata", "combine_data")
  file_format <- "fastx"

  expect_error(fastx_combine_files(files_dir = files_dir,
                                   output_file = NULL,
                                   file_ext = ".fa",
                                   file_format = file_format),
               "Invalid file_format. The files must be a fasta or fastq.")
})

test_that("error when directory does not exist", {

  files_dir <- "some_directory"

  expect_error(fastx_combine_files(files_dir = files_dir),
               paste("Directory does not exist:", files_dir))
})

test_that("error when no files are found in specified folder", {

  files_dir <- test_path("testdata", "empty_dir")
  file_ext <- ".fasta"
  file_format <- "fasta"

  expect_error(fastx_combine_files(files_dir = files_dir,
                                   file_ext = file_ext,
                                   file_format = file_format),
               paste("No", file_ext, "files found in the specified folder:", files_dir))
})

test_that("two fastq files can be combined, and written to fastq file", {

  files_dir <- test_path("testdata", "combine_data")
  output_file <- withr::local_tempfile()
  file_ext <- ".fastq"
  file_format <- "fastq"

  fastx_combine_files(files_dir = files_dir,
                      output_file = output_file,
                      file_ext = file_ext,
                      file_format = file_format)

  expect_equal(microseq::readFastq(output_file),
               microseq::readFastq(test_path("testdata", "output", "combine.fastq")))

})

test_that("two fastq files can be combined, and returned as fastq tibble", {

  files_dir <- test_path("testdata", "combine_data")
  output_file <- NULL
  file_ext <- ".fastq"
  file_format <- "fastq"

  expect_equal(fastx_combine_files(files_dir = files_dir,
                                   output_file = output_file,
                                   file_ext = file_ext,
                                   file_format = file_format),
               microseq::readFastq(test_path("testdata", "output", "combine.fastq")))

})

test_that("two fasta files can be combined, and written to fasta file", {

  files_dir <- test_path("testdata", "combine_data")
  output_file <- withr::local_tempfile()
  file_ext <- ".fasta"
  file_format <- "fasta"

  return_value <- fastx_combine_files(files_dir = files_dir,
                                      output_file = output_file,
                                      file_ext = file_ext,
                                      file_format = file_format)

  expect_null(return_value)

  expect_equal(microseq::readFasta(output_file),
               microseq::readFasta(test_path("testdata", "output", "combine.fasta")))

})

test_that("two fasta files can be combined, and returned as fasta tibble", {

  files_dir <- test_path("testdata", "combine_data")
  output_file <- NULL
  file_ext <- ".fasta"
  file_format <- "fasta"

  expect_equal(fastx_combine_files(files_dir = files_dir,
                                   output_file = output_file,
                                   file_ext = file_ext,
                                   file_format = file_format),
               microseq::readFasta(test_path("testdata", "output", "combine.fasta")))

})

test_that("existing output file is removed before combining", {

  files_dir <- test_path("testdata", "combine_data")
  output_file <- withr::local_tempfile()
  file_ext <- ".fastq"
  file_format <- "fastq"

  # Manually create the file so it exists before the function runs
  writeLines("This file will be overwritten", output_file)
  expect_true(file.exists(output_file))  # Confirm the file is there

  # Now run the function, which should remove the file
  fastx_combine_files(files_dir = files_dir,
                      output_file = output_file,
                      file_ext = file_ext,
                      file_format = file_format)

  # The function should overwrite the file with valid FASTQ data
  combined_data <- microseq::readFastq(output_file)
  expected_data <- microseq::readFastq(test_path("testdata", "output", "combine.fastq"))

  expect_equal(combined_data, expected_data)
})

