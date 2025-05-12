test_that("make_sintax_db produces a database", {

  input_df <- readRDS(test_path("testdata", "make_sintax_db.rds"))

  temp_out_db <- withr::local_tempfile()

  expect_invisible(make_sintax_db(taxonomy_table = input_df,
                                  outfile = temp_out_db))

  expect_true(file.exists(temp_out_db))
  expect_equal(
    microseq::readFasta(
      test_path("testdata", "output", "sintax_db_synthetic.fasta")),
    microseq::readFasta(temp_out_db))
})

test_that("make_sintax_db errors when required columns are missing", {
  input_df <- readRDS(test_path("testdata", "make_sintax_db.rds"))

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-Header), outfile),
    "must have a column named Header, with a unique text for each sequence"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-Sequence), outfile),
    "must have a column named Sequence, with the sequences"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-domain), outfile),
    "must have a column named domain"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-phylum), outfile),
    "must have a column named phylum"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-class), outfile),
    "must have a column named class"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-order), outfile),
    "must have a column named order"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-family), outfile),
    "must have a column named family"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-genus), outfile),
    "must have a column named genus"
  )

  expect_error(
    make_sintax_db(input_df |> dplyr::select(-species), outfile),
    "must have a column named species"
  )
})
