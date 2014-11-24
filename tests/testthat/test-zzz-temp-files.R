context("tmp file deletion")

test_that('all temporary files are deleted', {
  temp_files <- list.files(tempdir(), pattern = '^jaatha_[0-9]+_')
  for (temp_file in temp_files) {
    warning("temp file not deleted: ", temp_file)
  }
  expect_equal(length(temp_files), 0)
})