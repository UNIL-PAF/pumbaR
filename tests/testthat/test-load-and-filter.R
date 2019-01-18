context("Load and filter MQ data")

proteinGroups_path <- "Conde_9508_sub.txt"
pg <- load_MQ(proteinGroups_path)

test_that("load_MQ filters correctly", {
  expect_equal(nrow(pg), 404)
  expect_equal(ncol(pg), 57)
})

test_that("intensities are correctly extracted", {
  ints <- get_intensities(pg)

  expect_equal(nrow(ints), 404)
  expect_equal(ncol(ints), 45)
})

test_that("correct number of slices are extracted", {
  slice_numbers <- get_slice_numbers(pg)

  expect_equal(length(slice_numbers), 45)
})






