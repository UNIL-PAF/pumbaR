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

test_that("slices are ignored", {
  pg_2 <- load_MQ(proteinGroups_path, ignore_slices=c(24,25))
  slice_numbers <- get_slice_numbers(pg_2)
  expect_equal(sum(pg_2[(get_intensity_columns(pg_2)[c(24,25)])]), 0)
  expect_equal(sum(pg_2[(get_intensity_columns(pg_2)[c(1)])]), 215960050)
  expect_equal(length(slice_numbers), 45)
})






