context("Load and filter MQ data")

proteinGroups_path <- "Conde_9508_sub.txt"
pg <- load_MQ(proteinGroups_path)

test_that("load_MQ filters correctly", {
  expect_equal(nrow(pg), 474)
  expect_equal(ncol(pg), 56)
})

test_that("load_MQ filters correctly second file", {
  proteinGroups_path_2 <- "Conde_9508_sub_2.txt"
  pg_2 <- load_MQ(proteinGroups_path_2)

  expect_equal(nrow(pg_2), 500)
  expect_equal(ncol(pg_2), 56)
})


test_that("intensities are correctly extracted", {
  ints <- get_intensities(pg)

  expect_equal(nrow(ints), 474)
  expect_equal(ncol(ints), 45)
})

test_that("correct number of slices are extracted", {
  slice_numbers <- get_slice_numbers(pg)

  expect_equal(length(slice_numbers), 45)
})






