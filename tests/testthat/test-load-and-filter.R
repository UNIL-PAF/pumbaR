context("Load and filter MQ data")

test_that("load_MQ filters correctly", {
  proteinGroups_path <- "Conde_9508_sub.txt"
  pg <- load_MQ(proteinGroups_path)

  expect_equal(nrow(pg), 474)
  expect_equal(ncol(pg), 56)
})






