context("Normalize intensitites")

proteinGroups_path <- "Conde_9508_sub.txt"
pg <- load_MQ(proteinGroups_path)

test_that("Normalize intensities", {
  ints <- get_intensities(pg)
  norm_ints <- normalize_intensities(ints)

  expect_equal(dim(ints), dim(norm_ints))
  expect_false(ints[1,37] == norm_ints[1,37])
})

test_that("Get normalized table", {
  norm_pg <- get_normalized_table(pg)

  expect_equal(nrow(norm_pg), 474)
  expect_equal(ncol(norm_pg), 48)
  expect_equal(colnames(norm_pg)[4], "Intensity.Norm.1")
})

test_that("Get maximal normalized intensity", {
  norm_pg <- get_normalized_table(pg)

  max_norm_int <- get_max_norm_int(norm_pg)
  expect_equal(max_norm_int, 0.1116033)
})
