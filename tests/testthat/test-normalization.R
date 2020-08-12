context("Normalize intensitites")

proteinGroups_path <- "Conde_9508_sub.txt"
pg <- load_MQ(proteinGroups_path)

proteinGroups_path_merged <- "HEK293_11971_11973_12019_sub.txt"
sample_name <- "11971"
pg_merged <- load_MQ(proteinGroups_path_merged, sample_name = sample_name)

test_that("Normalize intensities", {
  ints <- get_intensities(pg)
  norm_ints <- normalize_intensities(ints)

  expect_equal(dim(ints), dim(norm_ints))
  expect_false(ints[1,35] == norm_ints[1,35])
})

test_that("Get normalized table", {
  norm_pg <- get_normalized_table(pg)

  expect_equal(nrow(norm_pg), 404)
  expect_equal(ncol(norm_pg), 50)
  expect_equal(colnames(norm_pg)[6], "Intensity.Norm.1")
})

test_that("Get maximal normalized intensity", {
  norm_pg <- get_normalized_table(pg)

  max_norm_int <- get_max_norm_int(norm_pg)
  expect_equal(max_norm_int, 0.05246229)
})

test_that("Get maximal normalized intensity from merged", {
  ints <- get_intensities(pg_merged, sample_name = sample_name)
  expect_equal(length(ints), 46)

  norm_pg <- get_normalized_table(pg_merged, sample_name = sample_name)

  max_norm_int <- get_max_norm_int(norm_pg)
  expect_equal(max_norm_int, 0.04129393)
})


