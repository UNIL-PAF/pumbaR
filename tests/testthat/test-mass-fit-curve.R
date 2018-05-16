context("Mass-fit curve")

# load the protein group data
proteinGroups_path <- "Conde_9508_sub.txt"
pg <- load_MQ(proteinGroups_path)

test_that("Repeated entries are filtered", {
  pg_flt <- filter_repeated_entries(pg)

  expect_equal(dim(pg_flt), c(453, 56))
})

test_that("Weak intensity entries are filtered", {
  pg_flt <- filter_weak_intensity(pg)

  # the dimensions should be the same, but more 0
  expect_equal(dim(pg_flt), c(474, 56))
  expect_true(sum(pg == 0) < sum(pg_flt == 0))
})

test_that("Low density entries are filtered", {
  pg_flt <- filter_weak_intensity(pg, 0.8)
  pg_flt_2 <- filter_low_densities(pg_flt, 2)

  expect_equal(dim(pg_flt), c(474, 56))
  expect_true(sum(is.na(pg_flt)) < sum(is.na(pg_flt_2)))
})

test_that("Curve is fitted", {
  mass_fit <- fit_curve(pg)
  expect_equal(mass_fit$coefficients[[1]], 2.807391)
})

