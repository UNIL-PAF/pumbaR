context("Mass-fit curve")

# load the protein group data
proteinGroups_path <- "Conde_9508_sub.txt"
pg <- load_MQ(proteinGroups_path)

test_that("Repeated entries are filtered", {
  pg_flt <- filter_repeated_entries(pg)

  expect_equal(dim(pg_flt), c(393, 57))
})

test_that("Weak intensity entries are filtered", {
  pg_flt <- filter_weak_intensity(pg)

  # the dimensions should be the same, but more 0
  expect_equal(dim(pg_flt), c(404, 57))
  expect_true(sum(pg == 0, na.rm=TRUE) < sum(pg_flt == 0, na.rm=TRUE))
})

test_that("Low density entries are filtered", {
  pg_flt <- filter_weak_intensity(pg, 0.8)
  pg_flt_2 <- filter_low_densities(pg_flt, 2)

  expect_equal(dim(pg_flt), c(404, 57))
  #expect_true(sum(is.na(pg_flt)) < sum(is.na(pg_flt_2)))
})

test_that("Curve is fitted", {
  mass_fit <- fit_curve(pg)
  test_tol <- 1e-06
  res_diff <- abs(mass_fit$coefficients[[1]] - 1.693976)
  expect_true(res_diff < test_tol)
})

test_that("Filters and fits", {
  mass_fit <- filter_and_fit(pg, low_density_threshold = 2)
  test_tol <- 1e-03
  res_diff <- abs(mass_fit$coefficients[[1]] - 2.905583)
  expect_true(res_diff < test_tol)
})

