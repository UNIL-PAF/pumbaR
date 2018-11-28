context("Merge proteins")

# load the protein group data
data <- list()
data[[1]] <- list()
data[[1]][['mass_fit_params']] <- c(2.92541101792906,-0.0856835305957989,
                                    0.00229499350478187,-2.85888515674011e-05)
data[[1]][['ints']] <- c(1930800,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,3145700,4036800,6036000,6696500000,
                         8679700000,63454000,0,0,0,0,0,0,0,0,0)
data[[2]] <- list()
data[[2]][['mass_fit_params']] <- c(2.93157752998003,-0.0911120501100777,
                                    0.00261614574546816,-3.38080473782906e-05)
data[[2]][['ints']] <- c(0,0,0,424060,0,0,26561000,0,0,0,8257400,0,0,0,0,0,0,
                         0,0,4896700,3668500,0,0,0,0,52009000,0,0,0,0,0,0,
                         5325700,2393800,2973300,4061200000,3175600000,
                         35474000,25662000,18546000,0,0,7599400,0,0,5936300,0)
data[[3]] <- list()
data[[3]][['mass_fit_params']] <- c(2.86255232915833,-0.0860407095427795,
                                    0.00249869975354482,-3.31808222647113e-05)
data[[3]][['ints']] <- c(1053800,0,0,0,0,0,0,0,3159900,1061300,0,4672100,0,0,
                         2268100,0,0,0,0,0,0,0,0,0,3289300,0,0,0,5302800,
                         4653300,0,7113600,24130000,6285400,460980000,
                         1.3405e+10,3747900000,17037000,15226000,0,0,
                         17088000,3872700,0,14897000,8817200)

merged_protein <- merge_proteins(data, cut_size=100, loess_span=0.05)

test_that("Dimensions are correct", {
  expect_equal(dim(merged_protein), c(4700, 2))
})

test_that("Intensity max is correct", {
  expect_equal(max(merged_protein$y), 8712848594)
})

data_2 <- list()
data_2[[1]] <- data[[1]]
merged_protein_2 <- merge_proteins(data_2, cut_size=100, loess_span=0.05)

test_that("Data with only one sample works", {
  expect_equal(dim(merged_protein_2), c(4700, 2))
})

test_that("Data with only one sample is correct", {
  expect_equal(max(merged_protein_2$y), 8996897187)
})

