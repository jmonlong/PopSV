context("Mean/sd bin counts in reference")

test_that("Output correct format and non-missing values", {
  x = rnorm(100,100,10)
  msd = mean.sd.outlierR(x)
  expect_true(is.list(msd))
  expect_true(!is.null(msd$m))
  expect_true(!is.na(msd$m))
  expect_true(!is.null(msd$sd))
  expect_true(!is.na(msd$sd))
})

test_that("Robust to NAs", {
  x = rnorm(100,100,10)
  x[1] = NA
  msd = mean.sd.outlierR(x)  
  expect_true(!is.null(msd$m))
  expect_true(!is.na(msd$m))
  expect_true(!is.null(msd$sd))
  expect_true(!is.na(msd$sd))
})

test_that("Works different distributions", {
  x = runif(100,50,150)
  msd = mean.sd.outlierR(x)  
  expect_true(!is.null(msd$m))
  expect_true(!is.na(msd$m))
  expect_true(!is.null(msd$sd))
  expect_true(!is.na(msd$sd))
  x = rpois(100,3)
  msd = mean.sd.outlierR(x)  
  expect_true(!is.null(msd$m))
  expect_true(!is.na(msd$m))
  expect_true(!is.null(msd$sd))
  expect_true(!is.na(msd$sd))
})

test_that("Works with mixture of gaussian", {
  x = c(rnorm(70,100,10), rnorm(30,150,10))
  msd = mean.sd.outlierR(x)
  expect_true(abs(msd$m-100)<10)
})

test_that("Stops if negative value inputed", {
  x = rnorm(100,100,10)
  x[1] = -1
  expect_error(mean.sd.outlierR(x), "negative value")
})
