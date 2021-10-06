context("Testing BEMD")

set.seed(1)

test_that("bogus arguments throw error",{
  expect_error(bemd("bogus"))
  expect_error(bemd(7))
  expect_error(bemd(complex(2), num_imfs = -2))
  expect_error(bemd(complex(2), num_shiftings = -2))
  expect_error(bemd(complex(2), directions = -2))
})

test_that("output of BEMD is of correct size and form",{
  N <- 64
  set.seed(1)
  x <- rnorm(N) + rnorm(N) * 1i
  imfs <- bemd(x)
  expect_identical(dim(imfs), c(64L, 6L))
  expect_identical(class(imfs), c("mts", "ts", "matrix"))
  expect_identical(mode(imfs), "complex")
  x <- ts(x, start = 2000, frequency = 12)
  imfs <- emd(x, num_imfs = 3)
  expect_identical(tsp(imfs), tsp(x))
})


test_that("BEMD returns the original series for short series",{
  for (i in 1:3) {
    x <- rnorm(i) + rnorm(i) * 1i
    imfs <- bemd(x)
    expect_identical(ncol(imfs), 1L)
    expect_identical(x,c(imfs))
  }
})

test_that("EMD returns the same value each time",{
  N <- 64
  for (i in 1:10) {
    x <- rnorm(N) + rnorm(N) * 1i
    expect_identical(bemd(x, num_siftings = i), bemd(x, num_siftings = i))
  }
})

test_that("sum of imfs equals to original series",{
  N <- 64
  set.seed(1)
  x <- rnorm(N) + rnorm(N) * 1i
  imfs <- bemd(x)
  expect_equal(rowSums(imfs), x)
})

test_that("subsets of IMFs are identical for different num_imfs",{
  N <- 64
  set.seed(1)
  x <- rnorm(N) + rnorm(N) * 1i
  imfs3 <- bemd(x, num_imfs = 3)
  imfs4 <- bemd(x, num_imfs = 4)
  expect_identical(imfs3[, 1:2], imfs4[, 1:2])
})

test_that("num_imfs = 1 returns residual which equals data",{
  N <- 64
  set.seed(1)
  x <- rnorm(N) + rnorm(N) * 1i
  imfs <- bemd(x, num_imfs = 1)
  expect_identical(c(imfs), x)
})
