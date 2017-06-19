

all_in_bounds <- function(object, lower, upper) {
  all(object >= lower & object <= upper)
}


test_that("Check if values generated using rtnorm stay in bounds", {
  
  N <- 5000
  mu <- 0
  sigma <- 1
  
  a <- -1
  b <- 1
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

  a <- 0
  b <- 1
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

  a <- -1
  b <- 0
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

  a <- 4
  b <- 5
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

  a <- -5
  b <- -4
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

  a <- -Inf
  b <- 5
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

  a <- -5
  b <- Inf
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

  a <- -Inf
  b <- -4
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

  a <- 4
  b <- Inf
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_true(all_in_bounds(X, a, b))

})

