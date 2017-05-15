

all_in_bounds <- function(object, lower, upper, warning = FALSE) {
  all(object >= lower & object <= upper)
}


test_that("Check if values generated using rtnorm stay in bounds", {
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- -1
  b <- 1
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- 0
  b <- 1
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- -1
  b <- 0
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- 4
  b <- 5
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- -5
  b <- -4
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- -Inf
  b <- 5
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- -5
  b <- Inf
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- -Inf
  b <- -4
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  N <- 5000
  mu <- 0
  sigma <- 1
  a <- 4
  b <- Inf
  
  expect_true(all_in_bounds(rtnorm(N, mu, sigma, a, b), a, b))
  
  
})

