

test_that("K-S test values generated using rtnorm", {
  
  skip_on_cran()
  
  N <- 2500
  mu <- 0
  sigma <- 1
  
  a <- -1
  b <- 1
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
  a <- 0
  b <- 1
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
  a <- -1
  b <- 0
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
  a <- 4
  b <- 5
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
  a <- -5
  b <- -4
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
  a <- -Inf
  b <- 5
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
  a <- -5
  b <- Inf
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
  a <- -Inf
  b <- -4
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
  a <- 4
  b <- Inf
  X <- rtnorm(N, mu, sigma, a, b)
  
  expect_gt(ks.test(X, "ptnorm", mu, sigma, a, b)$p.value, 0.01)
  
})

