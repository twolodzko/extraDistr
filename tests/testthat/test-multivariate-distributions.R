

# extraDistr::ddirichlet was compared to MCMCpack::ddirichlet

test_that("Testing multivariate distributions", {
  
  expect_true(all.equal(rowSums(rmnom(5000, 50, c(2/10, 5/10, 3/10))), rep(50, 5000)))
  
  expect_true(all.equal(rowSums(rdirichlet(5000, c(2/10, 5/10, 3/10))), rep(1.0, 5000)))
  
  xx <- expand.grid(0:20, 0:20, 0:20)
  expect_equal(sum(ddirmnom(xx[rowSums(xx) == 20,], 20, c(2/10, 5/10, 3/10))), 1)
  
})
