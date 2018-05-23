

# extraDistr::ddirichlet was compared to MCMCpack::ddirichlet

test_that("Testing multivariate distributions", {
  
  skip_on_cran()

  expect_true(all.equal(rowSums(rmnom(5000, 50, c(2/10, 5/10, 3/10))), rep(50, 5000)))
  expect_true(all.equal(rowSums(rdirichlet(5000, c(2/10, 5/10, 3/10))), rep(1.0, 5000)))

  xx <- expand.grid(0:20, 0:20, 0:20)
  expect_equal(sum(ddirmnom(xx[rowSums(xx) == 20,], 20, c(2, 5, 3))), 1)
  expect_equal(sum(dmnom(xx[rowSums(xx) == 20,], 20, c(2/10, 5/10, 3/10))), 1)
  expect_equal(sum(dmvhyper(xx[rowSums(xx) == 35,], c(20, 20, 20), 35)), 1)

  p <- c(4, 5, 1, 6, 2)

  expect_equal(prop.table(colSums(rmnom(1e5, 100, p/sum(p)))),
               p/sum(p),
               tolerance = 1e-2)

  expect_equal(as.numeric(prop.table(table(rcat(1e5, p/sum(p))))),
               p/sum(p),
               tolerance = 1e-2)

  expect_equal(prop.table(colSums(rdirichlet(1e5, p))),
               p/sum(p),
               tolerance = 1e-2)

  expect_equal(prop.table(colSums(rdirmnom(1e5, 100, p))),
               p/sum(p),
               tolerance = 1e-2)

  n <- c(11, 24, 43, 7, 56)

  expect_equal(prop.table(colSums(rmvhyper(1e5, n, 100))),
               n/sum(n),
               tolerance = 1e-2)

})


test_that("Evaluate wrong parameters first", {

  expect_warning(expect_true(is.nan(dbvpois(-1, -1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(ddirichlet(c(2, 2), c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(ddirmnom(c(-1, 1, 1), 1.5, c(1, 1, 1)))))
  expect_warning(expect_true(is.nan(dmnom(c(-1, 1, 1), 1.5, c(1/3, 1/3, 1/3)))))
  expect_warning(expect_true(is.nan(dmvhyper(c(-1, 2, 2), c(2,3,4), -5))))

})


test_that("Check if rmnom and rdirmnom deal with underflow (#3, #7)", {
  
  skip_on_cran()
  
  expect_false(anyNA(rmnom(5000, 100, c(0.504115095275327, 2.669522645838e-39, 0, 2.58539638831141, 0))))
  expect_false(anyNA(rdirmnom(5000, 100, c(1.480592e+00, 1.394943e-03, 4.529932e-06, 3.263573e+00, 4.554952e-06))))
  
  p <- c(0, 0, 1, 0, 2.77053929981958e-18)
  
  expect_false(anyNA(rmnom(5000, 100, p)))
  expect_false(anyNA(rdirmnom(5000, 100, p + 1e-5)))
  
})
