
test_that("Zero probabilities for non-integers", {
  
  expect_warning(expect_equal(0, ddlaplace(0.5, 0, 0.5)))
  expect_warning(expect_equal(0, ddnorm(0.5)))
  expect_warning(expect_equal(0, ddgamma(0.5, 9, 1)))
  expect_warning(expect_equal(0, ddweibull(0.5, 0.5, 1)))
  expect_warning(expect_equal(0, ddunif(0.5, 0, 5)))
  expect_warning(expect_equal(0, dcat(0.5, c(0.5, 0.5))))
  expect_warning(expect_equal(0, dmnom(c(0.5, 1), 5, c(0.5, 0.5))))
  expect_warning(expect_equal(0, dnhyper(0.5, 60, 35, 15)))
  expect_warning(expect_equal(0, dbern(0.5)))
  expect_warning(expect_equal(0, dbbinom(0.5, 1, 1, 1)))
  expect_warning(expect_equal(0, dbnbinom(0.5, 1, 1, 1)))
  expect_warning(expect_equal(0, dgpois(0.5, 1, 1)))
  expect_warning(expect_equal(0, dlgser(0.5, 0.5)))
  expect_warning(expect_equal(0, dmvhyper(c(1, 1, 0.5), c(2, 2, 2), 3)))
  expect_warning(expect_equal(0, dskellam(0.5, 1, 1)))
  expect_warning(expect_equal(0, dtpois(0.5, lambda = 25, a = 0)))
  expect_warning(expect_equal(0, dtbinom(0.5, 100, 0.56, a = 0)))
  expect_warning(expect_equal(0, dzip(0.5, 1, 0.5)))
  expect_warning(expect_equal(0, dzib(0.5, 1, 0.5, 0.5)))
  expect_warning(expect_equal(0, dzinb(0.5, 1, 0.5, 0.5)))
  expect_warning(expect_equal(0, dmixpois(0.5, c(1,2,3), c(1/3,1/3,1/3))))
  
})

test_that("cdf vs cumsum(pdf)", {
  
  xx <- seq(-500, 500, by = 1)
  epsilon <- 1e-4 # sqrt(.Machine$double.eps)
  
  expect_equal(cumsum(ddlaplace(xx, 0, 0.5)), pdlaplace(xx, 0, 0.5), tolerance = epsilon)
  expect_equal(cumsum(ddnorm(xx, 0, 15)), pdnorm(xx, 0, 15), tolerance = epsilon)
  
  xx <- seq(0, 200, by = 1)
  expect_equal(cumsum(ddweibull(xx, .32, 1)), pdweibull(xx, .32, 1), tolerance = epsilon)
  expect_equal(cumsum(ddunif(xx, 1, 199)), pdunif(xx, 1, 199), tolerance = epsilon)
  expect_equal(cumsum(ddgamma(xx, 9, 1)), pdgamma(xx, 9, 1), tolerance = epsilon)
  
  p <- rdirichlet(1, rep(1, 100))
  expect_equal(cumsum(dcat(xx, p)), pcat(xx, p), tolerance = epsilon)
  expect_equal(cumsum(dnhyper(xx, 60, 35, 15)), pnhyper(xx, 60, 35, 15), tolerance = epsilon)
  
  dnhyper2 <- function(x, n, m, r) ifelse(x<r | x>=(n+r), 0, choose(x-1, r-1)*choose(m+n-x, m-r)/choose(m+n, n))
  expect_equal(dnhyper2(xx, 60, 35, 15), dnhyper(xx, 60, 35, 15), tolerance = epsilon)
  expect_equal(cumsum(dnhyper2(xx, 60, 35, 15)), pnhyper(xx, 60, 35, 15), tolerance = epsilon)
  
  expect_equal(cumsum(dtbinom(xx, 200, 0.5, a = 100)), ptbinom(xx, 200, 0.5, a = 100), tolerance = epsilon)
  expect_equal(cumsum(dtbinom(xx, 200, 0.5, b = 100)), ptbinom(xx, 200, 0.5, b = 100), tolerance = epsilon)
  
  expect_equal(cumsum(dbbinom(xx, 200, 5, 13)), pbbinom(xx, 200, 5, 13), tolerance = epsilon)
  expect_equal(cumsum(dbnbinom(xx, 70, 5, 13)), pbnbinom(xx, 70, 5, 13), tolerance = epsilon)
  expect_equal(cumsum(dgpois(xx, 500, 16)), pgpois(xx, 500, 16), tolerance = epsilon)
  
  expect_equal(cumsum(dlgser(xx, 0.9)), plgser(xx, 0.9), tolerance = epsilon)
  expect_equal(cumsum(dtpois(xx, 100, a = 80, b = 150)), ptpois(xx, 100, a = 80, b = 150), tolerance = epsilon)
  expect_equal(cumsum(dtbinom(xx, 100, 0.5, a = 80, b = 150)), ptbinom(xx, 100, 0.5, a = 80, b = 150), tolerance = epsilon)
  
  expect_equal(cumsum(dzip(xx, 70, 0.5)), pzip(xx, 70, 0.5), tolerance = epsilon)
  expect_equal(cumsum(dzib(xx, 200, 0.5, 0.5)), pzib(xx, 200, 0.5, 0.5), tolerance = epsilon)
  expect_equal(cumsum(dzinb(xx, 70, 0.5, 0.5)), pzinb(xx, 70, 0.5, 0.5), tolerance = epsilon)
  expect_equal(cumsum(dmixpois(xx, c(40,50,70), c(1/3,1/3,1/3))), pmixpois(xx, c(40,50,70), c(1/3,1/3,1/3)), tolerance = epsilon)
  
})

