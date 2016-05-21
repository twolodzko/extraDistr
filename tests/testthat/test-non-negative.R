

context("Testing distributions with non-negative support")

test_that("Zero probabilities for values <0", {
  
  expect_equal(0, ddweibull(-1, 0.5, 1))
  expect_equal(0, ddunif(-1, 0, 5))
  expect_equal(0, dcat(-1, c(0.5, 0.5)))
  expect_equal(0, dmnom(c(-1, 1), 5, c(0.5, 0.5)))
  expect_warning(expect_equal(0, dbern(-1)))
  expect_equal(0, dbbinom(-1, 1, 1, 1))
  expect_equal(0, dbnbinom(-1, 1, 1, 1))
  expect_equal(0, dgpois(-1, 1, 1))
  expect_equal(0, dmvhyper(c(1, 1, -1), c(2, 2, 2), 3))
  expect_equal(0, dzip(-1, 1, 0.5))
  expect_equal(0, dzib(-1, 1, 0.5, 0.5))
  expect_equal(0, dzinb(-1, 1, 0.5, 0.5))
  
  expect_equal(0, dinvchisq(-1, 1))
  expect_equal(0, dinvchisq(-1, 1, 1))
  expect_equal(0, dinvgamma(-1, 1, 1))
  expect_equal(0, dgompertz(-1, 1, 1))
  expect_equal(0, dgpois(-1, 1, 1))
  expect_equal(0, dlomax(-1, 1, 1))
  expect_equal(0, dpower(-1, 1, 0.5))
  expect_equal(0, drayleigh(-1, 1))
  expect_equal(0, dwald(-1, 1, 1))
  
  expect_equal(0, dhcauchy(-1, 1))
  expect_equal(0, dhnorm(-1, 1))
  expect_equal(0, dht(-1, 5, 1))
  
})

test_that("Zero probabilities for values x < mean", {
  
  expect_equal(0, dfatigue(-1, 1, 1, 0))
  expect_equal(0, dfrechet(-1, 1, 0, 1))
  expect_equal(0, dgpd(-1, 0, 1, 1))
  
})

test_that("Zero probabilities for values < 1", {
  
  expect_equal(c(0, 0), dlgser(c(-1, 0), 0.5))
  expect_equal(c(0, 0), dpareto(c(-1, 0), 1, 1))
  
})

