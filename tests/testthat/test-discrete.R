

test_that("Zero probabilities for non-integers", {
  
  expect_warning(expect_equal(0, ddlaplace(0.5, 0.5)))
  expect_warning(expect_equal(0, ddnorm(0.5)))
  expect_warning(expect_equal(0, ddweibull(0.5, 0.5, 1)))
  expect_warning(expect_equal(0, ddunif(0.5, 0, 5)))
  expect_warning(expect_equal(0, dcat(0.5, c(0.5, 0.5))))
  expect_warning(expect_equal(0, dmnom(c(0.5, 1), 5, c(0.5, 0.5))))
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

