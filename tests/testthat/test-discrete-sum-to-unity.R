
test_that("Discrete probabilities sum to unity", {

  expect_equal(sum(ddunif(-2:10, -2, 10)), 1)
  
  expect_equal(sum(dbbinom(0:10, 10, 1, 1)), 1)
  expect_equal(sum(dbbinom(0:10, 10, 1, 10)), 1)
  expect_equal(sum(dbbinom(0:10, 10, 10, 1)), 1)
  expect_equal(sum(dbbinom(0:10, 10, 10, 10)), 1)
  
  expect_equal(sum(dbnbinom(0:1e5, 10, 1, 1)), 1, tolerance = 1e-3)
  expect_equal(sum(dbnbinom(0:1e5, 10, 1, 10)), 1, tolerance = 1e-3)
  expect_equal(sum(dbnbinom(0:1e5, 10, 10, 1)), 1, tolerance = 1e-3)
  expect_equal(sum(dbnbinom(0:1e5, 10, 10, 10)), 1, tolerance = 1e-3)
  
  expect_equal(sum(dbvpois(expand.grid(0:100, 0:100), a = 10, b = 10, c = 10)), 1, tolerance = 1e-3)
  
  divBySum <- function(x) x/sum(x, na.rm = TRUE)
  
  expect_equal(sum(dcat(1:10, divBySum(runif(10)))), 1)
  
  expect_equal(sum(ddlaplace(-50:50, 0, 0.5)), 1)
  expect_equal(sum(ddlaplace(-50:50, 0, 0.1)), 1)
  expect_equal(sum(ddlaplace(-50:50, 0, 0.7)), 1)
  
  expect_equal(sum(ddnorm(-50:50, sd = 5)), 1)
  expect_equal(sum(ddnorm(-50:50, sd = 7)), 1)
  
  
  expect_equal(sum(dnhyper(0:100, 60, 35, 15)), 1)
  
  expect_equal(sum(dgpois(0:100, 10, 1)), 1)
  
  expect_equal(sum(dmnom(expand.grid(0:30, 0:30, 0:30), 30, c(0.2, 0.1, 0.7))), 1)
  expect_equal(sum(ddirmnom(expand.grid(0:30, 0:30, 0:30), 30, c(2, 3, 8))), 1)
  expect_equal(sum(dmvhyper(expand.grid(0:30, 0:30, 0:30), c(20, 30, 28), 30)), 1)
  
  expect_equal(sum(dskellam(-100:100, 10, 20)), 1)
  
  expect_equal(sum(dtpois(0:100, 30, a = -Inf, b = Inf)), 1)
  expect_equal(sum(dtpois(0:100, 30, a = 0, b = Inf)), 1)
  expect_equal(sum(dtpois(0:100, 30, a = 25, b = Inf)), 1)
  expect_equal(sum(dtpois(0:100, 30, a = -Inf, b = 32)), 1)
  expect_equal(sum(dtpois(0:100, 30, a = 25, b = 32)), 1)
  
  expect_equal(sum(dtbinom(0:100, 100, 0.29, a = -Inf, b = Inf)), 1)
  expect_equal(sum(dtbinom(0:100, 100, 0.29, a = 0, b = Inf)), 1)
  expect_equal(sum(dtbinom(0:100, 100, 0.29, a = 25, b = Inf)), 1)
  expect_equal(sum(dtbinom(0:100, 100, 0.29, a = -Inf, b = 32)), 1)
  expect_equal(sum(dtbinom(0:100, 100, 0.29, a = 25, b = 32)), 1)
  
  expect_equal(sum(dtbinom(0:100, 100, 0.5, a = -Inf, b = Inf)), 1)
  expect_equal(sum(dtbinom(0:100, 100, 0.5, a = 0, b = Inf)), 1)
  expect_equal(sum(dtbinom(0:100, 100, 0.5, a = 25, b = Inf)), 1)
  expect_equal(sum(dtbinom(0:100, 100, 0.5, a = -Inf, b = 32)), 1)
  expect_equal(sum(dtbinom(0:100, 100, 0.5, a = 25, b = 32)), 1)
  
  expect_equal(sum(dlgser(0:100, 0.6)), 1)
  expect_equal(sum(dlgser(0:100, 0.1)), 1)
  expect_equal(sum(dlgser(0:100, 0.8)), 1)
  
  
  expect_equal(sum(dmixpois(0:100, c(10, 20, 5), c(0.2, 0.5, 0.3))), 1)
  
  expect_equal(sum(dzib(0:1000, 30, 0.33, 0.4)), 1)
  expect_equal(sum(dzinb(0:1000, 30, 0.33, 0.4)), 1)
  expect_equal(sum(dzip(0:1000, 30, 0.4)), 1)
  
})
