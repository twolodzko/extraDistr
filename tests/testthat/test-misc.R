

test_that("other tests", {
  
  expect_true(all(rsign(1e5) %in% c(-1, 1)))
  
  expect_silent(qcat(seq(0, 1, by = 0.1), c(1, 1, 1), labels = letters[1:3]))
  expect_warning(qcat(seq(0, 1, by = 0.1), c(1, 1, 1), labels = letters[1:2]))
  
  expect_silent(rcat(10, c(1, 1, 1), labels = letters[1:3]))
  expect_warning(rcat(10, c(1, 1, 1), labels = letters[1:2]))
  
  expect_silent(rcatlp(10, log(c(1, 1, 1)/3), labels = letters[1:3]))
  expect_warning(rcatlp(10, log(c(1, 1, 1)/3), labels = letters[1:2]))

  expect_error(dbvpois(1:10, a = 1, b = 1, c= 1))
  
  expect_silent(!anyNA(rtlambda(5000, -1)))
  expect_silent(!anyNA(rtlambda(5000, 0)))
  expect_silent(!anyNA(rtlambda(5000, 1)))
  
  expect_silent(dbvnorm(mtcars[, 1:2]))
  expect_error(dbvnorm(mtcars))
  
  xx <- seq(-6, 6, by = 0.01)
  expect_identical(dnorm(xx), dtnorm(xx))
  expect_identical(pnorm(xx), ptnorm(xx))
  
  pp <- seq(0, 1, by = 0.01)
  expect_identical(qnorm(pp), qtnorm(pp))
  
  x <- c(-Inf, -100, -10, -5, -1, -0.5, 0, 0.5, 1, 5, 10, 100, Inf)
  
  expect_equal(pmixnorm(x, c(1,2,3), c(1,2,3), c(1/3,1/3,1/3), lower.tail = TRUE),
               1 - pmixnorm(x, c(1,2,3), c(1,2,3), c(1/3,1/3,1/3), lower.tail = FALSE))
  expect_equal(suppressWarnings(pmixpois(x, c(1,2,3), c(1/3,1/3,1/3), lower.tail = TRUE)),
               1 - suppressWarnings(pmixpois(x, c(1,2,3), c(1/3,1/3,1/3), lower.tail = FALSE)))
  
})
