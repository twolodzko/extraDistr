

x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
       0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
n <- length(x)


test_that("Compare PDF's/PMF's to pure-R benchmarks", {
  
  expect_warning(expect_equal(dbbinom(x, 100, 1, 10), dbbinomR(x, 100, 1, 10)))
  expect_warning(expect_equal(dbnbinom(x[1:(n-2)], 100, 1, 10),
                              dbnbinomR(x[1:(n-2)], 100, 1, 10))) # numerical precission in R version
  expect_equal(dbetapr(x, 1, 1, 1), dbetaprR(x, 1, 1, 1))
  expect_equal(dfatigue(x, 1, 1, 0), dfatigueR(x, 1, 1, 0))
  expect_warning(expect_equal(ddlaplace(x, 0, 0.5), ddlaplaceR(x, 0, 0.5)))
  expect_warning(expect_equal(ddweibull(x, 0.5, 1), ddweibullR(x, 0.5, 1)))
  expect_equal(dfrechet(x, 1, 1, 1), dfrechetR(x, 1, 1, 1))
  expect_warning(expect_equal(dgpois(x[1:(n-1)], 1, 1),
                              dgpoisR(x[1:(n-1)], 1, 1))) # numerical precission in R version
  expect_equal(dgev(x, 1, 1, 1), dgevR(x, 1, 1, 1))
  expect_equal(dgompertz(x, 1, 1), dgompertzR(x, 1, 1))
  expect_equal(dgpd(x, 1, 1, 1), dgpdR(x, 1, 1, 1))
  expect_equal(dgumbel(x, 1, 1), dgumbelR(x, 1, 1))
  expect_equal(dinvgamma(x, 1, 1), dinvgammaR(x, 1, 1))
  expect_equal(dlaplace(x, -1, 5), dlaplaceR(x, -1, 5))
  expect_warning(expect_equal(dlgser(x, 0.5), dlgserR(x, 0.5)))
  expect_equal(dlomax(x, 1, 0.5), dlomaxR(x, 1, 0.5))
  expect_equal(dpareto(x, 1, 1), dparetoR(x, 1, 1))
  expect_equal(dpower(x, 1, 1), dpowerR(x, 1, 1))
  expect_equal(drayleigh(x, 1), drayleighR(x, 1))
  expect_equal(dsgomp(x, 0.5, 1), dsgompR(x, 0.5, 1))
  
})


test_that("Compare dhuber to hoa implementation", {
  
  skip_on_cran()
  skip_if_not_installed("hoa")

  expect_equal(dhuber(x), hoa::dHuber(x))
  expect_equal(phuber(x), hoa::pHuber(x))
  
})


test_that("Compare GEV and GPD to evd implementation", {
  
  skip_on_cran()
  skip_if_not_installed("evd")

  expect_equal(dgev(x), evd::dgev(x))
  expect_equal(pgev(x), evd::pgev(x))
  
  expect_equal(dgpd(x), evd::dgpd(x))
  expect_equal(pgpd(x), evd::pgpd(x))
  
})


test_that("Compare zero-inflated distributions to actuar implementation", {
  
  skip_on_cran()
  skip_if_not_installed("actuar")
  
  expect_equal(dzib(x, 45, 0.7, 0.2), actuar::dzmbinom(x, 45, 0.7, 0.2),
               tolerance = 1e-4)
  expect_equal(pzib(x, 45, 0.7, 0.2), actuar::pzmbinom(x, 45, 0.7, 0.2),
               tolerance = 1e-4)
  
  expect_equal(dzinb(x, 45, 0.7, 0.2), actuar::dzmnbinom(x, 45, 0.7, 0.2),
               tolerance = 1e-4)
  expect_equal(pzinb(x, 45, 0.7, 0.2), actuar::pzmnbinom(x, 45, 0.7, 0.2),
               tolerance = 1e-4)
  
  expect_equal(dzip(x, 7, 0.2), actuar::dzmpois(x, 7, 0.2),
               tolerance = 1e-3)
  expect_equal(pzip(x, 7, 0.2), actuar::pzmpois(x, 7, 0.2),
               tolerance = 1e-3)
  
})


test_that("Compare ddirichlet to Compositional implementation", {
  
  skip_on_cran()
  skip_if_not_installed("Compositional")
  
  alpha <- runif(5, 0, 3)
  x <- rdirichlet(5000, alpha)
  
  expect_equal(ddirichlet(x, alpha),
               Compositional::ddiri(x, alpha, logged = FALSE))
  
})


test_that("Compare dmnom to dmultinom from base R", {

  n <- 100
  p <- runif(5)
  p <- p/sum(p)
  
  x <- rmnom(5000, n, p)
  
  expect_equal(dmnom(x, n, p), apply(x, 1, dmultinom, n, p))
  
})


test_that("Bivariate Poisson distribution agrees with code by Karlis and Ntzoufras", {
  
  x <- rbvpois(1000, 7, 8, 5)

  expect_equal(
    dbvpois(x[,1], x[,2], 7, 8, 5),
    pbivpois(x[,1], x[,2], 7, 8, 5)
  )
  
})


test_that("Skellam distribution agrees with the implementation in skellam package", {
  
  skip_on_cran()
  skip_if_not_installed("skellam")
  
  x <- extraDistr::rskellam(5000, 7, 8)
  
  expect_equal(
    extraDistr::dskellam(x, 7, 8),
    skellam::dskellam(x, 7, 8)
  )
  
})

