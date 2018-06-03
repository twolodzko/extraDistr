

test_that("Compare PDF's/PMF's to pure-R benchmarks", {
  
  x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
         0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
  n <- length(x)
  
  expect_warning(expect_equal(dbbinom(x, 100, 1, 10), dbbinomR(x, 100, 1, 10)))
  expect_warning(expect_equal(dbbinom(x, 100, 1000, 0.001), dbbinomR(x, 100, 1000, 0.001)))
  expect_warning(expect_equal(dbnbinom(x[1:(n-2)], 100, 1, 10),
                              dbnbinomR(x[1:(n-2)], 100, 1, 10))) # numerical precission in R version
  expect_equal(dbetapr(x, 1, 1, 1), dbetaprR(x, 1, 1, 1))
  expect_equal(dbetapr(x, 0.0001, 0.0001, 0.0001), dbetaprR(x, 0.0001, 0.0001, 0.0001))
  expect_equal(dbetapr(x, 100, 100, 100), dbetaprR(x, 100, 100, 100))
  expect_equal(dfatigue(x, 1, 1, 0), dfatigueR(x, 1, 1, 0))
  expect_equal(dfatigue(x, 0.0001, 0.0001, -100), dfatigueR(x, 0.0001, 0.0001, -100))
  expect_warning(expect_equal(ddlaplace(x, 0, 0.5), ddlaplaceR(x, 0, 0.5)))
  expect_warning(expect_equal(ddlaplace(x, 0, 0.999), ddlaplaceR(x, 0, 0.999)))
  expect_warning(expect_equal(ddlaplace(x, 0, 0.0001), ddlaplaceR(x, 0, 0.0001)))
  expect_warning(expect_equal(ddweibull(x, 0.5, 1), ddweibullR(x, 0.5, 1)))
  expect_warning(expect_equal(ddweibull(x, 0.0001, 0.0001), ddweibullR(x, 0.0001, 0.0001)))
  expect_warning(expect_equal(ddweibull(x, 0.9999, 1000), ddweibullR(x, 0.9999, 1000)))
  expect_equal(dfrechet(x, 1, 1, 1), dfrechetR(x, 1, 1, 1))
  expect_warning(expect_equal(dgpois(x[1:(n-1)], 1, 1),
                              dgpoisR(x[1:(n-1)], 1, 1))) # numerical precission in R version
  expect_equal(dgev(x, 1, 1, 1), dgevR(x, 1, 1, 1))
  expect_equal(dgev(x[2:n], 1, 1, 0), dgevR(x[2:n], 1, 1, 0)) # numerical precission in R version
  expect_equal(dgev(x, 1, 1, -1), dgevR(x, 1, 1, -1))
  expect_equal(dgompertz(x, 1, 1), dgompertzR(x, 1, 1))
  expect_equal(dgpd(x, 1, 1, 1), dgpdR(x, 1, 1, 1))
  expect_equal(dgpd(x, 1, 1, 0), dgpdR(x, 1, 1, 0))
  expect_equal(dgpd(x, 1, 1, -1), dgpdR(x, 1, 1, -1))
  expect_equal(dgumbel(x, 1, 1), dgumbelR(x, 1, 1))
  expect_equal(dinvgamma(x, 1, 1), dinvgammaR(x, 1, 1))
  expect_equal(dinvgamma(x, 1.2, 0.9), dinvgammaR(x, 1.2, 0.9))
  expect_equal(dlaplace(x, -1, 5), dlaplaceR(x, -1, 5))
  expect_equal(dlaplace(x, 9999, 0.000001), dlaplaceR(x, 9999, 0.000001))
  expect_warning(expect_equal(dlgser(x, 0.5), dlgserR(x, 0.5)))
  expect_warning(expect_equal(dlgser(x, 0.00001), dlgserR(x, 0.00001)))
  expect_warning(expect_equal(dlgser(x, 0.9999), dlgserR(x, 0.9999)))
  expect_equal(dlomax(x, 1, 0.5), dlomaxR(x, 1, 0.5))
  expect_equal(dpareto(x, 1, 1), dparetoR(x, 1, 1))
  expect_equal(dpower(x, 1, 1), dpowerR(x, 1, 1))
  expect_equal(drayleigh(x, 1), drayleighR(x, 1))
  expect_equal(dsgomp(x, 0.5, 1), dsgompR(x, 0.5, 1))
  
})


test_that("Compare dinvgamma to actuar implementation", {
  
  skip_on_cran()
  skip_if_not_installed("actuar")
  
  x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
         0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
  
  expect_equal(dinvgamma(x, 1.2, 0.9), actuar::dinvgamma(x, 1.2, scale=0.9))
  expect_equal(pinvgamma(x, 1.2, 0.9), actuar::pinvgamma(x, 1.2, scale=0.9))
  
})


test_that("Compare dhuber to hoa implementation", {
  
  x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
         0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
  n <- length(x)
  
  skip_on_cran()
  skip_if_not_installed("hoa")

  expect_equal(dhuber(x), hoa::dHuber(x))
  expect_equal(phuber(x), hoa::pHuber(x))
  
})


test_that("Compare triangular to triangle implementation", {
  
  x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
         0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
  n <- length(x)
  
  skip_on_cran()
  skip_if_not_installed("triangle")
  
  expect_equal(dtriang(x, -1, 1), triangle::dtriangle(x, -1, 1))
  expect_equal(dtriang(x, 0, 1, 0.8), triangle::dtriangle(x, 0, 1, 0.8))
  
  expect_equal(ptriang(x, -1, 1), triangle::ptriangle(x, -1, 1))
  expect_equal(ptriang(x, 0, 1, 0.8), triangle::ptriangle(x, 0, 1, 0.8))
  
})


test_that("Compare GEV and GPD to evd implementation", {
  
  x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
         0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
  n <- length(x)
  
  skip_on_cran()
  skip_if_not_installed("evd")

  expect_equal(dgev(x, 0, 1, 0), evd::dgev(x, 0, 1, 0))
  expect_equal(dgev(x, 0, 1, 1), evd::dgev(x, 0, 1, 1))
  expect_equal(dgev(x, 0, 1, -1), evd::dgev(x, 0, 1, -1))
  
  expect_equal(pgev(x, 0, 1, 0), evd::pgev(x, 0, 1, 0))
  expect_equal(pgev(x, 0, 1, 1), evd::pgev(x, 0, 1, 1))
  expect_equal(pgev(x, 0, 1, -1), evd::pgev(x, 0, 1, -1))
  
  expect_equal(dgpd(x, 0, 1, 0), evd::dgpd(x, 0, 1, 0))
  expect_equal(dgpd(x, 0, 1, 1), evd::dgpd(x, 0, 1, 1))
  expect_equal(dgpd(x, 0, 1, -1), evd::dgpd(x, 0, 1, -1))
  
  expect_equal(pgpd(x, 0, 1, 0), evd::pgpd(x, 0, 1, 0))
  expect_equal(pgpd(x, 0, 1, 1), evd::pgpd(x, 0, 1, 1))
  expect_equal(pgpd(x, 0, 1, -1), evd::pgpd(x, 0, 1, -1))
  
})


test_that("Compare to distributions from VGAM package", {
  
  x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
         0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
  n <- length(x)
  
  skip_on_cran()
  skip_if_not_installed("VGAM")
  
  expect_warning(expect_equal(dzib(x, 45, 0.7, 0.2), VGAM::dzibinom(x, 45, 0.7, 0.2)))
  expect_warning(expect_equal(dzib(x, 45, 0.7, 0.0001), VGAM::dzibinom(x, 45, 0.7, 0.0001)))
  expect_warning(expect_equal(dzib(x, 45, 0.7, 0.9999), VGAM::dzibinom(x, 45, 0.7, 0.9999)))
  expect_equal(pzib(x, 45, 0.7, 0.2), VGAM::pzibinom(x, 45, 0.7, 0.2))
  expect_equal(pzib(x, 45, 0.7, 0.0001), VGAM::pzibinom(x, 45, 0.7, 0.0001))
  expect_equal(pzib(x, 45, 0.7, 0.9999), VGAM::pzibinom(x, 45, 0.7, 0.9999))
  
  expect_warning(expect_equal(dzinb(x, 45, 0.7, 0.2), VGAM::dzinegbin(x, 45, 0.7, NULL, 0.2)))
  expect_warning(expect_equal(dzinb(x, 45, 0.7, 0.0001), VGAM::dzinegbin(x, 45, 0.7, NULL, 0.0001)))
  expect_warning(expect_equal(dzinb(x, 45, 0.7, 0.9999), VGAM::dzinegbin(x, 45, 0.7, NULL, 0.9999)))
  expect_equal(pzinb(x, 45, 0.7, 0.2), VGAM::pzinegbin(x, 45, 0.7, NULL, 0.2))
  expect_equal(pzinb(x, 45, 0.7, 0.0001), VGAM::pzinegbin(x, 45, 0.7, NULL, 0.0001))
  expect_equal(pzinb(x, 45, 0.7, 0.9999), VGAM::pzinegbin(x, 45, 0.7, NULL, 0.9999))
  
  expect_warning(expect_equal(dzip(x, 7, 0.2), VGAM::dzipois(x, 7, 0.2)))
  expect_warning(expect_equal(dzip(x, 7, 0.0001), VGAM::dzipois(x, 7, 0.0001)))
  expect_warning(expect_equal(dzip(x, 7, 0.9999), VGAM::dzipois(x, 7, 0.9999)))
  expect_equal(pzip(x, 7, 0.2), VGAM::pzipois(x, 7, 0.2))
  expect_equal(pzip(x, 7, 0.0001), VGAM::pzipois(x, 7, 0.0001))
  expect_equal(pzip(x, 7, 0.9999), VGAM::pzipois(x, 7, 0.9999))
  
  expect_equal(dlaplace(x), VGAM::dlaplace(x))
  expect_equal(plaplace(x), VGAM::plaplace(x))
  
  expect_equal(dpareto(x, 1, 1), VGAM::dpareto(x, 1, 1))
  expect_equal(ppareto(x, 1, 1), suppressWarnings(VGAM::ppareto(x, 1, 1)))
  
  expect_warning(expect_equal(dbbinom(x, 100, 2, 7), VGAM::dbetabinom.ab(x, 100, 2, 7)))
  expect_equal(pbbinom(x, 100, 2, 7), VGAM::pbetabinom.ab(x, 100, 2, 7))
  
  expect_equal(dfrechet(x, 1, 0, 1), VGAM::dfrechet(x, 0, 1, 1))
  expect_equal(pfrechet(x, 1, 0, 1), VGAM::pfrechet(x, 0, 1, 1))
  
  expect_equal(dgumbel(x, 0, 1), VGAM::dgumbel(x, 0, 1))
  expect_equal(pgumbel(x, 0, 1), VGAM::pgumbel(x, 0, 1))
  
  expect_equal(dgompertz(x, 1, 1), VGAM::dgompertz(x, 1, 1))
  expect_equal(pgompertz(x, 1, 1), VGAM::pgompertz(x, 1, 1))  
  
  expect_equal(dkumar(x, 3, 7), VGAM::dkumar(x, 3, 7))
  expect_equal(pkumar(x, 3, 7), suppressWarnings(VGAM::pkumar(x, 3, 7)))
  
  expect_equal(dslash(x), VGAM::dslash(x))
  # expect_equal(pslash(x), VGAM::pslash(x))
  
  expect_equal(dhuber(x), VGAM::dhuber(x, 1.345))
  expect_equal(phuber(x), VGAM::phuber(x, 1.345))
  
  expect_warning(expect_equal(dlgser(x, 0.5), VGAM::dlog(x, 0.5)))
  expect_warning(expect_equal(dlgser(x, 0.0001), VGAM::dlog(x, 0.0001)))
  expect_warning(expect_equal(dlgser(x, 0.9999), VGAM::dlog(x, 0.9999)))
  
  expect_equal(drayleigh(x), VGAM::drayleigh(x))
  expect_equal(prayleigh(x), VGAM::prayleigh(x))
  
  expect_warning(expect_equal(
    extraDistr::dskellam(x, 7, 8),
    VGAM::dskellam(x, 7, 8)
  ))

})


test_that("Check against the parameter values tested in VGAM", {
  
  x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
         0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
  n <- length(x)
  
  skip_on_cran()
  skip_if_not_installed("VGAM")

  expect_warning(expect_equal(dbbinom(x, 10, 0.8, 1.2),
                              VGAM::dbetabinom.ab(x, 10, 0.8, 1.2)))
  
  expect_equal(dpareto(x, 1.9, 2.3), dparetoR(x, 1.9, 2.3))
  expect_equal(dpareto(x, 2.3, 1.9), VGAM::dpareto(x, 1.9, 2.3)) # reverse order of params!
  
  expect_equal(dlst(x, 3, -0.9, 2), dt((x+0.9)/2, df = 3)/2)
  
  expect_equal(dlaplace(x, -0.9, 2), VGAM::dlaplace(x, -0.9, 2))
  
  expect_warning(expect_equal(dbern(x, 0.3), dbinom(x, 1, 0.3)))

})


test_that("Compare with LaplacesDemon package implementations", {
  
  skip_on_cran()
  skip_if_not_installed("LaplacesDemon")
  
  alpha <- runif(5, 0, 3)
  x <- rdirichlet(5000, alpha)
  
  expect_equal(ddirichlet(x, alpha),
               LaplacesDemon::ddirichlet(x, alpha))
  
  alpha <- c(0.0001, 0.0001, 0.0001, 0.0001, 0.0001)
  
  expect_equal(ddirichlet(x, alpha),
               LaplacesDemon::ddirichlet(x, alpha))
  
  alpha <- c(1000, 1000, 1000, 1000, 1000)
  
  expect_equal(ddirichlet(x, alpha),
               LaplacesDemon::ddirichlet(x, alpha))
  
  alpha <- c(1e-4, 10000, 100, 1e-5, 1000)
  
  expect_equal(ddirichlet(x, alpha),
               LaplacesDemon::ddirichlet(x, alpha))
  
  # x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
  #        0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
  # 
  # expect_equal(dbern(x, 0.32), LaplacesDemon::dbern(x, 0.32))
  # expect_equal(dgpd(x, 0.2, 3.4, 1.4), LaplacesDemon::dgpd(x, 0.2, 3.4, 1.4))
  # expect_equal(dhcauchy(x, 2.6), LaplacesDemon::dhalfcauchy(x, 2.6))
  # expect_equal(dhnorm(x, 8.3), LaplacesDemon::dhalfnorm(x, 8.3))
  # expect_equal(dht(x, sigma=3.4, nu=7), LaplacesDemon::dhalft(x, scale=3.4, nu=7))
  # expect_equal(dinvchisq(x, 2.3, 5), LaplacesDemon::dinvchisq(x, 2.3, 5))
  # expect_equal(dinvgamma(x, 7.6, 3), LaplacesDemon::dinvgamma(x, 7.6, 3))
  # expect_equal(dwald(x, -2, 6), LaplacesDemon::dinvgaussian(x, -2, 6))
  # expect_equal(dpareto(x, 5, 6), LaplacesDemon::dpareto(x, 5, 6))
  # 
  # expect_equal(dcat(), LaplacesDemon::dcat())
  
})


test_that("Compare dmnom to dmultinom from base R", {

  n <- 100
  p <- runif(5)
  p <- p/sum(p)
  
  x <- rmnom(5000, n, p)
  
  expect_equal(dmnom(x, n, p), apply(x, 1, dmultinom, n, p))
  
  p <- c(1e-4, 10000, 100, 1e-5, 1000)
  p <- p/sum(p)
  
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

