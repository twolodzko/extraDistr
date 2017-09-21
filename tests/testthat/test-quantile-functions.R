

test_that("Zeros in quantile functions", {
  
  expect_true(!is.nan(qbetapr(0, 1, 1, 1)))
  expect_true(!is.nan(qfatigue(0, 1)))
  expect_true(!is.nan(qcat(0, c(0.5, 0.5))))
  expect_true(!is.nan(qdweibull(0, 0.5, 1)))  
  expect_true(!is.nan(qfrechet(0)))
  expect_true(!is.nan(qgev(0, 1, 1, 1)))
  expect_true(!is.nan(qgompertz(0, 1, 1)))
  expect_true(!is.nan(qgpd(0, 1, 1, 1)))
  expect_true(!is.nan(qgumbel(0)))
  expect_true(!is.nan(qhuber(0)))
  expect_true(!is.nan(qhcauchy(0, 1)))
  expect_true(!is.nan(qhnorm(0, 1)))
  expect_true(!is.nan(qht(0, 5, 1)))
  expect_true(!is.nan(qinvgamma(0, 1, 1)))
  expect_true(!is.nan(qlaplace(0)))
  expect_true(!is.nan(qlgser(0, 0.5)))
  expect_true(!is.nan(qlomax(0, 1, 1)))
  expect_true(!is.nan(qnhyper(0, 60, 35, 15)))
  expect_true(!is.nan(qlst(0, df = 2)))
  expect_true(!is.nan(qpareto(0)))
  expect_true(!is.nan(qpower(0, 1, 1)))
  expect_true(!is.nan(qprop(0, 10, 0.5)))
  expect_true(!is.nan(qrayleigh(0)))
  expect_true(!is.nan(qtlambda(0, 0.5)))
  expect_true(!is.nan(qtbinom(0, 100, 0.83, 76, 86)))
  expect_true(!is.nan(qzip(0, 1, 0.5)))
  expect_true(!is.nan(qzib(0, 1, 1, 0.5)))
  expect_true(!is.nan(qzinb(0, 1, 1, 0.5)))
  
  expect_true(!is.nan(qtpois(0, lambda = 5, a = 0)))
  expect_true(!is.nan(qtpois(0, lambda = 5, a = 6)))
  
  expect_true(!is.nan(pdgamma(0, 9, 1)))
  expect_true(!is.nan(pdnorm(0, 1, 2)))
  
})

test_that("Ones in quantile functions", {
  
  expect_true(!is.nan(qbetapr(1, 1, 1, 1)))
  expect_true(!is.nan(qfatigue(1, 1)))
  expect_true(!is.nan(qcat(1, c(0.5, 0.5))))
  expect_true(!is.nan(qdweibull(1, 0.5, 1)))  
  expect_true(!is.nan(qfrechet(1)))
  expect_true(!is.nan(qgev(1, 1, 1, 1)))
  expect_true(!is.nan(qgompertz(1, 1, 1)))
  expect_true(!is.nan(qgpd(1, 1, 1, 1)))
  expect_true(!is.nan(qgumbel(1)))
  expect_true(!is.nan(qhuber(1)))
  expect_true(!is.nan(qhcauchy(1, 1)))
  expect_true(!is.nan(qhnorm(1, 1)))
  expect_true(!is.nan(qht(1, 5, 1)))
  expect_true(!is.nan(qinvgamma(1, 1, 1)))
  expect_true(!is.nan(qlaplace(1)))
  expect_true(!is.nan(qlgser(1, 0.5)))
  expect_true(!is.nan(qlomax(1, 1, 1)))
  expect_true(!is.nan(qnhyper(1, 60, 35, 15)))
  expect_true(!is.nan(qlst(1, df = 2)))
  expect_true(!is.nan(qpareto(1)))
  expect_true(!is.nan(qpower(1, 1, 1)))
  expect_true(!is.nan(qprop(1, 10, 0.5)))
  expect_true(!is.nan(qrayleigh(1)))
  expect_true(!is.nan(qtlambda(1, 0.5)))
  expect_true(!is.nan(qtbinom(1, 100, 0.83, 76, 86)))
  expect_true(!is.nan(qzip(1, 1, 0.5)))
  expect_true(!is.nan(qzib(1, 1, 1, 0.5)))
  expect_true(!is.nan(qzinb(1, 1, 1, 0.5)))
  
  expect_true(!is.nan(qtpois(1, lambda = 5, a = 0)))
  expect_true(!is.nan(qtpois(1, lambda = 5, a = 6)))
  
})



test_that("Checking p = F(F^-1(p))", {
  
  pp <- seq(0, 1, by = 0.001)
  
  expect_equal(pp, pbetapr(qbetapr(pp, 1, 1, 1), 1, 1, 1))
  expect_equal(pp, pfatigue(qfatigue(pp, 1), 1))
  expect_equal(pp, pfrechet(qfrechet(pp)))
  expect_equal(pp, pgev(qgev(pp, 1, 1, 1), 1, 1, 1))
  expect_equal(pp, pgompertz(qgompertz(pp, 1, 1), 1, 1))
  expect_equal(pp, pgpd(qgpd(pp, 1, 1, 1), 1, 1, 1))
  expect_equal(pp, pgumbel(qgumbel(pp)))
  expect_equal(pp, phuber(qhuber(pp)))
  
  expect_equal(pp, phcauchy(qhcauchy(pp, 1), 1))
  expect_equal(pp, phnorm(qhnorm(pp, 1), 1))
  expect_equal(pp, pht(qht(pp, 5, 1), 5, 1))

  expect_equal(pp, pinvgamma(qinvgamma(pp, 1, 1), 1, 1))
  expect_equal(pp, plaplace(qlaplace(pp)))
  expect_equal(pp, plomax(qlomax(pp, 1, 1), 1, 1))
  expect_equal(pp, plst(qlst(pp, df = 2), df = 2))
  expect_equal(pp, ppareto(qpareto(pp)))
  expect_equal(pp, ppower(qpower(pp, 1, 1), 1, 1))
  expect_equal(pp, pprop(qprop(pp, 10, 0.5), 10, 0.5))
  expect_equal(pp, prayleigh(qrayleigh(pp)))

})



