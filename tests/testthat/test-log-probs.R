

test_that("Check if log-probabilities are logs of probabilities (PMF's and PDF's)", {
  
  x <- c(-Inf, -100, -10, -5, -1, -0.5, 0, 0.5, 1, 5, 10, 100, Inf)
  
  expect_equal(suppressWarnings(dbbinom(x, 1, 1, 1, log = TRUE)),
               log(suppressWarnings(dbbinom(x, 1, 1, 1))))
  expect_equal(suppressWarnings(dbern(x, 0.5, log = TRUE)),
               log(suppressWarnings(dbern(x, 0.5))))
  expect_equal(suppressWarnings(dbetapr(x, 1, 1, 1, log = TRUE)),
               log(suppressWarnings(dbetapr(x, 1, 1, 1))))
  expect_equal(dbhatt(x, sigma = 1, log = TRUE),
               log(dbhatt(x, sigma = 1)))
  expect_equal(suppressWarnings(dbnbinom(x, 1, 1, 1, log = TRUE)),
               log(suppressWarnings(dbnbinom(x, 1, 1, 1))))
  expect_equal(dbvnorm(x, x, sd1 = 1, log = TRUE),
               log(dbvnorm(x, x, sd1 = 1)))
  expect_equal(suppressWarnings(dbvpois(x, x, 1, 1, 1, log = TRUE)),
               log(suppressWarnings(dbvpois(x, x, 1, 1, 1))))
  expect_equal(suppressWarnings(dcat(x, c(0.5, 0.5), log = TRUE)),
               log(suppressWarnings(dcat(x, c(0.5, 0.5)))))
  expect_equal(ddirichlet(c(0.5, 0.5), c(1, 0.5), log = TRUE),
               log(ddirichlet(c(0.5, 0.5), c(1, 0.5))))
  expect_equal(suppressWarnings(ddlaplace(x, 0, scale = 0.5, log = TRUE)),
               log(suppressWarnings(ddlaplace(x, 0, scale = 0.5))))
  expect_equal(suppressWarnings(ddnorm(x, sd = 1, log = TRUE)),
               log(suppressWarnings(ddnorm(x, sd = 1))))
  expect_equal(suppressWarnings(ddgamma(x, 9, 1, log = TRUE)),
               log(suppressWarnings(ddgamma(x, 9, 1))))
  expect_equal(suppressWarnings(ddunif(x, min = 10, max = 100, log = TRUE)),
               log(suppressWarnings(ddunif(x, min = 10, max = 100))))
  expect_equal(suppressWarnings(ddweibull(x, 0.5, 1, log = TRUE)),
               log(suppressWarnings(ddweibull(x, 0.5, 1))))
  expect_equal(dfatigue(x, 1, 1, log = TRUE),
               log(dfatigue(x, 1, 1)))
  expect_equal(dfrechet(x, lambda = 1, log = TRUE),
               log(dfrechet(x, lambda = 1)))
  expect_equal(dgev(x, 1, 1, 1, log = TRUE),
               log(dgev(x, 1, 1, 1)))
  # expect_equal(dgompertz(x, 1, 1, log = TRUE),
  #              log(dgompertz(x, 1, 1)))
  expect_equal(dgpd(x, 1, 1, 1, log = TRUE),
               log(dgpd(x, 1, 1, 1)))
  expect_equal(suppressWarnings(dgpois(x, 1, 1, log = TRUE)),
               log(suppressWarnings(dgpois(x, 1, 1))))
  # expect_equal(dgumbel(x, sigma = 1, log = TRUE),
  #              log(dgumbel(x, sigma = 1)))
  expect_equal(dhcauchy(x, 1, log = TRUE),
               log(dhcauchy(x, 1)))
  # expect_equal(dhnorm(x, 1, log = TRUE),
  #              log(dhnorm(x, 1)))  # numerical precission
  expect_equal(dht(x, 5, 1, log = TRUE),
               log(dht(x, 5, 1)))
  expect_equal(dhuber(x, 0, 1, 1, log = TRUE),
               log(dhuber(x, 0, 1, 1)))
  expect_equal(dinvgamma(x, 1, 1, log = TRUE),
               log(dinvgamma(x, 1, 1)))
  expect_equal(dinvchisq(x, 1, 1, log = TRUE),
               log(dinvchisq(x, 1, 1)))
  expect_equal(dkumar(x, 1, 1, log = TRUE),
               log(dkumar(x, 1, 1)))
  expect_equal(dlaplace(x, 0, 1, log = TRUE),
               log(dlaplace(x, 0, 1)))
  expect_equal(suppressWarnings(dlgser(x, 0.5, log = TRUE)),
               log(suppressWarnings(dlgser(x, 0.5))))
  expect_equal(dlomax(x, 1, 1, log = TRUE),
               log(dlomax(x, 1, 1)))
  expect_equal(dmixnorm(x, c(1,2,3), c(1,2,3), c(1/3,1/3,1/3), log = TRUE),
               log(dmixnorm(x, c(1,2,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_equal(suppressWarnings(dmixpois(x, c(1,2,3), c(1/3,1/3,1/3), log = TRUE)),
               log(suppressWarnings(dmixpois(x, c(1,2,3), c(1/3,1/3,1/3)))))
  expect_equal(suppressWarnings(dnhyper(x, 60, 35, 15, log = TRUE)),
               log(suppressWarnings(dnhyper(x, 60, 35, 15))))
  expect_equal(ddirmnom(c(1, 1, 1), 2, c(1, 1, 1), log = TRUE),
               log(ddirmnom(c(1, 1, 1), 2, c(1, 1, 1))))
  expect_equal(dmnom(c(1, 1, 1), 2, c(1/3, 1/3, 1/3), log = TRUE),
               log(dmnom(c(1, 1, 1), 2, c(1/3, 1/3, 1/3))))
  expect_equal(dmvhyper(c(1, 2, 2), c(2,3,4), 5, log = TRUE),
               log(dmvhyper(c(1, 2, 2), c(2,3,4), 5)))
  expect_equal(dnsbeta(x, 1, 1, -2, 2, log = TRUE),
               log(dnsbeta(x, 1, 1, -2, 2)))
  expect_equal(dlst(x, 2, 0, 1, log = TRUE),
               log(dlst(x, 2, 0, 1)))
  expect_equal(dpareto(x, 1, 1, log = TRUE),
               log(dpareto(x, 1, 1)))
  expect_equal(dprop(x, 10, 0.5, log = TRUE),
               log(dprop(x, 10, 0.5)))
  # expect_equal(drayleigh(x, 1, log = TRUE),
  #              log(drayleigh(x, 1)))
  expect_equal(suppressWarnings(dskellam(x, 1, 1, log = TRUE)),
               log(suppressWarnings(dskellam(x, 1, 1))))
  expect_equal(dsgomp(x, 0.4, 1, log = TRUE),
               log(dsgomp(x, 0.4, 1)))
  expect_equal(dslash(x, sigma = 1, log = TRUE),
               log(dslash(x, sigma = 1)))
  expect_equal(dtnorm(x, 0, 1, 1, 2, log = TRUE), 
               log(dtnorm(x, 0, 1, 1, 2)))
  expect_equal(suppressWarnings(dtpois(x, lambda = 25, a = 0, log = TRUE)),
               log(suppressWarnings(dtpois(x, lambda = 25, a = 0))))
  expect_equal(suppressWarnings(dtbinom(x, 100, 0.67, a = 60, b = 70, log = TRUE)),
               log(suppressWarnings(dtbinom(x, 100, 0.67, a = 60, b = 70))))
  expect_equal(dtriang(x, 1, 2, 1.5, log = TRUE),
               log(dtriang(x, 1, 2, 1.5)))
  expect_equal(dwald(x, 1, 1, log = TRUE),
               log(dwald(x, 1, 1)))
  expect_equal(suppressWarnings(dzip(x, 1, 0.5, log = TRUE)),
               log(suppressWarnings(dzip(x, 1, 0.5))))
  expect_equal(suppressWarnings(dzib(x, 1, 0.5, 0.5, log = TRUE)),
               log(suppressWarnings(dzib(x, 1, 0.5, 0.5))))
  expect_equal(suppressWarnings(dzinb(x, 1, 0.5, 0.5, log = TRUE)),
               log(suppressWarnings(dzinb(x, 1, 0.5, 0.5))))
  
})


test_that("Check if log-probabilities are logs of probabilities (CDF's)", {
  
  x <- c(-Inf, -100, -10, -5, -1, -0.5, 0, 0.5, 1, 5, 10, 100, Inf)
  
  expect_equal(suppressWarnings(pbbinom(x, 1, 1, 1, log.p = TRUE)),
               log(suppressWarnings(pbbinom(x, 1, 1, 1))))
  expect_equal(suppressWarnings(pbern(x, 0.5, log.p = TRUE)),
               log(suppressWarnings(pbern(x, 0.5))))
  expect_equal(suppressWarnings(pbetapr(x, 1, 1, 1, log.p = TRUE)),
               log(suppressWarnings(pbetapr(x, 1, 1, 1))))
  expect_equal(pbhatt(x, sigma = 1, log.p = TRUE),
               log(pbhatt(x, sigma = 1)))
  expect_equal(suppressWarnings(pbnbinom(x, 1, 1, 1, log.p = TRUE)),
               log(suppressWarnings(pbnbinom(x, 1, 1, 1))))
  expect_equal(suppressWarnings(pcat(x, c(0.5, 0.5), log.p = TRUE)),
               log(suppressWarnings(pcat(x, c(0.5, 0.5)))))
  expect_equal(suppressWarnings(pdlaplace(x, 0, scale = 0.5, log.p = TRUE)),
               log(suppressWarnings(pdlaplace(x, 0, scale = 0.5))))
  # expect_equal(suppressWarnings(pdnorm(x, sd = 1, log.p = TRUE)),
  #              log(suppressWarnings(pdnorm(x, sd = 1))))
  expect_equal(suppressWarnings(pdgamma(x, 9, 1, log.p = TRUE)),
               log(suppressWarnings(pdgamma(x, 9, 1))))
  expect_equal(suppressWarnings(pdunif(x, min = 10, max = 100, log.p = TRUE)),
               log(suppressWarnings(pdunif(x, min = 10, max = 100))))
  expect_equal(suppressWarnings(pdweibull(x, 0.5, 1, log.p = TRUE)),
               log(suppressWarnings(pdweibull(x, 0.5, 1))))
  expect_equal(pfatigue(x, 1, 1, log.p = TRUE),
               log(pfatigue(x, 1, 1)))
  expect_equal(pfrechet(x, lambda = 1, log.p = TRUE),
               log(pfrechet(x, lambda = 1)))
  expect_equal(pgev(x, 1, 1, 1, log.p = TRUE),
               log(pgev(x, 1, 1, 1)))
  # expect_equal(pgompertz(x, 1, 1, log.p = TRUE),
  #              log(pgompertz(x, 1, 1)))
  expect_equal(pgpd(x, 1, 1, 1, log.p = TRUE),
               log(pgpd(x, 1, 1, 1)))
  expect_equal(suppressWarnings(pgpois(x, 1, 1, log.p = TRUE)),
               log(suppressWarnings(pgpois(x, 1, 1))))
  # expect_equal(pgumbel(x, sigma = 1, log.p = TRUE),
  #              log(pgumbel(x, sigma = 1)))
  expect_equal(phcauchy(x, 1, log.p = TRUE),
               log(phcauchy(x, 1)))
  expect_equal(phnorm(x, 1, log.p = TRUE),
               log(phnorm(x, 1)))
  expect_equal(pht(x, 5, 1, log.p = TRUE),
               log(pht(x, 5, 1)))
  expect_equal(phuber(x, 0, 1, 1, log.p = TRUE),
               log(phuber(x, 0, 1, 1)))
  expect_equal(pinvgamma(x, 1, 1, log.p = TRUE),
               log(pinvgamma(x, 1, 1)))
  expect_equal(pinvchisq(x, 1, 1, log.p = TRUE),
               log(pinvchisq(x, 1, 1)))
  expect_equal(pkumar(x, 1, 1, log.p = TRUE),
               log(pkumar(x, 1, 1)))
  expect_equal(plaplace(x, 0, 1, log.p = TRUE),
               log(plaplace(x, 0, 1)))
  expect_equal(suppressWarnings(plgser(x, 0.5, log.p = TRUE)),
               log(suppressWarnings(plgser(x, 0.5))))
  expect_equal(plomax(x, 1, 1, log.p = TRUE),
               log(plomax(x, 1, 1)))
  expect_equal(pmixnorm(x, c(1,2,3), c(1,2,3), c(1/3,1/3,1/3), log.p = TRUE),
               log(pmixnorm(x, c(1,2,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_equal(suppressWarnings(pmixpois(x, c(1,2,3), c(1/3,1/3,1/3), log.p = TRUE)),
               log(suppressWarnings(pmixpois(x, c(1,2,3), c(1/3,1/3,1/3)))))
  expect_equal(suppressWarnings(pnhyper(x, 60, 35, 15, log.p = TRUE)),
               log(suppressWarnings(pnhyper(x, 60, 35, 15))))
  expect_equal(pnsbeta(x, 1, 1, -2, 2, log.p = TRUE),
               log(pnsbeta(x, 1, 1, -2, 2)))
  expect_equal(plst(x, 2, 0, 1, log.p = TRUE),
               log(plst(x, 2, 0, 1)))
  expect_equal(ppareto(x, 1, 1, log.p = TRUE),
               log(ppareto(x, 1, 1)))
  expect_equal(pprop(x, 10, 0.5, log.p = TRUE),
               log(pprop(x, 10, 0.5)))
  # expect_equal(prayleigh(x, 1, log.p = TRUE),
  #              log(prayleigh(x, 1)))
  expect_equal(psgomp(x, 0.4, 1, log.p = TRUE),
               log(psgomp(x, 0.4, 1)))
  expect_equal(pslash(x, sigma = 1, log.p = TRUE),
               log(pslash(x, sigma = 1)))
  expect_equal(ptnorm(x, 0, 1, 1, 2, log.p = TRUE), 
               log(ptnorm(x, 0, 1, 1, 2)))
  expect_equal(suppressWarnings(ptpois(x, lambda = 25, a = 0, log.p = TRUE)),
               log(suppressWarnings(ptpois(x, lambda = 25, a = 0))))
  expect_equal(suppressWarnings(ptbinom(x, 100, 0.67, a = 60, b = 70, log.p = TRUE)),
               log(suppressWarnings(ptbinom(x, 100, 0.67, a = 60, b = 70))))
  expect_equal(ptriang(x, 1, 2, 1.5, log.p = TRUE),
               log(ptriang(x, 1, 2, 1.5)))
  expect_equal(pwald(x, 1, 1, log.p = TRUE),
               log(pwald(x, 1, 1)))
  expect_equal(suppressWarnings(pzip(x, 1, 0.5, log.p = TRUE)),
               log(suppressWarnings(pzip(x, 1, 0.5))))
  expect_equal(suppressWarnings(pzib(x, 1, 0.5, 0.5, log.p = TRUE)),
               log(suppressWarnings(pzib(x, 1, 0.5, 0.5))))
  expect_equal(suppressWarnings(pzinb(x, 1, 0.5, 0.5, log.p = TRUE)),
               log(suppressWarnings(pzinb(x, 1, 0.5, 0.5))))
  
})


