


test_that("Wrong parameter values in PDF and PMF functions", {
  
  expect_warning(expect_true(is.nan(dbbinom(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(dbbinom(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(dbbinom(1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dbern(1, -1))))
  expect_warning(expect_true(is.nan(dbern(1, 2))))
  
  expect_warning(expect_true(is.nan(dbhatt(1, sigma = -1))))
  expect_warning(expect_true(is.nan(dbhatt(1, a = -1))))
  
  expect_warning(expect_true(is.nan(dbnbinom(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(dbnbinom(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(dbnbinom(1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dbvnorm(1, 1, sd1 = -1))))
  expect_warning(expect_true(is.nan(dbvnorm(1, 1, sd2 = -1))))
  expect_warning(expect_true(is.nan(dbvnorm(1, 1, cor = -2))))
  expect_warning(expect_true(is.nan(dbvnorm(1, 1, cor = 2))))
  
  expect_warning(expect_true(is.nan(dbvpois(1, 1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(dbvpois(1, 1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(dbvpois(1, 1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dcat(1, c(0.5, 1)))))
  expect_warning(expect_true(is.nan(dcat(1, c(1, 0.5)))))
  expect_warning(expect_true(is.nan(dcat(1, c(0.5, 2)))))
  expect_warning(expect_true(is.nan(dcat(1, c(2, 0.5)))))
  expect_warning(expect_true(is.nan(dcat(1, c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(dcat(1, c(0.5, -1)))))
  
  expect_warning(expect_true(is.nan(ddirichlet(c(0.5, 0.5), c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(ddirichlet(c(0.5, 0.5), c(0.5, -1)))))
  
  expect_warning(expect_true(is.nan(ddnorm(1, sd = -1))))
  
  expect_warning(expect_true(is.nan(ddunif(1, min = 10, max = 1))))
  expect_warning(expect_true(is.nan(ddunif(1, min = 0, max = Inf))))
  expect_warning(expect_true(is.nan(ddunif(1, min = -Inf, max = Inf))))
  expect_warning(expect_true(is.nan(ddunif(1, min = Inf, max = -Inf))))
  
  expect_warning(expect_true(is.nan(ddweibull(1, -1, 1)))) 
  expect_warning(expect_true(is.nan(ddweibull(1, 2, 1)))) 
  expect_warning(expect_true(is.nan(ddweibull(1, 0.5, -1))))
  
  expect_warning(expect_true(is.nan(dfatigue(1, -1, 1))))
  expect_warning(expect_true(is.nan(dfatigue(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dfrechet(1, lambda = -1))))
  expect_warning(expect_true(is.nan(dfrechet(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(dgev(1, 1, -1, 1))))
  
  expect_warning(expect_true(is.nan(dgompertz(1, -1, 1))))
  expect_warning(expect_true(is.nan(dgompertz(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dgpd(1, 1, -1, 1))))
  
  expect_warning(expect_true(is.nan(dgpois(1, -1, 1))))
  expect_warning(expect_true(is.nan(dgpois(1, 1, -1))))
  expect_warning(expect_true(is.nan(dgpois(1, 1, scale = 0))))
  
  expect_warning(expect_true(is.nan(dgumbel(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(dhcauchy(1, -1))))
  
  expect_warning(expect_true(is.nan(dhnorm(1, -1))))
  
  expect_warning(expect_true(is.nan(dht(1, 5, -1))))
  
  expect_warning(expect_true(is.nan(dhuber(1, 0, -1, 1))))
  expect_warning(expect_true(is.nan(dhuber(1, 0, 1, -1))))
  
  expect_warning(expect_true(is.nan(dinvgamma(1, -1, 1))))
  expect_warning(expect_true(is.nan(dinvgamma(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dinvchisq(1, -1, 1))))
  expect_warning(expect_true(is.nan(dinvchisq(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dkumar(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(dkumar(0.5, 1, -1))))
  
  expect_warning(expect_true(is.nan(dlaplace(1, 0, -1))))
  
  expect_warning(expect_true(is.nan(dlgser(1, -1))))
  expect_warning(expect_true(is.nan(dlgser(1, 2))))
  
  expect_warning(expect_true(is.nan(dlomax(1, -1, 1))))
  expect_warning(expect_true(is.nan(dlomax(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,1)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/4,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/4,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,1/4)))))
  
  expect_warning(expect_true(is.nan(dmixpois(0, c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(1/3,1,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(1/3,1/3,1)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(1/4,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(1/3,1/4,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(1/3,1/3,1/4)))))
  
  expect_warning(expect_true(is.nan(dmnom(c(1, 1, 1), 3, c(1/2, 1/2, 1/2)))))
  expect_warning(expect_true(is.nan(dmnom(c(1, 1, 1), 3, c(2, 1/3, 1/3)))))
  expect_warning(expect_true(is.nan(dmnom(c(1, 1, 1), 3, c(1/3, 2, 1/3)))))
  expect_warning(expect_true(is.nan(dmnom(c(1, 1, 1), 3, c(1/3, 1/3, 2)))))

  expect_warning(expect_true(is.nan(dmvhyper(c(1, 1, 2), c(2,3,4), 99))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 1, 2), c(-2,3,4), 99))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 1, 2), c(2,-3,4), 99))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 1, 2), c(2,3,-4), 99))))
  
  expect_warning(expect_true(is.nan(dnsbeta(0.5, -1, 1, -2, 2))))
  expect_warning(expect_true(is.nan(dnsbeta(0.5, 1, -1, -2, 2))))
  expect_warning(expect_true(is.nan(dnsbeta(0.5, 1, 1, 2, -2))))
  
  expect_warning(expect_true(is.nan(dnst(1, -2, 0, 1))))
  expect_warning(expect_true(is.nan(dnst(1, 2, 0, -1))))
  
  expect_warning(expect_true(is.nan(dpareto(0, -1, 1))))
  expect_warning(expect_true(is.nan(dpareto(0, 1, -1))))
  
  # expect_warning(expect_true(is.nan(dpower(1, 1, 1)))) # no restrictions
  
  expect_warning(expect_true(is.nan(dprop(1, -10, 0.5))))
  expect_warning(expect_true(is.nan(dprop(1, 10, -1))))
  expect_warning(expect_true(is.nan(dprop(1, 10, 2))))
  
  expect_warning(expect_true(is.nan(drayleigh(0, -1))))
  
  expect_warning(expect_true(is.nan(dskellam(1, -1, 1))))
  expect_warning(expect_true(is.nan(dskellam(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dslash(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(dtnorm(1, 0, -1, -2, 2))))
  expect_warning(expect_true(is.nan(dtnorm(1, 0, 1, 2, -2))))
  expect_warning(expect_true(is.nan(dtnorm(1, 0, 1, 0, 0)))) 
  
  expect_warning(expect_true(is.nan(dtpois(1, lambda = -5, s = 0))))
  expect_warning(expect_true(is.nan(dtpois(1, lambda = -5, s = 6))))
  expect_warning(expect_true(is.nan(dtpois(1, lambda = 5, s = -1))))

  expect_warning(expect_true(is.nan(dtriang(1, 0, 0, 0))))
  expect_warning(expect_true(is.nan(dtriang(1, 1, -1, 0))))
  expect_warning(expect_true(is.nan(dtriang(1, -1, 1, 2))))
  expect_warning(expect_true(is.nan(dtriang(1, -1, 1, -2))))
  
  expect_warning(expect_true(is.nan(dwald(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dzip(1, -1, 0.5))))
  expect_warning(expect_true(is.nan(dzip(1, 1, -1))))
  expect_warning(expect_true(is.nan(dzip(1, 1, 2))))
  
  expect_warning(expect_true(is.nan(dzib(1, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.nan(dzib(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(dzib(1, 1, 0.5, 2))))
  expect_warning(expect_true(is.nan(dzib(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(dzib(1, 1, 0.5, 2))))
  
  expect_warning(expect_true(is.nan(dzinb(1, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.nan(dzinb(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(dzinb(1, 1, 0.5, 2))))
  expect_warning(expect_true(is.nan(dzinb(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(dzinb(1, 1, 0.5, 2))))
  
})





test_that("Wrong parameter values in CDF functions", {
  
  expect_warning(expect_true(is.nan(pbbinom(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(pbbinom(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(pbbinom(1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pbbinom(1, c(-1, 1), c(1, 1), c(1, 1))[1])))
  expect_warning(expect_true(is.nan(pbbinom(1, c(1, 1), c(-1, 1), c(1, 1))[1])))
  expect_warning(expect_true(is.nan(pbbinom(1, c(1, 1), c(1, 1), c(-1, 1))[1])))
  
  expect_warning(expect_true(is.nan(pbbinom(1, c(1, -1), c(1, 1), c(1, 1))[2])))
  expect_warning(expect_true(is.nan(pbbinom(1, c(1, 1), c(1, -1), c(1, 1))[2])))
  expect_warning(expect_true(is.nan(pbbinom(1, c(1, 1), c(1, 1), c(1, -1))[2])))

  expect_warning(expect_true(is.nan(pbern(1, -1))))
  expect_warning(expect_true(is.nan(pbern(1, 2))))
  
  expect_warning(expect_true(is.nan(pbhatt(1, sigma = -1))))
  expect_warning(expect_true(is.nan(pbhatt(1, a = -1))))
  
  expect_warning(expect_true(is.nan(pbnbinom(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(pbnbinom(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(pbnbinom(1, 1, 1, -1))))
  
  # expect_warning(expect_true(is.nan(pbvnorm(1, 1, sd1 = -1))))
  # expect_warning(expect_true(is.nan(pbvnorm(1, 1, sd2 = -1))))
  # expect_warning(expect_true(is.nan(pbvnorm(1, 1, cor = -2))))
  # expect_warning(expect_true(is.nan(pbvnorm(1, 1, cor = 2))))
  
  # expect_warning(expect_true(is.nan(pbvpois(1, 1, -1, 1, 1))))
  # expect_warning(expect_true(is.nan(pbvpois(1, 1, 1, -1, 1))))
  # expect_warning(expect_true(is.nan(pbvpois(1, 1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pcat(1, c(0.5, 1)))))
  expect_warning(expect_true(is.nan(pcat(1, c(1, 0.5)))))
  expect_warning(expect_true(is.nan(pcat(1, c(0.5, 2)))))
  expect_warning(expect_true(is.nan(pcat(1, c(2, 0.5)))))
  expect_warning(expect_true(is.nan(pcat(1, c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(pcat(1, c(0.5, -1)))))
  
  # expect_warning(expect_true(is.nan(pdirichlet(c(0.5, 0.5), c(-1, 0.5)))))
  # expect_warning(expect_true(is.nan(pdirichlet(c(0.5, 0.5), c(0.5, -1)))))
  
  expect_warning(expect_true(is.nan(pdnorm(1, sd = -1))))
  
  expect_warning(expect_true(is.nan(pdunif(1, min = 10, max = 1))))
  expect_warning(expect_true(is.nan(pdunif(1, min = 0, max = Inf))))
  expect_warning(expect_true(is.nan(pdunif(1, min = -Inf, max = Inf))))
  expect_warning(expect_true(is.nan(pdunif(1, min = Inf, max = -Inf))))
  
  expect_warning(expect_true(is.nan(pdweibull(1, -1, 1)))) 
  expect_warning(expect_true(is.nan(pdweibull(1, 2, 1)))) 
  expect_warning(expect_true(is.nan(pdweibull(1, 0.5, -1))))
  
  expect_warning(expect_true(is.nan(pfatigue(1, -1, 1))))
  expect_warning(expect_true(is.nan(pfatigue(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pfrechet(1, lambda = -1))))
  expect_warning(expect_true(is.nan(pfrechet(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(pgev(1, 1, -1, 1))))
  
  expect_warning(expect_true(is.nan(pgompertz(1, -1, 1))))
  expect_warning(expect_true(is.nan(pgompertz(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pgpd(1, 1, -1, 1))))
  
  expect_warning(expect_true(is.nan(pgpois(1, -1, 1))))
  expect_warning(expect_true(is.nan(pgpois(1, 1, -1))))
  expect_warning(expect_true(is.nan(pgpois(1, 1, scale = 0))))
  
  expect_warning(expect_true(is.nan(pgumbel(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(phcauchy(1, -1))))
  
  expect_warning(expect_true(is.nan(phnorm(1, -1))))
  
  expect_warning(expect_true(is.nan(pht(1, 5, -1))))
  
  expect_warning(expect_true(is.nan(phuber(1, 0, -1, 1))))
  expect_warning(expect_true(is.nan(phuber(1, 0, 1, -1))))
  
  expect_warning(expect_true(is.nan(pinvgamma(1, -1, 1))))
  expect_warning(expect_true(is.nan(pinvgamma(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pinvchisq(1, -1, 1))))
  expect_warning(expect_true(is.nan(pinvchisq(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pkumar(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(pkumar(0.5, 1, -1))))
  
  expect_warning(expect_true(is.nan(plaplace(1, 0, -1))))
  
  expect_warning(expect_true(is.nan(plgser(1, -1))))
  expect_warning(expect_true(is.nan(plgser(1, 2))))
  
  expect_warning(expect_true(is.nan(plomax(1, -1, 1))))
  expect_warning(expect_true(is.nan(plomax(1, 1, -1))))
  
  # expect_warning(expect_true(is.nan(pmnom(c(1, 1, 1), 3, c(1/2, 1/2, 1/2)))))
  # expect_warning(expect_true(is.nan(pmnom(c(1, 1, 1), 3, c(2, 1/3, 1/3)))))
  # expect_warning(expect_true(is.nan(pmnom(c(1, 1, 1), 3, c(1/3, 2, 1/3)))))
  # expect_warning(expect_true(is.nan(pmnom(c(1, 1, 1), 3, c(1/3, 1/3, 2)))))
  
  # expect_warning(expect_true(is.nan(pmvhyper(c(1, 1, 2), c(2,3,4), 99))))
  
  
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,1)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/4,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/4,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,1/4)))))
  
  expect_warning(expect_true(is.nan(pmixpois(0, c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(1/3,1,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(1/3,1/3,1)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(1/4,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(1/3,1/4,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(1/3,1/3,1/4)))))
  
  expect_warning(expect_true(is.nan(pnsbeta(0.5, -1, 1, -2, 2))))
  expect_warning(expect_true(is.nan(pnsbeta(0.5, 1, -1, -2, 2))))
  expect_warning(expect_true(is.nan(pnsbeta(0.5, 1, 1, 2, -2))))
  
  expect_warning(expect_true(is.nan(pnst(1, -2, 0, 1))))
  expect_warning(expect_true(is.nan(pnst(1, 2, 0, -1))))
  
  expect_warning(expect_true(is.nan(ppareto(0, -1, 1))))
  expect_warning(expect_true(is.nan(ppareto(0, 1, -1))))
  
  # expect_warning(expect_true(is.nan(ppower(1, 1, 1)))) # no restrictions
  
  expect_warning(expect_true(is.nan(pprop(1, -10, 0.5))))
  expect_warning(expect_true(is.nan(pprop(1, 10, -1))))
  expect_warning(expect_true(is.nan(pprop(1, 10, 2))))
  
  expect_warning(expect_true(is.nan(prayleigh(0, -1))))
  
  # expect_warning(expect_true(is.nan(pskellam(1, -1, 1))))
  # expect_warning(expect_true(is.nan(pskellam(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pslash(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(ptnorm(1, 0, -1, -2, 2))))
  expect_warning(expect_true(is.nan(ptnorm(1, 0, 1, 2, -2))))
  expect_warning(expect_true(is.nan(ptnorm(1, 0, 1, 0, 0))))
  
  expect_warning(expect_true(is.nan(ptpois(1, lambda = -5, s = 0))))
  expect_warning(expect_true(is.nan(ptpois(1, lambda = -5, s = 6))))
  expect_warning(expect_true(is.nan(ptpois(1, lambda = 5, s = -1))))
  
  expect_warning(expect_true(is.nan(ptriang(1, 0, 0, 0))))
  expect_warning(expect_true(is.nan(ptriang(1, 1, -1, 0))))
  expect_warning(expect_true(is.nan(ptriang(1, -1, 1, 2))))
  expect_warning(expect_true(is.nan(ptriang(1, -1, 1, -2))))
  
  expect_warning(expect_true(is.nan(pwald(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pzip(1, -1, 0.5))))
  expect_warning(expect_true(is.nan(pzip(1, 1, -1))))
  expect_warning(expect_true(is.nan(pzip(1, 1, 2))))
  
  expect_warning(expect_true(is.nan(pzib(1, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.nan(pzib(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(pzib(1, 1, 0.5, 2))))
  expect_warning(expect_true(is.nan(pzib(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(pzib(1, 1, 0.5, 2))))
  
  expect_warning(expect_true(is.nan(pzinb(1, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.nan(pzinb(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(pzinb(1, 1, 0.5, 2))))
  expect_warning(expect_true(is.nan(pzinb(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(pzinb(1, 1, 0.5, 2))))
  
})




test_that("Wrong parameter values in inverse CDF functions", {
  
  # expect_warning(expect_true(is.nan(qbbinom(0.5, -1, 1, 1))))
  # expect_warning(expect_true(is.nan(qbbinom(0.5, 1, -1, 1))))
  # expect_warning(expect_true(is.nan(qbbinom(0.5, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(qbern(0.5, -1))))
  expect_warning(expect_true(is.nan(qbern(0.5, 2))))
  
  # expect_warning(expect_true(is.nan(qbhatt(0.5, sigma = -1))))
  # expect_warning(expect_true(is.nan(qbhatt(0.5, a = -1))))
  
  # expect_warning(expect_true(is.nan(qbnbinom(0.5, -1, 1, 1))))
  # expect_warning(expect_true(is.nan(qbnbinom(0.5, 1, -1, 1))))
  # expect_warning(expect_true(is.nan(qbnbinom(0.5, 1, 1, -1))))
  
  # expect_warning(expect_true(is.nan(qbvnorm(0.5, 1, sd1 = -1))))
  # expect_warning(expect_true(is.nan(qbvnorm(0.5, 1, sd2 = -1))))
  # expect_warning(expect_true(is.nan(qbvnorm(0.5, 1, cor = -2))))
  # expect_warning(expect_true(is.nan(qbvnorm(0.5, 1, cor = 2))))
  
  # expect_warning(expect_true(is.nan(qbvpois(0.5, 1, -1, 1, 1))))
  # expect_warning(expect_true(is.nan(qbvpois(0.5, 1, 1, -1, 1))))
  # expect_warning(expect_true(is.nan(qbvpois(0.5, 1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(qcat(0.5, c(0.5, 1)))))
  expect_warning(expect_true(is.nan(qcat(0.5, c(1, 0.5)))))
  expect_warning(expect_true(is.nan(qcat(0.5, c(0.5, 2)))))
  expect_warning(expect_true(is.nan(qcat(0.5, c(2, 0.5)))))
  expect_warning(expect_true(is.nan(qcat(0.5, c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(qcat(0.5, c(0.5, -1)))))
  
  # expect_warning(expect_true(is.nan(qdirichlet(c(0.5, 0.5), c(-1, 0.5)))))
  # expect_warning(expect_true(is.nan(qdirichlet(c(0.5, 0.5), c(0.5, -1)))))
  
  expect_warning(expect_true(is.nan(qdnorm(0.5, sd = -1))))
  
  expect_warning(expect_true(is.nan(qdunif(0.5, min = 10, max = 1))))
  expect_warning(expect_true(is.nan(qdunif(0.5, min = 0, max = Inf))))
  expect_warning(expect_true(is.nan(qdunif(0.5, min = -Inf, max = Inf))))
  expect_warning(expect_true(is.nan(qdunif(0.5, min = Inf, max = -Inf))))
  
  expect_warning(expect_true(is.nan(qdweibull(0.5, -1, 1)))) 
  expect_warning(expect_true(is.nan(qdweibull(0.5, 2, 1)))) 
  expect_warning(expect_true(is.nan(qdweibull(0.5, 0.5, -1))))
  
  expect_warning(expect_true(is.nan(qfatigue(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(qfatigue(0.5, 1, -1))))
  
  expect_warning(expect_true(is.nan(qfrechet(0.5, lambda = -1))))
  expect_warning(expect_true(is.nan(qfrechet(0.5, sigma = -1))))
  
  expect_warning(expect_true(is.nan(qgev(0.5, 1, -1, 1))))
  
  expect_warning(expect_true(is.nan(qgompertz(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(qgompertz(0.5, 1, -1))))
  
  expect_warning(expect_true(is.nan(qgpd(0.5, 1, -1, 1))))
  
  # expect_warning(expect_true(is.nan(qgpois(0.5, -1, 1))))
  # expect_warning(expect_true(is.nan(qgpois(0.5, 1, -1))))
  # expect_warning(expect_true(is.nan(qgpois(0.5, 1, scale = 0))))
  
  expect_warning(expect_true(is.nan(qgumbel(0.5, sigma = -1))))
  
  expect_warning(expect_true(is.nan(qhcauchy(0.5, -1))))
  
  expect_warning(expect_true(is.nan(qhnorm(0.5, -1))))
  
  expect_warning(expect_true(is.nan(qht(0.5, 5, -1))))
  
  expect_warning(expect_true(is.nan(qhuber(0.5, 0, -1, 1))))
  expect_warning(expect_true(is.nan(qhuber(0.5, 0, 1, -1))))
  
  expect_warning(expect_true(is.nan(qinvgamma(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(qinvgamma(0.5, 1, -1))))
  
  expect_warning(expect_true(is.nan(qinvchisq(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(qinvchisq(0.5, 1, -1))))
  
  expect_warning(expect_true(is.nan(qkumar(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(qkumar(0.5, 1, -1))))
  
  expect_warning(expect_true(is.nan(qlaplace(0.5, 0, -1))))
  
  expect_warning(expect_true(is.nan(qlgser(0.5, -1))))
  expect_warning(expect_true(is.nan(qlgser(0.5, 2))))
  
  expect_warning(expect_true(is.nan(qlomax(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(qlomax(0.5, 1, -1))))
  
  # expect_warning(expect_true(is.nan(qmnom(c(0.5, 1, 1), 3, c(1/2, 1/2, 1/2)))))
  # expect_warning(expect_true(is.nan(qmnom(c(0.5, 1, 1), 3, c(2, 1/3, 1/3)))))
  # expect_warning(expect_true(is.nan(qmnom(c(0.5, 1, 1), 3, c(1/3, 2, 1/3)))))
  # expect_warning(expect_true(is.nan(qmnom(c(0.5, 1, 1), 3, c(1/3, 1/3, 2)))))
  
  # expect_warning(expect_true(is.nan(qmvhyper(c(0.5, 1, 2), c(2,3,4), 99))))
  
  expect_warning(expect_true(is.nan(qnsbeta(0.5, -1, 1, -2, 2))))
  expect_warning(expect_true(is.nan(qnsbeta(0.5, 1, -1, -2, 2))))
  expect_warning(expect_true(is.nan(qnsbeta(0.5, 1, 1, 2, -2))))
  
  expect_warning(expect_true(is.nan(qnst(0.5, -2, 0, 1))))
  expect_warning(expect_true(is.nan(qnst(0.5, 2, 0, -1))))
  
  expect_warning(expect_true(is.nan(qpareto(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(qpareto(0.5, 1, -1))))
  
  # expect_warning(expect_true(is.nan(qpower(0.5, 1, 1)))) # no restrictions
  
  expect_warning(expect_true(is.nan(qprop(0.5, -10, 0.5))))
  expect_warning(expect_true(is.nan(qprop(0.5, 10, -1))))
  expect_warning(expect_true(is.nan(qprop(0.5, 10, 2))))
  
  expect_warning(expect_true(is.nan(qrayleigh(0, -1))))
  
  # expect_warning(expect_true(is.nan(qskellam(0.5, -1, 1))))
  # expect_warning(expect_true(is.nan(qskellam(0.5, 1, -1))))
  
  # expect_warning(expect_true(is.nan(qslash(0.5, sigma = -1))))
  
  expect_warning(expect_true(is.nan(qtnorm(0.5, 0, -1, -2, 2))))
  expect_warning(expect_true(is.nan(qtnorm(0.5, 0, 1, 2, -2))))
  expect_warning(expect_true(is.nan(qtnorm(0.5, 0, 1, 0, 0))))
  
  expect_warning(expect_true(is.nan(qtpois(0.5, lambda = -5, s = 0))))
  expect_warning(expect_true(is.nan(qtpois(0.5, lambda = -5, s = 6))))
  expect_warning(expect_true(is.nan(qtpois(0.5, lambda = 5, s = -1))))
  
  expect_warning(expect_true(is.nan(qtriang(0.5, 0, 0, 0))))
  expect_warning(expect_true(is.nan(qtriang(0.5, 1, -1, 0))))
  expect_warning(expect_true(is.nan(qtriang(0.5, -1, 1, 2))))
  expect_warning(expect_true(is.nan(qtriang(0.5, -1, 1, -2))))
  
  # expect_warning(expect_true(is.nan(qwald(0.5, 1, -1))))
  
  expect_warning(expect_true(is.nan(qzip(0.5, -1, 0.5))))
  expect_warning(expect_true(is.nan(qzip(0.5, 1, -1))))
  expect_warning(expect_true(is.nan(qzip(0.5, 1, 2))))
  
  expect_warning(expect_true(is.nan(qzib(0.5, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.nan(qzib(0.5, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(qzib(0.5, 1, 0.5, 2))))
  expect_warning(expect_true(is.nan(qzib(0.5, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(qzib(0.5, 1, 0.5, 2))))
  
  expect_warning(expect_true(is.nan(qzinb(0.5, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.nan(qzinb(0.5, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(qzinb(0.5, 1, 0.5, 2))))
  expect_warning(expect_true(is.nan(qzinb(0.5, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(qzinb(0.5, 1, 0.5, 2))))
  
})




test_that("Wrong parameter values in RNG functions", {
  
  expect_warning(expect_true(is.nan(rbbinom(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(rbbinom(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(rbbinom(1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(rbern(1, -1))))
  expect_warning(expect_true(is.nan(rbern(1, 2))))
  
  expect_warning(expect_true(is.nan(rbhatt(1, sigma = -1))))
  expect_warning(expect_true(is.nan(rbhatt(1, a = -1))))
  
  expect_warning(expect_true(all(is.nan(rbnbinom(1, -1, 1, 1)))))
  expect_warning(expect_true(all(is.nan(rbnbinom(1, 1, -1, 1)))))
  expect_warning(expect_true(all(is.nan(rbnbinom(1, 1, 1, -1)))))
  
  expect_warning(expect_true(all(is.nan(rbvnorm(1, sd1 = -1)))))
  expect_warning(expect_true(all(is.nan(rbvnorm(1, sd2 = -1)))))
  expect_warning(expect_true(all(is.nan(rbvnorm(1, cor = -2)))))
  expect_warning(expect_true(all(is.nan(rbvnorm(1, cor = 2)))))
  
  expect_warning(expect_true(all(is.nan(rbvpois(1, -1, 1, 1)))))
  expect_warning(expect_true(all(is.nan(rbvpois(1, 1, -1, 1)))))
  expect_warning(expect_true(all(is.nan(rbvpois(1, 1, 1, -1)))))
  
  expect_warning(expect_true(is.nan(rcat(1, c(0.5, 1)))))
  expect_warning(expect_true(is.nan(rcat(1, c(1, 0.5)))))
  expect_warning(expect_true(is.nan(rcat(1, c(0.5, 2)))))
  expect_warning(expect_true(is.nan(rcat(1, c(2, 0.5)))))
  expect_warning(expect_true(is.nan(rcat(1, c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(rcat(1, c(0.5, -1)))))
  
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, 1, 0.5, 0.5), byrow = T, ncol = 2))[1])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(1, 0.5, 0.5, 0.5), byrow = T, ncol = 2))[1])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, 2, 0.5, 0.5), byrow = T, ncol = 2))[1])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(2, 0.5, 0.5, 0.5), byrow = T, ncol = 2))[1])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(-1, 0.5, 0.5, 0.5), byrow = T, ncol = 2))[1])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, -1, 0.5, 0.5), byrow = T, ncol = 2))[1])))
  
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, 0.5, 0.5, 1), byrow = T, ncol = 2))[2])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, 0.5, 1, 0.5), byrow = T, ncol = 2))[2])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, 0.5, 0.5, 2), byrow = T, ncol = 2))[2])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, 0.5, 2, 0.5), byrow = T, ncol = 2))[2])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, 0.5, -1, 0.5), byrow = T, ncol = 2))[2])))
  expect_warning(expect_true(is.nan(rcat(2, matrix(c(0.5, 0.5, 0.5, -1), byrow = T, ncol = 2))[2])))
  
  expect_warning(expect_true(all(is.nan(rdirichlet(1, c(-1, 0.5))))))
  expect_warning(expect_true(all(is.nan(rdirichlet(1, c(0.5, -1))))))
  
  expect_warning(expect_true(is.nan(rdnorm(1, sd = -1))))
  
  expect_warning(expect_true(is.nan(rdunif(1, min = 10, max = 1))))
  expect_warning(expect_true(is.nan(rdunif(1, min = 0, max = Inf))))
  expect_warning(expect_true(is.nan(rdunif(1, min = -Inf, max = Inf))))
  expect_warning(expect_true(is.nan(rdunif(1, min = Inf, max = -Inf))))
  
  expect_warning(expect_true(is.nan(rdweibull(1, -1, 1)))) 
  expect_warning(expect_true(is.nan(rdweibull(1, 2, 1)))) 
  expect_warning(expect_true(is.nan(rdweibull(1, 0.5, -1))))
  
  expect_warning(expect_true(is.nan(rfatigue(1, -1, 1))))
  expect_warning(expect_true(is.nan(rfatigue(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(rfrechet(1, lambda = -1))))
  expect_warning(expect_true(is.nan(rfrechet(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(rgev(1, 1, -1, 1))))
  
  expect_warning(expect_true(is.nan(rgompertz(1, -1, 1))))
  expect_warning(expect_true(is.nan(rgompertz(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(rgpd(1, 1, -1, 1))))
  
  expect_warning(expect_true(is.nan(rgpois(1, -1, 1))))
  expect_warning(expect_true(is.nan(rgpois(1, 1, -1))))
  expect_warning(expect_true(is.nan(rgpois(1, 1, scale = 0))))
  
  expect_warning(expect_true(is.nan(rgumbel(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(rhcauchy(1, -1))))
  
  expect_warning(expect_true(is.nan(rhnorm(1, -1))))
  
  expect_warning(expect_true(is.nan(rht(1, 5, -1))))
  
  expect_warning(expect_true(is.nan(rhuber(1, 0, -1, 1))))
  expect_warning(expect_true(is.nan(rhuber(1, 0, 1, -1))))
  
  expect_warning(expect_true(is.nan(rinvgamma(1, -1, 1))))
  expect_warning(expect_true(is.nan(rinvgamma(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(rinvchisq(1, -1, 1))))
  expect_warning(expect_true(is.nan(rinvchisq(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(rkumar(1, -1, 1))))
  expect_warning(expect_true(is.nan(rkumar(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(rlaplace(1, 0, -1))))
  
  expect_warning(expect_true(is.nan(rlgser(1, -1))))
  expect_warning(expect_true(is.nan(rlgser(1, 2))))
  
  expect_warning(expect_true(is.nan(rlomax(1, -1, 1))))
  expect_warning(expect_true(is.nan(rlomax(1, 1, -1))))
  
  
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(1,2,3), c(1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,1,1/3)))))
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,1/3,1)))))
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/4,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,1/4,1/3)))))
  expect_warning(expect_true(is.nan(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,1/3,1/4)))))
  
  expect_warning(expect_true(is.nan(rmixpois(1, c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixpois(1, c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixpois(1, c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixpois(1, c(1,2,3), c(1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixpois(1, c(1,2,3), c(1/3,1,1/3)))))
  expect_warning(expect_true(is.nan(rmixpois(1, c(1,2,3), c(1/3,1/3,1)))))
  expect_warning(expect_true(is.nan(rmixpois(1, c(1,2,3), c(1/4,1/3,1/3)))))
  expect_warning(expect_true(is.nan(rmixpois(1, c(1,2,3), c(1/3,1/4,1/3)))))
  expect_warning(expect_true(is.nan(rmixpois(1, c(1,2,3), c(1/3,1/3,1/4)))))
  
  expect_warning(expect_true(all(is.nan(rmnom(c(1, 1, 1), 3, c(1/2, 1/2, 1/2))))))
  expect_warning(expect_true(all(is.nan(rmnom(c(1, 1, 1), 3, c(2, 1/3, 1/3))))))
  expect_warning(expect_true(all(is.nan(rmnom(c(1, 1, 1), 3, c(1/3, 2, 1/3))))))
  expect_warning(expect_true(all(is.nan(rmnom(c(1, 1, 1), 3, c(1/3, 1/3, 2))))))
  
  expect_warning(expect_true(all(is.nan(rmvhyper(1, c(2,3,4), 99)))))
  expect_warning(expect_true(all(is.nan(rmvhyper(1, c(-2,3,4), 99)))))
  expect_warning(expect_true(all(is.nan(rmvhyper(1, c(2,-3,4), 99)))))
  expect_warning(expect_true(all(is.nan(rmvhyper(1, c(2,3,-4), 99)))))
  
  expect_warning(expect_true(is.nan(rnsbeta(1, -1, 1, -2, 2))))
  expect_warning(expect_true(is.nan(rnsbeta(1, 1, -1, -2, 2))))
  expect_warning(expect_true(is.nan(rnsbeta(1, 1, 1, 2, -2))))
  
  expect_warning(expect_true(is.nan(rnst(1, -2, 0, 1))))
  expect_warning(expect_true(is.nan(rnst(1, 2, 0, -1))))
  
  expect_warning(expect_true(is.nan(rpareto(1, -1, 1))))
  expect_warning(expect_true(is.nan(rpareto(1, 1, -1))))
  
  # expect_warning(expect_true(is.nan(rpower(1, 1, 1)))) # no restrictions
  
  expect_warning(expect_true(is.nan(rprop(1, -10, 0.5))))
  expect_warning(expect_true(is.nan(rprop(1, 10, -1))))
  expect_warning(expect_true(is.nan(rprop(1, 10, 2))))
  
  expect_warning(expect_true(is.nan(rrayleigh(1, -1))))
  
  expect_warning(expect_true(is.nan(rskellam(1, -1, 1))))
  expect_warning(expect_true(is.nan(rskellam(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(rslash(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(rtnorm(1, 0, -1, -2, 2))))
  expect_warning(expect_true(is.nan(rtnorm(1, 0, 1, 2, -2))))
  expect_warning(expect_true(is.nan(rtnorm(1, 0, 1, 0, 0))))
  
  expect_warning(expect_true(is.nan(rtpois(1, lambda = -5, s = 0))))
  expect_warning(expect_true(is.nan(rtpois(1, lambda = -5, s = 6))))
  expect_warning(expect_true(is.nan(rtpois(1, lambda = 5, s = -1))))
  
  expect_warning(expect_true(is.nan(rtriang(1, 0, 0, 0))))
  expect_warning(expect_true(is.nan(rtriang(1, 1, -1, 0))))
  expect_warning(expect_true(is.nan(rtriang(1, -1, 1, 2))))
  expect_warning(expect_true(is.nan(rtriang(1, -1, 1, -2))))
  
  expect_warning(expect_true(is.nan(rwald(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(rzip(1, -1, 0.5))))
  expect_warning(expect_true(is.nan(rzip(1, 1, -1))))
  expect_warning(expect_true(is.nan(rzip(1, 1, 2))))
  
  expect_warning(expect_true(is.nan(rzib(1, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.nan(rzib(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(rzib(1, 1, 0.5, 2))))
  expect_warning(expect_true(is.nan(rzib(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(rzib(1, 1, 0.5, 2))))
  
  expect_warning(expect_true(is.nan(rzinb(1, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.nan(rzinb(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(rzinb(1, 1, 0.5, 2))))
  expect_warning(expect_true(is.nan(rzinb(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.nan(rzinb(1, 1, 0.5, 2))))
  
})