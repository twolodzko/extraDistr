


test_that("Wrong parameter values in PDF and PMF functions", {
  
  expect_warning(expect_true(is.nan(dbbinom(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(dbbinom(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(dbbinom(1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dbern(1, -1))))
  expect_warning(expect_true(is.nan(dbern(1, 2))))
  
  expect_warning(expect_true(is.nan(dbetapr(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(dbetapr(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(dbetapr(1, 1, 1, -1))))
  
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

  expect_warning(expect_true(is.nan(dcat(1, c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(dcat(1, c(0.5, -1)))))
  
  expect_warning(expect_true(is.nan(ddirichlet(c(0.5, 0.5), c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(ddirichlet(c(0.5, 0.5), c(0.5, -1)))))
  
  expect_warning(expect_true(is.nan(ddlaplace(1, 0, scale = -1))))
  expect_warning(expect_true(is.nan(ddlaplace(1, 0, scale = 2))))
  
  expect_warning(expect_true(is.nan(ddnorm(1, sd = -1))))
  
  expect_warning(expect_true(is.nan(ddgamma(1, -9, 1))))
  expect_warning(expect_true(is.nan(ddgamma(1, 9, -1))))
  
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
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(-1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,-1,1/3)))))
  expect_warning(expect_true(is.nan(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,-1)))))
  
  expect_warning(expect_true(is.nan(dmixpois(0, c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(-1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(1/3,-1,1/3)))))
  expect_warning(expect_true(is.nan(dmixpois(0, c(1,2,3), c(1/3,1/3,-1)))))
  
  expect_warning(expect_true(is.nan(dnhyper(1, 60.5, 35, 15))))
  expect_warning(expect_true(is.nan(dnhyper(1, 60, 35.5, 15))))
  expect_warning(expect_true(is.nan(dnhyper(1, 60, 35, 15.5))))
  expect_warning(expect_true(is.nan(dnhyper(1, -60, 35, 15))))
  expect_warning(expect_true(is.nan(dnhyper(1, 60, -35, 15))))
  expect_warning(expect_true(is.nan(dnhyper(1, 60, 35, -15))))
  expect_warning(expect_true(is.nan(dnhyper(1, 60, 35, 40))))
  
  expect_warning(expect_true(is.nan(dnhyper(1, 60.5, 35, 15))))
  expect_warning(expect_true(is.nan(dnhyper(1, 60.5, 35, 15))))
  expect_warning(expect_true(is.nan(dnhyper(1, 60.5, 35, 15))))
  expect_warning(expect_true(is.nan(ddirmnom(c(1, 1, 1), 1.5, c(1, 1, 1)))))
  expect_warning(expect_true(is.nan(ddirmnom(c(1, 1, 1), -3, c(1, 1, 1)))))
  expect_warning(expect_true(is.nan(ddirmnom(c(1, 1, 1), 3, c(-1, 1, 1)))))
  expect_warning(expect_true(is.nan(ddirmnom(c(1, 1, 1), 3, c(1, -1, 1)))))
  expect_warning(expect_true(is.nan(ddirmnom(c(1, 1, 1), 3, c(1, 1, -1)))))
  
  expect_warning(expect_true(is.nan(dmnom(c(1, 1, 1), 1.5, c(1/3, 1/3, 1/3)))))
  expect_warning(expect_true(is.nan(dmnom(c(1, 1, 1), 3, c(-1, 1/3, 1/3)))))
  expect_warning(expect_true(is.nan(dmnom(c(1, 1, 1), 3, c(1/3, -1, 1/3)))))
  expect_warning(expect_true(is.nan(dmnom(c(1, 1, 1), 3, c(1/3, 1/3, -1)))))

  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(2.5,3,4), 5))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(2,3.5,4), 5))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(2,3,4.5), 5))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(2,3,4), 5.5))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(-2,3,4), 5))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(2,-3,4), 5))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(2,3,-4), 5))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(2,3,4), -5))))
  expect_warning(expect_true(is.nan(dmvhyper(c(1, 2, 2), c(2,3,4), 85))))
  
  expect_warning(expect_true(is.nan(dnsbeta(0.5, -1, 1, -2, 2))))
  expect_warning(expect_true(is.nan(dnsbeta(0.5, 1, -1, -2, 2))))
  expect_warning(expect_true(is.nan(dnsbeta(0.5, 1, 1, 2, -2))))
  
  expect_warning(expect_true(is.nan(dlst(1, -2, 0, 1))))
  expect_warning(expect_true(is.nan(dlst(1, 2, 0, -1))))
  
  expect_warning(expect_true(is.nan(dpareto(0, -1, 1))))
  expect_warning(expect_true(is.nan(dpareto(0, 1, -1))))

  expect_warning(expect_true(is.nan(dprop(1, -10, 0.5))))
  expect_warning(expect_true(is.nan(dprop(1, 10, -1))))
  expect_warning(expect_true(is.nan(dprop(1, 10, 2))))
  
  expect_warning(expect_true(is.nan(drayleigh(0, -1))))
  
  expect_warning(expect_true(is.nan(dsgomp(1, -0.4, 1))))
  expect_warning(expect_true(is.nan(dsgomp(1, 0.4, -1))))
  
  expect_warning(expect_true(is.nan(dskellam(1, -1, 1))))
  expect_warning(expect_true(is.nan(dskellam(1, 1, -1))))
  
  expect_warning(expect_true(is.nan(dslash(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(dtnorm(1, 0, -1, -2, 2))))
  expect_warning(expect_true(is.nan(dtnorm(1, 0, 1, 2, -2))))
  expect_warning(expect_true(is.nan(dtnorm(1, 0, 1, 0, 0)))) 
  
  expect_warning(expect_true(is.nan(dtpois(1, lambda = -5, a = 0))))
  expect_warning(expect_true(is.nan(dtpois(1, lambda = -5, a = 6))))
  expect_warning(expect_true(is.nan(dtpois(1, lambda = -5, a = 6, b = 5))))

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
  
  expect_warning(expect_true(is.nan(pbetapr(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(pbetapr(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(pbetapr(1, 1, 1, -1))))

  expect_warning(expect_true(is.nan(pbern(1, -1))))
  expect_warning(expect_true(is.nan(pbern(1, 2))))
  
  expect_warning(expect_true(is.nan(pbhatt(1, sigma = -1))))
  expect_warning(expect_true(is.nan(pbhatt(1, a = -1))))
  
  expect_warning(expect_true(is.nan(pbnbinom(1, -1, 1, 1))))
  expect_warning(expect_true(is.nan(pbnbinom(1, 1, -1, 1))))
  expect_warning(expect_true(is.nan(pbnbinom(1, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(pcat(1, c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(pcat(1, c(0.5, -1)))))

  expect_warning(expect_true(is.nan(pdlaplace(1, 0, scale = -1))))
  expect_warning(expect_true(is.nan(pdlaplace(1, 0, scale = 2))))
  
  expect_warning(expect_true(is.nan(pdnorm(1, sd = -1))))
  
  expect_warning(expect_true(is.nan(pdgamma(1, -9, 1))))
  expect_warning(expect_true(is.nan(pdgamma(1, 9, -1))))
  
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
  
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(-1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,-1,1/3)))))
  expect_warning(expect_true(is.nan(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,-1)))))
  
  expect_warning(expect_true(is.nan(pmixpois(0, c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(-1,1/3,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(1/3,-1,1/3)))))
  expect_warning(expect_true(is.nan(pmixpois(0, c(1,2,3), c(1/3,1/3,-1)))))
  
  expect_warning(expect_true(is.nan(pnhyper(1, 60.5, 35, 15))))
  expect_warning(expect_true(is.nan(pnhyper(1, 60, 35.5, 15))))
  expect_warning(expect_true(is.nan(pnhyper(1, 60, 35, 15.5))))
  expect_warning(expect_true(is.nan(pnhyper(1, -60, 35, 15))))
  expect_warning(expect_true(is.nan(pnhyper(1, 60, -35, 15))))
  expect_warning(expect_true(is.nan(pnhyper(1, 60, 35, -15))))
  expect_warning(expect_true(is.nan(pnhyper(1, 60, 35, 40))))
  
  expect_warning(expect_true(is.nan(pnsbeta(0.5, -1, 1, -2, 2))))
  expect_warning(expect_true(is.nan(pnsbeta(0.5, 1, -1, -2, 2))))
  expect_warning(expect_true(is.nan(pnsbeta(0.5, 1, 1, 2, -2))))
  
  expect_warning(expect_true(is.nan(plst(1, -2, 0, 1))))
  expect_warning(expect_true(is.nan(plst(1, 2, 0, -1))))
  
  expect_warning(expect_true(is.nan(ppareto(0, -1, 1))))
  expect_warning(expect_true(is.nan(ppareto(0, 1, -1))))

  expect_warning(expect_true(is.nan(pprop(1, -10, 0.5))))
  expect_warning(expect_true(is.nan(pprop(1, 10, -1))))
  expect_warning(expect_true(is.nan(pprop(1, 10, 2))))
  
  expect_warning(expect_true(is.nan(prayleigh(0, -1))))
  
  expect_warning(expect_true(is.nan(psgomp(1, -0.4, 1))))
  expect_warning(expect_true(is.nan(psgomp(1, 0.4, -1))))

  expect_warning(expect_true(is.nan(pslash(1, sigma = -1))))
  
  expect_warning(expect_true(is.nan(ptnorm(1, 0, -1, -2, 2))))
  expect_warning(expect_true(is.nan(ptnorm(1, 0, 1, 2, -2))))
  expect_warning(expect_true(is.nan(ptnorm(1, 0, 1, 0, 0))))
  
  expect_warning(expect_true(is.nan(ptpois(1, lambda = -5, a = 0))))
  expect_warning(expect_true(is.nan(ptpois(1, lambda = -5, a = 6))))
  expect_warning(expect_true(is.nan(ptpois(1, lambda = -5, a = 6, b = 5))))

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




test_that("Wrong parameter values in quantile functions", {

  expect_warning(expect_true(is.nan(qbetapr(0.5, -1, 1, 1))))
  expect_warning(expect_true(is.nan(qbetapr(0.5, 1, -1, 1))))
  expect_warning(expect_true(is.nan(qbetapr(0.5, 1, 1, -1))))
  
  expect_warning(expect_true(is.nan(qbern(0.5, -1))))
  expect_warning(expect_true(is.nan(qbern(0.5, 2))))

  expect_warning(expect_true(is.nan(qcat(0.5, c(-1, 0.5)))))
  expect_warning(expect_true(is.nan(qcat(0.5, c(0.5, -1)))))
  
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
  
  expect_warning(expect_true(is.nan(qnhyper(0.5, 60.5, 35, 15))))
  expect_warning(expect_true(is.nan(qnhyper(0.5, 60, 35.5, 15))))
  expect_warning(expect_true(is.nan(qnhyper(0.5, 60, 35, 15.5))))
  expect_warning(expect_true(is.nan(qnhyper(0.5, -60, 35, 15))))
  expect_warning(expect_true(is.nan(qnhyper(0.5, 60, -35, 15))))
  expect_warning(expect_true(is.nan(qnhyper(0.5, 60, 35, -15))))
  expect_warning(expect_true(is.nan(qnhyper(0.5, 60, 35, 40))))

  expect_warning(expect_true(is.nan(qnsbeta(0.5, -1, 1, -2, 2))))
  expect_warning(expect_true(is.nan(qnsbeta(0.5, 1, -1, -2, 2))))
  expect_warning(expect_true(is.nan(qnsbeta(0.5, 1, 1, 2, -2))))
  
  expect_warning(expect_true(is.nan(qlst(0.5, -2, 0, 1))))
  expect_warning(expect_true(is.nan(qlst(0.5, 2, 0, -1))))
  
  expect_warning(expect_true(is.nan(qpareto(0.5, -1, 1))))
  expect_warning(expect_true(is.nan(qpareto(0.5, 1, -1))))

  expect_warning(expect_true(is.nan(qprop(0.5, -10, 0.5))))
  expect_warning(expect_true(is.nan(qprop(0.5, 10, -1))))
  expect_warning(expect_true(is.nan(qprop(0.5, 10, 2))))
  
  expect_warning(expect_true(is.nan(qrayleigh(0, -1))))

  expect_warning(expect_true(is.nan(qtnorm(0.5, 0, -1, -2, 2))))
  expect_warning(expect_true(is.nan(qtnorm(0.5, 0, 1, 2, -2))))
  expect_warning(expect_true(is.nan(qtnorm(0.5, 0, 1, 0, 0))))
  
  expect_warning(expect_true(is.nan(qtpois(0.5, lambda = -5, a = 0))))
  expect_warning(expect_true(is.nan(qtpois(0.5, lambda = -5, a = 6))))
  expect_warning(expect_true(is.nan(qtpois(0.5, lambda = -5, a = 6, b = 5))))

  expect_warning(expect_true(is.nan(qtriang(0.5, 0, 0, 0))))
  expect_warning(expect_true(is.nan(qtriang(0.5, 1, -1, 0))))
  expect_warning(expect_true(is.nan(qtriang(0.5, -1, 1, 2))))
  expect_warning(expect_true(is.nan(qtriang(0.5, -1, 1, -2))))

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
  
  expect_warning(expect_true(is.na(rbbinom(1, -1, 1, 1))))
  expect_warning(expect_true(is.na(rbbinom(1, 1, -1, 1))))
  expect_warning(expect_true(is.na(rbbinom(1, 1, 1, -1))))
  
  expect_warning(expect_true(is.na(rbetapr(1, -1, 1, 1))))
  expect_warning(expect_true(is.na(rbetapr(1, 1, -1, 1))))
  expect_warning(expect_true(is.na(rbetapr(1, 1, 1, -1))))
  
  expect_warning(expect_true(is.na(rbern(1, -1))))
  expect_warning(expect_true(is.na(rbern(1, 2))))
  
  expect_warning(expect_true(is.na(rbhatt(1, sigma = -1))))
  expect_warning(expect_true(is.na(rbhatt(1, a = -1))))
  
  expect_warning(expect_true(all(is.na(rbnbinom(1, -1, 1, 1)))))
  expect_warning(expect_true(all(is.na(rbnbinom(1, 1, -1, 1)))))
  expect_warning(expect_true(all(is.na(rbnbinom(1, 1, 1, -1)))))
  
  expect_warning(expect_true(all(is.na(rbvnorm(1, sd1 = -1)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, sd2 = -1)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, cor = -2)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, cor = 2)))))
  
  expect_warning(expect_true(all(is.na(rbvpois(1, -1, 1, 1)))))
  expect_warning(expect_true(all(is.na(rbvpois(1, 1, -1, 1)))))
  expect_warning(expect_true(all(is.na(rbvpois(1, 1, 1, -1)))))
  
  expect_warning(expect_true(is.na(rcat(1, c(-1, 0.5)))))
  expect_warning(expect_true(is.na(rcat(1, c(0.5, -1)))))
  expect_warning(expect_true(is.na(rcat(2, matrix(c(-1, 0.5, 0.5, 0.5), byrow = T, ncol = 2))[1])))
  expect_warning(expect_true(is.na(rcat(2, matrix(c(0.5, -1, 0.5, 0.5), byrow = T, ncol = 2))[1])))
  expect_warning(expect_true(is.na(rcat(2, matrix(c(0.5, 0.5, -1, 0.5), byrow = T, ncol = 2))[2])))
  expect_warning(expect_true(is.na(rcat(2, matrix(c(0.5, 0.5, 0.5, -1), byrow = T, ncol = 2))[2])))
  
  expect_warning(expect_true(all(is.na(rdirichlet(1, c(-1, 0.5))))))
  expect_warning(expect_true(all(is.na(rdirichlet(1, c(0.5, -1))))))
  
  expect_warning(expect_true(is.na(rdnorm(1, sd = -1))))
  
  expect_warning(expect_true(is.nan(rdgamma(1, -9, 1))))
  expect_warning(expect_true(is.nan(rdgamma(1, 9, -1))))
  
  expect_warning(expect_true(is.na(rdunif(1, min = 10, max = 1))))
  expect_warning(expect_true(is.na(rdunif(1, min = 0, max = Inf))))
  expect_warning(expect_true(is.na(rdunif(1, min = -Inf, max = Inf))))
  expect_warning(expect_true(is.na(rdunif(1, min = Inf, max = -Inf))))
  
  expect_warning(expect_true(is.na(rdweibull(1, -1, 1)))) 
  expect_warning(expect_true(is.na(rdweibull(1, 2, 1)))) 
  expect_warning(expect_true(is.na(rdweibull(1, 0.5, -1))))
  
  expect_warning(expect_true(is.na(rfatigue(1, -1, 1))))
  expect_warning(expect_true(is.na(rfatigue(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rfrechet(1, lambda = -1))))
  expect_warning(expect_true(is.na(rfrechet(1, sigma = -1))))
  
  expect_warning(expect_true(is.na(rgev(1, 1, -1, 1))))
  
  expect_warning(expect_true(is.na(rgompertz(1, -1, 1))))
  expect_warning(expect_true(is.na(rgompertz(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rgpd(1, 1, -1, 1))))
  
  expect_warning(expect_true(is.na(rgpois(1, -1, 1))))
  expect_warning(expect_true(is.na(rgpois(1, 1, -1))))
  expect_warning(expect_true(is.na(rgpois(1, 1, scale = 0))))
  
  expect_warning(expect_true(is.na(rgumbel(1, sigma = -1))))
  
  expect_warning(expect_true(is.na(rhcauchy(1, -1))))
  
  expect_warning(expect_true(is.na(rhnorm(1, -1))))
  
  expect_warning(expect_true(is.na(rht(1, 5, -1))))
  
  expect_warning(expect_true(is.na(rhuber(1, 0, -1, 1))))
  expect_warning(expect_true(is.na(rhuber(1, 0, 1, -1))))
  
  expect_warning(expect_true(is.na(rinvgamma(1, -1, 1))))
  expect_warning(expect_true(is.na(rinvgamma(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rinvchisq(1, -1, 1))))
  expect_warning(expect_true(is.na(rinvchisq(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rkumar(1, -1, 1))))
  expect_warning(expect_true(is.na(rkumar(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rlaplace(1, 0, -1))))
  
  expect_warning(expect_true(is.na(rlgser(1, -1))))
  expect_warning(expect_true(is.na(rlgser(1, 2))))
  
  expect_warning(expect_true(is.na(rlomax(1, -1, 1))))
  expect_warning(expect_true(is.na(rlomax(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(-1,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,-1,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,1/3,-1)))))
  
  expect_warning(expect_true(is.na(rmixpois(1, c(-1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,-2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,-3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,3), c(-1,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,3), c(1/3,-1,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,3), c(1/3,1/3,-1)))))
  
  expect_warning(expect_true(is.na(rnhyper(1, 60.5, 35, 15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, 35.5, 15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, 35, 15.5))))
  expect_warning(expect_true(is.na(rnhyper(1, -60, 35, 15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, -35, 15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, 35, -15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, 35, 40))))
  
  expect_warning(expect_true(all(is.na(rdirmnom(1, 1.5, c(1, 1, 1))))))
  expect_warning(expect_true(all(is.na(rdirmnom(1, -3, c(1, 1, 1))))))
  expect_warning(expect_true(all(is.na(rdirmnom(1, 3, c(-1, 1, 1))))))
  expect_warning(expect_true(all(is.na(rdirmnom(1, 3, c(1, -1, 1))))))
  expect_warning(expect_true(all(is.na(rdirmnom(1, 3, c(1, 1, -1))))))
  
  expect_warning(expect_true(all(is.na(rmnom(1, 3, c(-1, 1/3, 1/3))))))
  expect_warning(expect_true(all(is.na(rmnom(1, 3, c(1/3, -1, 1/3))))))
  expect_warning(expect_true(all(is.na(rmnom(1, 3, c(1/3, 1/3, -1))))))
  
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(2,3,4), 99)))))
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(-2,3,4), 99)))))
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(2,-3,4), 99)))))
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(2,3,-4), 99)))))
  
  expect_warning(expect_true(is.na(rnsbeta(1, -1, 1, -2, 2))))
  expect_warning(expect_true(is.na(rnsbeta(1, 1, -1, -2, 2))))
  expect_warning(expect_true(is.na(rnsbeta(1, 1, 1, 2, -2))))
  
  expect_warning(expect_true(is.na(rlst(1, -2, 0, 1))))
  expect_warning(expect_true(is.na(rlst(1, 2, 0, -1))))
  
  expect_warning(expect_true(is.na(rpareto(1, -1, 1))))
  expect_warning(expect_true(is.na(rpareto(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rprop(1, -10, 0.5))))
  expect_warning(expect_true(is.na(rprop(1, 10, -1))))
  expect_warning(expect_true(is.na(rprop(1, 10, 2))))
  
  expect_warning(expect_true(is.na(rrayleigh(1, -1))))
  
  expect_warning(expect_true(is.na(rsgomp(1, -0.4, 1))))
  expect_warning(expect_true(is.na(rsgomp(1, 0.4, -1))))
  
  expect_warning(expect_true(is.na(rskellam(1, -1, 1))))
  expect_warning(expect_true(is.na(rskellam(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rslash(1, sigma = -1))))
  
  expect_warning(expect_true(is.na(rtnorm(1, 0, -1, -2, 2))))
  expect_warning(expect_true(is.na(rtnorm(1, 0, 1, 2, -2))))
  expect_warning(expect_true(is.na(rtnorm(1, 0, 1, 0, 0))))
  
  expect_warning(expect_true(is.na(rtpois(1, lambda = -5, a = 0))))
  expect_warning(expect_true(is.na(rtpois(1, lambda = -5, a = 6))))
  expect_warning(expect_true(is.na(rtpois(1, lambda = -5, a = 6, b = 5))))

  expect_warning(expect_true(is.na(rtriang(1, 0, 0, 0))))
  expect_warning(expect_true(is.na(rtriang(1, 1, -1, 0))))
  expect_warning(expect_true(is.na(rtriang(1, -1, 1, 2))))
  expect_warning(expect_true(is.na(rtriang(1, -1, 1, -2))))
  
  expect_warning(expect_true(is.na(rwald(1, 1, -1))))
  
  expect_warning(expect_true(is.na(rzip(1, -1, 0.5))))
  expect_warning(expect_true(is.na(rzip(1, 1, -1))))
  expect_warning(expect_true(is.na(rzip(1, 1, 2))))
  
  expect_warning(expect_true(is.na(rzib(1, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.na(rzib(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.na(rzib(1, 1, 0.5, 2))))
  expect_warning(expect_true(is.na(rzib(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.na(rzib(1, 1, 0.5, 2))))
  
  expect_warning(expect_true(is.na(rzinb(1, -1, 0.5, 0.5))))
  expect_warning(expect_true(is.na(rzinb(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.na(rzinb(1, 1, 0.5, 2))))
  expect_warning(expect_true(is.na(rzinb(1, 1, -1, 0.5))))
  expect_warning(expect_true(is.na(rzinb(1, 1, 0.5, 2))))
  
})