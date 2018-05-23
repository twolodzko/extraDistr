

test_that("Missing values in PDF and PMF functions", {
  
  expect_true(is.na(dbbinom(NA, 1, 1, 1)))
  expect_true(is.na(dbbinom(1, NA, 1, 1)))
  expect_true(is.na(dbbinom(1, 1, NA, 1)))
  expect_true(is.na(dbbinom(1, 1, 1, NA)))
  
  expect_true(is.na(dbern(NA, 0.5)))
  expect_true(is.na(dbern(1, NA)))
  
  expect_true(is.na(dbetapr(NA, 1, 1, 1)))
  expect_true(is.na(dbetapr(1, NA, 1, 1)))
  expect_true(is.na(dbetapr(1, 1, NA, 1)))
  expect_true(is.na(dbetapr(1, 1, 1, NA)))
  
  expect_true(is.na(dbhatt(NA, 1, 1, 1)))
  expect_true(is.na(dbhatt(1, NA, 1, 1)))
  expect_true(is.na(dbhatt(1, 1, NA, 1)))
  expect_true(is.na(dbhatt(1, 1, 1, NA)))
  
  expect_true(is.na(dbnbinom(NA, 1, 1, 1)))
  expect_true(is.na(dbnbinom(1, NA, 1, 1)))
  expect_true(is.na(dbnbinom(1, 1, NA, 1)))
  expect_true(is.na(dbnbinom(1, 1, 1, NA)))
  
  expect_true(is.na(dbvnorm(NA, 1, 1, 1, 1, 1, 0.5)))
  expect_true(is.na(dbvnorm(1, NA, 1, 1, 1, 1, 0.5)))
  expect_true(is.na(dbvnorm(1, 1, NA, 1, 1, 1, 0.5)))
  expect_true(is.na(dbvnorm(1, 1, 1, NA, 1, 1, 0.5)))
  expect_true(is.na(dbvnorm(1, 1, 1, 1, NA, 1, 0.5)))
  expect_true(is.na(dbvnorm(1, 1, 1, 1, 1, NA, 0.5)))
  expect_true(is.na(dbvnorm(1, 1, 1, 1, 1, 1, NA)))
  
  expect_true(is.na(dbvpois(NA, 1, 1, 1, 1)))
  expect_true(is.na(dbvpois(1, NA, 1, 1, 1)))
  expect_true(is.na(dbvpois(1, 1, NA, 1, 1)))
  expect_true(is.na(dbvpois(1, 1, 1, NA, 1)))
  expect_true(is.na(dbvpois(1, 1, 1, 1, NA)))
  
  expect_true(is.na(dcat(NA, c(0.5, 0.5))))
  expect_true(is.na(dcat(1, c(NA, 0.5))))
  expect_true(is.na(dcat(1, c(0.5, NA))))
  
  expect_true(is.na(ddirichlet(c(NA, 0.5), c(0.5, 0.5))))
  expect_true(is.na(ddirichlet(c(0.5, NA), c(0.5, 0.5))))
  expect_true(is.na(ddirichlet(c(0.5, 0.5), c(NA, 0.5))))
  expect_true(is.na(ddirichlet(c(0.5, 0.5), c(0.5, NA))))
  
  expect_true(is.na(ddlaplace(NA, 1, 1)))
  expect_true(is.na(ddlaplace(1, NA, 1)))
  expect_true(is.na(ddlaplace(1, 1, NA)))
  
  expect_true(is.na(ddnorm(NA, 1, 1)))
  expect_true(is.na(ddnorm(1, NA, 1)))
  expect_true(is.na(ddnorm(1, 1, NA)))
  
  expect_true(is.na(ddgamma(NA, 9, 1)))
  expect_true(is.na(ddgamma(1, NA, 1)))
  expect_true(is.na(ddgamma(1, 9, NA)))
  
  expect_true(is.na(ddunif(NA, 1, 10)))
  expect_true(is.na(ddunif(1, NA, 10)))
  expect_true(is.na(ddunif(1, 1, NA)))
  
  expect_true(is.na(ddweibull(NA, 1, 1)))
  expect_true(is.na(ddweibull(1, NA, 1))) 
  expect_true(is.na(ddweibull(1, 1, NA)))
  
  expect_true(is.na(dfatigue(NA, 1, 1)))
  expect_true(is.na(dfatigue(1, NA, 1)))
  expect_true(is.na(dfatigue(1, 1, NA)))
  
  expect_true(is.na(dfrechet(NA, 1, 1, 1)))
  expect_true(is.na(dfrechet(1, NA, 1, 1)))
  expect_true(is.na(dfrechet(1, 1, NA, 1)))
  expect_true(is.na(dfrechet(1, 1, 1, NA)))

  expect_true(is.na(dgev(NA, 1, 1, 1)))
  expect_true(is.na(dgev(1, NA, 1, 1)))
  expect_true(is.na(dgev(1, 1, NA, 1)))
  expect_true(is.na(dgev(1, 1, 1, NA)))
  
  expect_true(is.na(dgompertz(NA, 1, 1)))
  expect_true(is.na(dgompertz(1, NA, 1)))
  expect_true(is.na(dgompertz(1, 1, NA)))
  
  expect_true(is.na(dgpd(NA, 1, 1, 1)))
  expect_true(is.na(dgpd(1, NA, 1, 1)))
  expect_true(is.na(dgpd(1, 1, NA, 1)))
  expect_true(is.na(dgpd(1, 1, 1, NA)))
  
  expect_true(is.na(dgpois(NA, 1, 1)))
  expect_true(is.na(dgpois(1, NA, 1)))
  expect_true(is.na(dgpois(1, 1, NA)))
  
  expect_true(is.na(dgumbel(NA, 1, 1)))
  expect_true(is.na(dgumbel(1, NA, 1)))
  expect_true(is.na(dgumbel(1, 1, NA)))
  
  expect_true(is.na(dhcauchy(NA, 1)))
  expect_true(is.na(dhcauchy(1, NA)))
  
  expect_true(is.na(dhnorm(NA, 1)))
  expect_true(is.na(dhnorm(1, NA)))
  
  expect_true(is.na(dht(NA, 5, 1)))
  expect_true(is.na(dht(1, NA, 1)))
  expect_true(is.na(dht(1, 5, NA)))
  
  expect_true(is.na(dhuber(NA, 0, 1, 1)))
  expect_true(is.na(dhuber(1, NA, 1, 1)))
  expect_true(is.na(dhuber(1, 0, NA, 1)))
  expect_true(is.na(dhuber(1, 0, 1, NA)))
  
  expect_true(is.na(dinvgamma(NA, 1, 1)))
  expect_true(is.na(dinvgamma(1, NA, 1)))
  expect_true(is.na(dinvgamma(1, 1, NA)))
  
  expect_true(is.na(dinvchisq(NA, 1, 1)))
  expect_true(is.na(dinvchisq(1, NA, 1)))
  expect_true(is.na(dinvchisq(1, 1, NA)))
  
  expect_true(is.na(dkumar(NA, 1, 1)))
  expect_true(is.na(dkumar(0.5, NA, 1)))
  expect_true(is.na(dkumar(0.5, 1, NA)))
  
  expect_true(is.na(dlaplace(NA, 0, 1)))
  expect_true(is.na(dlaplace(1, NA, 1)))
  expect_true(is.na(dlaplace(1, 0, NA)))
  
  expect_true(is.na(dlgser(NA, 0.5)))
  expect_true(is.na(dlgser(1, NA)))

  expect_true(is.na(dlomax(NA, 1, 1)))
  expect_true(is.na(dlomax(1, NA, 1)))
  expect_true(is.na(dlomax(1, 1, NA)))

  expect_true(is.na(dmixnorm(NA, c(1,2,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixnorm(0, c(NA,2,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixnorm(0, c(1,NA,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixnorm(0, c(1,2,NA), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixnorm(0, c(1,2,3), c(NA,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,NA,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,2,3), c(NA,1/3,1/3))))
  expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,NA,1/3))))
  expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,NA))))

  expect_true(is.na(dmixpois(NA, c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixpois(0, c(NA,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixpois(0, c(1,NA,3), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixpois(0, c(1,2,NA), c(1/3,1/3,1/3))))
  expect_true(is.na(dmixpois(0, c(1,2,3), c(NA,1/3,1/3))))
  expect_true(is.na(dmixpois(0, c(1,2,3), c(1/3,NA,1/3))))
  expect_true(is.na(dmixpois(0, c(1,2,3), c(1/3,1/3,NA))))

  expect_true(is.na(dnhyper(NA, 60, 35, 15)))
  expect_true(is.na(dnhyper(1, NA, 35, 15)))
  expect_true(is.na(dnhyper(1, 60, NA, 15)))
  expect_true(is.na(dnhyper(1, 60, 35, NA)))
  
  expect_true(is.na(ddirmnom(c(NA, 1, 1), 3, c(1, 1, 1))))
  expect_true(is.na(ddirmnom(c(1, NA, 1), 3, c(1, 1, 1))))
  expect_true(is.na(ddirmnom(c(1, 1, NA), 3, c(1, 1, 1))))
  expect_true(is.na(ddirmnom(c(1, 1, 1), NA, c(1, 1, 1))))
  expect_true(is.na(ddirmnom(c(1, 1, 1), 3, c(NA, 1, 1))))
  expect_true(is.na(ddirmnom(c(1, 1, 1), 3, c(1, NA, 1))))
  expect_true(is.na(ddirmnom(c(1, 1, 1), 3, c(1, 1, NA))))

  expect_true(is.na(dmnom(c(NA, 1, 1), 3, c(1/3, 1/3, 1/3))))
  expect_true(is.na(dmnom(c(1, NA, 1), 3, c(1/3, 1/3, 1/3))))
  expect_true(is.na(dmnom(c(1, 1, NA), 3, c(1/3, 1/3, 1/3))))
  expect_true(is.na(dmnom(c(1, 1, 1), NA, c(1/3, 1/3, 1/3))))
  expect_true(is.na(dmnom(c(1, 1, 1), 3, c(NA, 1/3, 1/3))))
  expect_true(is.na(dmnom(c(1, 1, 1), 3, c(1/3, NA, 1/3))))
  expect_true(is.na(dmnom(c(1, 1, 1), 3, c(1/3, 1/3, NA))))

  expect_true(is.na(dmvhyper(c(NA, 2, 2), c(2,3,4), 5)))
  expect_true(is.na(dmvhyper(c(1, NA, 2), c(2,3,4), 5)))
  expect_true(is.na(dmvhyper(c(1, 2, NA), c(2,3,4), 5)))
  expect_true(is.na(dmvhyper(c(1, 2, 2), c(NA,3,4), 5)))
  expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,NA,4), 5)))
  expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,3,NA), 5)))
  expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,3,4), NA)))

  expect_true(is.na(dnsbeta(NA, 1, 1, -2, 2)))
  expect_true(is.na(dnsbeta(0.5, NA, 1, -2, 2)))
  expect_true(is.na(dnsbeta(0.5, 1, NA, -2, 2)))
  expect_true(is.na(dnsbeta(0.5, 1, 1, NA, 2)))
  expect_true(is.na(dnsbeta(0.5, 1, 1, -2, NA)))

  expect_true(is.na(dlst(NA, 2, 0, 1)))
  expect_true(is.na(dlst(1, NA, 0, 1)))
  expect_true(is.na(dlst(1, 2, NA, 1)))
  expect_true(is.na(dlst(1, 2, 0, NA)))

  expect_true(is.na(dpareto(NA, 1, 1)))
  expect_true(is.na(dpareto(1, NA, 1)))
  expect_true(is.na(dpareto(1, 1, NA)))
  
  expect_true(is.na(dpower(NA, 1, 1)))
  expect_true(is.na(dpower(1, NA, 1)))
  expect_true(is.na(dpower(1, 1, NA)))
  
  expect_true(is.na(dprop(NA, 10, 0.5)))
  expect_true(is.na(dprop(1, NA, 0.5)))
  expect_true(is.na(dprop(1, 10, NA)))

  expect_true(is.na(drayleigh(NA, 1)))
  expect_true(is.na(drayleigh(1, NA)))

  expect_true(is.na(dskellam(NA, 1, 1)))
  expect_true(is.na(dskellam(1, NA, 1)))
  expect_true(is.na(dskellam(1, 1, NA)))
  
  expect_true(is.na(dsgomp(NA, 0.4, 1)))
  expect_true(is.na(dsgomp(1, NA, 1)))
  expect_true(is.na(dsgomp(1, 0.4, NA)))

  expect_true(is.na(dslash(NA, 1, 1)))
  expect_true(is.na(dslash(1, NA, 1)))
  expect_true(is.na(dslash(1, 1, NA)))

  expect_true(is.na(dtnorm(NA, 0, 1, -2, 2)))
  expect_true(is.na(dtnorm(1, NA, 1, -2, 2)))
  expect_true(is.na(dtnorm(1, 0, NA, -2, 2)))
  expect_true(is.na(dtnorm(1, 0, 1, NA, 2)))
  expect_true(is.na(dtnorm(1, 0, 1, -2, NA)))

  expect_true(is.na(dtpois(NA, 5, 0)))
  expect_true(is.na(dtpois(1, NA, 0)))
  expect_true(is.na(dtpois(1, 5, NA)))

  expect_true(is.na(dtriang(NA, 0, 1, 0.5)))
  expect_true(is.na(dtriang(0.5, NA, 1, 0.5)))
  expect_true(is.na(dtriang(0.5, 0, NA, 0.5)))
  expect_true(is.na(dtriang(0.5, 0, 1, NA)))

  expect_true(is.na(dwald(NA, 1, 1)))
  expect_true(is.na(dwald(1, NA, 1)))
  expect_true(is.na(dwald(1, 1, NA)))

  expect_true(is.na(dzip(NA, 1, 0.5)))
  expect_true(is.na(dzip(1, NA, 0.5)))
  expect_true(is.na(dzip(1, 1, NA)))

  expect_true(is.na(dzib(NA, 1, 0.5, 0.5)))
  expect_true(is.na(dzib(1, NA, 0.5, 0.5)))
  expect_true(is.na(dzib(1, 1, NA, 0.5)))
  expect_true(is.na(dzib(1, 1, 0.5, NA)))

  expect_true(is.na(dzinb(NA, 1, 0.5, 0.5)))
  expect_true(is.na(dzinb(1, NA, 0.5, 0.5)))
  expect_true(is.na(dzinb(1, 1, NA, 0.5)))
  expect_true(is.na(dzinb(1, 1, 0.5, NA)))

})





test_that("Wrong parameter values in CDF functions", {
  
  expect_true(is.na(pbbinom(NA, 1, 1, 1)))
  expect_true(is.na(pbbinom(1, NA, 1, 1)))
  expect_true(is.na(pbbinom(1, 1, NA, 1)))
  expect_true(is.na(pbbinom(1, 1, 1, NA)))
  
  expect_true(is.na(pbern(NA, 0.5)))
  expect_true(is.na(pbern(1, NA)))
  
  expect_true(is.na(pbetapr(NA, 1, 1, 1)))
  expect_true(is.na(pbetapr(1, NA, 1, 1)))
  expect_true(is.na(pbetapr(1, 1, NA, 1)))
  expect_true(is.na(pbetapr(1, 1, 1, NA)))
  
  expect_true(is.na(pbhatt(NA, 1, 1, 1)))
  expect_true(is.na(pbhatt(1, NA, 1, 1)))
  expect_true(is.na(pbhatt(1, 1, NA, 1)))
  expect_true(is.na(pbhatt(1, 1, 1, NA)))
  
  expect_true(is.na(pbnbinom(NA, 1, 1, 1)))
  expect_true(is.na(pbnbinom(1, NA, 1, 1)))
  expect_true(is.na(pbnbinom(1, 1, NA, 1)))
  expect_true(is.na(pbnbinom(1, 1, 1, NA)))
  
  expect_true(is.na(pcat(NA, c(0.5, 0.5))))
  expect_true(is.na(pcat(1, c(NA, 0.5))))
  expect_true(is.na(pcat(1, c(0.5, NA))))

  expect_true(is.na(pdlaplace(NA, 1, 1)))
  expect_true(is.na(pdlaplace(1, NA, 1)))
  expect_true(is.na(pdlaplace(1, 1, NA)))
  
  expect_true(is.na(pdnorm(NA, 1, 1)))
  expect_true(is.na(pdnorm(1, NA, 1)))
  expect_true(is.na(pdnorm(1, 1, NA)))
  
  expect_true(is.na(pdgamma(NA, 9, 1)))
  expect_true(is.na(pdgamma(1, NA, 1)))
  expect_true(is.na(pdgamma(1, 9, NA)))
  
  expect_true(is.na(pdunif(NA, 1, 10)))
  expect_true(is.na(pdunif(1, NA, 10)))
  expect_true(is.na(pdunif(1, 1, NA)))
  
  expect_true(is.na(pdweibull(NA, 1, 1)))
  expect_true(is.na(pdweibull(1, NA, 1))) 
  expect_true(is.na(pdweibull(1, 1, NA)))
  
  expect_true(is.na(pfatigue(NA, 1, 1)))
  expect_true(is.na(pfatigue(1, NA, 1)))
  expect_true(is.na(pfatigue(1, 1, NA)))
  
  expect_true(is.na(pfrechet(NA, 1, 1, 1)))
  expect_true(is.na(pfrechet(1, NA, 1, 1)))
  expect_true(is.na(pfrechet(1, 1, NA, 1)))
  expect_true(is.na(pfrechet(1, 1, 1, NA)))
  
  expect_true(is.na(pgev(NA, 1, 1, 1)))
  expect_true(is.na(pgev(1, NA, 1, 1)))
  expect_true(is.na(pgev(1, 1, NA, 1)))
  expect_true(is.na(pgev(1, 1, 1, NA)))
  
  expect_true(is.na(pgompertz(NA, 1, 1)))
  expect_true(is.na(pgompertz(1, NA, 1)))
  expect_true(is.na(pgompertz(1, 1, NA)))
  
  expect_true(is.na(pgpd(NA, 1, 1, 1)))
  expect_true(is.na(pgpd(1, NA, 1, 1)))
  expect_true(is.na(pgpd(1, 1, NA, 1)))
  expect_true(is.na(pgpd(1, 1, 1, NA)))
  
  expect_true(is.na(pgpois(NA, 1, 1)))
  expect_true(is.na(pgpois(1, NA, 1)))
  expect_true(is.na(pgpois(1, 1, NA)))
  
  expect_true(is.na(pgumbel(NA, 1, 1)))
  expect_true(is.na(pgumbel(1, NA, 1)))
  expect_true(is.na(pgumbel(1, 1, NA)))
  
  expect_true(is.na(phcauchy(NA, 1)))
  expect_true(is.na(phcauchy(1, NA)))
  
  expect_true(is.na(phnorm(NA, 1)))
  expect_true(is.na(phnorm(1, NA)))
  
  expect_true(is.na(pht(NA, 5, 1)))
  expect_true(is.na(pht(1, NA, 1)))
  expect_true(is.na(pht(1, 5, NA)))
  
  expect_true(is.na(phuber(NA, 0, 1, 1)))
  expect_true(is.na(phuber(1, NA, 1, 1)))
  expect_true(is.na(phuber(1, 0, NA, 1)))
  expect_true(is.na(phuber(1, 0, 1, NA)))
  
  expect_true(is.na(pinvgamma(NA, 1, 1)))
  expect_true(is.na(pinvgamma(1, NA, 1)))
  expect_true(is.na(pinvgamma(1, 1, NA)))
  
  expect_true(is.na(pinvchisq(NA, 1, 1)))
  expect_true(is.na(pinvchisq(1, NA, 1)))
  expect_true(is.na(pinvchisq(1, 1, NA)))
  
  expect_true(is.na(pkumar(NA, 1, 1)))
  expect_true(is.na(pkumar(0.5, NA, 1)))
  expect_true(is.na(pkumar(0.5, 1, NA)))
  
  expect_true(is.na(plaplace(NA, 0, 1)))
  expect_true(is.na(plaplace(1, NA, 1)))
  expect_true(is.na(plaplace(1, 0, NA)))
  
  expect_true(is.na(plgser(NA, 0.5)))
  expect_true(is.na(plgser(1, NA)))
  
  expect_true(is.na(plomax(NA, 1, 1)))
  expect_true(is.na(plomax(1, NA, 1)))
  expect_true(is.na(plomax(1, 1, NA)))
  
  expect_true(is.na(pmixnorm(NA, c(1,2,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixnorm(0, c(NA,2,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixnorm(0, c(1,NA,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixnorm(0, c(1,2,NA), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixnorm(0, c(1,2,3), c(NA,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,NA,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,2,3), c(NA,1/3,1/3))))
  expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,NA,1/3))))
  expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,NA))))
  
  expect_true(is.na(pmixpois(NA, c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixpois(0, c(NA,2,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixpois(0, c(1,NA,3), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixpois(0, c(1,2,NA), c(1/3,1/3,1/3))))
  expect_true(is.na(pmixpois(0, c(1,2,3), c(NA,1/3,1/3))))
  expect_true(is.na(pmixpois(0, c(1,2,3), c(1/3,NA,1/3))))
  expect_true(is.na(pmixpois(0, c(1,2,3), c(1/3,1/3,NA))))

  expect_true(is.na(pnhyper(NA, 60, 35, 15)))
  expect_true(is.na(pnhyper(1, NA, 35, 15)))
  expect_true(is.na(pnhyper(1, 60, NA, 15)))
  expect_true(is.na(pnhyper(1, 60, 35, NA)))
  
  expect_true(is.na(pnsbeta(NA, 1, 1, -2, 2)))
  expect_true(is.na(pnsbeta(0.5, NA, 1, -2, 2)))
  expect_true(is.na(pnsbeta(0.5, 1, NA, -2, 2)))
  expect_true(is.na(pnsbeta(0.5, 1, 1, NA, 2)))
  expect_true(is.na(pnsbeta(0.5, 1, 1, -2, NA)))
  
  expect_true(is.na(plst(NA, 2, 0, 1)))
  expect_true(is.na(plst(1, NA, 0, 1)))
  expect_true(is.na(plst(1, 2, NA, 1)))
  expect_true(is.na(plst(1, 2, 0, NA)))
  
  expect_true(is.na(ppareto(NA, 1, 1)))
  expect_true(is.na(ppareto(1, NA, 1)))
  expect_true(is.na(ppareto(1, 1, NA)))
  
  expect_true(is.na(ppower(NA, 1, 1)))
  expect_true(is.na(ppower(1, NA, 1)))
  expect_true(is.na(ppower(1, 1, NA)))
  
  expect_true(is.na(pprop(NA, 10, 0.5)))
  expect_true(is.na(pprop(1, NA, 0.5)))
  expect_true(is.na(pprop(1, 10, NA)))
  
  expect_true(is.na(prayleigh(NA, 1)))
  expect_true(is.na(prayleigh(1, NA)))
  
  expect_true(is.na(psgomp(NA, 0.4, 1)))
  expect_true(is.na(psgomp(1, NA, 1)))
  expect_true(is.na(psgomp(1, 0.4, NA)))

  expect_true(is.na(pslash(NA, 1, 1)))
  expect_true(is.na(pslash(1, NA, 1)))
  expect_true(is.na(pslash(1, 1, NA)))
  
  expect_true(is.na(ptnorm(NA, 0, 1, -2, 2)))
  expect_true(is.na(ptnorm(1, NA, 1, -2, 2)))
  expect_true(is.na(ptnorm(1, 0, NA, -2, 2)))
  expect_true(is.na(ptnorm(1, 0, 1, NA, 2)))
  expect_true(is.na(ptnorm(1, 0, 1, -2, NA)))
  
  expect_true(is.na(ptpois(NA, 5, 0)))
  expect_true(is.na(ptpois(1, NA, 0)))
  expect_true(is.na(ptpois(1, 5, NA)))
  
  expect_true(is.na(ptriang(NA, 0, 1, 0.5)))
  expect_true(is.na(ptriang(0.5, NA, 1, 0.5)))
  expect_true(is.na(ptriang(0.5, 0, NA, 0.5)))
  expect_true(is.na(ptriang(0.5, 0, 1, NA)))
  
  expect_true(is.na(pwald(NA, 1, 1)))
  expect_true(is.na(pwald(1, NA, 1)))
  expect_true(is.na(pwald(1, 1, NA)))
  
  expect_true(is.na(pzip(NA, 1, 0.5)))
  expect_true(is.na(pzip(1, NA, 0.5)))
  expect_true(is.na(pzip(1, 1, NA)))
  
  expect_true(is.na(pzib(NA, 1, 0.5, 0.5)))
  expect_true(is.na(pzib(1, NA, 0.5, 0.5)))
  expect_true(is.na(pzib(1, 1, NA, 0.5)))
  expect_true(is.na(pzib(1, 1, 0.5, NA)))
  
  expect_true(is.na(pzinb(NA, 1, 0.5, 0.5)))
  expect_true(is.na(pzinb(1, NA, 0.5, 0.5)))
  expect_true(is.na(pzinb(1, 1, NA, 0.5)))
  expect_true(is.na(pzinb(1, 1, 0.5, NA)))
  
})




test_that("Wrong parameter values in inverse CDF functions", {

  expect_true(is.na(qbern(NA, 0.5)))
  expect_true(is.na(qbern(0.5, NA)))
  
  expect_true(is.na(qbetapr(NA, 1, 1, 1)))
  expect_true(is.na(qbetapr(0.5, NA, 1, 1)))
  expect_true(is.na(qbetapr(0.5, 1, NA, 1)))
  expect_true(is.na(qbetapr(0.5, 1, 1, NA)))
  
  expect_true(is.na(qcat(NA, c(0.5, 0.5))))
  expect_true(is.na(qcat(0.5, c(NA, 0.5))))
  expect_true(is.na(qcat(0.5, c(0.5, NA))))

  expect_true(is.na(qdunif(NA, 1, 10)))
  expect_true(is.na(qdunif(0.5, NA, 10)))
  expect_true(is.na(qdunif(0.5, 1, NA)))
  
  expect_true(is.na(qdweibull(NA, 1, 1)))
  expect_true(is.na(qdweibull(0.5, NA, 1))) 
  expect_true(is.na(qdweibull(0.5, 1, NA)))
  
  expect_true(is.na(qfatigue(NA, 1, 1)))
  expect_true(is.na(qfatigue(0.5, NA, 1)))
  expect_true(is.na(qfatigue(0.5, 1, NA)))
  
  expect_true(is.na(qfrechet(NA, 1, 1, 1)))
  expect_true(is.na(qfrechet(0.5, NA, 1, 1)))
  expect_true(is.na(qfrechet(0.5, 1, NA, 1)))
  expect_true(is.na(qfrechet(0.5, 1, 1, NA)))
  
  expect_true(is.na(qgev(NA, 1, 1, 1)))
  expect_true(is.na(qgev(0.5, NA, 1, 1)))
  expect_true(is.na(qgev(0.5, 1, NA, 1)))
  expect_true(is.na(qgev(0.5, 1, 1, NA)))
  
  expect_true(is.na(qgompertz(NA, 1, 1)))
  expect_true(is.na(qgompertz(0.5, NA, 1)))
  expect_true(is.na(qgompertz(0.5, 1, NA)))
  
  expect_true(is.na(qgpd(NA, 1, 1, 1)))
  expect_true(is.na(qgpd(0.5, NA, 1, 1)))
  expect_true(is.na(qgpd(0.5, 1, NA, 1)))
  expect_true(is.na(qgpd(0.5, 1, 1, NA)))
  
  expect_true(is.na(qgumbel(NA, 1, 1)))
  expect_true(is.na(qgumbel(0.5, NA, 1)))
  expect_true(is.na(qgumbel(0.5, 1, NA)))
  
  expect_true(is.na(qhcauchy(NA, 1)))
  expect_true(is.na(qhcauchy(0.5, NA)))
  
  expect_true(is.na(qhnorm(NA, 1)))
  expect_true(is.na(qhnorm(0.5, NA)))
  
  expect_true(is.na(qht(NA, 5, 1)))
  expect_true(is.na(qht(0.5, NA, 1)))
  expect_true(is.na(qht(0.5, 5, NA)))
  
  expect_true(is.na(qhuber(NA, 0, 1, 1)))
  expect_true(is.na(qhuber(0.5, NA, 1, 1)))
  expect_true(is.na(qhuber(0.5, 0, NA, 1)))
  expect_true(is.na(qhuber(0.5, 0, 1, NA)))
  
  expect_true(is.na(qinvgamma(NA, 1, 1)))
  expect_true(is.na(qinvgamma(0.5, NA, 1)))
  expect_true(is.na(qinvgamma(0.5, 1, NA)))
  
  expect_true(is.na(qinvchisq(NA, 1, 1)))
  expect_true(is.na(qinvchisq(0.5, NA, 1)))
  expect_true(is.na(qinvchisq(0.5, 1, NA)))
  
  expect_true(is.na(qkumar(NA, 1, 1)))
  expect_true(is.na(qkumar(0.5, NA, 1)))
  expect_true(is.na(qkumar(0.5, 1, NA)))
  
  expect_true(is.na(qlaplace(NA, 0, 1)))
  expect_true(is.na(qlaplace(0.5, NA, 1)))
  expect_true(is.na(qlaplace(0.5, 0, NA)))
  
  expect_true(is.na(qlgser(NA, 0.5)))
  expect_true(is.na(qlgser(0.5, NA)))
  
  expect_true(is.na(qlomax(NA, 1, 1)))
  expect_true(is.na(qlomax(0.5, NA, 1)))
  expect_true(is.na(qlomax(0.5, 1, NA)))

  expect_true(is.na(qnhyper(NA, 60, 35, 15)))
  expect_true(is.na(qnhyper(0.5, NA, 35, 15)))
  expect_true(is.na(qnhyper(0.5, 60, NA, 15)))
  expect_true(is.na(qnhyper(0.5, 60, 35, NA)))
  
  expect_true(is.na(qnsbeta(NA, 1, 1, -2, 2)))
  expect_true(is.na(qnsbeta(0.5, NA, 1, -2, 2)))
  expect_true(is.na(qnsbeta(0.5, 1, NA, -2, 2)))
  expect_true(is.na(qnsbeta(0.5, 1, 1, NA, 2)))
  expect_true(is.na(qnsbeta(0.5, 1, 1, -2, NA)))
  
  expect_true(is.na(qlst(NA, 2, 0, 1)))
  expect_true(is.na(qlst(0.5, NA, 0, 1)))
  expect_true(is.na(qlst(0.5, 2, NA, 1)))
  expect_true(is.na(qlst(0.5, 2, 0, NA)))
  
  expect_true(is.na(qpareto(NA, 1, 1)))
  expect_true(is.na(qpareto(0.5, NA, 1)))
  expect_true(is.na(qpareto(0.5, 1, NA)))
  
  expect_true(is.na(dpower(NA, 1, 1)))
  expect_true(is.na(dpower(0.5, NA, 1)))
  expect_true(is.na(dpower(0.5, 1, NA)))
  
  expect_true(is.na(qprop(NA, 10, 0.5)))
  expect_true(is.na(qprop(0.5, NA, 0.5)))
  expect_true(is.na(qprop(0.5, 10, NA)))
  
  expect_true(is.na(qrayleigh(NA, 1)))
  expect_true(is.na(qrayleigh(0.5, NA)))
  
  expect_true(is.na(qtlambda(NA, 0.5)))
  expect_true(is.na(qtlambda(0, NA)))

  expect_true(is.na(qtnorm(NA, 0, 1, -2, 2)))
  expect_true(is.na(qtnorm(0.5, NA, 1, -2, 2)))
  expect_true(is.na(qtnorm(0.5, 0, NA, -2, 2)))
  expect_true(is.na(qtnorm(0.5, 0, 1, NA, 2)))
  expect_true(is.na(qtnorm(0.5, 0, 1, -2, NA)))
  
  expect_true(is.na(qtpois(NA, 5, 0)))
  expect_true(is.na(qtpois(0.5, NA, 0)))
  expect_true(is.na(qtpois(0.5, 5, NA)))
  
  expect_true(is.na(qtriang(NA, 0, 1, 0.5)))
  expect_true(is.na(qtriang(0.5, NA, 1, 0.5)))
  expect_true(is.na(qtriang(0.5, 0, NA, 0.5)))
  expect_true(is.na(qtriang(0.5, 0, 1, NA)))
  
  expect_true(is.na(qzip(NA, 1, 0.5)))
  expect_true(is.na(qzip(0.5, NA, 0.5)))
  expect_true(is.na(qzip(0.5, 1, NA)))
  
  expect_true(is.na(qzib(NA, 1, 0.5, 0.5)))
  expect_true(is.na(qzib(0.5, NA, 0.5, 0.5)))
  expect_true(is.na(qzib(0.5, 1, NA, 0.5)))
  expect_true(is.na(qzib(0.5, 1, 0.5, NA)))
  
  expect_true(is.na(qzinb(NA, 1, 0.5, 0.5)))
  expect_true(is.na(qzinb(0.5, NA, 0.5, 0.5)))
  expect_true(is.na(qzinb(0.5, 1, NA, 0.5)))
  expect_true(is.na(qzinb(0.5, 1, 0.5, NA)))

})

 
 

test_that("Wrong parameter values in RNG functions", {

  expect_warning(expect_true(is.na(rbbinom(1, NA, 1, 1))))
  expect_warning(expect_true(is.na(rbbinom(1, 1, NA, 1))))
  expect_warning(expect_true(is.na(rbbinom(1, 1, 1, NA))))

  expect_warning(expect_true(is.na(rbern(1, NA))))
  
  expect_warning(expect_true(is.na(rbetapr(1, NA, 1, 1))))
  expect_warning(expect_true(is.na(rbetapr(1, 1, NA, 1))))
  expect_warning(expect_true(is.na(rbetapr(1, 1, 1, NA))))

  expect_warning(expect_true(is.na(rbhatt(1, NA, 1, 1))))
  expect_warning(expect_true(is.na(rbhatt(1, 1, NA, 1))))
  expect_warning(expect_true(is.na(rbhatt(1, 1, 1, NA))))

  expect_warning(expect_true(is.na(rbnbinom(1, NA, 1, 1))))
  expect_warning(expect_true(is.na(rbnbinom(1, 1, NA, 1))))
  expect_warning(expect_true(is.na(rbnbinom(1, 1, 1, NA))))

  expect_warning(expect_true(all(is.na(rbvnorm(1, NA, 1, 1, 1, 0.5)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, 1, NA, 1, 1, 0.5)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, 1, 1, NA, 1, 0.5)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, 1, 1, 1, NA, 0.5)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, 1, 1, 1, 1, NA)))))
  
  expect_warning(expect_true(all(is.na(rbvpois(1, NA, 1, 1)))))
  expect_warning(expect_true(all(is.na(rbvpois(1, 1, NA, 1)))))
  expect_warning(expect_true(all(is.na(rbvpois(1, 1, 1, NA)))))

  expect_warning(expect_true(is.na(rcat(1, c(NA, 0.5)))))
  expect_warning(expect_true(is.na(rcat(1, c(0.5, NA)))))
  
  expect_warning(expect_true(is.na(rcatlp(1, c(NA, 0.5)))))
  expect_warning(expect_true(is.na(rcatlp(1, c(0.5, NA)))))

  expect_warning(expect_true(all(is.na(rdirichlet(1, c(NA, 0.5))))))
  expect_warning(expect_true(all(is.na(rdirichlet(1, c(0.5, NA))))))
  
  expect_warning(expect_true(is.na(rdlaplace(1, NA, 1))))
  expect_warning(expect_true(is.na(rdlaplace(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rdnorm(1, NA, 1))))
  expect_warning(expect_true(is.na(rdnorm(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rdgamma(1, NA, 1))))
  expect_warning(expect_true(is.na(rdgamma(1, 9, NA))))
  
  expect_warning(expect_true(is.na(rdunif(1, NA, 10))))
  expect_warning(expect_true(is.na(rdunif(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rdweibull(1, NA, 1)))) 
  expect_warning(expect_true(is.na(rdweibull(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rfatigue(1, NA, 1))))
  expect_warning(expect_true(is.na(rfatigue(1, 1, NA))))

  expect_warning(expect_true(is.na(rfrechet(1, NA, 1, 1))))
  expect_warning(expect_true(is.na(rfrechet(1, 1, NA, 1))))
  expect_warning(expect_true(is.na(rfrechet(1, 1, 1, NA))))
  
  expect_warning(expect_true(is.na(rgev(1, NA, 1, 1))))
  expect_warning(expect_true(is.na(rgev(1, 1, NA, 1))))
  expect_warning(expect_true(is.na(rgev(1, 1, 1, NA))))

  expect_warning(expect_true(is.na(rgompertz(1, NA, 1))))
  expect_warning(expect_true(is.na(rgompertz(1, 1, NA))))

  expect_warning(expect_true(is.na(rgpd(1, NA, 1, 1))))
  expect_warning(expect_true(is.na(rgpd(1, 1, NA, 1))))
  expect_warning(expect_true(is.na(rgpd(1, 1, 1, NA))))

  expect_warning(expect_true(is.na(rgpois(1, NA, 1))))
  expect_warning(expect_true(is.na(rgpois(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rgumbel(1, NA, 1))))
  expect_warning(expect_true(is.na(rgumbel(1, 1, NA))))

  expect_warning(expect_true(is.na(rhcauchy(1, NA))))

  expect_warning(expect_true(is.na(rhnorm(1, NA))))
  
  expect_warning(expect_true(is.na(rht(1, NA, 1))))
  expect_warning(expect_true(is.na(rht(1, 5, NA))))
  
  expect_warning(expect_true(is.na(rhuber(1, NA, 1, 1))))
  expect_warning(expect_true(is.na(rhuber(1, 0, NA, 1))))
  expect_warning(expect_true(is.na(rhuber(1, 0, 1, NA))))
  
  # expect_warning(expect_true(is.na(rinvgamma(1, NA, 1))))
  # expect_warning(expect_true(is.na(rinvgamma(1, 1, NA))))
  # 
  # expect_warning(expect_true(is.na(rinvchisq(1, NA, 1))))
  # expect_warning(expect_true(is.na(rinvchisq(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rkumar(1, NA, 1))))
  expect_warning(expect_true(is.na(rkumar(1, 1, NA))))

  expect_warning(expect_true(is.na(rlaplace(1, NA, 1))))
  expect_warning(expect_true(is.na(rlaplace(1, 0, NA))))
  
  expect_warning(expect_true(is.na(rlgser(1, NA))))
  
  expect_warning(expect_true(is.na(rlomax(1, NA, 1))))
  expect_warning(expect_true(is.na(rlomax(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rmixnorm(1, c(NA,2,3), c(1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,NA,3), c(1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,NA), c(1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(NA,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,NA,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(NA,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,NA,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,1/3,NA)))))

  expect_warning(expect_true(is.na(rmixpois(1, c(NA,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,NA,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,NA), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,3), c(NA,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,3), c(1/3,NA,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,3), c(1/3,1/3,NA)))))
  
  expect_warning(expect_true(is.na(rnhyper(1, NA, 35, 15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, NA, 15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, 35, NA))))
  
  expect_warning(expect_true(all(is.na(rdirmnom(1, NA, c(1, 1, 1))))))
  expect_warning(expect_true(all(is.na(rdirmnom(1, 3, c(NA, 1, 1))))))
  expect_warning(expect_true(all(is.na(rdirmnom(1, 3, c(1, NA, 1))))))
  expect_warning(expect_true(all(is.na(rdirmnom(1, 3, c(1, 1, NA))))))
  
  expect_warning(expect_true(all(is.na(rmnom(1, NA, c(1/3, 1/3, 1/3))))))
  expect_warning(expect_true(all(is.na(rmnom(1, 3, c(NA, 1/3, 1/3))))))
  expect_warning(expect_true(all(is.na(rmnom(1, 3, c(1/3, NA, 1/3))))))
  expect_warning(expect_true(all(is.na(rmnom(1, 3, c(1/3, 1/3, NA))))))
  
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(NA,3,4), 5)))))
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(2,NA,4), 5)))))
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(2,3,NA), 5)))))
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(2,3,4), NA)))))
  
  expect_warning(expect_true(is.na(rnsbeta(1, NA, 1, -2, 2))))
  expect_warning(expect_true(is.na(rnsbeta(1, 1, NA, -2, 2))))
  expect_warning(expect_true(is.na(rnsbeta(1, 1, 1, NA, 2))))
  expect_warning(expect_true(is.na(rnsbeta(1, 1, 1, -2, NA))))

  expect_warning(expect_true(is.na(rlst(1, NA, 0, 1))))
  expect_warning(expect_true(is.na(rlst(1, 2, NA, 1))))
  expect_warning(expect_true(is.na(rlst(1, 2, 0, NA))))
  
  expect_warning(expect_true(is.na(rpareto(1, NA, 1))))
  expect_warning(expect_true(is.na(rpareto(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rpower(1, NA, 1))))
  expect_warning(expect_true(is.na(rpower(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rprop(1, NA, 0.5))))
  expect_warning(expect_true(is.na(rprop(1, 10, NA))))
  
  expect_warning(expect_true(is.na(rrayleigh(1, NA))))
  
  expect_warning(expect_true(is.na(rtlambda(1, NA))))
  
  expect_warning(expect_true(is.na(rsgomp(1, NA, 1))))
  expect_warning(expect_true(is.na(rsgomp(1, 0.4, NA))))
  
  expect_warning(expect_true(is.na(rskellam(1, NA, 1))))
  expect_warning(expect_true(is.na(rskellam(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rslash(1, NA, 1))))
  expect_warning(expect_true(is.na(rslash(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rtnorm(1, NA, 1, -2, 2))))
  expect_warning(expect_true(is.na(rtnorm(1, 0, NA, -2, 2))))
  expect_warning(expect_true(is.na(rtnorm(1, 0, 1, NA, 2))))
  expect_warning(expect_true(is.na(rtnorm(1, 0, 1, -2, NA))))

  expect_warning(expect_true(is.na(rtpois(1, NA, 0))))
  expect_warning(expect_true(is.na(rtpois(1, 5, NA))))
  
  expect_warning(expect_true(is.na(rtriang(1, NA, 1, 0.5))))
  expect_warning(expect_true(is.na(rtriang(1, 0, NA, 0.5))))
  expect_warning(expect_true(is.na(rtriang(1, 0, 1, NA))))
  
  expect_warning(expect_true(is.na(rwald(1, NA, 1))))
  expect_warning(expect_true(is.na(rwald(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rzip(1, NA, 0.5))))
  expect_warning(expect_true(is.na(rzip(1, 1, NA))))
  
  expect_warning(expect_true(is.na(rzib(1, NA, 0.5, 0.5))))
  expect_warning(expect_true(is.na(rzib(1, 1, NA, 0.5))))
  expect_warning(expect_true(is.na(rzib(1, 1, 0.5, NA))))
  
  expect_warning(expect_true(is.na(rzinb(1, NA, 0.5, 0.5))))
  expect_warning(expect_true(is.na(rzinb(1, 1, NA, 0.5))))
  expect_warning(expect_true(is.na(rzinb(1, 1, 0.5, NA))))

})