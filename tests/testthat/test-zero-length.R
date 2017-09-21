

is_zero_length <- function(x) length(x) == 0L

test_that("Zero-length in PDF and PMF functions", {
  
  expect_true(is_zero_length(dbbinom(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(dbbinom(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dbbinom(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(dbbinom(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dbern(numeric(0), 0.5)))
  expect_true(is_zero_length(dbern(1, numeric(0))))
  
  expect_true(is_zero_length(dbetapr(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(dbetapr(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dbetapr(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(dbetapr(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dbhatt(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(dbhatt(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dbhatt(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(dbhatt(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dbnbinom(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(dbnbinom(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dbnbinom(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(dbnbinom(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dbvnorm(numeric(0), 1, 1, 1, 1, 1, 0.5)))
  expect_true(is_zero_length(dbvnorm(1, numeric(0), 1, 1, 1, 1, 0.5)))
  expect_true(is_zero_length(dbvnorm(1, 1, numeric(0), 1, 1, 1, 0.5)))
  expect_true(is_zero_length(dbvnorm(1, 1, 1, numeric(0), 1, 1, 0.5)))
  expect_true(is_zero_length(dbvnorm(1, 1, 1, 1, numeric(0), 1, 0.5)))
  expect_true(is_zero_length(dbvnorm(1, 1, 1, 1, 1, numeric(0), 0.5)))
  expect_true(is_zero_length(dbvnorm(1, 1, 1, 1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dbvpois(numeric(0), 1, 1, 1, 1)))
  expect_true(is_zero_length(dbvpois(1, numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(dbvpois(1, 1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dbvpois(1, 1, 1, numeric(0), 1)))
  expect_true(is_zero_length(dbvpois(1, 1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dcat(numeric(0), c(0.5, 0.5))))
  expect_true(is_zero_length(dcat(1, numeric(0))))
  expect_true(is_zero_length(dcat(1, matrix(1, 0, 0))))
  
  expect_true(is_zero_length(ddirichlet(numeric(0), c(0.5, 0.5))))
  expect_true(is_zero_length(ddirichlet(c(0.5, 0.5), numeric(0))))
  
  expect_true(is_zero_length(ddlaplace(numeric(0), 1, 1)))
  expect_true(is_zero_length(ddlaplace(1, numeric(0), 1)))
  expect_true(is_zero_length(ddlaplace(1, 1, numeric(0))))
  
  expect_true(is_zero_length(ddnorm(numeric(0), 1, 1)))
  expect_true(is_zero_length(ddnorm(1, numeric(0), 1)))
  expect_true(is_zero_length(ddnorm(1, 1, numeric(0))))
  
  expect_true(is_zero_length(ddgamma(numeric(0), 9, 1)))
  expect_true(is_zero_length(ddgamma(1, numeric(0), 1)))
  expect_true(is_zero_length(ddgamma(1, 9, numeric(0))))
  
  expect_true(is_zero_length(ddunif(numeric(0), 1, 10)))
  expect_true(is_zero_length(ddunif(1, numeric(0), 10)))
  expect_true(is_zero_length(ddunif(1, 1, numeric(0))))
  
  expect_true(is_zero_length(ddweibull(numeric(0), 1, 1)))
  expect_true(is_zero_length(ddweibull(1, numeric(0), 1))) 
  expect_true(is_zero_length(ddweibull(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dfatigue(numeric(0), 1, 1)))
  expect_true(is_zero_length(dfatigue(1, numeric(0), 1)))
  expect_true(is_zero_length(dfatigue(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dfrechet(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(dfrechet(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dfrechet(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(dfrechet(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dgev(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(dgev(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dgev(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(dgev(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dgompertz(numeric(0), 1, 1)))
  expect_true(is_zero_length(dgompertz(1, numeric(0), 1)))
  expect_true(is_zero_length(dgompertz(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dgpd(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(dgpd(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dgpd(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(dgpd(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(dgpois(numeric(0), 1, 1)))
  expect_true(is_zero_length(dgpois(1, numeric(0), 1)))
  expect_true(is_zero_length(dgpois(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dgumbel(numeric(0), 1, 1)))
  expect_true(is_zero_length(dgumbel(1, numeric(0), 1)))
  expect_true(is_zero_length(dgumbel(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dhcauchy(numeric(0), 1)))
  expect_true(is_zero_length(dhcauchy(1, numeric(0))))
  
  expect_true(is_zero_length(dhnorm(numeric(0), 1)))
  expect_true(is_zero_length(dhnorm(1, numeric(0))))
  
  expect_true(is_zero_length(dht(numeric(0), 5, 1)))
  expect_true(is_zero_length(dht(1, numeric(0), 1)))
  expect_true(is_zero_length(dht(1, 5, numeric(0))))
  
  expect_true(is_zero_length(dhuber(numeric(0), 0, 1, 1)))
  expect_true(is_zero_length(dhuber(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(dhuber(1, 0, numeric(0), 1)))
  expect_true(is_zero_length(dhuber(1, 0, 1, numeric(0))))
  
  expect_true(is_zero_length(dinvgamma(numeric(0), 1, 1)))
  expect_true(is_zero_length(dinvgamma(1, numeric(0), 1)))
  expect_true(is_zero_length(dinvgamma(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dinvchisq(numeric(0), 1, 1)))
  expect_true(is_zero_length(dinvchisq(1, numeric(0), 1)))
  expect_true(is_zero_length(dinvchisq(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dkumar(numeric(0), 1, 1)))
  expect_true(is_zero_length(dkumar(0.5, numeric(0), 1)))
  expect_true(is_zero_length(dkumar(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(dlaplace(numeric(0), 0, 1)))
  expect_true(is_zero_length(dlaplace(1, numeric(0), 1)))
  expect_true(is_zero_length(dlaplace(1, 0, numeric(0))))
  
  expect_true(is_zero_length(dlgser(numeric(0), 0.5)))
  expect_true(is_zero_length(dlgser(1, numeric(0))))
  
  expect_true(is_zero_length(dlomax(numeric(0), 1, 1)))
  expect_true(is_zero_length(dlomax(1, numeric(0), 1)))
  expect_true(is_zero_length(dlomax(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dmixnorm(numeric(0), c(1,2,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(dmixnorm(0, numeric(0), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(dmixnorm(0, c(1,2,3), numeric(0), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(dmixnorm(0, c(1,2,3), c(1,2,3), numeric(0))))
  
  expect_true(is_zero_length(dmixpois(numeric(0), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(dmixpois(0, numeric(0), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(dmixpois(0, c(1,2,3), numeric(0))))
  
  expect_true(is_zero_length(dnhyper(numeric(0), 60, 35, 15)))
  expect_true(is_zero_length(dnhyper(1, numeric(0), 35, 15)))
  expect_true(is_zero_length(dnhyper(1, 60, numeric(0), 15)))
  expect_true(is_zero_length(dnhyper(1, 60, 35, numeric(0))))
  
  expect_true(is_zero_length(ddirmnom(numeric(0), 3, c(1, 1, 1))))
  expect_true(is_zero_length(ddirmnom(c(1, 1, 1), numeric(0), c(1, 1, 1))))
  expect_true(is_zero_length(ddirmnom(c(1, 1, 1), 3, numeric(0))))
  
  expect_true(is_zero_length(dmnom(numeric(0), 3, c(1/3, 1/3, 1/3))))
  expect_true(is_zero_length(dmnom(c(1, 1, 1), numeric(0), c(1/3, 1/3, 1/3))))
  expect_true(is_zero_length(dmnom(c(1, 1, 1), 3, numeric(0))))
  
  expect_true(is_zero_length(dmvhyper(numeric(0), c(2,3,4), 5)))
  expect_true(is_zero_length(dmvhyper(c(1, 2, 2), numeric(0), 5)))
  expect_true(is_zero_length(dmvhyper(c(1, 2, 2), c(2,3,4), numeric(0))))
  
  expect_true(is_zero_length(dnsbeta(numeric(0), 1, 1, -2, 2)))
  expect_true(is_zero_length(dnsbeta(0.5, numeric(0), 1, -2, 2)))
  expect_true(is_zero_length(dnsbeta(0.5, 1, numeric(0), -2, 2)))
  expect_true(is_zero_length(dnsbeta(0.5, 1, 1, numeric(0), 2)))
  expect_true(is_zero_length(dnsbeta(0.5, 1, 1, -2, numeric(0))))
  
  expect_true(is_zero_length(dlst(numeric(0), 2, 0, 1)))
  expect_true(is_zero_length(dlst(1, numeric(0), 0, 1)))
  expect_true(is_zero_length(dlst(1, 2, numeric(0), 1)))
  expect_true(is_zero_length(dlst(1, 2, 0, numeric(0))))
  
  expect_true(is_zero_length(dpareto(numeric(0), 1, 1)))
  expect_true(is_zero_length(dpareto(1, numeric(0), 1)))
  expect_true(is_zero_length(dpareto(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dprop(numeric(0), 10, 0.5)))
  expect_true(is_zero_length(dprop(1, numeric(0), 0.5)))
  expect_true(is_zero_length(dprop(1, 10, numeric(0))))
  
  expect_true(is_zero_length(drayleigh(numeric(0), 1)))
  expect_true(is_zero_length(drayleigh(1, numeric(0))))
  
  expect_true(is_zero_length(dskellam(numeric(0), 1, 1)))
  expect_true(is_zero_length(dskellam(1, numeric(0), 1)))
  expect_true(is_zero_length(dskellam(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dsgomp(numeric(0), 0.4, 1)))
  expect_true(is_zero_length(dsgomp(1, numeric(0), 1)))
  expect_true(is_zero_length(dsgomp(1, 0.4, numeric(0))))
  
  expect_true(is_zero_length(dslash(numeric(0), 1, 1)))
  expect_true(is_zero_length(dslash(1, numeric(0), 1)))
  expect_true(is_zero_length(dslash(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dtnorm(numeric(0), 0, 1, -2, 2)))
  expect_true(is_zero_length(dtnorm(1, numeric(0), 1, -2, 2)))
  expect_true(is_zero_length(dtnorm(1, 0, numeric(0), -2, 2)))
  expect_true(is_zero_length(dtnorm(1, 0, 1, numeric(0), 2)))
  expect_true(is_zero_length(dtnorm(1, 0, 1, -2, numeric(0))))
  
  expect_true(is_zero_length(dtpois(numeric(0), 5, 0)))
  expect_true(is_zero_length(dtpois(1, numeric(0), 0)))
  expect_true(is_zero_length(dtpois(1, 5, numeric(0))))
  
  expect_true(is_zero_length(dtriang(numeric(0), 0, 1, 0.5)))
  expect_true(is_zero_length(dtriang(0.5, numeric(0), 1, 0.5)))
  expect_true(is_zero_length(dtriang(0.5, 0, numeric(0), 0.5)))
  expect_true(is_zero_length(dtriang(0.5, 0, 1, numeric(0))))
  
  expect_true(is_zero_length(dwald(numeric(0), 1, 1)))
  expect_true(is_zero_length(dwald(1, numeric(0), 1)))
  expect_true(is_zero_length(dwald(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dzip(numeric(0), 1, 0.5)))
  expect_true(is_zero_length(dzip(1, numeric(0), 0.5)))
  expect_true(is_zero_length(dzip(1, 1, numeric(0))))
  
  expect_true(is_zero_length(dzib(numeric(0), 1, 0.5, 0.5)))
  expect_true(is_zero_length(dzib(1, numeric(0), 0.5, 0.5)))
  expect_true(is_zero_length(dzib(1, 1, numeric(0), 0.5)))
  expect_true(is_zero_length(dzib(1, 1, 0.5, numeric(0))))
  
  expect_true(is_zero_length(dzinb(numeric(0), 1, 0.5, 0.5)))
  expect_true(is_zero_length(dzinb(1, numeric(0), 0.5, 0.5)))
  expect_true(is_zero_length(dzinb(1, 1, numeric(0), 0.5)))
  expect_true(is_zero_length(dzinb(1, 1, 0.5, numeric(0))))
  
})





test_that("Zero-length in CDF functions", {
  
  expect_true(is_zero_length(pbbinom(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(pbbinom(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(pbbinom(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(pbbinom(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(pbern(numeric(0), 0.5)))
  expect_true(is_zero_length(pbern(1, numeric(0))))
  
  expect_true(is_zero_length(pbetapr(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(pbetapr(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(pbetapr(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(pbetapr(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(pbhatt(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(pbhatt(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(pbhatt(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(pbhatt(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(pbnbinom(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(pbnbinom(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(pbnbinom(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(pbnbinom(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(pcat(numeric(0), c(0.5, 0.5))))
  expect_true(is_zero_length(pcat(1, numeric(0))))
  expect_true(is_zero_length(pcat(1, matrix(1, 0, 0))))
  
  expect_true(is_zero_length(pdlaplace(numeric(0), 1, 1)))
  expect_true(is_zero_length(pdlaplace(1, numeric(0), 1)))
  expect_true(is_zero_length(pdlaplace(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pdnorm(numeric(0), 1, 1)))
  expect_true(is_zero_length(pdnorm(1, numeric(0), 1)))
  expect_true(is_zero_length(pdnorm(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pdgamma(numeric(0), 9, 1)))
  expect_true(is_zero_length(pdgamma(1, numeric(0), 1)))
  expect_true(is_zero_length(pdgamma(1, 9, numeric(0))))
  
  expect_true(is_zero_length(pdunif(numeric(0), 1, 10)))
  expect_true(is_zero_length(pdunif(1, numeric(0), 10)))
  expect_true(is_zero_length(pdunif(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pdweibull(numeric(0), 1, 1)))
  expect_true(is_zero_length(pdweibull(1, numeric(0), 1))) 
  expect_true(is_zero_length(pdweibull(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pfatigue(numeric(0), 1, 1)))
  expect_true(is_zero_length(pfatigue(1, numeric(0), 1)))
  expect_true(is_zero_length(pfatigue(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pfrechet(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(pfrechet(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(pfrechet(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(pfrechet(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(pgev(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(pgev(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(pgev(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(pgev(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(pgompertz(numeric(0), 1, 1)))
  expect_true(is_zero_length(pgompertz(1, numeric(0), 1)))
  expect_true(is_zero_length(pgompertz(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pgpd(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(pgpd(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(pgpd(1, 1, numeric(0), 1)))
  expect_true(is_zero_length(pgpd(1, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(pgpois(numeric(0), 1, 1)))
  expect_true(is_zero_length(pgpois(1, numeric(0), 1)))
  expect_true(is_zero_length(pgpois(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pgumbel(numeric(0), 1, 1)))
  expect_true(is_zero_length(pgumbel(1, numeric(0), 1)))
  expect_true(is_zero_length(pgumbel(1, 1, numeric(0))))
  
  expect_true(is_zero_length(phcauchy(numeric(0), 1)))
  expect_true(is_zero_length(phcauchy(1, numeric(0))))
  
  expect_true(is_zero_length(phnorm(numeric(0), 1)))
  expect_true(is_zero_length(phnorm(1, numeric(0))))
  
  expect_true(is_zero_length(pht(numeric(0), 5, 1)))
  expect_true(is_zero_length(pht(1, numeric(0), 1)))
  expect_true(is_zero_length(pht(1, 5, numeric(0))))
  
  expect_true(is_zero_length(phuber(numeric(0), 0, 1, 1)))
  expect_true(is_zero_length(phuber(1, numeric(0), 1, 1)))
  expect_true(is_zero_length(phuber(1, 0, numeric(0), 1)))
  expect_true(is_zero_length(phuber(1, 0, 1, numeric(0))))
  
  expect_true(is_zero_length(pinvgamma(numeric(0), 1, 1)))
  expect_true(is_zero_length(pinvgamma(1, numeric(0), 1)))
  expect_true(is_zero_length(pinvgamma(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pinvchisq(numeric(0), 1, 1)))
  expect_true(is_zero_length(pinvchisq(1, numeric(0), 1)))
  expect_true(is_zero_length(pinvchisq(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pkumar(numeric(0), 1, 1)))
  expect_true(is_zero_length(pkumar(0.5, numeric(0), 1)))
  expect_true(is_zero_length(pkumar(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(plaplace(numeric(0), 0, 1)))
  expect_true(is_zero_length(plaplace(1, numeric(0), 1)))
  expect_true(is_zero_length(plaplace(1, 0, numeric(0))))
  
  expect_true(is_zero_length(plgser(numeric(0), 0.5)))
  expect_true(is_zero_length(plgser(1, numeric(0))))
  
  expect_true(is_zero_length(plomax(numeric(0), 1, 1)))
  expect_true(is_zero_length(plomax(1, numeric(0), 1)))
  expect_true(is_zero_length(plomax(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pmixnorm(numeric(0), c(1,2,3), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(pmixnorm(0, numeric(0), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(pmixnorm(0, c(1,2,3), numeric(0), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(pmixnorm(0, c(1,2,3), c(1,2,3), numeric(0))))
  
  expect_true(is_zero_length(pmixpois(numeric(0), c(1,2,3), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(pmixpois(0, numeric(0), c(1/3,1/3,1/3))))
  expect_true(is_zero_length(pmixpois(0, c(1,2,3), numeric(0))))
  
  expect_true(is_zero_length(pnhyper(numeric(0), 60, 35, 15)))
  expect_true(is_zero_length(pnhyper(1, numeric(0), 35, 15)))
  expect_true(is_zero_length(pnhyper(1, 60, numeric(0), 15)))
  expect_true(is_zero_length(pnhyper(1, 60, 35, numeric(0))))
  
  expect_true(is_zero_length(pnsbeta(numeric(0), 1, 1, -2, 2)))
  expect_true(is_zero_length(pnsbeta(0.5, numeric(0), 1, -2, 2)))
  expect_true(is_zero_length(pnsbeta(0.5, 1, numeric(0), -2, 2)))
  expect_true(is_zero_length(pnsbeta(0.5, 1, 1, numeric(0), 2)))
  expect_true(is_zero_length(pnsbeta(0.5, 1, 1, -2, numeric(0))))
  
  expect_true(is_zero_length(plst(numeric(0), 2, 0, 1)))
  expect_true(is_zero_length(plst(1, numeric(0), 0, 1)))
  expect_true(is_zero_length(plst(1, 2, numeric(0), 1)))
  expect_true(is_zero_length(plst(1, 2, 0, numeric(0))))
  
  expect_true(is_zero_length(ppareto(numeric(0), 1, 1)))
  expect_true(is_zero_length(ppareto(1, numeric(0), 1)))
  expect_true(is_zero_length(ppareto(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pprop(numeric(0), 10, 0.5)))
  expect_true(is_zero_length(pprop(1, numeric(0), 0.5)))
  expect_true(is_zero_length(pprop(1, 10, numeric(0))))
  
  expect_true(is_zero_length(prayleigh(numeric(0), 1)))
  expect_true(is_zero_length(prayleigh(1, numeric(0))))
  
  expect_true(is_zero_length(psgomp(numeric(0), 0.4, 1)))
  expect_true(is_zero_length(psgomp(1, numeric(0), 1)))
  expect_true(is_zero_length(psgomp(1, 0.4, numeric(0))))
  
  expect_true(is_zero_length(pslash(numeric(0), 1, 1)))
  expect_true(is_zero_length(pslash(1, numeric(0), 1)))
  expect_true(is_zero_length(pslash(1, 1, numeric(0))))
  
  expect_true(is_zero_length(ptnorm(numeric(0), 0, 1, -2, 2)))
  expect_true(is_zero_length(ptnorm(1, numeric(0), 1, -2, 2)))
  expect_true(is_zero_length(ptnorm(1, 0, numeric(0), -2, 2)))
  expect_true(is_zero_length(ptnorm(1, 0, 1, numeric(0), 2)))
  expect_true(is_zero_length(ptnorm(1, 0, 1, -2, numeric(0))))
  
  expect_true(is_zero_length(ptpois(numeric(0), 5, 0)))
  expect_true(is_zero_length(ptpois(1, numeric(0), 0)))
  expect_true(is_zero_length(ptpois(1, 5, numeric(0))))
  
  expect_true(is_zero_length(ptriang(numeric(0), 0, 1, 0.5)))
  expect_true(is_zero_length(ptriang(0.5, numeric(0), 1, 0.5)))
  expect_true(is_zero_length(ptriang(0.5, 0, numeric(0), 0.5)))
  expect_true(is_zero_length(ptriang(0.5, 0, 1, numeric(0))))
  
  expect_true(is_zero_length(pwald(numeric(0), 1, 1)))
  expect_true(is_zero_length(pwald(1, numeric(0), 1)))
  expect_true(is_zero_length(pwald(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pzip(numeric(0), 1, 0.5)))
  expect_true(is_zero_length(pzip(1, numeric(0), 0.5)))
  expect_true(is_zero_length(pzip(1, 1, numeric(0))))
  
  expect_true(is_zero_length(pzib(numeric(0), 1, 0.5, 0.5)))
  expect_true(is_zero_length(pzib(1, numeric(0), 0.5, 0.5)))
  expect_true(is_zero_length(pzib(1, 1, numeric(0), 0.5)))
  expect_true(is_zero_length(pzib(1, 1, 0.5, numeric(0))))
  
  expect_true(is_zero_length(pzinb(numeric(0), 1, 0.5, 0.5)))
  expect_true(is_zero_length(pzinb(1, numeric(0), 0.5, 0.5)))
  expect_true(is_zero_length(pzinb(1, 1, numeric(0), 0.5)))
  expect_true(is_zero_length(pzinb(1, 1, 0.5, numeric(0))))
  
})




test_that("Zero-length in inverse CDF functions", {
  
  expect_true(is_zero_length(qbern(numeric(0), 0.5)))
  expect_true(is_zero_length(qbern(0.5, numeric(0))))
  
  expect_true(is_zero_length(qbetapr(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(qbetapr(0.5, numeric(0), 1, 1)))
  expect_true(is_zero_length(qbetapr(0.5, 1, numeric(0), 1)))
  expect_true(is_zero_length(qbetapr(0.5, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(qcat(numeric(0), c(0.5, 0.5))))
  expect_true(is_zero_length(qcat(0.5, numeric(0))))
  expect_true(is_zero_length(qcat(0.5, matrix(1, 0, 0))))
  
  expect_true(is_zero_length(qdunif(numeric(0), 1, 10)))
  expect_true(is_zero_length(qdunif(0.5, numeric(0), 10)))
  expect_true(is_zero_length(qdunif(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qdweibull(numeric(0), 1, 1)))
  expect_true(is_zero_length(qdweibull(0.5, numeric(0), 1))) 
  expect_true(is_zero_length(qdweibull(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qfatigue(numeric(0), 1, 1)))
  expect_true(is_zero_length(qfatigue(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qfatigue(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qfrechet(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(qfrechet(0.5, numeric(0), 1, 1)))
  expect_true(is_zero_length(qfrechet(0.5, 1, numeric(0), 1)))
  expect_true(is_zero_length(qfrechet(0.5, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(qgev(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(qgev(0.5, numeric(0), 1, 1)))
  expect_true(is_zero_length(qgev(0.5, 1, numeric(0), 1)))
  expect_true(is_zero_length(qgev(0.5, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(qgompertz(numeric(0), 1, 1)))
  expect_true(is_zero_length(qgompertz(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qgompertz(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qgpd(numeric(0), 1, 1, 1)))
  expect_true(is_zero_length(qgpd(0.5, numeric(0), 1, 1)))
  expect_true(is_zero_length(qgpd(0.5, 1, numeric(0), 1)))
  expect_true(is_zero_length(qgpd(0.5, 1, 1, numeric(0))))
  
  expect_true(is_zero_length(qgumbel(numeric(0), 1, 1)))
  expect_true(is_zero_length(qgumbel(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qgumbel(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qhcauchy(numeric(0), 1)))
  expect_true(is_zero_length(qhcauchy(0.5, numeric(0))))
  
  expect_true(is_zero_length(qhnorm(numeric(0), 1)))
  expect_true(is_zero_length(qhnorm(0.5, numeric(0))))
  
  expect_true(is_zero_length(qht(numeric(0), 5, 1)))
  expect_true(is_zero_length(qht(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qht(0.5, 5, numeric(0))))
  
  expect_true(is_zero_length(qhuber(numeric(0), 0, 1, 1)))
  expect_true(is_zero_length(qhuber(0.5, numeric(0), 1, 1)))
  expect_true(is_zero_length(qhuber(0.5, 0, numeric(0), 1)))
  expect_true(is_zero_length(qhuber(0.5, 0, 1, numeric(0))))
  
  expect_true(is_zero_length(qinvgamma(numeric(0), 1, 1)))
  expect_true(is_zero_length(qinvgamma(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qinvgamma(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qinvchisq(numeric(0), 1, 1)))
  expect_true(is_zero_length(qinvchisq(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qinvchisq(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qkumar(numeric(0), 1, 1)))
  expect_true(is_zero_length(qkumar(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qkumar(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qlaplace(numeric(0), 0, 1)))
  expect_true(is_zero_length(qlaplace(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qlaplace(0.5, 0, numeric(0))))
  
  expect_true(is_zero_length(qlgser(numeric(0), 0.5)))
  expect_true(is_zero_length(qlgser(0.5, numeric(0))))
  
  expect_true(is_zero_length(qlomax(numeric(0), 1, 1)))
  expect_true(is_zero_length(qlomax(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qlomax(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qnhyper(numeric(0), 60, 35, 15)))
  expect_true(is_zero_length(qnhyper(0.5, numeric(0), 35, 15)))
  expect_true(is_zero_length(qnhyper(0.5, 60, numeric(0), 15)))
  expect_true(is_zero_length(qnhyper(0.5, 60, 35, numeric(0))))
  
  expect_true(is_zero_length(qnsbeta(numeric(0), 1, 1, -2, 2)))
  expect_true(is_zero_length(qnsbeta(0.5, numeric(0), 1, -2, 2)))
  expect_true(is_zero_length(qnsbeta(0.5, 1, numeric(0), -2, 2)))
  expect_true(is_zero_length(qnsbeta(0.5, 1, 1, numeric(0), 2)))
  expect_true(is_zero_length(qnsbeta(0.5, 1, 1, -2, numeric(0))))
  
  expect_true(is_zero_length(qlst(numeric(0), 2, 0, 1)))
  expect_true(is_zero_length(qlst(0.5, numeric(0), 0, 1)))
  expect_true(is_zero_length(qlst(0.5, 2, numeric(0), 1)))
  expect_true(is_zero_length(qlst(0.5, 2, 0, numeric(0))))
  
  expect_true(is_zero_length(qpareto(numeric(0), 1, 1)))
  expect_true(is_zero_length(qpareto(0.5, numeric(0), 1)))
  expect_true(is_zero_length(qpareto(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qprop(numeric(0), 10, 0.5)))
  expect_true(is_zero_length(qprop(0.5, numeric(0), 0.5)))
  expect_true(is_zero_length(qprop(0.5, 10, numeric(0))))
  
  expect_true(is_zero_length(qrayleigh(numeric(0), 1)))
  expect_true(is_zero_length(qrayleigh(0.5, numeric(0))))
  
  expect_true(is_zero_length(qtlambda(numeric(0), 0.5)))
  expect_true(is_zero_length(qtlambda(0, numeric(0))))
  
  expect_true(is_zero_length(qtnorm(numeric(0), 0, 1, -2, 2)))
  expect_true(is_zero_length(qtnorm(0.5, numeric(0), 1, -2, 2)))
  expect_true(is_zero_length(qtnorm(0.5, 0, numeric(0), -2, 2)))
  expect_true(is_zero_length(qtnorm(0.5, 0, 1, numeric(0), 2)))
  expect_true(is_zero_length(qtnorm(0.5, 0, 1, -2, numeric(0))))
  
  expect_true(is_zero_length(qtpois(numeric(0), 5, 0)))
  expect_true(is_zero_length(qtpois(0.5, numeric(0), 0)))
  expect_true(is_zero_length(qtpois(0.5, 5, numeric(0))))
  
  expect_true(is_zero_length(qtriang(numeric(0), 0, 1, 0.5)))
  expect_true(is_zero_length(qtriang(0.5, numeric(0), 1, 0.5)))
  expect_true(is_zero_length(qtriang(0.5, 0, numeric(0), 0.5)))
  expect_true(is_zero_length(qtriang(0.5, 0, 1, numeric(0))))
  
  expect_true(is_zero_length(qzip(numeric(0), 1, 0.5)))
  expect_true(is_zero_length(qzip(0.5, numeric(0), 0.5)))
  expect_true(is_zero_length(qzip(0.5, 1, numeric(0))))
  
  expect_true(is_zero_length(qzib(numeric(0), 1, 0.5, 0.5)))
  expect_true(is_zero_length(qzib(0.5, numeric(0), 0.5, 0.5)))
  expect_true(is_zero_length(qzib(0.5, 1, numeric(0), 0.5)))
  expect_true(is_zero_length(qzib(0.5, 1, 0.5, numeric(0))))
  
  expect_true(is_zero_length(qzinb(numeric(0), 1, 0.5, 0.5)))
  expect_true(is_zero_length(qzinb(0.5, numeric(0), 0.5, 0.5)))
  expect_true(is_zero_length(qzinb(0.5, 1, numeric(0), 0.5)))
  expect_true(is_zero_length(qzinb(0.5, 1, 0.5, numeric(0))))
  
})




test_that("Zero-length in RNG functions", {
  
  expect_warning(expect_true(is.na(rbbinom(1, numeric(0), 1, 1))))
  expect_warning(expect_true(is.na(rbbinom(1, 1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rbbinom(1, 1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rbern(1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rbetapr(1, numeric(0), 1, 1))))
  expect_warning(expect_true(is.na(rbetapr(1, 1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rbetapr(1, 1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rbhatt(1, numeric(0), 1, 1))))
  expect_warning(expect_true(is.na(rbhatt(1, 1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rbhatt(1, 1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rbnbinom(1, numeric(0), 1, 1))))
  expect_warning(expect_true(is.na(rbnbinom(1, 1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rbnbinom(1, 1, 1, numeric(0)))))
  
  expect_warning(expect_true(all(is.na(rbvnorm(1, numeric(0), 1, 1, 1, 0.5)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, 1, numeric(0), 1, 1, 0.5)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, 1, 1, numeric(0), 1, 0.5)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, 1, 1, 1, numeric(0), 0.5)))))
  expect_warning(expect_true(all(is.na(rbvnorm(1, 1, 1, 1, 1, numeric(0))))))
  
  expect_warning(expect_true(all(is.na(rbvpois(1, numeric(0), 1, 1)))))
  expect_warning(expect_true(all(is.na(rbvpois(1, 1, numeric(0), 1)))))
  expect_warning(expect_true(all(is.na(rbvpois(1, 1, 1, numeric(0))))))
  
  expect_warning(expect_true(is.na(rcat(1, numeric(0)))))
  expect_warning(expect_true(is.na(rcat(1, matrix(1, 0, 0)))))
  
  expect_warning(expect_true(is.na(rcatlp(1, numeric(0)))))
  expect_warning(expect_true(is.na(rcatlp(1, matrix(1, 0, 0)))))
  
  expect_warning(expect_true(all(is.na(rdirichlet(1, numeric(0))))))
  
  expect_warning(expect_true(is.na(rdlaplace(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rdlaplace(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rdnorm(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rdnorm(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rdgamma(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rdgamma(1, 9, numeric(0)))))
  
  expect_warning(expect_true(is.na(rdunif(1, numeric(0), 10))))
  expect_warning(expect_true(is.na(rdunif(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rdweibull(1, numeric(0), 1)))) 
  expect_warning(expect_true(is.na(rdweibull(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rfatigue(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rfatigue(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rfrechet(1, numeric(0), 1, 1))))
  expect_warning(expect_true(is.na(rfrechet(1, 1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rfrechet(1, 1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rgev(1, numeric(0), 1, 1))))
  expect_warning(expect_true(is.na(rgev(1, 1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rgev(1, 1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rgompertz(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rgompertz(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rgpd(1, numeric(0), 1, 1))))
  expect_warning(expect_true(is.na(rgpd(1, 1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rgpd(1, 1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rgpois(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rgpois(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rgumbel(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rgumbel(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rhcauchy(1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rhnorm(1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rht(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rht(1, 5, numeric(0)))))
  
  expect_warning(expect_true(is.na(rhuber(1, numeric(0), 1, 1))))
  expect_warning(expect_true(is.na(rhuber(1, 0, numeric(0), 1))))
  expect_warning(expect_true(is.na(rhuber(1, 0, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rinvgamma(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rinvgamma(1, 1, numeric(0)))))

  expect_warning(expect_true(is.na(rinvchisq(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rinvchisq(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rkumar(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rkumar(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rlaplace(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rlaplace(1, 0, numeric(0)))))
  
  expect_warning(expect_true(is.na(rlgser(1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rlomax(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rlomax(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rmixnorm(1, numeric(0), c(1,2,3), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), numeric(0), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), numeric(0)))))
  
  expect_warning(expect_true(is.na(rmixpois(1, numeric(0), c(1/3,1/3,1/3)))))
  expect_warning(expect_true(is.na(rmixpois(1, c(1,2,3), numeric(0)))))
  
  expect_warning(expect_true(is.na(rnhyper(1, numeric(0), 35, 15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, numeric(0), 15))))
  expect_warning(expect_true(is.na(rnhyper(1, 60, 35, numeric(0)))))
  
  expect_warning(expect_true(all(is.na(rdirmnom(1, numeric(0), c(1, 1, 1))))))
  expect_warning(expect_true(all(is.na(rdirmnom(1, 3, numeric(0))))))
  
  expect_warning(expect_true(all(is.na(rmnom(1, numeric(0), c(1/3, 1/3, 1/3))))))
  expect_warning(expect_true(all(is.na(rmnom(1, 3, numeric(0))))))
  
  expect_warning(expect_true(all(is.na(rmvhyper(1, numeric(0), 5)))))
  expect_warning(expect_true(all(is.na(rmvhyper(1, c(2,3,4), numeric(0))))))
  
  expect_warning(expect_true(is.na(rnsbeta(1, numeric(0), 1, -2, 2))))
  expect_warning(expect_true(is.na(rnsbeta(1, 1, numeric(0), -2, 2))))
  expect_warning(expect_true(is.na(rnsbeta(1, 1, 1, numeric(0), 2))))
  expect_warning(expect_true(is.na(rnsbeta(1, 1, 1, -2, numeric(0)))))
  
  expect_warning(expect_true(is.na(rlst(1, numeric(0), 0, 1))))
  expect_warning(expect_true(is.na(rlst(1, 2, numeric(0), 1))))
  expect_warning(expect_true(is.na(rlst(1, 2, 0, numeric(0)))))
  
  expect_warning(expect_true(is.na(rpareto(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rpareto(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rprop(1, numeric(0), 0.5))))
  expect_warning(expect_true(is.na(rprop(1, 10, numeric(0)))))
  
  expect_warning(expect_true(is.na(rrayleigh(1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rtlambda(1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rsgomp(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rsgomp(1, 0.4, numeric(0)))))
  
  expect_warning(expect_true(is.na(rskellam(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rskellam(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rslash(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rslash(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rtnorm(1, numeric(0), 1, -2, 2))))
  expect_warning(expect_true(is.na(rtnorm(1, 0, numeric(0), -2, 2))))
  expect_warning(expect_true(is.na(rtnorm(1, 0, 1, numeric(0), 2))))
  expect_warning(expect_true(is.na(rtnorm(1, 0, 1, -2, numeric(0)))))
  
  expect_warning(expect_true(is.na(rtpois(1, numeric(0), 0))))
  expect_warning(expect_true(is.na(rtpois(1, 5, numeric(0)))))
  
  expect_warning(expect_true(is.na(rtriang(1, numeric(0), 1, 0.5))))
  expect_warning(expect_true(is.na(rtriang(1, 0, numeric(0), 0.5))))
  expect_warning(expect_true(is.na(rtriang(1, 0, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rwald(1, numeric(0), 1))))
  expect_warning(expect_true(is.na(rwald(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rzip(1, numeric(0), 0.5))))
  expect_warning(expect_true(is.na(rzip(1, 1, numeric(0)))))
  
  expect_warning(expect_true(is.na(rzib(1, numeric(0), 0.5, 0.5))))
  expect_warning(expect_true(is.na(rzib(1, 1, numeric(0), 0.5))))
  expect_warning(expect_true(is.na(rzib(1, 1, 0.5, numeric(0)))))
  
  expect_warning(expect_true(is.na(rzinb(1, numeric(0), 0.5, 0.5))))
  expect_warning(expect_true(is.na(rzinb(1, 1, numeric(0), 0.5))))
  expect_warning(expect_true(is.na(rzinb(1, 1, 0.5, numeric(0)))))
  
})