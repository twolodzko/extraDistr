

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
  
  # expect_true(is.na(dlgser(1, -1)))
  # expect_true(is.na(dlgser(1, 2)))
  # 
  # expect_true(is.na(dlomax(1, -1, 1)))
  # expect_true(is.na(dlomax(1, 1, -1)))
  # 
  # expect_true(is.na(dmixnorm(0, c(1,2,3), c(-1,2,3), c(1/3,1/3,1/3))))
  # expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,-2,3), c(1/3,1/3,1/3))))
  # expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,2,-3), c(1/3,1/3,1/3))))
  # expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,2,3), c(-1,1/3,1/3))))
  # expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,-1,1/3))))
  # expect_true(is.na(dmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,-1))))
  # 
  # expect_true(is.na(dmixpois(0, c(-1,2,3), c(1/3,1/3,1/3))))
  # expect_true(is.na(dmixpois(0, c(1,-2,3), c(1/3,1/3,1/3))))
  # expect_true(is.na(dmixpois(0, c(1,2,-3), c(1/3,1/3,1/3))))
  # expect_true(is.na(dmixpois(0, c(1,2,3), c(-1,1/3,1/3))))
  # expect_true(is.na(dmixpois(0, c(1,2,3), c(1/3,-1,1/3))))
  # expect_true(is.na(dmixpois(0, c(1,2,3), c(1/3,1/3,-1))))
  # 
  # expect_true(is.na(ddirmnom(c(1, 1, 1), 1.5, c(1, 1, 1))))
  # expect_true(is.na(ddirmnom(c(1, 1, 1), -3, c(1, 1, 1))))
  # expect_true(is.na(ddirmnom(c(1, 1, 1), 3, c(-1, 1, 1))))
  # expect_true(is.na(ddirmnom(c(1, 1, 1), 3, c(1, -1, 1))))
  # expect_true(is.na(ddirmnom(c(1, 1, 1), 3, c(1, 1, -1))))
  # 
  # expect_true(is.na(dmnom(c(1, 1, 1), 1.5, c(1/3, 1/3, 1/3))))
  # expect_true(is.na(dmnom(c(1, 1, 1), 3, c(-1, 1/3, 1/3))))
  # expect_true(is.na(dmnom(c(1, 1, 1), 3, c(1/3, -1, 1/3))))
  # expect_true(is.na(dmnom(c(1, 1, 1), 3, c(1/3, 1/3, -1))))
  # 
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(2.5,3,4), 5)))
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,3.5,4), 5)))
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,3,4.5), 5)))
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,3,4), 5.5)))
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(-2,3,4), 5)))
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,-3,4), 5)))
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,3,-4), 5)))
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,3,4), -5)))
  # expect_true(is.na(dmvhyper(c(1, 2, 2), c(2,3,4), 85)))
  # 
  # expect_true(is.na(dnsbeta(0.5, -1, 1, -2, 2)))
  # expect_true(is.na(dnsbeta(0.5, 1, -1, -2, 2)))
  # expect_true(is.na(dnsbeta(0.5, 1, 1, 2, -2)))
  # 
  # expect_true(is.na(dnst(1, -2, 0, 1)))
  # expect_true(is.na(dnst(1, 2, 0, -1)))
  # 
  # expect_true(is.na(dpareto(0, -1, 1)))
  # expect_true(is.na(dpareto(0, 1, -1)))
  # 
  # expect_true(is.na(dprop(1, -10, 0.5)))
  # expect_true(is.na(dprop(1, 10, -1)))
  # expect_true(is.na(dprop(1, 10, 2)))
  # 
  # expect_true(is.na(drayleigh(0, -1)))
  # 
  # expect_true(is.na(dskellam(1, -1, 1)))
  # expect_true(is.na(dskellam(1, 1, -1)))
  # 
  # expect_true(is.na(dslash(1, sigma = -1)))
  # 
  # expect_true(is.na(dtnorm(1, 0, -1, -2, 2)))
  # expect_true(is.na(dtnorm(1, 0, 1, 2, -2)))
  # expect_true(is.na(dtnorm(1, 0, 1, 0, 0)))
  # 
  # expect_true(is.na(dtpois(1, lambda = -5, s = 0)))
  # expect_true(is.na(dtpois(1, lambda = -5, s = 6)))
  # expect_true(is.na(dtpois(1, lambda = 5, s = -1)))
  # 
  # expect_true(is.na(dtriang(1, 0, 0, 0)))
  # expect_true(is.na(dtriang(1, 1, -1, 0)))
  # expect_true(is.na(dtriang(1, -1, 1, 2)))
  # expect_true(is.na(dtriang(1, -1, 1, -2)))
  # 
  # expect_true(is.na(dwald(1, 1, -1)))
  # 
  # expect_true(is.na(dzip(1, -1, 0.5)))
  # expect_true(is.na(dzip(1, 1, -1)))
  # expect_true(is.na(dzip(1, 1, 2)))
  # 
  # expect_true(is.na(dzib(1, -1, 0.5, 0.5)))
  # expect_true(is.na(dzib(1, 1, -1, 0.5)))
  # expect_true(is.na(dzib(1, 1, 0.5, 2)))
  # expect_true(is.na(dzib(1, 1, -1, 0.5)))
  # expect_true(is.na(dzib(1, 1, 0.5, 2)))
  # 
  # expect_true(is.na(dzinb(1, -1, 0.5, 0.5)))
  # expect_true(is.na(dzinb(1, 1, -1, 0.5)))
  # expect_true(is.na(dzinb(1, 1, 0.5, 2)))
  # expect_true(is.na(dzinb(1, 1, -1, 0.5)))
  # expect_true(is.na(dzinb(1, 1, 0.5, 2)))
  
})


# 
# 
# 
# test_that("Wrong parameter values in CDF functions", {
# 
#   expect_true(is.na(pbbinom(1, -1, 1, 1)))
#   expect_true(is.na(pbbinom(1, 1, -1, 1)))
#   expect_true(is.na(pbbinom(1, 1, 1, -1)))
# 
#   expect_true(is.na(pbbinom(1, c(-1, 1), c(1, 1), c(1, 1))[1]))
#   expect_true(is.na(pbbinom(1, c(1, 1), c(-1, 1), c(1, 1))[1]))
#   expect_true(is.na(pbbinom(1, c(1, 1), c(1, 1), c(-1, 1))[1]))
# 
#   expect_true(is.na(pbbinom(1, c(1, -1), c(1, 1), c(1, 1))[2]))
#   expect_true(is.na(pbbinom(1, c(1, 1), c(1, -1), c(1, 1))[2]))
#   expect_true(is.na(pbbinom(1, c(1, 1), c(1, 1), c(1, -1))[2]))
# 
#   expect_true(is.na(pbetapr(1, -1, 1, 1)))
#   expect_true(is.na(pbetapr(1, 1, -1, 1)))
#   expect_true(is.na(pbetapr(1, 1, 1, -1)))
# 
#   expect_true(is.na(pbern(1, -1)))
#   expect_true(is.na(pbern(1, 2)))
# 
#   expect_true(is.na(pbhatt(1, sigma = -1)))
#   expect_true(is.na(pbhatt(1, a = -1)))
# 
#   expect_true(is.na(pbnbinom(1, -1, 1, 1)))
#   expect_true(is.na(pbnbinom(1, 1, -1, 1)))
#   expect_true(is.na(pbnbinom(1, 1, 1, -1)))
# 
#   expect_true(is.na(pcat(1, c(-1, 0.5))))
#   expect_true(is.na(pcat(1, c(0.5, -1))))
# 
#   expect_true(is.na(pdlaplace(1, scale = -1)))
#   expect_true(is.na(pdlaplace(1, scale = 2)))
# 
#   expect_true(is.na(pdnorm(1, sd = -1)))
# 
#   expect_true(is.na(pdunif(1, min = 10, max = 1)))
#   expect_true(is.na(pdunif(1, min = 0, max = Inf)))
#   expect_true(is.na(pdunif(1, min = -Inf, max = Inf)))
#   expect_true(is.na(pdunif(1, min = Inf, max = -Inf)))
# 
#   expect_true(is.na(pdweibull(1, -1, 1)))
#   expect_true(is.na(pdweibull(1, 2, 1)))
#   expect_true(is.na(pdweibull(1, 0.5, -1)))
# 
#   expect_true(is.na(pfatigue(1, -1, 1)))
#   expect_true(is.na(pfatigue(1, 1, -1)))
# 
#   expect_true(is.na(pfrechet(1, lambda = -1)))
#   expect_true(is.na(pfrechet(1, sigma = -1)))
# 
#   expect_true(is.na(pgev(1, 1, -1, 1)))
# 
#   expect_true(is.na(pgompertz(1, -1, 1)))
#   expect_true(is.na(pgompertz(1, 1, -1)))
# 
#   expect_true(is.na(pgpd(1, 1, -1, 1)))
# 
#   expect_true(is.na(pgpois(1, -1, 1)))
#   expect_true(is.na(pgpois(1, 1, -1)))
#   expect_true(is.na(pgpois(1, 1, scale = 0)))
# 
#   expect_true(is.na(pgumbel(1, sigma = -1)))
# 
#   expect_true(is.na(phcauchy(1, -1)))
# 
#   expect_true(is.na(phnorm(1, -1)))
# 
#   expect_true(is.na(pht(1, 5, -1)))
# 
#   expect_true(is.na(phuber(1, 0, -1, 1)))
#   expect_true(is.na(phuber(1, 0, 1, -1)))
# 
#   expect_true(is.na(pinvgamma(1, -1, 1)))
#   expect_true(is.na(pinvgamma(1, 1, -1)))
# 
#   expect_true(is.na(pinvchisq(1, -1, 1)))
#   expect_true(is.na(pinvchisq(1, 1, -1)))
# 
#   expect_true(is.na(pkumar(0.5, -1, 1)))
#   expect_true(is.na(pkumar(0.5, 1, -1)))
# 
#   expect_true(is.na(plaplace(1, 0, -1)))
# 
#   expect_true(is.na(plgser(1, -1)))
#   expect_true(is.na(plgser(1, 2)))
# 
#   expect_true(is.na(plomax(1, -1, 1)))
#   expect_true(is.na(plomax(1, 1, -1)))
# 
#   expect_true(is.na(pmixnorm(0, c(1,2,3), c(-1,2,3), c(1/3,1/3,1/3))))
#   expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,-2,3), c(1/3,1/3,1/3))))
#   expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,2,-3), c(1/3,1/3,1/3))))
#   expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,2,3), c(-1,1/3,1/3))))
#   expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,-1,1/3))))
#   expect_true(is.na(pmixnorm(0, c(1,2,3), c(1,2,3), c(1/3,1/3,-1))))
# 
#   expect_true(is.na(pmixpois(0, c(-1,2,3), c(1/3,1/3,1/3))))
#   expect_true(is.na(pmixpois(0, c(1,-2,3), c(1/3,1/3,1/3))))
#   expect_true(is.na(pmixpois(0, c(1,2,-3), c(1/3,1/3,1/3))))
#   expect_true(is.na(pmixpois(0, c(1,2,3), c(-1,1/3,1/3))))
#   expect_true(is.na(pmixpois(0, c(1,2,3), c(1/3,-1,1/3))))
#   expect_true(is.na(pmixpois(0, c(1,2,3), c(1/3,1/3,-1))))
# 
#   expect_true(is.na(pnsbeta(0.5, -1, 1, -2, 2)))
#   expect_true(is.na(pnsbeta(0.5, 1, -1, -2, 2)))
#   expect_true(is.na(pnsbeta(0.5, 1, 1, 2, -2)))
# 
#   expect_true(is.na(pnst(1, -2, 0, 1)))
#   expect_true(is.na(pnst(1, 2, 0, -1)))
# 
#   expect_true(is.na(ppareto(0, -1, 1)))
#   expect_true(is.na(ppareto(0, 1, -1)))
# 
#   expect_true(is.na(pprop(1, -10, 0.5)))
#   expect_true(is.na(pprop(1, 10, -1)))
#   expect_true(is.na(pprop(1, 10, 2)))
# 
#   expect_true(is.na(prayleigh(0, -1)))
# 
#   expect_true(is.na(pslash(1, sigma = -1)))
# 
#   expect_true(is.na(ptnorm(1, 0, -1, -2, 2)))
#   expect_true(is.na(ptnorm(1, 0, 1, 2, -2)))
#   expect_true(is.na(ptnorm(1, 0, 1, 0, 0)))
# 
#   expect_true(is.na(ptpois(1, lambda = -5, s = 0)))
#   expect_true(is.na(ptpois(1, lambda = -5, s = 6)))
#   expect_true(is.na(ptpois(1, lambda = 5, s = -1)))
# 
#   expect_true(is.na(ptriang(1, 0, 0, 0)))
#   expect_true(is.na(ptriang(1, 1, -1, 0)))
#   expect_true(is.na(ptriang(1, -1, 1, 2)))
#   expect_true(is.na(ptriang(1, -1, 1, -2)))
# 
#   expect_true(is.na(pwald(1, 1, -1)))
# 
#   expect_true(is.na(pzip(1, -1, 0.5)))
#   expect_true(is.na(pzip(1, 1, -1)))
#   expect_true(is.na(pzip(1, 1, 2)))
# 
#   expect_true(is.na(pzib(1, -1, 0.5, 0.5)))
#   expect_true(is.na(pzib(1, 1, -1, 0.5)))
#   expect_true(is.na(pzib(1, 1, 0.5, 2)))
#   expect_true(is.na(pzib(1, 1, -1, 0.5)))
#   expect_true(is.na(pzib(1, 1, 0.5, 2)))
# 
#   expect_true(is.na(pzinb(1, -1, 0.5, 0.5)))
#   expect_true(is.na(pzinb(1, 1, -1, 0.5)))
#   expect_true(is.na(pzinb(1, 1, 0.5, 2)))
#   expect_true(is.na(pzinb(1, 1, -1, 0.5)))
#   expect_true(is.na(pzinb(1, 1, 0.5, 2)))
# 
# })
# 
# 
# 
# 
# test_that("Wrong parameter values in inverse CDF functions", {
# 
#   expect_true(is.na(qbetapr(0.5, -1, 1, 1)))
#   expect_true(is.na(qbetapr(0.5, 1, -1, 1)))
#   expect_true(is.na(qbetapr(0.5, 1, 1, -1)))
# 
#   expect_true(is.na(qbern(0.5, -1)))
#   expect_true(is.na(qbern(0.5, 2)))
# 
#   expect_true(is.na(qcat(0.5, c(-1, 0.5))))
#   expect_true(is.na(qcat(0.5, c(0.5, -1))))
# 
#   expect_true(is.na(qdnorm(0.5, sd = -1)))
# 
#   expect_true(is.na(qdunif(0.5, min = 10, max = 1)))
#   expect_true(is.na(qdunif(0.5, min = 0, max = Inf)))
#   expect_true(is.na(qdunif(0.5, min = -Inf, max = Inf)))
#   expect_true(is.na(qdunif(0.5, min = Inf, max = -Inf)))
# 
#   expect_true(is.na(qdweibull(0.5, -1, 1)))
#   expect_true(is.na(qdweibull(0.5, 2, 1)))
#   expect_true(is.na(qdweibull(0.5, 0.5, -1)))
# 
#   expect_true(is.na(qfatigue(0.5, -1, 1)))
#   expect_true(is.na(qfatigue(0.5, 1, -1)))
# 
#   expect_true(is.na(qfrechet(0.5, lambda = -1)))
#   expect_true(is.na(qfrechet(0.5, sigma = -1)))
# 
#   expect_true(is.na(qgev(0.5, 1, -1, 1)))
# 
#   expect_true(is.na(qgompertz(0.5, -1, 1)))
#   expect_true(is.na(qgompertz(0.5, 1, -1)))
# 
#   expect_true(is.na(qgpd(0.5, 1, -1, 1)))
# 
#   expect_true(is.na(qgumbel(0.5, sigma = -1)))
# 
#   expect_true(is.na(qhcauchy(0.5, -1)))
# 
#   expect_true(is.na(qhnorm(0.5, -1)))
# 
#   expect_true(is.na(qht(0.5, 5, -1)))
# 
#   expect_true(is.na(qhuber(0.5, 0, -1, 1)))
#   expect_true(is.na(qhuber(0.5, 0, 1, -1)))
# 
#   expect_true(is.na(qinvgamma(0.5, -1, 1)))
#   expect_true(is.na(qinvgamma(0.5, 1, -1)))
# 
#   expect_true(is.na(qinvchisq(0.5, -1, 1)))
#   expect_true(is.na(qinvchisq(0.5, 1, -1)))
# 
#   expect_true(is.na(qkumar(0.5, -1, 1)))
#   expect_true(is.na(qkumar(0.5, 1, -1)))
# 
#   expect_true(is.na(qlaplace(0.5, 0, -1)))
# 
#   expect_true(is.na(qlgser(0.5, -1)))
#   expect_true(is.na(qlgser(0.5, 2)))
# 
#   expect_true(is.na(qlomax(0.5, -1, 1)))
#   expect_true(is.na(qlomax(0.5, 1, -1)))
# 
#   expect_true(is.na(qnsbeta(0.5, -1, 1, -2, 2)))
#   expect_true(is.na(qnsbeta(0.5, 1, -1, -2, 2)))
#   expect_true(is.na(qnsbeta(0.5, 1, 1, 2, -2)))
# 
#   expect_true(is.na(qnst(0.5, -2, 0, 1)))
#   expect_true(is.na(qnst(0.5, 2, 0, -1)))
# 
#   expect_true(is.na(qpareto(0.5, -1, 1)))
#   expect_true(is.na(qpareto(0.5, 1, -1)))
# 
#   expect_true(is.na(qprop(0.5, -10, 0.5)))
#   expect_true(is.na(qprop(0.5, 10, -1)))
#   expect_true(is.na(qprop(0.5, 10, 2)))
# 
#   expect_true(is.na(qrayleigh(0, -1)))
# 
#   expect_true(is.na(qtnorm(0.5, 0, -1, -2, 2)))
#   expect_true(is.na(qtnorm(0.5, 0, 1, 2, -2)))
#   expect_true(is.na(qtnorm(0.5, 0, 1, 0, 0)))
# 
#   expect_true(is.na(qtpois(0.5, lambda = -5, s = 0)))
#   expect_true(is.na(qtpois(0.5, lambda = -5, s = 6)))
#   expect_true(is.na(qtpois(0.5, lambda = 5, s = -1)))
# 
#   expect_true(is.na(qtriang(0.5, 0, 0, 0)))
#   expect_true(is.na(qtriang(0.5, 1, -1, 0)))
#   expect_true(is.na(qtriang(0.5, -1, 1, 2)))
#   expect_true(is.na(qtriang(0.5, -1, 1, -2)))
# 
#   expect_true(is.na(qzip(0.5, -1, 0.5)))
#   expect_true(is.na(qzip(0.5, 1, -1)))
#   expect_true(is.na(qzip(0.5, 1, 2)))
# 
#   expect_true(is.na(qzib(0.5, -1, 0.5, 0.5)))
#   expect_true(is.na(qzib(0.5, 1, -1, 0.5)))
#   expect_true(is.na(qzib(0.5, 1, 0.5, 2)))
#   expect_true(is.na(qzib(0.5, 1, -1, 0.5)))
#   expect_true(is.na(qzib(0.5, 1, 0.5, 2)))
# 
#   expect_true(is.na(qzinb(0.5, -1, 0.5, 0.5)))
#   expect_true(is.na(qzinb(0.5, 1, -1, 0.5)))
#   expect_true(is.na(qzinb(0.5, 1, 0.5, 2)))
#   expect_true(is.na(qzinb(0.5, 1, -1, 0.5)))
#   expect_true(is.na(qzinb(0.5, 1, 0.5, 2)))
# 
# })
# 
# 
# 
# 
# test_that("Wrong parameter values in RNG functions", {
# 
#   expect_true(is.na(rbbinom(1, -1, 1, 1)))
#   expect_true(is.na(rbbinom(1, 1, -1, 1)))
#   expect_true(is.na(rbbinom(1, 1, 1, -1)))
# 
#   expect_true(is.na(rbetapr(1, -1, 1, 1)))
#   expect_true(is.na(rbetapr(1, 1, -1, 1)))
#   expect_true(is.na(rbetapr(1, 1, 1, -1)))
# 
#   expect_true(is.na(rbern(1, -1)))
#   expect_true(is.na(rbern(1, 2)))
# 
#   expect_true(is.na(rbhatt(1, sigma = -1)))
#   expect_true(is.na(rbhatt(1, a = -1)))
# 
#   expect_true(all(is.na(rbnbinom(1, -1, 1, 1))))
#   expect_true(all(is.na(rbnbinom(1, 1, -1, 1))))
#   expect_true(all(is.na(rbnbinom(1, 1, 1, -1))))
# 
#   expect_true(all(is.na(rbvnorm(1, sd1 = -1))))
#   expect_true(all(is.na(rbvnorm(1, sd2 = -1))))
#   expect_true(all(is.na(rbvnorm(1, cor = -2))))
#   expect_true(all(is.na(rbvnorm(1, cor = 2))))
# 
#   expect_true(all(is.na(rbvpois(1, -1, 1, 1))))
#   expect_true(all(is.na(rbvpois(1, 1, -1, 1))))
#   expect_true(all(is.na(rbvpois(1, 1, 1, -1))))
# 
#   expect_true(is.na(rcat(1, c(-1, 0.5))))
#   expect_true(is.na(rcat(1, c(0.5, -1))))
#   expect_true(is.na(rcat(2, matrix(c(-1, 0.5, 0.5, 0.5), byrow = T, ncol = 2))[1]))
#   expect_true(is.na(rcat(2, matrix(c(0.5, -1, 0.5, 0.5), byrow = T, ncol = 2))[1]))
#   expect_true(is.na(rcat(2, matrix(c(0.5, 0.5, -1, 0.5), byrow = T, ncol = 2))[2]))
#   expect_true(is.na(rcat(2, matrix(c(0.5, 0.5, 0.5, -1), byrow = T, ncol = 2))[2]))
# 
#   expect_true(all(is.na(rdirichlet(1, c(-1, 0.5)))))
#   expect_true(all(is.na(rdirichlet(1, c(0.5, -1)))))
# 
#   expect_true(is.na(rdnorm(1, sd = -1)))
# 
#   expect_true(is.na(rdunif(1, min = 10, max = 1)))
#   expect_true(is.na(rdunif(1, min = 0, max = Inf)))
#   expect_true(is.na(rdunif(1, min = -Inf, max = Inf)))
#   expect_true(is.na(rdunif(1, min = Inf, max = -Inf)))
# 
#   expect_true(is.na(rdweibull(1, -1, 1)))
#   expect_true(is.na(rdweibull(1, 2, 1)))
#   expect_true(is.na(rdweibull(1, 0.5, -1)))
# 
#   expect_true(is.na(rfatigue(1, -1, 1)))
#   expect_true(is.na(rfatigue(1, 1, -1)))
# 
#   expect_true(is.na(rfrechet(1, lambda = -1)))
#   expect_true(is.na(rfrechet(1, sigma = -1)))
# 
#   expect_true(is.na(rgev(1, 1, -1, 1)))
# 
#   expect_true(is.na(rgompertz(1, -1, 1)))
#   expect_true(is.na(rgompertz(1, 1, -1)))
# 
#   expect_true(is.na(rgpd(1, 1, -1, 1)))
# 
#   expect_true(is.na(rgpois(1, -1, 1)))
#   expect_true(is.na(rgpois(1, 1, -1)))
#   expect_true(is.na(rgpois(1, 1, scale = 0)))
# 
#   expect_true(is.na(rgumbel(1, sigma = -1)))
# 
#   expect_true(is.na(rhcauchy(1, -1)))
# 
#   expect_true(is.na(rhnorm(1, -1)))
# 
#   expect_true(is.na(rht(1, 5, -1)))
# 
#   expect_true(is.na(rhuber(1, 0, -1, 1)))
#   expect_true(is.na(rhuber(1, 0, 1, -1)))
# 
#   expect_true(is.na(rinvgamma(1, -1, 1)))
#   expect_true(is.na(rinvgamma(1, 1, -1)))
# 
#   expect_true(is.na(rinvchisq(1, -1, 1)))
#   expect_true(is.na(rinvchisq(1, 1, -1)))
# 
#   expect_true(is.na(rkumar(1, -1, 1)))
#   expect_true(is.na(rkumar(1, 1, -1)))
# 
#   expect_true(is.na(rlaplace(1, 0, -1)))
# 
#   expect_true(is.na(rlgser(1, -1)))
#   expect_true(is.na(rlgser(1, 2)))
# 
#   expect_true(is.na(rlomax(1, -1, 1)))
#   expect_true(is.na(rlomax(1, 1, -1)))
# 
#   expect_true(is.na(rmixnorm(1, c(1,2,3), c(-1,2,3), c(1/3,1/3,1/3))))
#   expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,-2,3), c(1/3,1/3,1/3))))
#   expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,-3), c(1/3,1/3,1/3))))
#   expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(-1,1/3,1/3))))
#   expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,-1,1/3))))
#   expect_true(is.na(rmixnorm(1, c(1,2,3), c(1,2,3), c(1/3,1/3,-1))))
# 
#   expect_true(is.na(rmixpois(1, c(-1,2,3), c(1/3,1/3,1/3))))
#   expect_true(is.na(rmixpois(1, c(1,-2,3), c(1/3,1/3,1/3))))
#   expect_true(is.na(rmixpois(1, c(1,2,-3), c(1/3,1/3,1/3))))
#   expect_true(is.na(rmixpois(1, c(1,2,3), c(-1,1/3,1/3))))
#   expect_true(is.na(rmixpois(1, c(1,2,3), c(1/3,-1,1/3))))
#   expect_true(is.na(rmixpois(1, c(1,2,3), c(1/3,1/3,-1))))
# 
#   expect_true(all(is.na(rdirmnom(1, 1.5, c(1, 1, 1)))))
#   expect_true(all(is.na(rdirmnom(1, -3, c(1, 1, 1)))))
#   expect_true(all(is.na(rdirmnom(1, 3, c(-1, 1, 1)))))
#   expect_true(all(is.na(rdirmnom(1, 3, c(1, -1, 1)))))
#   expect_true(all(is.na(rdirmnom(1, 3, c(1, 1, -1)))))
# 
#   expect_true(all(is.na(rmnom(1, 3, c(-1, 1/3, 1/3)))))
#   expect_true(all(is.na(rmnom(1, 3, c(1/3, -1, 1/3)))))
#   expect_true(all(is.na(rmnom(1, 3, c(1/3, 1/3, -1)))))
# 
#   expect_true(all(is.na(rmvhyper(1, c(2,3,4), 99))))
#   expect_true(all(is.na(rmvhyper(1, c(-2,3,4), 99))))
#   expect_true(all(is.na(rmvhyper(1, c(2,-3,4), 99))))
#   expect_true(all(is.na(rmvhyper(1, c(2,3,-4), 99))))
# 
#   expect_true(is.na(rnsbeta(1, -1, 1, -2, 2)))
#   expect_true(is.na(rnsbeta(1, 1, -1, -2, 2)))
#   expect_true(is.na(rnsbeta(1, 1, 1, 2, -2)))
# 
#   expect_true(is.na(rnst(1, -2, 0, 1)))
#   expect_true(is.na(rnst(1, 2, 0, -1)))
# 
#   expect_true(is.na(rpareto(1, -1, 1)))
#   expect_true(is.na(rpareto(1, 1, -1)))
# 
#   expect_true(is.na(rprop(1, -10, 0.5)))
#   expect_true(is.na(rprop(1, 10, -1)))
#   expect_true(is.na(rprop(1, 10, 2)))
# 
#   expect_true(is.na(rrayleigh(1, -1)))
# 
#   expect_true(is.na(rskellam(1, -1, 1)))
#   expect_true(is.na(rskellam(1, 1, -1)))
# 
#   expect_true(is.na(rslash(1, sigma = -1)))
# 
#   expect_true(is.na(rtnorm(1, 0, -1, -2, 2)))
#   expect_true(is.na(rtnorm(1, 0, 1, 2, -2)))
#   expect_true(is.na(rtnorm(1, 0, 1, 0, 0)))
# 
#   expect_true(is.na(rtpois(1, lambda = -5, s = 0)))
#   expect_true(is.na(rtpois(1, lambda = -5, s = 6)))
#   expect_true(is.na(rtpois(1, lambda = 5, s = -1)))
# 
#   expect_true(is.na(rtriang(1, 0, 0, 0)))
#   expect_true(is.na(rtriang(1, 1, -1, 0)))
#   expect_true(is.na(rtriang(1, -1, 1, 2)))
#   expect_true(is.na(rtriang(1, -1, 1, -2)))
# 
#   expect_true(is.na(rwald(1, 1, -1)))
# 
#   expect_true(is.na(rzip(1, -1, 0.5)))
#   expect_true(is.na(rzip(1, 1, -1)))
#   expect_true(is.na(rzip(1, 1, 2)))
# 
#   expect_true(is.na(rzib(1, -1, 0.5, 0.5)))
#   expect_true(is.na(rzib(1, 1, -1, 0.5)))
#   expect_true(is.na(rzib(1, 1, 0.5, 2)))
#   expect_true(is.na(rzib(1, 1, -1, 0.5)))
#   expect_true(is.na(rzib(1, 1, 0.5, 2)))
# 
#   expect_true(is.na(rzinb(1, -1, 0.5, 0.5)))
#   expect_true(is.na(rzinb(1, 1, -1, 0.5)))
#   expect_true(is.na(rzinb(1, 1, 0.5, 2)))
#   expect_true(is.na(rzinb(1, 1, -1, 0.5)))
#   expect_true(is.na(rzinb(1, 1, 0.5, 2)))
# 
# })