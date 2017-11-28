
source("testthat/helper_pure_r_implementations.R")

x <- c(-1e5, -100, -10, -5.5, -5, -1.01, -1, -0.5, 0.001, 0,
       0.001, 0.5, 1, 1.01, 5, 5.5, 10, 100, 1e5)
n <- length(x)

nsim <- 50000


if (requireNamespace("rbenchmark", quietly = TRUE)) {

  print(benchmark(dbbinom(x, 100, 1, 10), dbbinomR(x, 100, 1, 10), 
            replications = nsim))
  print(benchmark(dbnbinom(x[1:(n-2)], 100, 1, 10),
            dbnbinomR(x[1:(n-2)], 100, 1, 10), 
            replications = nsim))
  print(benchmark(dbetapr(x, 1, 1, 1), dbetaprR(x, 1, 1, 1), 
            replications = nsim))
  print(benchmark(dfatigue(x, 1, 1, 0), dfatigueR(x, 1, 1, 0), 
            replications = nsim))
  print(benchmark(ddlaplace(x, 0, 0.5), ddlaplaceR(x, 0, 0.5), 
            replications = nsim))
  print(benchmark(ddweibull(x, 0.5, 1), ddweibullR(x, 0.5, 1), 
            replications = nsim))
  print(benchmark(dfrechet(x, 1, 1, 1), dfrechetR(x, 1, 1, 1), 
            replications = nsim))
  print(benchmark(dgpois(x[1:(n-1)], 1, 1),
               dgpoisR(x[1:(n-1)], 1, 1), 
            replications = nsim))
  print(benchmark(dgev(x, 1, 1, 1), dgevR(x, 1, 1, 1), 
            replications = nsim))
  print(benchmark(dgompertz(x, 1, 1), dgompertzR(x, 1, 1), 
            replications = nsim))
  print(benchmark(dgpd(x, 1, 1, 1), dgpdR(x, 1, 1, 1), 
            replications = nsim))
  print(benchmark(dgumbel(x, 1, 1), dgumbelR(x, 1, 1), 
            replications = nsim))
  print(benchmark(dinvgamma(x, 1, 1), dinvgammaR(x, 1, 1), 
            replications = nsim))
  print(benchmark(dlaplace(x, -1, 5), dlaplaceR(x, -1, 5), 
            replications = nsim))
  print(benchmark(dlgser(x, 0.5), dlgserR(x, 0.5), 
            replications = nsim))
  print(benchmark(dlomax(x, 1, 0.5), dlomaxR(x, 1, 0.5), 
            replications = nsim))
  print(benchmark(dpareto(x, 1, 1), dparetoR(x, 1, 1), 
            replications = nsim))
  print(benchmark(dpower(x, 1, 1), dpowerR(x, 1, 1), 
            replications = nsim))
  print(benchmark(drayleigh(x, 1), drayleighR(x, 1), 
            replications = nsim))
  print(benchmark(dsgomp(x, 0.5, 1), dsgompR(x, 0.5, 1), 
            replications = nsim))
  
  
  if (requireNamespace("hoa", quietly = TRUE)) {
    
    print(benchmark(dhuber(x), hoa::dHuber(x), 
          replications = nsim))
    print(benchmark(phuber(x), hoa::pHuber(x), 
          replications = nsim))
    
  }
  
  
  if (requireNamespace("evd", quietly = TRUE)) {
    
    print(benchmark(dgev(x), evd::dgev(x), 
          replications = nsim))
    print(benchmark(pgev(x), evd::pgev(x), 
          replications = nsim))
    
    print(benchmark(dgpd(x), evd::dgpd(x), 
          replications = nsim))
    print(benchmark(pgpd(x), evd::pgpd(x), 
          replications = nsim))
    
  }
  
  
  if (requireNamespace("Compositional", quietly = TRUE)) {
    
    alpha <- runif(5, 0, 3)
    x <- rdirichlet(5000, alpha)
    
    print(benchmark(ddirichlet(x, alpha),
                    Compositional::ddiri(x, alpha, logged = FALSE), 
          replications = 500))
    
  }
  
  
  if (requireNamespace("actuar", quietly = TRUE)) {

    print(benchmark(dzib(x, 45, 0.7, 0.2), actuar::dzmbinom(x, 45, 0.7, 0.2), 
                    replications = nsim))
    print(benchmark(pzib(x, 45, 0.7, 0.2), actuar::pzmbinom(x, 45, 0.7, 0.2), 
                    replications = nsim))
    
    print(benchmark(dzinb(x, 45, 0.7, 0.2), actuar::dzmnbinom(x, 45, 0.7, 0.2), 
                    replications = nsim))
    print(benchmark(pzinb(x, 45, 0.7, 0.2), actuar::pzmnbinom(x, 45, 0.7, 0.2), 
                    replications = nsim))
    
    print(benchmark(dzip(x, 7, 0.2), actuar::dzmpois(x, 7, 0.2), 
                    replications = nsim))
    print(benchmark(pzip(x, 7, 0.2), actuar::pzmpois(x, 7, 0.2), 
                    replications = nsim))
    
  }
  
  
  if (requireNamespace("skellam", quietly = TRUE)) {

    x <- extraDistr::rskellam(5000, 7, 8)
    
    print(benchmark(
      extraDistr::dskellam(x, 7, 8),
      skellam::dskellam(x, 7, 8),
      replications = 5000
    ))
    
  }
  
  
  n <- 100
  p <- runif(5)
  p <- p/sum(p)
  
  x <- rmnom(5000, n, p)
  
  print(benchmark(dmnom(x, n, p), apply(x, 1, dmultinom, n, p), 
        replications = 500))
  
  
  x <- rbvpois(1000, 7, 8, 5)
  
  print(benchmark(
    dbvpois(x[,1], x[,2], 7, 8, 5),
    pbivpois(x[,1], x[,2], 7, 8, 5), 
    replications = 500
  ))

}

  