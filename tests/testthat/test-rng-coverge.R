

probCoverage <- function(stub, ..., n = 5000L) {
  rfoo <- eval(as.name(paste0("r", stub)))
  pfoo <- eval(as.name(paste0("p", stub)))
  diff(pfoo(range(rfoo(n, ...)), ...))
}


test_that("Coverage of RNG's", {
  
  skip_on_cran()

  expect_gte(probCoverage("betapr", 1, 1, 1), 0.99)
  
  expect_gte(probCoverage("bhatt", 1, 1, 1), 0.99)
  
  expect_gte(probCoverage("fatigue", 1, 1), 0.99)
  
  expect_gte(probCoverage("frechet", 1, 1, 1), 0.99)
  
  expect_gte(probCoverage("gev", 1, 1, 1), 0.99)
  
  expect_gte(probCoverage("gompertz", 1, 1), 0.99)
  
  expect_gte(probCoverage("gpd", 1, 1, 1), 0.99)
  
  expect_gte(probCoverage("gumbel", 1, 1), 0.99)
  
  expect_gte(probCoverage("hcauchy", 1), 0.99)
  
  expect_gte(probCoverage("hnorm", 1), 0.99)
  
  expect_gte(probCoverage("ht", 5, 1), 0.99)
  
  expect_gte(probCoverage("huber", 0, 1, 1), 0.99)
  
  expect_gte(probCoverage("invgamma", 1, 1), 0.99)
  
  expect_gte(probCoverage("invchisq", 1, 1), 0.99)
  
  expect_gte(probCoverage("kumar", 1, 1), 0.99)
  expect_gte(probCoverage("kumar", 100, 1), 0.99)
  expect_gte(probCoverage("kumar", 1, 100), 0.99)
  expect_gte(probCoverage("kumar", 100, 100), 0.99)
  
  expect_gte(probCoverage("laplace", 0, 1), 0.99)
  expect_gte(probCoverage("laplace", 0, 1000), 0.99)
  
  expect_gte(probCoverage("lomax", 1, 0.001), 0.99)
  expect_gte(probCoverage("lomax", 1, 0.5), 0.99)
  expect_gte(probCoverage("lomax", 1, 0.999), 0.99)
  
  expect_gte(probCoverage("mixnorm", c(1,2,3), c(1,2,3), c(1/3,1/3,1/3)), 0.99)
  
  expect_gte(probCoverage("nsbeta", 1, 1, -2, 2), 0.99)
  
  expect_gte(probCoverage("lst", 2, 0, 1), 0.99)
  
  expect_gte(probCoverage("pareto", 1, 1), 0.99)
  
  expect_gte(probCoverage("prop", 10, 0.5), 0.99)
  expect_gte(probCoverage("prop", 100, 0.5), 0.99)
  expect_gte(probCoverage("prop", 1000, 0.5), 0.99)
  expect_gte(probCoverage("prop", 10, 0.01), 0.99)
  expect_gte(probCoverage("prop", 10, 0.5), 0.99)
  expect_gte(probCoverage("prop", 10, 0.99), 0.99)
  
  expect_gte(probCoverage("rayleigh", 1), 0.99)
  
  expect_gte(probCoverage("sgomp", 0.4, 1), 0.99)
  
  expect_gte(probCoverage("slash", 1, 1), 0.99)
  
  expect_gte(probCoverage("tnorm", 0, 1, -1, 1), 0.99)
  expect_gte(probCoverage("tnorm", 0, 1, -2, 2), 0.99)
  expect_gte(probCoverage("tnorm", 0, 1, 2, Inf), 0.99)
  expect_gte(probCoverage("tnorm", 0, 1, 4, Inf), 0.99)
  expect_gte(probCoverage("tnorm", 0, 1, -Inf, -2), 0.99)
  expect_gte(probCoverage("tnorm", 0, 1, -Inf, -4), 0.99)
  expect_gte(probCoverage("tnorm", 0, 1, -6, -4), 0.99)
  expect_gte(probCoverage("tnorm", 0, 1, 4, 6), 0.99)
  
  expect_gte(probCoverage("triang"), 0.99)
  expect_gte(probCoverage("triang", 0, 1, 0.5), 0.99)
  
  expect_gte(probCoverage("wald", 1, 1), 0.99)
  
})

