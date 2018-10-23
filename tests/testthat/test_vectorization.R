
test_that("Fix for #16 works", {
  
  test_dat <- data.frame(
      q = c(1, 10, 10, 100),
      size = c(10, 100, 1000, 1000)
    )
  
  # This failed in #16
  with(
    test_dat,
    pbbinom(q, size, 1, 1)
  )
  
  with(
    test_dat,
    pbnbinom(q, size, 1, 1)
  )
  
})