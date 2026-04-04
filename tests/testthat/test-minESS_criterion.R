test_that("minESS_criterion returns a numeric scalar", {
  toy_stan_dat <- list(
    d_A = 4,
    d_B = 2,
    d_C = 1,
    m = 3
  )

  result <- minESS_criterion(toy_stan_dat)

  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(is.finite(result))
})

test_that("minESS_criterion matches direct mcmcse::minESS calculation", {
  toy_stan_dat <- list(
    d_A = 4,
    d_B = 2,
    d_C = 1,
    m = 3
  )

  num_param <- get_num_param(toy_stan_dat)
  expected <- as.numeric(mcmcse::minESS(num_param, alpha = 0.05, eps = 0.05))

  result <- minESS_criterion(toy_stan_dat, alpha = 0.05, eps = 0.05)

  expect_equal(result, expected)
})
