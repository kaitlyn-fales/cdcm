test_that("get_num_param returns correct number of parameters", {
  toy_stan_dat <- list(
    d_A = 4,
    d_B = 2,
    d_C = 1,
    m = 3
  )

  result <- get_num_param(toy_stan_dat)

  # expected = 4 + 2 + 1 + 3*3 = 7 + 9 = 16
  expect_equal(result, 16)
})

test_that("get_num_param returns a scalar numeric value", {
  toy_stan_dat <- list(
    d_A = 1,
    d_B = 1,
    d_C = 1,
    m = 1
  )

  result <- get_num_param(toy_stan_dat)

  expect_length(result, 1)
  expect_true(is.numeric(result))
})

test_that("get_num_param errors with informative message when fields missing", {
  bad_dat <- list(d_A = 1, d_B = 1, m = 2)

  expect_error(
    get_num_param(bad_dat),
    "missing required component"
  )
})
