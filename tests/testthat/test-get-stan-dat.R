toy_dat <- list(
  times = c(2, 4, 6, 8),
  u = matrix(c(
    0,
    1,
    0,
    1
  ), ncol = 1),
  y_obs = matrix(c(
    0.1, 0.2,
    0.0, 0.1,
    0.2, 0.3,
    0.1, 0.4
  ), ncol = 2, byrow = TRUE)
)

toy_idxs <- list(
  A_idxs = matrix(c(
    1, 1,
    2, 1,
    1, 2,
    2, 2
  ), ncol = 2, byrow = TRUE),
  B_idxs = matrix(c(
    2, 1, 1
  ), ncol = 3, byrow = TRUE),
  C_idxs = matrix(c(
    1, 1
  ), ncol = 2, byrow = TRUE)
)

test_that("get_stan_dat returns a list and warns for weak toy design", {
  expect_warning(
    stan_dat <- get_stan_dat(toy_dat, toy_idxs),
    "Experimental design may not satisfy sufficient identifiability conditions"
  )

  expect_type(stan_dat, "list")
  expect_true(length(stan_dat) > 0)
})

test_that("get_stan_dat errors for malformed idxs", {
  bad_idxs <- list(
    A_idxs = matrix(c(1, 1), ncol = 2),
    B_idxs = NULL,
    C_idxs = matrix(c(1, 1), ncol = 2)
  )

  expect_error(
    get_stan_dat(toy_dat, bad_idxs),
    "`hypothesis_idxs\\$B_idxs` must be a matrix."
  )
})

test_that("get_stan_dat errors when u has wrong number of rows", {
  bad_dat <- toy_dat
  bad_dat$u <- bad_dat$u[-1, , drop = FALSE]

  expect_error(
    get_stan_dat(bad_dat, toy_idxs),
    "same number of rows"
  )
})

test_that("get_stan_dat errors when times are not strictly increasing", {
  bad_dat <- toy_dat
  bad_dat$times <- c(1, 2, 2, 3)

  expect_error(
    get_stan_dat(bad_dat, toy_idxs),
    "strictly increasing"
  )
})

test_that("get_stan_dat errors when times are not equally spaced", {
  bad_dat <- toy_dat
  bad_dat$times <- c(1, 2, 4, 5)

  expect_error(
    get_stan_dat(bad_dat, toy_idxs),
    "equally spaced"
  )
})

test_that("get_stan_dat errors when times are not positive", {
  bad_dat <- toy_dat
  bad_dat$times <- c(0, 1, 2, 3)

  expect_error(
    get_stan_dat(bad_dat, toy_idxs),
    "positive values"
  )
})
