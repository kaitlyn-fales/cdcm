test_that("get_changepoints returns terminal row when there are no changes", {
  input <- matrix(0, nrow = 4, ncol = 2)

  result <- get_changepoints(input)

  expect_equal(result, 4)
})

test_that("get_changepoints detects a single change and appends terminal row", {
  input <- matrix(
    c(0, 0,
      0, 0,
      1, 0,
      1, 0),
    ncol = 2,
    byrow = TRUE
  )

  result <- get_changepoints(input)

  expect_equal(result, c(3, 4))
})

test_that("get_changepoints deduplicates simultaneous changes across columns", {
  input <- matrix(
    c(0, 0,
      0, 0,
      1, 1,
      1, 1),
    ncol = 2,
    byrow = TRUE
  )

  result <- get_changepoints(input)

  expect_equal(result, c(3, 4))
})

test_that("get_changepoints returns all change rows plus terminal row", {
  input <- matrix(
    c(0, 0,
      1, 0,
      1, 0,
      1, 1,
      0, 1),
    ncol = 2,
    byrow = TRUE
  )

  result <- get_changepoints(input)

  expect_equal(result, c(2, 4, 5, 5))
})

test_that("get_changepoints errors when input is not a matrix", {
  input <- c(0, 0, 1, 1)

  expect_error(
    get_changepoints(input),
    "`input` must be a matrix."
  )
})

test_that("get_changepoints errors when input has no rows", {
  input <- matrix(numeric(0), nrow = 0, ncol = 2)

  expect_error(
    get_changepoints(input),
    "at least one row"
  )
})

test_that("get_changepoints errors when input is not numeric", {
  input <- matrix(c("a", "a", "b", "b"), ncol = 1)

  expect_error(
    get_changepoints(input),
    "`input` must be numeric."
  )
})

test_that("get_last_draws_for_init extracts last draw by parameter prefix for one chain", {
  draws <- setNames(
    data.frame(
      c(1, 1),
      c(1, 2),
      c(1, 2),
      c(-10, -9),
      c(1.0, 1.1),
      c(1.5, 1.6),
      c(0.1, 0.2),
      c(0.5, 0.6),
      check.names = FALSE
    ),
    c(".chain", ".iteration", ".draw", "lp__", "sigma[1]", "sigma[2]", "beta[1]", "nu_A[1]")
  )

  result <- get_last_draws_for_init(draws)

  expect_true(is.list(result))
  expect_named(result, c("sigma", "beta", "nu_A"))
  expect_equal(result$sigma, c(1.1, 1.6))
  expect_equal(result$beta, 0.2)
  expect_equal(result$nu_A, 0.6)
})

test_that("get_pathfinder_init extracts last draw by parameter prefix for one chain", {
  draws <- setNames(
    data.frame(
      c(1, 1),
      c(1, 2),
      c(1, 2),
      c(-10, -9),
      c(-11, -10),
      c(1.0, 1.1),
      c(1.5, 1.6),
      c(0.1, 0.2),
      c(0.5, 0.6),
      check.names = FALSE
    ),
    c(".chain", ".iteration", ".draw", "lp__", "lp_approx__", "sigma[1]", "sigma[2]", "beta[1]", "nu_A[1]")
  )

  result <- get_pathfinder_init(draws)

  expect_true(is.list(result))
  expect_named(result, c("sigma", "beta", "nu_A"))
  expect_equal(result$sigma, c(1.1, 1.6))
  expect_equal(result$beta, 0.2)
  expect_equal(result$nu_A, 0.6)
})

