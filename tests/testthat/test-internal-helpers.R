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

test_that("get_AnalyticSol handles an input change at the final row", {
  parms <- list(
    A = matrix(c(-0.5, 0.2, 0, -0.5), 2, 2),
    B = list(matrix(0, 2, 2)),
    C = matrix(c(1, 0), 2, 1)
  )
  times <- 0:6

  # Input switches on at the last row, so the terminal row is listed twice
  input_end <- matrix(c(0, 1, 1, 0, 0, 0, 1), ncol = 1)
  ch_end <- get_changepoints(input_end)
  expect_equal(ch_end, c(2, 4, 7, 7))

  # Same design without the final switch; the last-row input never acts on z,
  # so the two solutions must agree
  input_ref <- input_end
  input_ref[7, 1] <- 0

  sol_end <- get_AnalyticSol(input_end, ch_end, times, parms, m = 2)
  sol_ref <- get_AnalyticSol(input_ref, get_changepoints(input_ref), times, parms, m = 2)

  expect_true(all(is.finite(sol_end)))
  expect_equal(sol_end, sol_ref)
})

