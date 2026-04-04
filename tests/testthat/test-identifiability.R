test_that("get_block_lengths returns full length when there are no changepoints", {
  input <- matrix(0, nrow = 5, ncol = 2)

  result <- get_block_lengths(input)

  expect_equal(result, 5)
})

test_that("get_block_lengths returns correct lengths for one changepoint", {
  input <- matrix(
    c(0, 0,
      0, 0,
      1, 0,
      1, 0,
      1, 0),
    ncol = 2,
    byrow = TRUE
  )

  result <- get_block_lengths(input)

  expect_equal(result, c(2, 3))
})

test_that("get_block_lengths returns correct lengths for multiple changepoints", {
  input <- matrix(
    c(0, 0,
      0, 0,
      1, 0,
      1, 0,
      1, 1,
      1, 1,
      0, 1),
    ncol = 2,
    byrow = TRUE
  )

  result <- get_block_lengths(input)

  expect_equal(result, c(2, 2, 2, 1))
})

test_that("get_block_lengths errors when input is not a matrix", {
  input <- c(0, 0, 1, 1)

  expect_error(
    get_block_lengths(input),
    "`input` must be a matrix."
  )
})

test_that("get_block_lengths errors when input has no rows", {
  input <- matrix(numeric(0), nrow = 0, ncol = 2)

  expect_error(
    get_block_lengths(input),
    "at least one row"
  )
})

test_that("check_design_identifiability returns expected structure for a strong design", {
  toy_data <- list(
    times = 1:12,
    u = matrix(
      c(0, 0,
        0, 0,
        0, 0,
        1, 0,
        1, 0,
        1, 0,
        0, 1,
        0, 1,
        0, 1,
        1, 1,
        1, 1,
        1, 1),
      ncol = 2,
      byrow = TRUE
    ),
    y_obs = matrix(rnorm(12 * 2), ncol = 2)
  )

  result <- check_design_identifiability(toy_data, warn = FALSE)

  expect_true(is.list(result))
  expect_named(
    result,
    c(
      "T", "m", "n_u", "changepoints", "block_lengths", "first_block_zero",
      "n_distinct_block_inputs", "cond_A1_first", "cond_A1_later",
      "cond_A1", "cond_A2", "cond_A3", "all_ok"
    )
  )

  expect_equal(result$T, 12)
  expect_equal(result$m, 2)
  expect_equal(result$n_u, 2)
  expect_true(result$cond_A1)
  expect_true(result$cond_A2)
  expect_true(result$cond_A3)
  expect_true(result$all_ok)
})

test_that("check_design_identifiability detects short blocks", {
  toy_data <- list(
    times = 1:8,
    u = matrix(
      c(0, 0,
        0, 0,
        1, 0,
        1, 0,
        0, 1,
        0, 1,
        1, 1,
        1, 1),
      ncol = 2,
      byrow = TRUE
    ),
    y_obs = matrix(rnorm(8 * 2), ncol = 2)
  )

  result <- check_design_identifiability(toy_data, warn = FALSE)

  expect_false(result$cond_A1)
  expect_false(result$all_ok)
})

test_that("check_design_identifiability detects insufficient number of observations", {
  toy_data <- list(
    times = 1:9,
    u = matrix(
      c(0, 0,
        0, 0,
        0, 0,
        1, 0,
        1, 0,
        1, 0,
        0, 1,
        0, 1,
        0, 1),
      ncol = 2,
      byrow = TRUE
    ),
    y_obs = matrix(rnorm(9 * 2), ncol = 2)
  )

  result <- check_design_identifiability(toy_data, warn = FALSE)

  expect_false(result$cond_A2)
  expect_false(result$all_ok)
})

test_that("check_design_identifiability detects too few distinct blockwise inputs", {
  toy_data <- list(
    times = 1:12,
    u = matrix(
      c(0, 0,
        0, 0,
        0, 0,
        1, 0,
        1, 0,
        1, 0,
        0, 0,
        0, 0,
        0, 0,
        1, 0,
        1, 0,
        1, 0),
      ncol = 2,
      byrow = TRUE
    ),
    y_obs = matrix(rnorm(12 * 2), ncol = 2)
  )

  result <- check_design_identifiability(toy_data, warn = FALSE)

  expect_false(result$cond_A3)
  expect_false(result$all_ok)
})

test_that("check_design_identifiability warns for short blocks", {
  toy_data <- list(
    times = 1:12,
    u = matrix(
      c(0, 0,
        0, 0,
        1, 0,
        1, 0,
        0, 1,
        0, 1,
        1, 1,
        1, 1,
        0, 0,
        0, 0,
        1, 0,
        1, 0),
      ncol = 2,
      byrow = TRUE
    ),
    y_obs = matrix(rnorm(12 * 2), ncol = 2)
  )

  expect_warning(
    check_design_identifiability(toy_data, warn = TRUE),
    "Block-length condition failed"
  )
})

test_that("check_design_identifiability errors when required data components are missing", {
  bad_data <- list(
    times = 1:10,
    u = matrix(0, nrow = 10, ncol = 2)
  )

  expect_error(
    check_design_identifiability(bad_data, warn = FALSE),
    "missing required component"
  )
})


