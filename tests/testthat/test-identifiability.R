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
      "cond_A1",
      "cond_A2",
      "all_ok",
      "equally_spaced",
      "n_eligible_blocks",
      "rank_U_aug"
    )
  )

  expect_true(result$equally_spaced)
  expect_equal(result$n_eligible_blocks, 4)
  expect_equal(result$rank_U_aug, 3)
  expect_true(result$cond_A1)
  expect_true(result$cond_A2)
  expect_true(result$all_ok)
})

test_that("check_design_identifiability detects too few eligible blocks", {
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

  expect_equal(result$n_eligible_blocks, 0)
  expect_false(result$cond_A1)
  expect_false(result$all_ok)
})

test_that("check_design_identifiability detects non-full-rank block inputs", {
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

  expect_true(result$cond_A1)
  expect_lt(result$rank_U_aug, 3)
  expect_false(result$cond_A2)
  expect_false(result$all_ok)
})

test_that("check_design_identifiability allows extra unusable blocks if enough usable blocks exist", {
  toy_data <- list(
    times = 1:14,
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
        1, 1,
        0, 0,
        1, 0),
      ncol = 2,
      byrow = TRUE
    ),
    y_obs = matrix(rnorm(14 * 2), ncol = 2)
  )

  result <- check_design_identifiability(toy_data, warn = FALSE)

  expect_equal(result$n_eligible_blocks, 4)
  expect_true(result$cond_A1)
  expect_true(result$cond_A2)
  expect_true(result$all_ok)
})


test_that("check_latent_identifiability passes for a simple identifiable case", {
  A <- matrix(c(-0.5, 0,
                0, -1.0), nrow = 2, byrow = TRUE)

  B <- list(matrix(0, nrow = 2, ncol = 2))

  u <- matrix(
    c(0,
      0,
      0,
      1,
      1,
      1),
    ncol = 1
  )

  z <- matrix(
    c(0, 0,
      1, 0,
      0, 1,
      1, 1,
      2, 1,
      1, 2),
    ncol = 2,
    byrow = TRUE
  )

  result <- check_latent_identifiability(A, B, u, z, warn = FALSE)

  expect_true(result$cond_A3)
  expect_true(result$cond_A4)
  expect_true(result$all_ok)
  expect_true(is.data.frame(result$block_results))
})

test_that("check_latent_identifiability detects repeated eigenvalues", {
  A <- diag(c(-0.5, -0.5))
  B <- list(matrix(0, nrow = 2, ncol = 2))

  u <- matrix(rep(0, 3), ncol = 1)

  z <- matrix(
    c(0, 0,
      1, 0,
      0, 1),
    ncol = 2,


    byrow = TRUE
  )

  result <- check_latent_identifiability(A, B, u, z, warn = FALSE)

  expect_false(result$cond_A3)
  expect_true(result$cond_A4)
  expect_false(result$all_ok)
})

test_that("check_latent_identifiability detects non-invertible latent data matrix", {
  A <- matrix(c(-0.5, 0,
                0, -1.0), nrow = 2, byrow = TRUE)

  B <- list(matrix(0, nrow = 2, ncol = 2))

  u <- matrix(rep(0, 3), ncol = 1)

  z <- matrix(
    c(1, 1,
      1, 1,
      1, 1),
    ncol = 2,
    byrow = TRUE
  )

  result <- check_latent_identifiability(A, B, u, z, warn = FALSE)

  expect_true(result$cond_A3)
  expect_false(result$cond_A4)
  expect_false(result$all_ok)
})

test_that("simulate_cdcm optionally returns identifiability diagnostics", {
  nu <- list(
    nu_A = c(0.4, 0.3, -0.1, 0.15),
    nu_B = c(-0.2),
    nu_C = c(0.7)
  )

  hypothesis_idxs <- list(
    A_idxs = matrix(c(2,1,
                      1,2,
                      1,1,
                      2,2), byrow = TRUE, ncol = 2),
    B_idxs = matrix(c(1,2,1), byrow = TRUE, ncol = 3),
    C_idxs = matrix(c(1,1), byrow = TRUE, ncol = 2)
  )

  U <- matrix(c(0,0,0,0,1,1,1,1,0,0), ncol = 1)
  sim <- simulate_cdcm(
    nu = nu,
    hypothesis_idxs = hypothesis_idxs,
    nscan = 10,
    m = 2,
    n_u = 1,
    TR = 2,
    U = U,
    SNR = 5,
    seed = 1,
    check_solution = FALSE,
    check_identifiability = TRUE
  )

  expect_named(sim, c("simulated_data", "true_signals", "identifiability"))
  expect_named(sim$identifiability, c("design", "latent"))
  expect_true(is.list(sim$identifiability$design))
  expect_true(is.list(sim$identifiability$latent))
})


