test_that("simulate_cdcm returns objects with expected structure", {
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

  out <- simulate_cdcm(
    nu = nu,
    hypothesis_idxs = hypothesis_idxs,
    nscan = 10,
    m = 2,
    n_u = 1,
    TR = 2,
    U = U,
    SNR = 5,
    check_solution = FALSE,
    seed = 1
  )

  expect_type(out, "list")
  expect_named(out, c("simulated_data", "true_signals"))

  expect_equal(length(out$simulated_data$times), 10)
  expect_equal(dim(out$simulated_data$u), c(10, 1))
  expect_equal(dim(out$simulated_data$y_obs), c(10, 2))

  expect_equal(dim(out$true_signals$z_signal), c(10, 2))
  expect_equal(dim(out$true_signals$y_signal), c(10, 2))
})

test_that("simulate_cdcm is reproducible with a fixed seed", {
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

  out1 <- simulate_cdcm(
    nu = nu,
    hypothesis_idxs = hypothesis_idxs,
    nscan = 10,
    m = 2,
    n_u = 1,
    TR = 2,
    U = U,
    SNR = 5,
    check_solution = FALSE,
    seed = 123
  )

  out2 <- simulate_cdcm(
    nu = nu,
    hypothesis_idxs = hypothesis_idxs,
    nscan = 10,
    m = 2,
    n_u = 1,
    TR = 2,
    U = U,
    SNR = 5,
    check_solution = FALSE,
    seed = 123
  )

  expect_equal(out1$simulated_data$y_obs, out2$simulated_data$y_obs)
})

test_that("simulate_cdcm errors for non-binary U", {
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

  U <- matrix(c(0, 1, 0, 2, 0, 1, 0, 1, 0, 1), ncol = 1)

  expect_error(
    simulate_cdcm(
      nu = nu,
      hypothesis_idxs = hypothesis_idxs,
      nscan = 10,
      m = 2,
      n_u = 1,
      TR = 2,
      U = U,
      SNR = 5,
      check_solution = FALSE
    ),
    "All entries of `U` must be in \\{0, 1\\}"
  )
})

test_that("handle_seed works as expected", {
  expect_equal(handle_seed(NULL, 3), 1:3)
  expect_equal(handle_seed(5, 3), rep(5, 3))
  expect_equal(handle_seed(c(1, 2, 3), 3), c(1, 2, 3))
  expect_error(handle_seed(c(1, 2), 3))
})

test_that("simulate_cdcm accepts U with nscan rows", {
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

  expect_message(
    out <- simulate_cdcm(
      nu = nu,
      hypothesis_idxs = hypothesis_idxs,
      nscan = 10,
      m = 2,
      n_u = 1,
      TR = 2,
      U = U,
      SNR = 5,
      check_solution = FALSE,
      seed = 1
    ),
    "prepended"
  )

  expect_equal(dim(out$true_signals$z_signal), c(10, 2))
})

test_that("simulate_cdcm accepts U with nscan plus one rows", {
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

  U_scan <- matrix(c(0,0,0,0,1,1,1,1,0,0), ncol = 1)
  U_aug <- rbind(U_scan[0, , drop = FALSE], U_scan)

  out <- simulate_cdcm(
    nu = nu,
    hypothesis_idxs = hypothesis_idxs,
    nscan = 10,
    m = 2,
    n_u = 1,
    TR = 2,
    U = U_aug,
    SNR = 5,
    check_solution = FALSE,
    seed = 1
  )

  expect_equal(dim(out$simulated_data$u), c(10, 1))
})

test_that("simulate_cdcm errors for bad U dimensions", {
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

  U <- matrix(0, nrow = 8, ncol = 1)

  expect_error(
    simulate_cdcm(
      nu = nu,
      hypothesis_idxs = hypothesis_idxs,
      nscan = 10,
      m = 2,
      n_u = 1,
      TR = 2,
      U = U,
      SNR = 5,
      check_solution = FALSE
    ),
    "must have either `nscan` rows or `nscan \\+ 1` rows"
  )
})

test_that("simulate_cdcm errors when nu and hypothesis indices do not match", {
  nu_bad <- list(
    nu_A = c(0.4, 0.3),
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

  U <- matrix(c(0,1,0,1,0,1,0,1,0,1), ncol = 1)

  expect_error(
    simulate_cdcm(
      nu = nu_bad,
      hypothesis_idxs = hypothesis_idxs,
      nscan = 10,
      m = 2,
      n_u = 1,
      TR = 2,
      U = U,
      SNR = 5,
      check_solution = FALSE
    ),
    "nu\\$nu_A"
  )
})

test_that("simulate_cdcm runs with check_solution enabled", {
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

  expect_no_warning(
    simulate_cdcm(
      nu = nu,
      hypothesis_idxs = hypothesis_idxs,
      nscan = 10,
      m = 2,
      n_u = 1,
      TR = 2,
      U = U,
      SNR = 5,
      check_solution = TRUE,
      seed = 1
    )
  )
})

test_that("simulation_format_check returns standardized inputs", {
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

  U <- matrix(c(0,1,0,1,0,1,0,1,0,1), ncol = 1)

  out <- simulation_format_check(
    nu = nu,
    hypothesis_idxs = hypothesis_idxs,
    nscan = 10,
    m = 2,
    n_u = 1,
    TR = 2,
    U = U
  )

  expect_equal(length(out$times), 11)
  expect_equal(dim(out$input), c(11, 1))
  expect_equal(out$input[1, 1], U[1, 1])
})

test_that("struct_paramMats places parameters correctly", {
  idxs <- list(
    A_idxs = matrix(c(2,1,
                      1,2,
                      1,1,
                      2,2), byrow = TRUE, ncol = 2),
    B_idxs = matrix(c(1,2,1), byrow = TRUE, ncol = 3),
    C_idxs = matrix(c(1,1), byrow = TRUE, ncol = 2)
  )

  nu <- list(
    nu_A = c(0.4, 0.3, -0.1, 0.15),
    nu_B = c(-0.2),
    nu_C = c(0.7)
  )

  out <- struct_paramMats(m = 2, n_u = 1, idxs = idxs, nu = nu, reparam = FALSE)

  expect_equal(out$A[2, 1], 0.4)
  expect_equal(out$A[1, 2], 0.3)
  expect_equal(out$A[1, 1], -0.1)
  expect_equal(out$A[2, 2], 0.15)

  expect_equal(out$B[[1]][2, 1], -0.2)
  expect_equal(out$C[1, 1], 0.7)
})

test_that("struct_paramMats reparameterizes diagonal of A when requested", {
  idxs <- list(
    A_idxs = matrix(c(1,1,
                      2,2), byrow = TRUE, ncol = 2),
    B_idxs = matrix(c(1,1,1), byrow = TRUE, ncol = 3),
    C_idxs = matrix(c(1,1), byrow = TRUE, ncol = 2)
  )

  nu <- list(
    nu_A = c(log(2), log(4)),
    nu_B = c(0),
    nu_C = c(0)
  )

  out <- struct_paramMats(m = 2, n_u = 1, idxs = idxs, nu = nu, reparam = TRUE)

  expect_equal(diag(out$A), c(-1, -2))
})

test_that("add_noise is reproducible and returns expected dimensions", {
  y_conv <- cbind(1:10, 11:20)
  times <- 0:10

  out1 <- add_noise(y_conv, times, SNR = 5, seed = 123)
  out2 <- add_noise(y_conv, times, SNR = 5, seed = 123)

  expect_equal(dim(out1), c(10, 2))
  expect_equal(out1, out2)
})

test_that("add_noise with rho = 0 reproduces independent Gaussian noise", {
  y_conv <- cbind(1:10, 11:20)
  times <- 0:10

  out <- add_noise(y_conv, times, SNR = 5, seed = 123, rho = 0)

  expected <- y_conv
  for (j in 1:2) {
    set.seed(123)
    expected[, j] <- y_conv[, j] +
      stats::rnorm(10, mean = 0, sd = stats::sd(y_conv[, j]) / 5)
  }

  expect_equal(out, expected)
  expect_equal(out, add_noise(y_conv, times, SNR = 5, seed = 123))
})

test_that("add_noise with rho > 0 gives AR(1) noise with the SNR-implied sd", {
  n <- 20000
  y_conv <- cbind(sin(seq_len(n) / 50), cos(seq_len(n) / 50))
  times <- 0:n

  out <- add_noise(y_conv, times, SNR = 2, seed = 1, rho = c(0.6, 0))
  e <- out - y_conv

  acf1 <- function(x) stats::acf(x, lag.max = 1, plot = FALSE)$acf[2]

  expect_lt(abs(acf1(e[, 1]) - 0.6), 0.03)
  expect_lt(abs(acf1(e[, 2])), 0.03)
  # Marginal sd is set by SNR regardless of rho
  expect_lt(abs(stats::sd(e[, 1]) / (stats::sd(y_conv[, 1]) / 2) - 1), 0.05)
})

test_that("simulate_cdcm errors for invalid rho", {
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

  U <- matrix(c(rep(0, 5), rep(1, 5)), ncol = 1)

  for (bad_rho in list(1, c(0.2, 0.3, 0.4), NA_real_)) {
    expect_error(
      simulate_cdcm(nu = nu, hypothesis_idxs = hypothesis_idxs, nscan = 10,
                    m = 2, n_u = 1, TR = 2, U = U, SNR = 5,
                    check_solution = FALSE, rho = bad_rho),
      "`rho` must be"
    )
  }
})
