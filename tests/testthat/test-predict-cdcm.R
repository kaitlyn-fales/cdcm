test_that("predict_cdcm returns expected structure and dimensions", {
  hypothesis_idxs <- list(
    A_idxs = matrix(c(
      1, 1,
      2, 2,
      2, 1
    ), ncol = 2, byrow = TRUE),
    B_idxs = matrix(c(
      1, 2, 1
    ), ncol = 3, byrow = TRUE),
    C_idxs = matrix(c(
      1, 1
    ), ncol = 2, byrow = TRUE)
  )

  nscan <- 12
  m <- 2
  n_u <- 1
  TR <- 2

  U <- matrix(
    c(rep(0, 6), rep(1, 6)),
    ncol = 1
  )

  n_draws <- 20

  draws_df <- data.frame(
    lp__ = rnorm(n_draws),

    # noise sd, length m
    sigma.1 = rep(0.1, n_draws),
    sigma.2 = rep(0.1, n_draws),

    # nu_A, length 3
    nu_A.1 = rep(log(1), n_draws),   # transformed to -0.5
    nu_A.2 = rep(log(2), n_draws),   # transformed to -1.0
    nu_A.3 = rep(0.2, n_draws),

    # nu_B, length 1
    nu_B.1 = rep(0.1, n_draws),

    # nu_C, length 1
    nu_C.1 = rep(0.5, n_draws),

    # z0, length m
    z0.1 = rep(0, n_draws),
    z0.2 = rep(0, n_draws),

    # beta, length m
    beta.1 = rep(1, n_draws),
    beta.2 = rep(2, n_draws),

    .chain = rep(1, n_draws),
    .iteration = seq_len(n_draws),
    .draw = seq_len(n_draws)
  )

  out <- predict_cdcm(
    draws_df = posterior::as_draws_df(draws_df),
    hypothesis_idxs = hypothesis_idxs,
    nscan = nscan,
    m = m,
    n_u = n_u,
    TR = TR,
    U = U,
    check_identifiability = FALSE
  )

  expect_true(is.list(out))
  expect_named(
    out,
    c("times","z_pred", "y_pred", "posterior_summary")
  )

  expect_equal(length(out$times), nscan)
  expect_equal(dim(out$z_pred), c(nscan, m))
  expect_equal(dim(out$y_pred), c(nscan, m))

})

test_that("predict_cdcm returns identifiability diagnostics and detects A4 failure", {

  hypothesis_idxs <- list(
    A_idxs = matrix(c(
      1, 1,
      2, 2
    ), ncol = 2, byrow = TRUE),
    B_idxs = matrix(c(
      1, 2, 1
    ), ncol = 3, byrow = TRUE),
    C_idxs = matrix(c(
      1, 1
    ), ncol = 2, byrow = TRUE)
  )

  n_draws <- 10
  m <- 2
  n_u <- 1
  nscan <- 12

  draws_df <- data.frame(
    lp__ = rnorm(n_draws),

    sigma.1 = rep(0.1, n_draws),
    sigma.2 = rep(0.1, n_draws),

    nu_A.1 = rep(log(1), n_draws),
    nu_A.2 = rep(log(2), n_draws),

    nu_B.1 = rep(0.1, n_draws),

    nu_C.1 = rep(0.5, n_draws),

    # zero initial condition encourages A4 failure
    z0.1 = rep(0, n_draws),
    z0.2 = rep(0, n_draws),

    beta.1 = rep(0, n_draws),
    beta.2 = rep(0, n_draws),

    .chain = rep(1, n_draws),
    .iteration = seq_len(n_draws),
    .draw = seq_len(n_draws)
  )

  U <- matrix(
    c(rep(0, 3),
      rep(1, 3),
      rep(0, 3),
      rep(1, 3)),
    ncol = 1
  )

  expect_warning(
    out <- predict_cdcm(
      draws_df = posterior::as_draws_df(draws_df),
      hypothesis_idxs = hypothesis_idxs,
      nscan = nscan,
      m = m,
      n_u = n_u,
      TR = 2,
      U = U,
      check_identifiability = TRUE,
      identifiability_tol = 1e-8
    ),
    "A4 failed"
  )

  expect_true("identifiability" %in% names(out))

  expect_named(
    out$identifiability,
    c("design", "latent")
  )

  expect_true(out$identifiability$design$all_ok)

  expect_true(out$identifiability$latent$cond_A3)
  expect_false(out$identifiability$latent$cond_A4)
  expect_false(out$identifiability$latent$all_ok)
})

test_that("predict_cdcm errors when draws_df has wrong number of parameter columns", {
  hypothesis_idxs <- list(
    A_idxs = matrix(c(1, 1), ncol = 2),
    B_idxs = matrix(numeric(0), ncol = 3),
    C_idxs = matrix(c(1, 1), ncol = 2)
  )

  bad_draws_df <- data.frame(
    lp__ = rnorm(5),
    sigma.1 = rep(0.1, 5),
    .chain = rep(1, 5),
    .iteration = seq_len(5),
    .draw = seq_len(5)
  )

  expect_error(
    predict_cdcm(
      draws_df = posterior::as_draws_df(bad_draws_df),
      hypothesis_idxs = hypothesis_idxs,
      nscan = 8,
      m = 1,
      n_u = 1,
      TR = 2,
      U = matrix(rep(0, 8), ncol = 1)
    ),
    "does not match the expected"
  )
})

test_that("predict_cdcm returns identifiability diagnostics when all checks pass", {

  hypothesis_idxs <- list(
    A_idxs = matrix(c(
      1, 1,
      2, 2
    ), ncol = 2, byrow = TRUE),
    B_idxs = matrix(c(
      1, 2, 1
    ), ncol = 3, byrow = TRUE),
    C_idxs = matrix(c(
      1, 1,
      2, 1
    ), ncol = 2, byrow = TRUE)
  )

  n_draws <- 10
  m <- 2
  n_u <- 1
  nscan <- 12

  draws_df <- data.frame(
    lp__ = rnorm(n_draws),

    sigma.1 = rep(0.1, n_draws),
    sigma.2 = rep(0.1, n_draws),

    # diagonal A entries, transformed inside predict_cdcm
    nu_A.1 = rep(log(1), n_draws),  # becomes -0.5
    nu_A.2 = rep(log(2), n_draws),  # becomes -1.0

    # small modulatory effect
    nu_B.1 = rep(0.05, n_draws),

    # input enters both regions
    nu_C.1 = rep(0.5, n_draws),
    nu_C.2 = rep(0.3, n_draws),

    # nonzero initial condition to avoid trivial/collinear latent states
    z0.1 = rep(0.1, n_draws),
    z0.2 = rep(0.2, n_draws),

    beta.1 = rep(0, n_draws),
    beta.2 = rep(0, n_draws),

    .chain = rep(1, n_draws),
    .iteration = seq_len(n_draws),
    .draw = seq_len(n_draws)
  )

  U <- matrix(
    c(rep(0, 3),
      rep(1, 3),
      rep(0, 3),
      rep(1, 3)),
    ncol = 1
  )

  expect_no_warning(
    out <- predict_cdcm(
      draws_df = posterior::as_draws_df(draws_df),
      hypothesis_idxs = hypothesis_idxs,
      nscan = nscan,
      m = m,
      n_u = n_u,
      TR = 2,
      U = U,
      check_identifiability = TRUE,
      identifiability_tol = 1e-8
    )
  )

  expect_true("identifiability" %in% names(out))
  expect_named(out$identifiability, c("design", "latent"))

  expect_true(out$identifiability$design$cond_A1)
  expect_true(out$identifiability$design$cond_A2)
  expect_true(out$identifiability$design$all_ok)

  expect_true(out$identifiability$latent$cond_A3)
  expect_true(out$identifiability$latent$cond_A4)
  expect_true(out$identifiability$latent$all_ok)
})
