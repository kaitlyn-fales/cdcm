#' Predict CDCM latent and BOLD signals from posterior draws
#'
#' @param draws_df Posterior draws data frame. The first column is assumed to be
#'   the log-posterior column, the last three columns are assumed to be draw
#'   metadata, and the parameter columns are assumed to be ordered as:
#'   noise standard deviations, `nu_A`, `nu_B`, `nu_C`, and `z0`.
#' @param hypothesis_idxs List of model index matrices containing `A_idxs`,
#'   `B_idxs`, and `C_idxs`.
#' @param nscan Number of scan time points.
#' @param m Number of brain regions/ROIs, denoted d in the manuscript.
#' @param n_u Number of stimulus inputs, denoted m in the manuscript.
#' @param TR Repetition time.
#' @param U Input stimulus matrix.
#' @param check_identifiability Logical. If TRUE, runs post-fit plug-in
#'   identifiability diagnostics.
#' @param identifiability_tol Numeric tolerance for identifiability diagnostics.
#'
#' @return A list containing predicted latent states and predicted BOLD signal,
#' and optional identifiability check.
#'
#' @export
predict_cdcm <- function(draws_df, hypothesis_idxs, nscan, m, n_u, TR, U,
                         check_identifiability = FALSE,
                         identifiability_tol = 1e-8) {

  if (!is.data.frame(draws_df)) {
    stop("`draws_df` must be a data frame.", call. = FALSE)
  }

  n_A <- nrow(hypothesis_idxs$A_idxs)
  n_B <- nrow(hypothesis_idxs$B_idxs)
  n_C <- nrow(hypothesis_idxs$C_idxs)

  n_noise <- m
  n_z0 <- m

  n_param <- n_noise + n_A + n_B + n_C + n_z0

  param_cols <- seq(2, ncol(draws_df) - 3)

  if (length(param_cols) != n_param) {
    stop(
      "The number of parameter columns in `draws_df` does not match the expected ",
      "number based on `m` and `hypothesis_idxs`.",
      call. = FALSE
    )
  }

  draws <- suppressWarnings(posterior::as_draws_df(
    draws_df[, param_cols, drop = FALSE]
  ))

  nu_A_cols <- seq(n_noise + 1, length.out = n_A)
  nu_B_cols <- seq(max(nu_A_cols) + 1, length.out = n_B)
  nu_C_cols <- seq(max(nu_B_cols) + 1, length.out = n_C)
  z0_cols <- seq(max(nu_C_cols) + 1, length.out = n_z0)

  self_loop_idx <- which(hypothesis_idxs$A_idxs[, 1] == hypothesis_idxs$A_idxs[, 2])

  if (length(self_loop_idx) > 0) {
    self_loop_cols <- nu_A_cols[self_loop_idx]
    draws[, self_loop_cols] <- suppressWarnings(-0.5 * exp(as.matrix(draws[, self_loop_cols, drop = FALSE])))
  }

  post_summary <- posterior::summarise_draws(draws)
  post_mean <- post_summary$mean

  nu_A <- post_mean[nu_A_cols]
  nu_B <- post_mean[nu_B_cols]
  nu_C <- post_mean[nu_C_cols]
  z0 <- post_mean[z0_cols]

  nu <- list(
    nu_A = nu_A,
    nu_B = nu_B,
    nu_C = nu_C
  )

  sim <- suppressMessages(
    simulate_cdcm(
      nu = nu,
      hypothesis_idxs = hypothesis_idxs,
      nscan = nscan,
      m = m,
      n_u = n_u,
      TR = TR,
      U = U,
      SNR = 1,
      reparam = FALSE,
      z0 = z0,
      check_solution = FALSE,
      seed = 1,
      check_identifiability = check_identifiability,
      identifiability_tol = identifiability_tol
    )
  )

  out <- list(
    times = sim$simulated_data$times,
    z_pred = sim$true_signals$z_signal,
    y_pred = sim$true_signals$y_signal,
    posterior_summary = post_summary
  )

  if (check_identifiability) {
    out$identifiability <- sim$identifiability
  }

  return(out)
}
