#' Simulate fMRI time series under CDCM
#'
#' Simulates neural state trajectories and observed BOLD signals under a
#' CDCM with user-specified connectivity parameters and experimental inputs.
#'
#' This function generates:
#' \itemize{
#'   \item latent neural states via an analytic solution to a bilinear state equation
#'   \item observed BOLD signals via convolution with a canonical hemodynamic response function (HRF)
#'   \item additive Gaussian noise controlled by a specified signal-to-noise ratio (SNR)
#' }
#'
#' The user supplies the connectivity structure through `nu` and `hypothesis_idxs`,
#' along with a stimulus design matrix `U`. Internally, the function augments the
#' input matrix to include the initial time point required for solving the system.
#'
#' @param nu A list containing true parameter values:
#' \describe{
#'   \item{nu_A}{numeric vector of A (baseline connectivity) parameters}
#'   \item{nu_B}{numeric vector of B (modulatory) parameters}
#'   \item{nu_C}{numeric vector of C (input) parameters}
#' }
#'
#' @param hypothesis_idxs A list specifying the locations of nonzero parameters:
#' \describe{
#'   \item{A_idxs}{matrix with 2 columns (row, column indices)}
#'   \item{B_idxs}{matrix with 3 columns (input, row, column indices)}
#'   \item{C_idxs}{matrix with 2 columns (row, input indices)}
#' }
#'
#' @param nscan Number of scan time points.
#' @param m Number of brain regions (state dimension).
#' @param n_u Number of experimental inputs.
#' @param TR Repetition time (seconds).
#' @param U Input (stimulus) matrix of dimension `nscan x n_u` or `(nscan + 1) x n_u`.
#'   Values must be binary (0/1).
#' @param SNR Numeric. Signal-to-noise ratio for additive Gaussian noise,
#'   defined as the ratio of the standard deviation of the signal to that of
#'   the noise (i.e., \eqn{\mathrm{SNR} = \mathrm{sd}(\text{signal}) / \mathrm{sd}(\text{noise})}).
#'   Noise is generated such that
#'   \eqn{\mathrm{sd}(\text{noise}) = \mathrm{sd}(\text{signal}) / \mathrm{SNR}},
#'   implying a noise variance of
#'   \eqn{\sigma^2 = \mathrm{Var}(\text{signal}) / \mathrm{SNR}^2}.
#'   Must be positive.
#'
#' @param reparam Logical. If `TRUE`, reparameterizes diagonal elements of the A matrix
#'   to enforce stability.
#' @param z0 Optional initial state vector. If `NULL`, a default is used.
#' @param check_solution Logical. If `TRUE`, compares analytic and numeric ODE solutions
#'   to assess numerical stability.
#' @param tol Numeric tolerance for the analytic vs numeric solution comparison.
#' @param seed Optional seed (scalar or length `m`) for reproducible noise generation.
#'
#' @return A list with two components:
#' \describe{
#'   \item{simulated_data}{
#'     \describe{
#'       \item{times}{Vector of scan times (length `nscan`)}
#'       \item{u}{Input matrix at scan times (`nscan x n_u`)}
#'       \item{y_obs}{Observed BOLD signals (`nscan x m`)}
#'     }
#'   }
#'   \item{true_signals}{
#'     \describe{
#'       \item{z_signal}{Latent neural states (`nscan x m`)}
#'       \item{y_signal}{Noise-free BOLD signals (`nscan x m`)}
#'     }
#'   }
#' }
#'
#' @details
#' The neural dynamics follow the bilinear state space DCM. Within each block of
#' constant input, the system admits a closed-form solution, which is used here
#' for efficient simulation. The observed BOLD signal is obtained by convolving
#' the neural states with the canonical HRF (linear combination of two gamma distributions).
#'
#' @examples
#' nu <- list(
#'   nu_A = c(0.4, 0.3, -0.1, 0.15),
#'   nu_B = c(-0.2),
#'   nu_C = c(0.7)
#' )
#'
#' hypothesis_idxs <- list(
#'   A_idxs = matrix(c(2,1,
#'                     1,2,
#'                     1,1,
#'                     2,2), byrow = TRUE, ncol = 2),
#'   B_idxs = matrix(c(1,2,1), byrow = TRUE, ncol = 3),
#'   C_idxs = matrix(c(1,1), byrow = TRUE, ncol = 2)
#' )
#'
#' U <- matrix(c(0,rep(c(rep(0,5),rep(1,5)), 2)), ncol = 1)
#'
#' sim <- simulate_cdcm(
#'   nu = nu,
#'   hypothesis_idxs = hypothesis_idxs,
#'   nscan = 20,
#'   m = 2,
#'   n_u = 1,
#'   TR = 2,
#'   U = U,
#'   SNR = 5,
#'   seed = 1
#' )
#'
#' @export
simulate_cdcm <- function(nu, hypothesis_idxs, nscan, m, n_u, TR, U, SNR,
                          reparam = TRUE, z0 = NULL, check_solution = TRUE, tol = 1e-4,
                          seed = NULL) {

  # Check inputs and modify as needed for simulation (creation of times, checking U, etc.)
  checked <- simulation_format_check(nu, hypothesis_idxs, nscan, m, n_u, TR, U)

  nu <- checked$nu
  hypothesis_idxs <- checked$hypothesis_idxs
  times <- checked$times
  input <- checked$input

  if (length(SNR) != 1 || !is.numeric(SNR) || is.na(SNR) || SNR <= 0) {
    stop("`SNR` must be a single positive number.", call. = FALSE)
  }

  ### Now that all inputs are checked, proceed with simulating data ###

  # Structure parameters, reparametrize diag(A)
  paramMats <- struct_paramMats(m = m,
                                n_u = n_u,
                                idxs = hypothesis_idxs,
                                nu = nu,
                                reparam = reparam)

  # Running functions to get analytic solution
  # Obtain dataframe of changes based on input
  changes <- get_changepoints(input)

  # Running final function
  analytic_sol <- get_AnalyticSol(input, changes, times, parms = paramMats, m = m, v0 = z0)

  # Optional check between analytic and numeric solution for stability
  if (check_solution) {

    input_u = function(t){
      sapply(1:ncol(input), function(j) stats::approxfun(x = times,y = input[,j], rule=2, method="const",f=0)(t))
    }

    # Quickly check with forward solving ODE numerically
    numeric_sol <- deSolve::ode(y = as.numeric(analytic_sol[1, ]),
                                times = times,
                                func = linear,
                                parms= with(paramMats, list(A=A, B=B, C=C)),
                                input = input_u,
                                atol = 1e-6,
                                rtol = 1e-6
    )

    num_sol_aligned <- numeric_sol[, -1, drop = FALSE]
    max_diff <- max(abs(analytic_sol - num_sol_aligned))

    if (max_diff > tol) {
      warning(
        "Analytic and numeric solutions differ beyond tolerance (max diff = ",
        signif(max_diff, 3),
        "). This may indicate numerical instability or unsuitable parameter values.",
        call. = FALSE
      )
    }
  }

  # Convolve analytic solution with canonical HRF
  y_conv <- vapply(
    X = seq_len(m),
    FUN = function(i) HRF_mu(analytic_sol[-1, i], times[-1]),
    FUN.VALUE = numeric(length(times) - 1)
  )

  # Add Gaussian measurement noise according to specified SNR
  y_obs <- add_noise(y_conv, times, SNR, seed = seed)

  # Make simulated export objects
  simulated_data <- list(times = times[-1], u = input[-1, , drop = FALSE], y_obs = y_obs)
  true_signals <- list(z_signal = analytic_sol[-1,], y_signal = y_conv)

  # Return list of both objects
  return(list(simulated_data = simulated_data, true_signals = true_signals))
}

