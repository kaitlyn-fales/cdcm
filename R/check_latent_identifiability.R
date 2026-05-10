#' Check parameter and latent-state identifiability diagnostics
#'
#' Evaluates simulation- or plug-in-based diagnostics for the CDCM
#' identifiability assumptions involving the effective connectivity matrices
#' and the latent neural states.
#'
#' @param A A \eqn{d \times d} numeric matrix giving the baseline effective
#'   connectivity matrix.
#' @param B A list of \eqn{d \times d} numeric matrices giving the modulatory
#'   connectivity matrices. The length of \code{B} should match the number of
#'   stimulus inputs.
#' @param u A \eqn{T \times m} numeric matrix encoding the stimulus inputs at
#'   each observation time.
#' @param z A \eqn{T \times d} numeric matrix of latent neural states evaluated
#'   at the observation times. Rows correspond to time points and columns
#'   correspond to ROIs.
#' @param block_indices Optional integer vector giving the block indices to
#'   check. If \code{NULL}, all eligible blocks with at least \eqn{d + 1}
#'   observations are checked.
#' @param tol Numeric tolerance used for eigenvalue, rank, and determinant
#'   checks. Defaults to \code{1e-8}.
#' @param warn Logical; if \code{TRUE}, warnings are issued when diagnostics
#'   fail. Defaults to \code{TRUE}.
#'
#' @details
#' This function checks the CDCM assumptions involving the true parameters and
#' latent discrete observations. Unlike the design assumptions, these conditions
#' generally cannot be verified directly from real data because the true
#' parameters and latent states are unknown. Therefore, this function is intended
#' primarily for simulation studies, where \code{A}, \code{B}, and \code{z} are
#' known, or as a post-fit plug-in diagnostic using posterior summaries.
#'
#' For each checked block \eqn{b}, the function forms
#' \deqn{
#'   \tilde A^{(b)} = A + \sum_{i = 1}^m \tilde u_i^{(b)} B_i,
#' }
#' where \eqn{\tilde u^{(b)}} is the block-level stimulus vector, obtained from
#' the first row of \code{u} within block \eqn{b}. Assumption A3 is checked by
#' verifying that \eqn{\tilde A^{(b)}} has real and distinct eigenvalues.
#'
#' Assumption A4 is checked by constructing
#' \deqn{
#'   X^{(b)} = [x_0^{(b)}, \ldots, x_d^{(b)}],
#'   \qquad
#'   x_j^{(b)} = [z(t^{(b)} + t_j)^T, 1]^T,
#' }
#' from the first \eqn{d + 1} latent observations in block \eqn{b}, and checking
#' whether \eqn{X^{(b)}} has full rank.
#'
#' When used with posterior means from real data, these checks should be
#' interpreted only as plug-in diagnostics conditional on the quality and
#' convergence of the posterior fit.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{cond_A3}}{Logical indicating whether A3 is satisfied for all
#'   checked blocks.}
#'   \item{\code{cond_A4}}{Logical indicating whether A4 is satisfied for all
#'   checked blocks.}
#'   \item{\code{all_ok}}{Logical indicating whether both diagnostics are
#'   satisfied.}
#'   \item{\code{block_results}}{A data frame containing block-level diagnostic
#'   results.}
#' }
#'
#' @export
check_latent_identifiability <- function(A, B, u, z,
                                         block_indices = NULL,
                                         tol = 1e-8,
                                         warn = TRUE) {

  if (!is.matrix(A) || !is.numeric(A)) {
    stop("`A` must be a numeric matrix.", call. = FALSE)
  }

  if (nrow(A) != ncol(A)) {
    stop("`A` must be square.", call. = FALSE)
  }

  if (!is.list(B)) {
    stop("`B` must be a list of numeric matrices.", call. = FALSE)
  }

  if (!is.matrix(u) || !is.numeric(u)) {
    stop("`u` must be a numeric matrix.", call. = FALSE)
  }

  if (!is.matrix(z) || !is.numeric(z)) {
    stop("`z` must be a numeric matrix.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("`tol` must be a positive numeric scalar.", call. = FALSE)
  }

  if (!is.logical(warn) || length(warn) != 1) {
    stop("`warn` must be TRUE or FALSE.", call. = FALSE)
  }

  d <- nrow(A)
  m <- ncol(u)

  if (ncol(A) != d) {
    stop("`A` must be square.", call. = FALSE)
  }

  if (length(B) != m) {
    stop("`B` must have length equal to `ncol(u)`.", call. = FALSE)
  }

  for (i in seq_along(B)) {
    if (!is.matrix(B[[i]]) || !is.numeric(B[[i]])) {
      stop("Each element of `B` must be a numeric matrix.", call. = FALSE)
    }

    if (!all(dim(B[[i]]) == c(d, d))) {
      stop("Each element of `B` must have the same dimensions as `A`.", call. = FALSE)
    }
  }

  if (nrow(z) != nrow(u)) {
    stop("`z` and `u` must have the same number of rows.", call. = FALSE)
  }

  if (ncol(z) != d) {
    stop("`z` must have `nrow(A)` columns.", call. = FALSE)
  }

  changes <- get_changepoints(u)
  block_lengths <- get_block_lengths(u)

  block_starts <- if (length(changes) == 1) {
    1
  } else {
    c(1, changes[-length(changes)])
  }

  n_blocks <- length(block_lengths)

  if (is.null(block_indices)) {
    block_indices <- which(block_lengths >= (d + 1))
  } else {
    if (!is.numeric(block_indices) || any(block_indices < 1) || any(block_indices > n_blocks)) {
      stop("`block_indices` must contain valid block indices.", call. = FALSE)
    }

    block_indices <- as.integer(block_indices)
  }

  if (length(block_indices) == 0) {
    stop("No checked blocks contain at least `d + 1` observations.", call. = FALSE)
  }

  block_results <- lapply(block_indices, function(b) {

    u_b <- u[block_starts[b], ]

    A_tilde <- A
    for (i in seq_len(m)) {
      A_tilde <- A_tilde + u_b[i] * B[[i]]
    }

    eig <- eigen(A_tilde, only.values = TRUE)$values

    eig_real <- all(abs(Im(eig)) <= tol)
    eig_values_real_part <- Re(eig)

    eig_diffs <- abs(
      outer(eig_values_real_part, eig_values_real_part, "-")
    )

    min_eig_diff <- if (length(eig_values_real_part) > 1) {
      min(eig_diffs[upper.tri(eig_diffs)])
    } else {
      Inf
    }

    eig_distinct <- min_eig_diff > tol
    cond_A3_block <- eig_real && eig_distinct

    z_block <- z[
      block_starts[b]:(block_starts[b] + d),
      ,
      drop = FALSE
    ]

    X_b <- rbind(t(z_block), rep(1, d + 1))
    rank_X_b <- qr(X_b, tol = tol)$rank
    det_X_b <- det(X_b)

    cond_A4_block <- rank_X_b == (d + 1)

    data.frame(
      block = b,
      block_start = block_starts[b],
      block_length = block_lengths[b],
      eig_real = eig_real,
      eig_distinct = eig_distinct,
      min_eig_diff = min_eig_diff,
      rank_X = rank_X_b,
      det_X = det_X_b,
      cond_A3 = cond_A3_block,
      cond_A4 = cond_A4_block
    )
  })

  block_results <- do.call(rbind, block_results)

  cond_A3 <- all(block_results$cond_A3)
  cond_A4 <- all(block_results$cond_A4)
  all_ok <- cond_A3 && cond_A4

  if (warn) {
    if (!cond_A3) {
      warning(
        "A3 failed: at least one checked block has an effective connectivity matrix without distinct real eigenvalues.",
        call. = FALSE
      )
    }

    if (!cond_A4) {
      warning(
        "A4 failed: at least one checked block has a non-invertible latent-state data matrix.",
        call. = FALSE
      )
    }
  }

  return(list(
    cond_A3 = cond_A3,
    cond_A4 = cond_A4,
    all_ok = all_ok,
    block_results = block_results
  ))
}
