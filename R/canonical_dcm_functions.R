# Function to get changepoints in u when u(t) is any dimension
get_changepoints <- function(input){

  if (!is.matrix(input)) {
    stop("`input` must be a matrix.", call. = FALSE)
  }

  if (nrow(input) < 1) {
    stop("`input` must have at least one row.", call. = FALSE)
  }

  if (!is.numeric(input)) {
    stop("`input` must be numeric.", call. = FALSE)
  }

  # Find indices where changes occur
  change_idx <- which(abs(diff(input)) >= 1, arr.ind = TRUE)

  # If there are no changes, return terminal row only
  if (length(change_idx) == 0) {
    return(nrow(input))
  }

  # Get dataframe of indices when changes occur
  changes <- as.data.frame(change_idx)
  changes$row <- changes$row+1 # add 1 to row index to indicate the row the change starts
  changes <- changes[order(changes$row),] # ascending order

  # Make duplicates their own category to reflect when switching from input1 -> input2
  changes <- changes[!duplicated(changes$row),] # remove duplicates

  # Add in last row to reflect final block of times
  changes <- rbind(changes, c(nrow(input),0))$row

  # Return dataframe
  return(changes)
}

# Function to get block lengths for A1-A3 identifiability check
get_block_lengths <- function(input) {
  if (!is.matrix(input)) {
    stop("`input` must be a matrix.", call. = FALSE)
  }

  if (nrow(input) < 1) {
    stop("`input` must have at least one row.", call. = FALSE)
  }

  changes <- get_changepoints(input)

  # No internal changepoints; only terminal row
  if (length(changes) == 1) {
    return(nrow(input))
  }

  internal_changes <- changes[-length(changes)]  # true block starts after row 1
  terminal_row <- changes[length(changes)]

  starts <- c(1, internal_changes)
  ends <- c(internal_changes - 1, terminal_row)

  lengths <- ends - starts + 1
  return(lengths)
}

#' Check whether a block design satisfies sufficient identifiability diagnostics
#'
#' Evaluates design sufficient conditions for identifiability under a
#' block-design experimental input. These checks are intended as practical
#' design diagnostics for CDCM analyses and related simulation studies.
#'
#' @param data A list containing:
#'   \itemize{
#'     \item \code{times}: a numeric vector of observation times
#'     \item \code{u}: a \eqn{T \times n_u} matrix encoding the experimental
#'     stimulus inputs
#'     \item \code{y_obs}: a \eqn{T \times m} matrix of observed ROI-level BOLD
#'     signals
#'   }
#' @param warn Logical; if \code{TRUE}, warnings are issued when one or more
#'   sufficient conditions are not satisfied.
#'
#' @details
#' This function evaluates a set of practical design diagnostics motivated by
#' sufficient conditions for identifiability in the CDCM framework under
#' block-design experiments. In particular, it checks:
#' \itemize{
#'   \item whether block lengths are long enough relative to the number of ROIs
#'   \item whether the total number of observations is large enough relative to
#'   the number of ROIs and stimulus inputs
#'   \item whether enough distinct blockwise input configurations are observed
#'   when the number of stimulus dimensions is 2 or 3
#' }
#'
#' These checks are based on sufficient conditions and should be interpreted as
#' practical diagnostics rather than necessary conditions. Failure to satisfy
#' one or more of them does not necessarily imply non-identifiability, but it
#' may indicate that the design is weak for reliable estimation.
#'
#' The number of ROIs is taken to be \code{m = ncol(y_obs)}, the number of
#' stimulus dimensions is \code{n_u = ncol(u)}, and the number of observations
#' is \code{T = length(times)}.
#'
#' @return A list containing:
#' \describe{
#'   \item{T}{The number of observations.}
#'   \item{m}{The number of ROIs.}
#'   \item{n_u}{The number of stimulus dimensions.}
#'   \item{changepoints}{The changepoint vector returned internally by
#'   \code{get_changepoints()}.}
#'   \item{block_lengths}{The block lengths computed from the design matrix.}
#'   \item{first_block_zero}{Logical indicating whether the first input block is
#'   the all-zero input configuration.}
#'   \item{n_distinct_block_inputs}{The number of distinct blockwise input
#'   configurations observed.}
#'   \item{cond_A1_first}{Logical indicating whether the first block satisfies
#'   the minimum block-length requirement.}
#'   \item{cond_A1_later}{Logical indicating whether all later blocks satisfy
#'   the minimum block-length requirement.}
#'   \item{cond_A1}{Logical indicating whether the block-length condition is
#'   satisfied overall.}
#'   \item{cond_A2}{Logical indicating whether the observation-count condition is
#'   satisfied.}
#'   \item{cond_A3}{Logical or \code{NA}. For designs with \code{n_u = 2} or
#'   \code{3}, indicates whether enough distinct blockwise input settings are
#'   observed. For other values of \code{n_u}, this is \code{NA}.}
#'   \item{all_ok}{Logical indicating whether all evaluated sufficient-condition
#'   diagnostics are satisfied.}
#' }
#'
#' @export
check_design_identifiability <- function(data, warn = TRUE) {

  if (!is.list(data)) {
    stop("`data` must be a list.", call. = FALSE)
  }

  required_data_names <- c("times", "u", "y_obs")
  missing_data_names <- setdiff(required_data_names, names(data))
  if (length(missing_data_names) > 0) {
    stop(
      "`data` is missing required component(s): ",
      paste(missing_data_names, collapse = ", "),
      call. = FALSE
    )
  }

  times <- data$times
  u <- data$u
  y_obs <- data$y_obs

  if (!is.numeric(times) || !is.vector(times)) {
    stop("`data$times` must be a numeric vector.", call. = FALSE)
  }

  if (!is.matrix(u)) {
    stop("`data$u` must be a matrix.", call. = FALSE)
  }

  if (!is.matrix(y_obs)) {
    stop("`data$y_obs` must be a matrix.", call. = FALSE)
  }

  T_obs <- length(times)
  m <- ncol(y_obs)   # number of ROIs
  n_u <- ncol(u)     # number of stimuli

  if (nrow(u) != T_obs) {
    stop("`data$u` must have the same number of rows as length(data$times).", call. = FALSE)
  }

  if (nrow(y_obs) != T_obs) {
    stop("`data$y_obs` must have the same number of rows as length(data$times).", call. = FALSE)
  }

  changes <- get_changepoints(u)
  block_lengths <- get_block_lengths(u)

  first_block_zero <- all(u[1, ] == 0)

  # A1
  first_block_min <- if (first_block_zero) (m + 1) else (m + 2)

  cond_A1_first <- block_lengths[1] >= first_block_min
  cond_A1_later <- if (length(block_lengths) > 1) {
    all(block_lengths[-1] >= (m + 1))
  } else {
    TRUE
  }
  cond_A1 <- cond_A1_first && cond_A1_later

  # A2
  cond_A2 <- T_obs > (m + 1) * (n_u + 1)

  # A3
  cond_A3 <- NA
  n_distinct_block_inputs <- NA_integer_

  if (length(changes) == 1) {
    block_starts <- 1
  } else {
    block_starts <- c(1, changes[-length(changes)])
  }

  block_inputs <- u[block_starts, , drop = FALSE]
  n_distinct_block_inputs <- nrow(unique(block_inputs))

  if (n_u %in% c(2, 3)) {
    cond_A3 <- n_distinct_block_inputs >= (n_u + 1)
  }

  all_ok <- cond_A1 && cond_A2 && (is.na(cond_A3) || cond_A3)

  if (warn) {
    if (!cond_A1) {
      warning(
        "Block-length condition failed: the design may not satisfy the minimum block-length requirements for identifiability.",
        call. = FALSE
      )
    }
    if (!cond_A2) {
      warning(
        "Observation-count condition failed: `length(times)` is not greater than `(m + 1) * (n_u + 1)`.",
        call. = FALSE
      )
    }
    if (!is.na(cond_A3) && !cond_A3) {
      warning(
        "Distinct-input condition failed: fewer than `n_u + 1` distinct blockwise input settings were observed.",
        call. = FALSE
      )
    }
  }

  return(list(
    T = T_obs,
    m = m,
    n_u = n_u,
    changepoints = changes,
    block_lengths = block_lengths,
    first_block_zero = first_block_zero,
    n_distinct_block_inputs = n_distinct_block_inputs,
    cond_A1_first = cond_A1_first,
    cond_A1_later = cond_A1_later,
    cond_A1 = cond_A1,
    cond_A2 = cond_A2,
    cond_A3 = cond_A3,
    all_ok = all_ok
  ))
}

# Function to extract draws from pathfinder for init
get_pathfinder_init <- function(draws) {
  last_draws <- draws %>%
    dplyr::group_by(.data$.chain) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::starts_with("."))

  # Remove non-parameter columns retained after dropping metadata
  param_cols <- setdiff(names(last_draws), c("lp__", "lp_approx__"))
  prefixes <- unique(sub("\\[.*", "", param_cols))

  grouped_list <- list()

  for (prefix in prefixes) {
    cols_to_extract <- grep(paste0("^", prefix, "\\["), param_cols, value = TRUE)
    sublist_df <- last_draws[, cols_to_extract, drop = FALSE]
    grouped_list[[prefix]] <- as.numeric(sublist_df)
  }

  return(grouped_list)
}

# Function to extract last draw for next init to keep chain going
get_last_draws_for_init <- function(draws) {
  last_draws <- draws %>%
    dplyr::group_by(.data$.chain) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::starts_with("."))

  # Remove non-parameter column retained after dropping metadata
  param_cols <- setdiff(names(last_draws), "lp__")
  prefixes <- unique(sub("\\[.*", "", param_cols))

  grouped_list <- list()

  for (prefix in prefixes) {
    cols_to_extract <- grep(paste0("^", prefix, "\\["), param_cols, value = TRUE)
    sublist_df <- last_draws[, cols_to_extract, drop = FALSE]
    grouped_list[[prefix]] <- as.numeric(sublist_df)
  }

  return(grouped_list)
}

#' Construct a Stan data list for the CDCM model
#'
#' Prepares the input data and hypothesis indexing structures required by the
#' CDCM Stan program.
#'
#' @param data A list containing the observed time series and experimental
#'   design. It must contain:
#'   \itemize{
#'     \item \code{times}: a numeric vector of observation times
#'     \item \code{u}: a \eqn{T \times n_u} matrix encoding the experimental
#'     stimulus inputs
#'     \item \code{y_obs}: a \eqn{T \times m} matrix of observed ROI-level BOLD
#'     signals
#'   }
#' @param hypothesis_idxs A list specifying the hypothesized nonzero entries of
#'   the connectivity matrices. It must contain:
#'   \itemize{
#'     \item \code{A_idxs}: a 2-column matrix of \eqn{(i,j)} indices for the
#'     baseline connectivity matrix \eqn{A}
#'     \item \code{B_idxs}: a 3-column matrix of \eqn{(k,i,j)} indices for the
#'     modulatory matrices \eqn{B_k}, where \eqn{k} indexes the stimulus
#'     dimension
#'     \item \code{C_idxs}: a 2-column matrix of \eqn{(i,j)} indices for the
#'     driving-input matrix \eqn{C}
#'   }
#' @param ode_solver_type Integer code specifying the ODE solver used in the
#'   Stan program (default = 1; piecewise analytic solution, faster than CKRK numeric integration = 0).
#' @param tol Numeric tolerance used for both relative and absolute ODE solver
#'   tolerances (default = 1e-5; used for CKRK numeric integration).
#' @param sigma_nu Prior scale for off-diagonal neural connectivity parameters (default = 1).
#' @param sigma_nu_self Prior scale for diagonal self-connectivity parameters (default = 0.125; tighter prior is needed).
#' @param rate_sigma Rate parameter for the prior on observation noise.
#' @param sigma_z0 Prior scale for the initial neural state.
#' @param conv Integer indicator controlling whether convolution with the HRF is
#'   performed in the Stan model (default = 1; leave as-is unless you want latent z(t) instead of y(t)).
#' @param max_num_steps Maximum number of ODE solver steps allowed (default = 10^6; CKRK numeric solver only).
#'
#' @details
#' This function validates the structure of the supplied data and hypothesis
#' index objects, computes changepoints in the experimental design, and returns
#' a named list in the format required by the CDCM Stan model.
#'
#' In addition, the function evaluates a set of design
#' diagnostics via \code{\link{check_design_identifiability}}. If those
#' sufficient conditions are not satisfied, a warning is issued. This warning is
#' intended as a practical diagnostic only: failure to satisfy these conditions
#' does not necessarily imply non-identifiability, but it may indicate that the
#' design is weak for reliable estimation.
#'
#' The returned Stan data list can be used for either single-chain or
#' multi-chain sampling. Its structure does not depend on the number of chains.
#'
#' @return A named list containing the data, dimensions, changepoint structure,
#'   hypothesis index matrices, prior hyperparameters, and solver settings
#'   required by the CDCM Stan program.
#'
#' @export
get_stan_dat <- function(data, hypothesis_idxs, ode_solver_type = 1, tol = 1e-5,
                         sigma_nu = 1, sigma_nu_self = 0.125, rate_sigma = 0.5,
                         sigma_z0 = 0.3, conv = 1, max_num_steps = 10^6){
  # Assumes input data object is list with following structure:
  #   times: T x 1 vector with the time index, separated by TR (seconds)
  #   u: T x n_u matrix encoding the timing of experimental stimuli
  #   y_obs: T x m matrix encoding the scaled BOLD signal for each ROI
  #   scale: real number for scale factor of BOLD signal

  # Assumes input hypothesis_idxs is list with following structure:
  #   A_idxs: 2-column matrix where each row contains the i,j coord for A matrix hypothesis
  #   B_idxs: 3-column matrix where each row contains the k,i,j coord for B matrix hypothesis (k is stimulus)
  #   C_idxs: 2-column matrix where each row contains the i,j coord for C matrix hypothesis

  # Checks for data structure
  if (!is.list(data)) {
    stop("`data` must be a list.", call. = FALSE)
  }

  required_data_names <- c("times", "u", "y_obs")
  missing_data_names <- setdiff(required_data_names, names(data))
  if (length(missing_data_names) > 0) {
    stop(
      "`data` is missing required component(s): ",
      paste(missing_data_names, collapse = ", "),
      call. = FALSE
    )
  }

  # Extract components for readability
  times <- data$times
  u <- data$u
  y_obs <- data$y_obs

  # Check types
  if (!is.numeric(times) || !is.vector(times)) {
    stop("`data$times` must be a numeric vector.", call. = FALSE)
  }

  if (!is.matrix(u)) {
    stop("`data$u` must be a matrix.", call. = FALSE)
  }

  if (!is.matrix(y_obs)) {
    stop("`data$y_obs` must be a matrix.", call. = FALSE)
  }

  # Check matching number of rows
  n_time <- length(times)

  if (nrow(u) != n_time) {
    stop("`data$u` must have the same number of rows as length(data$times).", call. = FALSE)
  }

  if (nrow(y_obs) != n_time) {
    stop("`data$y_obs` must have the same number of rows as length(data$times).", call. = FALSE)
  }

  # Check times are positive
  if (any(times <= 0)) {
    stop("`data$times` must contain only positive values.", call. = FALSE)
  }

  # Check equally spaced
  if (length(times) > 1) {
    diffs <- diff(times)

    if (any(diffs <= 0)) {
      stop("`data$times` must be strictly increasing.", call. = FALSE)
    }

    if (any(abs(diffs - diffs[1]) > 1e-8)) {
      stop("`data$times` must be equally spaced.", call. = FALSE)
    }
  }

  # Format checks for hypothesis_idxs
  if (!is.list(hypothesis_idxs)) {
    stop("`hypothesis_idxs` must be a list.", call. = FALSE)
  }

  required_idx_names <- c("A_idxs", "B_idxs", "C_idxs")
  missing_idx_names <- setdiff(required_idx_names, names(hypothesis_idxs))
  if (length(missing_idx_names) > 0) {
    stop(
      "`hypothesis_idxs` is missing required component(s): ",
      paste(missing_idx_names, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.matrix(hypothesis_idxs$A_idxs)) {
    stop("`hypothesis_idxs$A_idxs` must be a matrix.", call. = FALSE)
  }
  if (!is.matrix(hypothesis_idxs$B_idxs)) {
    stop("`hypothesis_idxs$B_idxs` must be a matrix.", call. = FALSE)
  }
  if (!is.matrix(hypothesis_idxs$C_idxs)) {
    stop("`hypothesis_idxs$C_idxs` must be a matrix.", call. = FALSE)
  }

  if (ncol(hypothesis_idxs$A_idxs) != 2) {
    stop("`hypothesis_idxs$A_idxs` must have exactly 2 columns.", call. = FALSE)
  }
  if (ncol(hypothesis_idxs$B_idxs) != 3) {
    stop("`hypothesis_idxs$B_idxs` must have exactly 3 columns.", call. = FALSE)
  }
  if (ncol(hypothesis_idxs$C_idxs) != 2) {
    stop("`hypothesis_idxs$C_idxs` must have exactly 2 columns.", call. = FALSE)
  }

  # Require at least one row in each matrix (one par to be est in each of A, B, C)
  if (nrow(hypothesis_idxs$A_idxs) < 1) {
    stop("`hypothesis_idxs$A_idxs` must have at least one row.", call. = FALSE)
  }
  if (nrow(hypothesis_idxs$B_idxs) < 1) {
    stop("`hypothesis_idxs$B_idxs` must have at least one row.", call. = FALSE)
  }
  if (nrow(hypothesis_idxs$C_idxs) < 1) {
    stop("`hypothesis_idxs$C_idxs` must have at least one row.", call. = FALSE)
  }

  # Check indices are positive integers
  idx_mats <- list(
    A_idxs = hypothesis_idxs$A_idxs,
    B_idxs = hypothesis_idxs$B_idxs,
    C_idxs = hypothesis_idxs$C_idxs
  )

  for (nm in names(idx_mats)) {
    mat <- idx_mats[[nm]]
    if (!is.numeric(mat)) {
      stop("`hypothesis_idxs$", nm, "` must be numeric.", call. = FALSE)
    }
    if (any(mat <= 0) || any(mat %% 1 != 0)) {
      stop("`hypothesis_idxs$", nm, "` must contain positive integer indices.", call. = FALSE)
    }
  }

  # Light design identifiability check - still run just print warning
  design_check <- check_design_identifiability(data, warn = FALSE)

  if (!design_check$all_ok) {
    warning(
      "Experimental design may not satisfy sufficient identifiability conditions. ",
      "Run `check_design_identifiability()` for details.",
      call. = FALSE
    )
  }

  # Obtain dataframe of changes based on input
  input <- rbind(data$u[1,], data$u)
  changes <- get_changepoints(input)

  # Getting data in proper list structure for stan program
  stan_dat = with(data,
                  list(T = length(times),
                       m = ncol(y_obs),
                       n_u = ncol(u),
                       n_changes = length(changes),
                       change_pts = changes,
                       d_A = nrow(hypothesis_idxs$A_idxs),
                       d_B = nrow(hypothesis_idxs$B_idxs),
                       d_C = nrow(hypothesis_idxs$C_idxs),
                       sigma_nu = sigma_nu,
                       sigma_nu_self = sigma_nu_self,
                       rate_sigma = rate_sigma,
                       sigma_z0 = sigma_z0,
                       A_idxs = hypothesis_idxs$A_idxs,
                       B_idxs = hypothesis_idxs$B_idxs,
                       C_idxs = hypothesis_idxs$C_idxs,
                       tp= times,
                       u = u,
                       conv = conv,
                       y_obs = y_obs,
                       rel_tol = tol,
                       abs_tol = tol,
                       max_num_steps = max_num_steps,
                       ode_solver_type = ode_solver_type))

  return(stan_dat)
}


get_num_param <- function(data){

  required <- c("d_A", "d_B", "d_C", "m")
  missing <- setdiff(required, names(data))

  if (length(missing) > 0) {
    stop(
      "`data` is missing required component(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  num_nu <- data$d_A + data$d_B + data$d_C
  num_sigma_z0_beta <- 3*data$m
  num_param <- num_nu + num_sigma_z0_beta

  return(num_param)
}


#' Compute a minimum ESS threshold for CDCM sampling
#'
#' Computes the minimum effective sample size (ESS) threshold corresponding to a
#' specified confidence level and relative precision, using the number of model
#' parameters implied by a CDCM Stan data list.
#'
#' @param data A named list of data in the format required by the Stan model,
#'   typically produced by \code{\link{get_stan_dat}}.
#' @param alpha Significance level used in the \code{mcmcse::minESS()}
#'   calculation. Defaults to \code{0.05}.
#' @param eps Relative precision target used in the
#'   \code{mcmcse::minESS()} calculation. Defaults to \code{0.05}.
#'
#' @details
#' This function first determines the number of model parameters implied by the
#' supplied Stan data list, then uses \code{mcmcse::minESS()} to compute the
#' minimum effective sample size required to achieve the specified level
#' \code{alpha} and relative precision \code{eps}.
#'
#' This can be used to define the \code{ess_check} argument in
#' \code{\link{dcm_sample}}.
#'
#' @return A numeric scalar giving the minimum ESS threshold.
#'
#' @examples
#' \dontrun{
#' stan_dat <- get_stan_dat(data, hypothesis_idxs)
#' ess_target <- minESS_criterion(stan_dat)
#' }
#'
#' @export
minESS_criterion <- function(data, alpha = 0.05, eps = 0.05){
  # Determine number of model parameters
  num_param <- get_num_param(data)

  # Use number of parameters to determine minESS from mcmcse for specified level alpha and relative precision epsilon
  ess_check <- mcmcse::minESS(num_param, alpha = alpha, eps = eps)

  # Verify numeric
  ess_check <- as.numeric(ess_check)
  return(ess_check)
}

#' Generate initial values for CDCM posterior sampling
#'
#' Produces a list of initial values for the CDCM Stan model, optionally using
#' CmdStanR Pathfinder to obtain data-informed initializations.
#'
#' @param mod A compiled CmdStanR model object, typically returned by
#'   \code{\link{compile_cdcm}}.
#' @param data A named list of data in the format required by the Stan model,
#'   typically produced by \code{\link{get_stan_dat}}.
#' @param pathfinder_init Optional named list of initial values to seed the
#'   Pathfinder run. If \code{NULL}, default values are constructed internally.
#'
#' @details
#' If \code{pathfinder_init} is not supplied, default initial values are used:
#' observation noise parameters are initialized at 1, neural connectivity
#' parameters are initialized at 0, initial neural states are initialized at
#' 0.1, and hemodynamic scaling parameters are initialized at 0.
#'
#' The function then attempts to run CmdStanR Pathfinder and extract the final
#' Pathfinder draw in a format suitable for posterior sampling. If Pathfinder
#' fails for any reason, a backup initialization list based on the default
#' values is returned instead, and a message is issued.
#'
#' The object returned by this function is naturally a **single-chain**
#' initialization list. For multi-chain sampling via \code{mod$sample()}, users
#' may either:
#' \itemize{
#'   \item run \code{get_initial_vals()} separately to generate one
#'   initialization list per chain, or
#'   \item reuse/replicate a single initialization list across chains if that is
#'   appropriate for the analysis.
#' }
#'
#' @return A list of initial values suitable for use in Stan sampling. When
#'   Pathfinder succeeds, this is based on the final Pathfinder draw. When
#'   Pathfinder fails, a backup initialization list based on default values is
#'   returned instead.
#'
#' @export
get_initial_vals <- function(mod, data, pathfinder_init = NULL){

  required <- c("m", "d_A", "d_B", "d_C")
  missing <- setdiff(required, names(data))

  if (length(missing) > 0) {
    stop(
      "`data` is missing required component(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(pathfinder_init)) {pathfinder_init <- list(sigma = rep(1,data$m),
                                                         nu_A = rep(0,data$d_A),
                                                         nu_B = rep(0,data$d_B),
                                                         nu_C = rep(0,data$d_C),
                                                         z0 = rep(0.1,data$m),
                                                         beta = rep(0,data$m))}

  # tryCatch loop for initial value specification
  inits_list <- tryCatch(
    { # Try pathfinder for initial values
      pathfinder <- mod$pathfinder(data = data,
                                         seed = 1234,
                                         init = list(pathfinder_init),
                                         num_paths = 1)
      draws <- posterior::as_draws_df(pathfinder$draws())
      result <- get_pathfinder_init(draws)
      result
    },
    error = function(e) {
      # Initialize with 0 except for sigma and z0 if pathfinder fails for any reason
      message("Pathfinder failed: ", e$message)

      result_backup <- list(list(sigma = rep(1,data$m),
                                 nu_A = rep(0,data$d_A),
                                 nu_B = rep(0,data$d_B),
                                 nu_C = rep(0,data$d_C),
                                 z0 = rep(0.1,data$m),
                                 beta = rep(0,data$m)))
      result_backup
    }
  )

  return(inits_list)
}

#' Sample from the CDCM posterior in checkpointed chunks
#'
#' Runs posterior sampling for a compiled CDCM model in sequential chunks,
#' saving intermediate results to disk after each chunk and stopping when a
#' target multivariate effective sample size (ESS) is reached or a maximum
#' iteration limit is exceeded.
#'
#' This function is designed for checkpointed **single-chain** sampling and is
#' especially useful for long runs on HPC systems, where jobs may need to stop
#' at walltime limits. For standard multi-chain sampling, use the compiled
#' CmdStanR model directly via \code{mod$sample()}.
#'
#' @param mod A compiled CmdStanR model object, typically returned by
#'   \code{\link{compile_cdcm}}.
#' @param data A named list of data in the format required by the Stan model,
#'   typically produced by \code{\link{get_stan_dat}}.
#' @param inits_list A named list of initial values for model parameters,
#'   typically produced by \code{\link{get_initial_vals}}.
#' @param output_dir Character string giving the directory where checkpointed
#'   posterior draws and diagnostics should be saved.
#' @param basename Character string used as the base filename for saved output
#'   files.
#' @param ess_check Numeric scalar giving the target multivariate ESS used for
#'   convergence assessment.
#' @param refresh Integer giving the sampling progress refresh rate passed to
#'   CmdStanR.
#' @param warmup_iter Integer giving the number of warmup iterations for the
#'   initial sampling phase.
#' @param n_iter_chunk Integer giving the number of post-warmup sampling
#'   iterations to run per chunk.
#' @param max_iter Integer giving the maximum total number of post-warmup draws
#'   allowed across all chunks.
#' @param adapt_delta Numeric target acceptance rate passed to CmdStanR.
#' @param seed Integer random seed for reproducibility.
#'
#' @details
#' Sampling begins with an initial warmup phase, followed by posterior sampling
#' in chunks. After the first chunk, the last draw is used as the initial value
#' for the next chunk, and the adapted step size and inverse metric from the
#' previous fit are reused so that sampling can continue efficiently without
#' repeating warmup.
#'
#' After each completed chunk, posterior draws and NUTS diagnostics are saved to
#' disk in \code{output_dir} using filenames based on \code{basename}. This
#' checkpointing allows intermediate results to be recovered if sampling is
#' interrupted, for example due to HPC walltime limits.
#'
#' Convergence is assessed using the multivariate effective sample size (ESS),
#' computed via \code{mcmcse}, on the accumulated single-chain posterior draws.
#' Sampling continues until the target ESS specified by \code{ess_check} is
#' reached or the maximum number of iterations is exceeded.
#'
#' This function is intentionally restricted to **single-chain** sampling.
#' Its checkpointing logic, reuse of the last draw as initialization, reuse of
#' the adapted inverse metric and step size, and ESS calculation are all based
#' on a single continuing chain. If multiple chains are desired, users should
#' call the underlying CmdStanR sampler directly via \code{mod$sample()} (or the
#' compiled model object's \code{$sample()} method), specifying their preferred
#' number of chains and parallelization settings. In that case, the checkpointed
#' chunkwise workflow implemented here is not used.
#'
#' @return Invisibly returns a list with elements:
#' \describe{
#'   \item{draws}{A data frame containing the combined posterior draws across
#'   all completed sampling chunks.}
#'   \item{diagnostics}{A list of NUTS diagnostics for each completed sampling
#'   chunk.}
#'   \item{converged}{Logical indicating whether the ESS-based convergence
#'   criterion was satisfied.}
#'   \item{total_draws}{Total number of post-warmup draws generated across all
#'   completed chunks.}
#' }
#'
#' If sampling is interrupted before the function completes, checkpointed draws
#' and diagnostics saved to disk remain available, but no final return object is
#' produced.
#'
#' @section Multi-chain sampling:
#' For standard multi-chain posterior sampling, call the compiled CmdStanR model
#' directly rather than using \code{dcm_sample()}. For example:
#'
#' \preformatted{
#' mod <- compile_cdcm()
#' stan_dat <- get_stan_dat(data, hypothesis_idxs)
#' init_list <- get_initial_vals(mod, stan_dat)
#'
#' fit <- mod$sample(
#'   data = stan_dat,
#'   init = init_list,
#'   chains = 4,
#'   parallel_chains = 4,
#'   iter_warmup = 1000,
#'   iter_sampling = 1000,
#'   seed = 1234
#' )
#' }
#'
#' This approach provides standard multi-chain CmdStanR sampling, but does not
#' include the checkpointing or chunkwise ESS-based stopping rule implemented in
#' \code{dcm_sample()}.
#'
#' @export
dcm_sample <- function(mod, data, inits_list, output_dir, basename, ess_check,
                       refresh = 100, warmup_iter = 5000, n_iter_chunk = 1000,
                       max_iter = 10^6, adapt_delta = 0.9, seed = 1234){

  if (!dir.exists(output_dir)) {
    stop("`output_dir` must be an existing directory.", call. = FALSE)
  }

  if (!is.character(basename) || length(basename) != 1 || basename == "") {
    stop("`basename` must be a non-empty character string.", call. = FALSE)
  }

  # Make sure ess_check is not empty
  if (missing(ess_check)) {
    stop("`ess_check` must be provided. It should reflect the required multivariate ESS for convergence.", call. = FALSE)
  }

  # Running stan program to sample from posterior - initial warmup and starting to sample
  fit = mod$sample(
    data = data,
    init = list(inits_list), # initialize from pathfinder values
    refresh = refresh, # output frequency
    iter_warmup = warmup_iter, # warm-up iterations
    iter_sampling = n_iter_chunk, # sampling iterations
    seed = seed, # seed for reproducibility
    chains = 1,
    adapt_delta = adapt_delta,
    save_warmup = TRUE,
    output_dir = output_dir,
    output_basename = paste0(basename,"_out"),
    metric = "dense_e")

  total_draws <- n_iter_chunk
  chunk <- 1
  converged <- FALSE

  # Get inits to start sampling by checkpoints
  draws <- posterior::as_draws_df(fit$draws())
  last_init <- get_last_draws_for_init(draws)

  # List to store draws and list to store diagnostics
  all_draws <- list()
  all_diagnostics <- list()

  # Add draws to list
  all_draws[[1]] <- draws

  # Extract diagnostics information from fit object
  all_diagnostics[[1]] <- bayesplot::nuts_params(fit)

  # Inv_metric and step_size from warm-up
  inv_metric <- fit$inv_metric(matrix = F)$`1`
  step_size <- fit$metadata()$step_size_adaptation

  ############ Run checkpoints of iterations ####
  while (!converged && total_draws < max_iter) {
    message("Sampling chunk ", chunk)

    fit <- mod$sample(
      data = data,
      init = list(last_init),
      refresh = refresh, # output frequency
      iter_warmup = 0, # no warmup, only onto sampling
      iter_sampling = n_iter_chunk, # sampling iterations
      seed = seed, # seed for reproducibility
      chains = 1,
      adapt_delta = adapt_delta,
      adapt_engaged = FALSE,
      metric = "dense_e",
      step_size = step_size,
      inv_metric = inv_metric)

    draws <- posterior::as_draws_df(fit$draws())
    all_draws[[chunk+1]] <- draws

    # Combine all draws so far to assess convergence - save cumulative draws at each checkpoint
    draws_df <- dplyr::bind_rows(all_draws)
    draws_df$.iteration <- c(1:nrow(draws_df))
    draws_df$.draw <- c(1:nrow(draws_df))
    save(draws_df, file = paste0(output_dir,"/",basename,"_draws.RData"))

    # Combine all diagnostics information to use for plots later
    all_diagnostics[[chunk+1]] <- bayesplot::nuts_params(fit)
    save(all_diagnostics, file = paste0(output_dir,"/",basename,"_diagnostics.RData"))

    # convergence check: multivariate ESS and asymptotic covariance from momentLS package
    param_draws <- suppressWarnings(as.matrix(draws_df[,2:(ncol(draws_df)-3)]))
    avar <- momentLS::mtvMLSE(param_draws)$cov
    multi_ess <- mcmcse::multiESS(param_draws, covmat = avar)
    ess_ok <- !is.na(multi_ess) && multi_ess > ess_check

    message("Multivariate ESS = ", round(multi_ess, digits = 2))

    if (ess_ok) {
      message("Converged after ", total_draws + n_iter_chunk, " iterations.")
      converged <- TRUE
      break
    }

    # extract last draws to use as init for next chunk, plus inv_metric and step size from warmup
    last_init <- get_last_draws_for_init(draws)
    inv_metric <- fit$inv_metric(matrix = F)$`1`
    step_size <- fit$metadata()$step_size_adaptation

    total_draws <- total_draws + n_iter_chunk
    chunk <- chunk + 1
  }

  # Return list of draws, diagnostics, if converged, and total draws (post warm-up)
  # Draws and diagnostics are also checkpointed along the way and saved to drive
  return(invisible(list(
    draws = draws_df,
    diagnostics = all_diagnostics,
    converged = converged,
    total_draws = total_draws + n_iter_chunk
  )))
}











