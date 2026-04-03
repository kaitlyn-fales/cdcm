# Function to get changepoints in u when u(t) is any dimension
get_changepoints <- function(input){

  # Get dataframe of indices when changes occur
  changes <- as.data.frame(which(abs(diff(input)) >= 1, arr.ind = T))
  changes$row <- changes$row+1 # add 1 to row index to indicate the row the change starts
  changes <- changes[order(changes$row),] # ascending order

  # Make duplicates their own category to reflect when switching from input1 -> input2
  changes <- changes[!duplicated(changes$row),] # remove duplicates

  # Add in last row to reflect final block of times
  changes <- rbind(changes, c(nrow(input),0))$row

  # Return dataframe
  return(changes)
}

# Function to extract draws from pathfinder for init
get_pathfinder_init <- function(draws) {
  last_draws <- draws %>%
    dplyr::group_by(.chain) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::starts_with("."))

  # Extracting prefixes
  col_names <- names(last_draws)[-c(1:2)]
  prefixes <- unique(sub("\\[.*", "", col_names))

  # Use split() to group the column names by their extracted prefix
  grouped_list <- list()

  # Extract columns
  for (prefix in prefixes) {
    # Find columns matching the current prefix
    cols_to_extract <- grep(paste0("^", prefix, "\\["), colnames(last_draws), value = TRUE)

    # Extract these columns and convert to a list
    sublist_df <- last_draws[, cols_to_extract]
    grouped_list[[prefix]] <- as.numeric(sublist_df)
  }

  return(grouped_list)
}

# Function to extract last draw for next init to keep chain going
get_last_draws_for_init <- function(draws) {
  last_draws <- draws %>%
    dplyr::group_by(.chain) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::starts_with("."))

  # Extracting prefixes
  col_names <- names(last_draws)[-1]
  prefixes <- unique(sub("\\[.*", "", col_names))

  # Use split() to group the column names by their extracted prefix
  grouped_list <- list()

  # Extract columns
  for (prefix in prefixes) {
    # Find columns matching the current prefix
    cols_to_extract <- grep(paste0("^", prefix, "\\["), colnames(last_draws), value = TRUE)

    # Extract these columns and convert to a list
    sublist_df <- last_draws[, cols_to_extract]
    grouped_list[[prefix]] <- as.numeric(sublist_df)
  }

  return(grouped_list)
}

#' Construct Stan data list for CDCM
#'
#' Prepares input data and model specification into the list format required
#' for the canonical DCM Stan program. This includes time series data, stimulus
#' inputs, connectivity structure indices, and hyperparameters.
#'
#' @param data A list containing:
#' \itemize{
#'   \item \code{times}: Numeric vector of time points (length T)
#'   \item \code{u}: Matrix of stimulus inputs (T x n_u)
#'   \item \code{y_obs}: Matrix of observed BOLD signals (T x m)
#' }
#' @param hypothesis_idxs A list specifying model structure:
#' \itemize{
#'   \item \code{A_idxs}: Matrix of indices for A matrix (connections)
#'   \item \code{B_idxs}: Matrix of indices for B matrix (modulatory effects)
#'   \item \code{C_idxs}: Matrix of indices for C matrix (inputs)
#' }
#' @param ode_solver_type Integer specifying ODE solver (default 1 is piecewise analytic solution)
#' @param tol Numeric tolerance for ODE solver (default 1e-5, if using numeric CKRK integration)
#' @param sigma_nu Prior scale for connectivity parameters
#' @param sigma_nu_self Prior scale for self-connections in A and B
#' @param rate_sigma Prior rate parameter for variance terms
#' @param sigma_z0 Prior scale for initial neural states
#' @param conv Integer flag for convolution option (default 1 convolves with canonical HRF - recommended to not change unless you want latent neural signal)
#' @param max_num_steps Maximum number of ODE solver steps (for CKRK solver)
#'
#' @return A named list suitable for input to a Stan model.
#' @export
get_stan_dat <- function(data, hypothesis_idxs, ode_solver_type = 1, tol = 10^{-5},
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

  # Obtain dataframe of changes based on input
  input <- rbind(0,data$u)
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

#' Compute number of model parameters
#'
#' Calculates the total number of parameters in the canonical DCM model,
#' based on the dimensions of connectivity components and variance terms.
#' This is primarily used to determine convergence thresholds (e.g., required ESS).
#'
#' @param data A list containing model dimensions:
#' \itemize{
#'   \item \code{d_A}: Number of A parameters
#'   \item \code{d_B}: Number of B parameters
#'   \item \code{d_C}: Number of C parameters
#'   \item \code{m}: Number of regions (ROIs)
#' }
#'
#' @return Integer representing total number of parameters.
#' @export
get_num_param <- function(data){
  num_nu <- data$d_A + data$d_B + data$d_C
  num_sigma_z0_beta <- 3*data$m
  num_param <- num_nu + num_sigma_z0_beta

  return(num_param)
}

#' Generate initial values using Pathfinder
#'
#' Uses CmdStanR's Pathfinder algorithm to obtain initial values for MCMC
#' sampling. If Pathfinder fails, a fallback initialization is used.
#'
#' @param mod A CmdStanR model object.
#' @param data A list of data formatted for Stan.
#' @param pathfinder_init Optional named list of initial values for Pathfinder.
#' @param num_paths Number of Pathfinder optimization paths (default 1).
#'
#' @return A list of initial values suitable for use in \code{mod$sample()}.
#' @export
get_initial_vals <- function(mod, data, pathfinder_init = NULL, num_paths = 1){

  if (is.null(pathfinder_init) == T) {pathfinder_init <- list(sigma = rep(1,data$m),
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
                                         num_paths = num_paths)
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

#' Run MCMC sampling for canonical DCM with convergence checks
#'
#' Runs Hamiltonian Monte Carlo (NUTS) sampling using CmdStanR in chunks,
#' continuing until a specified multivariate effective sample size (ESS)
#' threshold is reached or a maximum number of iterations is exceeded.
#'
#' @param mod A CmdStanR model object.
#' @param data A list of data formatted for Stan.
#' @param inits_list Initial values for sampling.
#' @param output_dir Directory for saving output files.
#' @param basename Base name for saved output files.
#' @param metric Metric type for HMC (\code{"dense_e"} or \code{"diag_e"}).
#' @param refresh Frequency of sampler output.
#' @param warmup_iter Number of warmup iterations.
#' @param n_iter_chunk Number of sampling iterations per chunk.
#' @param max_iter Maximum total number of iterations.
#' @param adapt_delta Target acceptance probability.
#' @param seed Random seed.
#' @param chains Number of chains.
#' @param ess_check Required multivariate effective sample size for convergence.
#'
#' @details
#' Sampling proceeds in chunks, reusing the last draw as initialization for the
#' next chunk. Convergence is assessed using multivariate ESS computed via
#' \code{mcmcse}.
#'
#' @return Invisibly returns NULL. Sampling results and diagnostics are saved to disk.
#' @export
dcm_sample <- function(mod, data, inits_list, output_dir, basename, ess_check, metric = c("dense_e","diag_e"),
                       refresh = 100, warmup_iter = 5000, n_iter_chunk = 1000,
                       max_iter = 10^6, adapt_delta = 0.9, seed = 1234, chains = 1){

  # Make sure ess_check is not empty
  if (missing(ess_check)) {
    stop("`ess_check` must be provided. It should reflect the required multivariate ESS for convergence.")
  }

  # Running stan program to sample from posterior - initial warmup and starting to sample
  fit = mod$sample(
    data = data,
    init = list(inits_list), # initialize from pathfinder values
    refresh = refresh, # output frequency
    iter_warmup = warmup_iter, # warm-up iterations
    iter_sampling = n_iter_chunk, # sampling iterations
    seed = seed, # seed for reproducibility
    chains = chains,
    adapt_delta = adapt_delta,
    save_warmup = TRUE,
    output_dir = output_dir,
    output_basename = paste0(basename,"_out"),
    metric = metric)

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

    fit = mod$sample(
      data = data,
      init = list(last_init),
      refresh = refresh, # output frequency
      iter_warmup = 0, # no warmup, only onto sampling
      iter_sampling = n_iter_chunk, # sampling iterations
      seed = seed, # seed for reproducibility
      chains = chains,
      adapt_delta = adapt_delta,
      adapt_engaged = FALSE,
      metric = metric,
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
    ess_ok <- ifelse(is.na(multi_ess), F, ifelse(multi_ess > ess_check, T, F))

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
}











