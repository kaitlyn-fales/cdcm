# Simulation input checks
simulation_format_check <- function(nu, hypothesis_idxs, nscan, m, n_u, TR, U){
  # Check nscan
  if (length(nscan) != 1 || !is.numeric(nscan) || is.na(nscan) ||
      nscan <= 0 || nscan != as.integer(nscan)) {
    stop("`nscan` must be a single positive integer.")
  }
  nscan <- as.integer(nscan)

  # Check n_u
  if (length(n_u) != 1 || !is.numeric(n_u) || is.na(n_u) ||
      n_u <= 0 || n_u != as.integer(n_u)) {
    stop("`n_u` must be a single positive integer.")
  }
  n_u <- as.integer(n_u)

  # Check TR
  if (length(TR) != 1 || !is.numeric(TR) || is.na(TR) || TR <= 0) {
    stop("`TR` must be a single positive number.")
  }

  # Construct times vector including initial time 0
  times <- seq(0, nscan * TR, by = TR)

  # Check that U is supplied
  if (missing(U) || is.null(U)) {
    stop("`U` must be supplied.")
  }

  # Allow vector input only when n_u = 1
  if (is.vector(U)) {
    if (n_u == 1) {
      U <- matrix(U, ncol = 1)
      message("`U` was supplied as a vector and converted to a one-column matrix.")
    } else {
      stop("`U` must be a matrix with `n_u` columns.")
    }
  }

  # Check that U is a matrix
  if (!is.matrix(U)) {
    stop("`U` must be a matrix.")
  }

  # Check type and missingness
  if (!is.numeric(U) || anyNA(U)) {
    stop("`U` must be a numeric matrix with no missing values.")
  }

  # Check number of columns
  if (ncol(U) != n_u) {
    stop("`U` must have `n_u` columns.")
  }

  # Users will often supply an nscan x n_u design matrix.
  # For simulation, prepend the first row to create the time-0 initial-condition row.
  if (nrow(U) == nscan) {
    input <- rbind(U[1, , drop = FALSE], U)
    message("`U` had `nscan` rows, so the first row was prepended to create the time-0 initial-condition row.")
  } else if (nrow(U) == (nscan + 1)) {
    input <- U
  } else {
    stop("`U` must have either `nscan` rows or `nscan + 1` rows.")
  }

  # Check that U is binary
  if (!all(input %in% c(0, 1))) {
    stop("All entries of `U` must be in {0, 1}.")
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

  # Require at least one hypothesized nonzero parameter in each of A, B, and C
  if (nrow(hypothesis_idxs$A_idxs) < 1) {
    stop("`hypothesis_idxs$A_idxs` must have at least one row.", call. = FALSE)
  }
  if (nrow(hypothesis_idxs$B_idxs) < 1) {
    stop("`hypothesis_idxs$B_idxs` must have at least one row.", call. = FALSE)
  }
  if (nrow(hypothesis_idxs$C_idxs) < 1) {
    stop("`hypothesis_idxs$C_idxs` must have at least one row.", call. = FALSE)
  }

  # Check indices are numeric positive integers
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
    if (anyNA(mat)) {
      stop("`hypothesis_idxs$", nm, "` cannot contain missing values.", call. = FALSE)
    }
    if (any(!is.finite(mat))) {
      stop("`hypothesis_idxs$", nm, "` must contain only finite values.", call. = FALSE)
    }
    if (any(mat <= 0) || any(mat %% 1 != 0)) {
      stop("`hypothesis_idxs$", nm, "` must contain positive integer indices.", call. = FALSE)
    }
  }

  # Check indices are within bounds of m and n_u
  # A_idxs: both columns in 1, ..., m
  if (any(hypothesis_idxs$A_idxs[, 1] > m) || any(hypothesis_idxs$A_idxs[, 2] > m)) {
    stop("Entries of `hypothesis_idxs$A_idxs` must be <= `m`.", call. = FALSE)
  }

  # B_idxs: first column in 1, ..., n_u; second and third columns in 1, ..., m
  if (any(hypothesis_idxs$B_idxs[, 1] > n_u)) {
    stop("The first column of `hypothesis_idxs$B_idxs` must be <= `n_u`.", call. = FALSE)
  }
  if (any(hypothesis_idxs$B_idxs[, 2] > m) || any(hypothesis_idxs$B_idxs[, 3] > m)) {
    stop("The second and third columns of `hypothesis_idxs$B_idxs` must be <= `m`.", call. = FALSE)
  }

  # C_idxs: first column in 1, ..., m; second column in 1, ..., n_u
  if (any(hypothesis_idxs$C_idxs[, 1] > m)) {
    stop("The first column of `hypothesis_idxs$C_idxs` must be <= `m`.", call. = FALSE)
  }
  if (any(hypothesis_idxs$C_idxs[, 2] > n_u)) {
    stop("The second column of `hypothesis_idxs$C_idxs` must be <= `n_u`.", call. = FALSE)
  }


  # Format checks for nu
  if (!is.list(nu)) {
    stop("`nu` must be a list.", call. = FALSE)
  }

  required_nu_names <- c("nu_A", "nu_B", "nu_C")
  missing_nu_names <- setdiff(required_nu_names, names(nu))
  if (length(missing_nu_names) > 0) {
    stop(
      "`nu` is missing required component(s): ",
      paste(missing_nu_names, collapse = ", "),
      call. = FALSE
    )
  }

  nu_vecs <- list(
    nu_A = nu$nu_A,
    nu_B = nu$nu_B,
    nu_C = nu$nu_C
  )

  for (nm in names(nu_vecs)) {
    vec <- nu_vecs[[nm]]

    if (!is.numeric(vec)) {
      stop("`nu$", nm, "` must be numeric.", call. = FALSE)
    }
    if (is.matrix(vec) || is.data.frame(vec) || length(dim(vec)) > 1) {
      stop("`nu$", nm, "` must be a vector, not a matrix or array.", call. = FALSE)
    }
    if (length(vec) < 1) {
      stop("`nu$", nm, "` must contain at least one value.", call. = FALSE)
    }
    if (anyNA(vec)) {
      stop("`nu$", nm, "` cannot contain missing values.", call. = FALSE)
    }
    if (any(!is.finite(vec))) {
      stop("`nu$", nm, "` must contain only finite values.", call. = FALSE)
    }
  }

  # Check correspondence between nu and hypothesis_idxs
  if (length(nu$nu_A) != nrow(hypothesis_idxs$A_idxs)) {
    stop(
      "`nu$nu_A` must have length equal to nrow(`hypothesis_idxs$A_idxs`).",
      call. = FALSE
    )
  }

  if (length(nu$nu_B) != nrow(hypothesis_idxs$B_idxs)) {
    stop(
      "`nu$nu_B` must have length equal to nrow(`hypothesis_idxs$B_idxs`).",
      call. = FALSE
    )
  }

  if (length(nu$nu_C) != nrow(hypothesis_idxs$C_idxs)) {
    stop(
      "`nu$nu_C` must have length equal to nrow(`hypothesis_idxs$C_idxs`).",
      call. = FALSE
    )
  }
  # Return checked/validated inputs
  return(list(
    nu = nu,
    hypothesis_idxs = hypothesis_idxs,
    times = times,
    input = input,
    nscan = nscan,
    m = m,
    n_u = n_u,
    TR = TR
  ))
}

# Scenario a: u(t) = 0
eval_u_a <- function(t, z0, A, prev_index, times){
  expm::expAtv(A = A, t = times[t]-times[prev_index], v = z0)$eAtv
}

# Scenario b: u(t) != 0
eval_u_b <- function(t, z0, A, B, C, input, prev_index, times){
  A0 = A + Reduce("+",lapply(1:ncol(input), function(i) input[t-1,i]*B[[i]]))
  b0 = as.numeric(C%*%input[t-1,])
  zstar = -solve(A0,b0)
  expm::expAtv(A = A0,t = times[t]-times[prev_index],v = z0-zstar)$eAtv + zstar
}

# Function to get piecewise analytic solution to ODE
get_AnalyticSol <- function(input, changes, times, parms, m, v0 = NULL){
  # Get initial values if not specified
  if (is.null(v0) == T) {v0 <- rep(0.1,m)}

  # Matrix to store solution
  sol <- matrix(NA, nrow = length(times), ncol = m)
  sol[1,] <- v0 # add initial values

  # Initialize at u(t)=0
  prev_index <- 1

  # For loop to go from one change to the next
  for (i in 1:length(changes)){

    # Denote index for which we will work up to
    index <- changes[i]

    # Update solution based on change_type
    if (sum(input[index-1,]) == 0){
      sol[(prev_index+1):index,] <- t(sapply((prev_index+1):index, eval_u_a,
                                             z0 = sol[prev_index,], A = parms$A,
                                             prev_index = prev_index, times = times))
    } else {
      if(sum(input[index-1,]) != 0){
        sol[(prev_index+1):index,] <- t(sapply((prev_index+1):index, eval_u_b,
                                               z0 = sol[prev_index,], A = parms$A,
                                               B = parms$B, C = parms$C,
                                               input = input, prev_index = prev_index,
                                               times = times))
      }
    }

    # Update indices
    prev_index <- index

  }
  return(sol)
}

# Canonical HRF
HRF = function(t){stats::dgamma(x = t, shape = 6,rate = 1) - (1/6)*stats::dgamma(x = t,shape = 16,rate = 1)}

# Function for convolving z with canonical HRF
HRF_mu = function(mu,tp){
  stats::convolve(mu, rev(HRF(tp)),type="open")[1:length(tp)]
}

# Internal functions for random seeds
handle_seed <- function(seed, m) {
  if (is.null(seed)) {
    return(c(1:m))
  }

  if (!is.numeric(seed)) {
    stop("`seed` must be NULL or numeric.")
  }

  if (length(seed) == 1) {
    return(rep(seed, m))
  }

  if (length(seed) == m) {
    return(seed)
  }

  stop("`seed` must be NULL, a scalar, or a vector of length m.")
}

# Function for adding Gaussian noise at specified SNR to noiseless y(t)
add_noise <- function(y_conv,times,SNR,seed = NULL){
  m <- ncol(y_conv)

  # Seeds for reproducibility
  seed <- handle_seed(seed,m)

  # Empty vector of variance
  variance <-  numeric()

  # Empty matrix for y_obs
  y_obs <- matrix(NA, nrow = length(times)-1, ncol = ncol(y_conv))

  # Add Gaussian noise
  for (j in 1:m){
    set.seed(seed[j])
    variance[j] <- (stats::var(y_conv[,j])+(mean(y_conv[,j]))^2)/SNR
    y_obs[,j] <- y_conv[,j] + stats::rnorm(length(times)-1, mean = 0, sd = sqrt(variance[j]))
  }
  return(y_obs)
}

# Function for structuring paramMats
struct_paramMats = function(m, n_u, idxs, nu, reparam = TRUE){

  if (reparam == FALSE) {
    message("diag(A) is not being reparametrized, assuming properly transformed input values")
  }

  A = matrix(data = 0,nrow = m, ncol = m)
  for(i in 1:nrow(idxs$A_idxs)){
    A[idxs$A_idxs[i,1], idxs$A_idxs[i,2]] = with(nu,nu_A[i])
  }
  # Reparametrize diag(A) internally - if not T, needs to be done separately
  if (reparam) {
    diag(A) <- -0.5 * exp(diag(A))
  }

  B = lapply(1:n_u, function(x){matrix(data = 0,nrow = m, ncol = m)})
  for(i in 1:nrow(idxs$B_idxs)){
    B[[idxs$B_idxs[i,1]]][idxs$B_idxs[i,2], idxs$B_idxs[i,3]] = with(nu,nu_B[i])
  }

  C = matrix(data=0, nrow = m, ncol = n_u)
  for(i in 1:nrow(idxs$C_idxs)){
    C[idxs$C_idxs[i,1], idxs$C_idxs[i,2]] = with(nu,nu_C[i])
  }

  return(list(A=A,B=B,C=C))

}

# ODE function for deSolve
linear <- function(t, z, params, input) {
  # t is the current time point in the integration,
  # z is the current estimate of the variables in the ODE system
  A = params[["A"]]
  B = params[["B"]]
  C = params[["C"]]
  u = matrix(input(t),ncol=1)
  B_all = Reduce("+", lapply(1:nrow(u), function(i) u[i,1]*B[[i]]))
  dz <- (A+B_all)%*%z + C%*%u
  return(list(dz))
}




