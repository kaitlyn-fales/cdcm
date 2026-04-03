# Construct Stan data list for CDCM

Prepares input data and model specification into the list format
required for the canonical DCM Stan program. This includes time series
data, stimulus inputs, connectivity structure indices, and
hyperparameters.

## Usage

``` r
get_stan_dat(
  data,
  hypothesis_idxs,
  ode_solver_type = 1,
  tol = 10^{
-5
 },
  sigma_nu = 1,
  sigma_nu_self = 0.125,
  rate_sigma = 0.5,
  sigma_z0 = 0.3,
  conv = 1,
  max_num_steps = 10^6
)
```

## Arguments

- data:

  A list containing:

  - `times`: Numeric vector of time points (length T)

  - `u`: Matrix of stimulus inputs (T x n_u)

  - `y_obs`: Matrix of observed BOLD signals (T x m)

- hypothesis_idxs:

  A list specifying model structure:

  - `A_idxs`: Matrix of indices for A matrix (connections)

  - `B_idxs`: Matrix of indices for B matrix (modulatory effects)

  - `C_idxs`: Matrix of indices for C matrix (inputs)

- ode_solver_type:

  Integer specifying ODE solver (default 1 is piecewise analytic
  solution)

- tol:

  Numeric tolerance for ODE solver (default 1e-5, if using numeric CKRK
  integration)

- sigma_nu:

  Prior scale for connectivity parameters

- sigma_nu_self:

  Prior scale for self-connections in A and B

- rate_sigma:

  Prior rate parameter for variance terms

- sigma_z0:

  Prior scale for initial neural states

- conv:

  Integer flag for convolution option (default 1 convolves with
  canonical HRF - recommended to not change unless you want latent
  neural signal)

- max_num_steps:

  Maximum number of ODE solver steps (for CKRK solver)

## Value

A named list suitable for input to a Stan model.
