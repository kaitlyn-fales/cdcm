# Simulate fMRI time series under CDCM

Simulates neural state trajectories and observed BOLD signals under a
CDCM with user-specified connectivity parameters and experimental
inputs.

## Usage

``` r
simulate_cdcm(
  nu,
  hypothesis_idxs,
  nscan,
  m,
  n_u,
  TR,
  U,
  SNR,
  reparam = TRUE,
  z0 = NULL,
  check_solution = TRUE,
  tol = 1e-04,
  seed = NULL
)
```

## Arguments

- nu:

  A list containing true parameter values:

  nu_A

  :   numeric vector of A (baseline connectivity) parameters

  nu_B

  :   numeric vector of B (modulatory) parameters

  nu_C

  :   numeric vector of C (input) parameters

- hypothesis_idxs:

  A list specifying the locations of nonzero parameters:

  A_idxs

  :   matrix with 2 columns (row, column indices)

  B_idxs

  :   matrix with 3 columns (input, row, column indices)

  C_idxs

  :   matrix with 2 columns (row, input indices)

- nscan:

  Number of scan time points.

- m:

  Number of brain regions (state dimension).

- n_u:

  Number of experimental inputs.

- TR:

  Repetition time (seconds).

- U:

  Input (stimulus) matrix of dimension `nscan x n_u` or
  `(nscan + 1) x n_u`. Values must be binary (0/1).

- SNR:

  Numeric. Signal-to-noise ratio for additive Gaussian noise, defined as
  the ratio of the standard deviation of the signal to that of the noise
  (i.e., \\\mathrm{SNR} = \mathrm{sd}(\text{signal}) /
  \mathrm{sd}(\text{noise})\\). Noise is generated such that
  \\\mathrm{sd}(\text{noise}) = \mathrm{sd}(\text{signal}) /
  \mathrm{SNR}\\, implying a noise variance of \\\sigma^2 =
  \mathrm{Var}(\text{signal}) / \mathrm{SNR}^2\\. Must be positive.

- reparam:

  Logical. If `TRUE`, reparameterizes diagonal elements of the A matrix
  to enforce stability.

- z0:

  Optional initial state vector. If `NULL`, a default is used.

- check_solution:

  Logical. If `TRUE`, compares analytic and numeric ODE solutions to
  assess numerical stability.

- tol:

  Numeric tolerance for the analytic vs numeric solution comparison.

- seed:

  Optional seed (scalar or length `m`) for reproducible noise
  generation.

## Value

A list with two components:

- simulated_data:

  times

  :   Vector of scan times (length `nscan`)

  u

  :   Input matrix at scan times (`nscan x n_u`)

  y_obs

  :   Observed BOLD signals (`nscan x m`)

- true_signals:

  z_signal

  :   Latent neural states (`nscan x m`)

  y_signal

  :   Noise-free BOLD signals (`nscan x m`)

## Details

This function generates:

- latent neural states via an analytic solution to a bilinear state
  equation

- observed BOLD signals via convolution with a canonical hemodynamic
  response function (HRF)

- additive Gaussian noise controlled by a specified signal-to-noise
  ratio (SNR)

The user supplies the connectivity structure through `nu` and
`hypothesis_idxs`, along with a stimulus design matrix `U`. Internally,
the function augments the input matrix to include the initial time point
required for solving the system.

The neural dynamics follow the bilinear state space DCM. Within each
block of constant input, the system admits a closed-form solution, which
is used here for efficient simulation. The observed BOLD signal is
obtained by convolving the neural states with the canonical HRF (linear
combination of two gamma distributions).

## Examples

``` r
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

U <- matrix(c(0,rep(c(rep(0,5),rep(1,5)), 2)), ncol = 1)

sim <- simulate_cdcm(
  nu = nu,
  hypothesis_idxs = hypothesis_idxs,
  nscan = 20,
  m = 2,
  n_u = 1,
  TR = 2,
  U = U,
  SNR = 5,
  seed = 1
)
```
