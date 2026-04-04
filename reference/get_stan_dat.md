# Construct a Stan data list for the CDCM model

Prepares the input data and hypothesis indexing structures required by
the CDCM Stan program.

## Usage

``` r
get_stan_dat(
  data,
  hypothesis_idxs,
  ode_solver_type = 1,
  tol = 1e-05,
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

  A list containing the observed time series and experimental design. It
  must contain:

  - `times`: a numeric vector of observation times

  - `u`: a \\T \times n_u\\ matrix encoding the experimental stimulus
    inputs

  - `y_obs`: a \\T \times m\\ matrix of observed ROI-level BOLD signals

- hypothesis_idxs:

  A list specifying the hypothesized nonzero entries of the connectivity
  matrices. It must contain:

  - `A_idxs`: a 2-column matrix of \\(i,j)\\ indices for the baseline
    connectivity matrix \\A\\

  - `B_idxs`: a 3-column matrix of \\(k,i,j)\\ indices for the
    modulatory matrices \\B_k\\, where \\k\\ indexes the stimulus
    dimension

  - `C_idxs`: a 2-column matrix of \\(i,j)\\ indices for the
    driving-input matrix \\C\\

- ode_solver_type:

  Integer code specifying the ODE solver used in the Stan program
  (default = 1; piecewise analytic solution, faster than CKRK numeric
  integration = 0).

- tol:

  Numeric tolerance used for both relative and absolute ODE solver
  tolerances (default = 1e-5; used for CKRK numeric integration).

- sigma_nu:

  Prior scale for off-diagonal neural connectivity parameters (default =
  1).

- sigma_nu_self:

  Prior scale for diagonal self-connectivity parameters (default =
  0.125; tighter prior is needed).

- rate_sigma:

  Rate parameter for the prior on observation noise.

- sigma_z0:

  Prior scale for the initial neural state.

- conv:

  Integer indicator controlling whether convolution with the HRF is
  performed in the Stan model (default = 1; leave as-is unless you want
  latent z(t) instead of y(t)).

- max_num_steps:

  Maximum number of ODE solver steps allowed (default = 10^6; CKRK
  numeric solver only).

## Value

A named list containing the data, dimensions, changepoint structure,
hypothesis index matrices, prior hyperparameters, and solver settings
required by the CDCM Stan program.

## Details

This function validates the structure of the supplied data and
hypothesis index objects, computes changepoints in the experimental
design, and returns a named list in the format required by the CDCM Stan
model.

In addition, the function evaluates a set of design diagnostics via
[`check_design_identifiability`](https://kaitlyn-fales.github.io/cdcm/reference/check_design_identifiability.md).
If those sufficient conditions are not satisfied, a warning is issued.
This warning is intended as a practical diagnostic only: failure to
satisfy these conditions does not necessarily imply non-identifiability,
but it may indicate that the design is weak for reliable estimation.

The returned Stan data list can be used for either single-chain or
multi-chain sampling. Its structure does not depend on the number of
chains.
