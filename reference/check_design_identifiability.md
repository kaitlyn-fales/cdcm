# Check whether a block design satisfies sufficient identifiability diagnostics

Evaluates design sufficient conditions for identifiability under a
block-design experimental input. These checks are intended as practical
design diagnostics for CDCM analyses and related simulation studies.

## Usage

``` r
check_design_identifiability(data, warn = TRUE)
```

## Arguments

- data:

  A list containing:

  - `times`: a numeric vector of observation times

  - `u`: a \\T \times n_u\\ matrix encoding the experimental stimulus
    inputs

  - `y_obs`: a \\T \times m\\ matrix of observed ROI-level BOLD signals

- warn:

  Logical; if `TRUE`, warnings are issued when one or more sufficient
  conditions are not satisfied.

## Value

A list containing:

- T:

  The number of observations.

- m:

  The number of ROIs.

- n_u:

  The number of stimulus dimensions.

- changepoints:

  The changepoint vector returned internally by `get_changepoints()`.

- block_lengths:

  The block lengths computed from the design matrix.

- first_block_zero:

  Logical indicating whether the first input block is the all-zero input
  configuration.

- n_distinct_block_inputs:

  The number of distinct blockwise input configurations observed.

- cond_A1_first:

  Logical indicating whether the first block satisfies the minimum
  block-length requirement.

- cond_A1_later:

  Logical indicating whether all later blocks satisfy the minimum
  block-length requirement.

- cond_A1:

  Logical indicating whether the block-length condition is satisfied
  overall.

- cond_A2:

  Logical indicating whether the observation-count condition is
  satisfied.

- cond_A3:

  Logical or `NA`. For designs with `n_u = 2` or `3`, indicates whether
  enough distinct blockwise input settings are observed. For other
  values of `n_u`, this is `NA`.

- all_ok:

  Logical indicating whether all evaluated sufficient-condition
  diagnostics are satisfied.

## Details

This function evaluates a set of practical design diagnostics motivated
by sufficient conditions for identifiability in the CDCM framework under
block-design experiments. In particular, it checks:

- whether block lengths are long enough relative to the number of ROIs

- whether the total number of observations is large enough relative to
  the number of ROIs and stimulus inputs

- whether enough distinct blockwise input configurations are observed
  when the number of stimulus dimensions is 2 or 3

These checks are based on sufficient conditions and should be
interpreted as practical diagnostics rather than necessary conditions.
Failure to satisfy one or more of them does not necessarily imply
non-identifiability, but it may indicate that the design is weak for
reliable estimation.

The number of ROIs is taken to be `m = ncol(y_obs)`, the number of
stimulus dimensions is `n_u = ncol(u)`, and the number of observations
is `T = length(times)`.
