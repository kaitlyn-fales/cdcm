# Check whether a block design satisfies CDCM identifiability assumptions

Evaluates practical sufficient-condition diagnostics for the CDCM
identifiability result under a block-design experimental input. The
function checks whether the full scan contains enough usable blocks and
whether the corresponding block-level stimulus combinations are
sufficiently unconfounded.

## Usage

``` r
check_design_identifiability(data, tol = 1e-08, warn = TRUE)
```

## Arguments

- data:

  A list containing:

  - `times`: a numeric vector of equally spaced observation times.

  - `u`: a \\T \times m\\ numeric matrix encoding the stimulus inputs at
    each observation time.

  - `y_obs`: a \\T \times d\\ numeric matrix of observed ROI-level BOLD
    signals.

- tol:

  Numeric tolerance used for the equal-spacing and matrix-rank checks.
  Defaults to `1e-8`.

- warn:

  Logical; if `TRUE`, warnings are issued when one or more diagnostics
  are not satisfied. Defaults to `TRUE`.

## Value

A list containing:

- `cond_A1`:

  Logical indicating whether the A1 diagnostic is satisfied.

- `cond_A2`:

  Logical indicating whether the A2 diagnostic is satisfied.

- `all_ok`:

  Logical indicating whether both diagnostics are satisfied.

- `equally_spaced`:

  Logical indicating whether the observation times are equally spaced up
  to `tol`.

- `n_eligible_blocks`:

  Number of blocks containing at least \\d + 1\\ observations.

- `rank_U_aug`:

  Rank of the augmented block-level stimulus matrix formed from eligible
  blocks.

## Details

This function checks the experimental design assumptions used in the
CDCM identifiability result. The input matrix `u` is provided at the
observation level, with each row corresponding to \\u(t_j)\\. Since the
design is assumed to be blockwise constant, the block-level stimulus
vector \\\tilde u^{(b)}\\ is obtained by taking the first row of `u`
within each block.

Assumption (A1) is checked by verifying that the observation times are
equally spaced and that at least \\m + 1\\ blocks contain at least \\d +
1\\ observations, where \\d\\ is the number of ROIs and \\m\\ is the
number of stimulus inputs.

Assumption (A2) is checked by forming the augmented block-level stimulus
matrix \$\$ U^\* = \left\[ \begin{array}{cc} 1 & \tilde{u}^{(b_1)T} \\
\vdots & \vdots \\ 1 & \tilde{u}^{(b_q)T} \end{array} \right\] \$\$ from
the blocks satisfying A1's block-length requirement, and checking
whether this matrix has rank \\m + 1\\. This is equivalent to checking
whether there exists a subset of \\m + 1\\ usable blocks whose augmented
stimulus matrix is invertible.

These diagnostics are sufficient-condition checks. Failure to satisfy
one or more diagnostics does not necessarily imply that the model is not
identifiable, but it may indicate that the design is weak for reliable
estimation.
