# Check parameter and latent-state identifiability diagnostics

Evaluates simulation- or plug-in-based diagnostics for the CDCM
identifiability assumptions involving the effective connectivity
matrices and the latent neural states.

## Usage

``` r
check_latent_identifiability(
  A,
  B,
  u,
  z,
  block_indices = NULL,
  tol = 1e-08,
  warn = TRUE
)
```

## Arguments

- A:

  A \\d \times d\\ numeric matrix giving the baseline effective
  connectivity matrix.

- B:

  A list of \\d \times d\\ numeric matrices giving the modulatory
  connectivity matrices. The length of `B` should match the number of
  stimulus inputs.

- u:

  A \\T \times m\\ numeric matrix encoding the stimulus inputs at each
  observation time.

- z:

  A \\T \times d\\ numeric matrix of latent neural states evaluated at
  the observation times. Rows correspond to time points and columns
  correspond to ROIs.

- block_indices:

  Optional integer vector giving the block indices to check. If `NULL`,
  all eligible blocks with at least \\d + 1\\ observations are checked.

- tol:

  Numeric tolerance used for eigenvalue, rank, and determinant checks.
  Defaults to `1e-8`.

- warn:

  Logical; if `TRUE`, warnings are issued when diagnostics fail.
  Defaults to `TRUE`.

## Value

A list containing:

- `cond_A3`:

  Logical indicating whether A3 is satisfied for all checked blocks.

- `cond_A4`:

  Logical indicating whether A4 is satisfied for all checked blocks.

- `all_ok`:

  Logical indicating whether both diagnostics are satisfied.

- `block_results`:

  A data frame containing block-level diagnostic results.

## Details

This function checks the CDCM assumptions involving the true parameters
and latent discrete observations. Unlike the design assumptions, these
conditions generally cannot be verified directly from real data because
the true parameters and latent states are unknown. Therefore, this
function is intended primarily for simulation studies, where `A`, `B`,
and `z` are known, or as a post-fit plug-in diagnostic using posterior
summaries.

For each checked block \\b\\, the function forms \$\$ \tilde A^{(b)} =
A + \sum\_{i = 1}^m \tilde u_i^{(b)} B_i, \$\$ where \\\tilde u^{(b)}\\
is the block-level stimulus vector, obtained from the first row of `u`
within block \\b\\. Assumption A3 is checked by verifying that \\\tilde
A^{(b)}\\ has real and distinct eigenvalues.

Assumption A4 is checked by constructing \$\$ X^{(b)} = \[x_0^{(b)},
\ldots, x_d^{(b)}\], \qquad x_j^{(b)} = \[z(t^{(b)} + t_j)^T, 1\]^T,
\$\$ from the first \\d + 1\\ latent observations in block \\b\\, and
checking whether \\X^{(b)}\\ has full rank.

When used with posterior means from real data, these checks should be
interpreted only as plug-in diagnostics conditional on the quality and
convergence of the posterior fit.
