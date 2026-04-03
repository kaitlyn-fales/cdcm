# Run MCMC sampling for canonical DCM with convergence checks

Runs Hamiltonian Monte Carlo (NUTS) sampling using CmdStanR in chunks,
continuing until a specified multivariate effective sample size (ESS)
threshold is reached or a maximum number of iterations is exceeded.

## Usage

``` r
dcm_sample(
  mod,
  data,
  inits_list,
  output_dir,
  basename,
  ess_check,
  metric = c("dense_e", "diag_e"),
  refresh = 100,
  warmup_iter = 5000,
  n_iter_chunk = 1000,
  max_iter = 10^6,
  adapt_delta = 0.9,
  seed = 1234,
  chains = 1
)
```

## Arguments

- mod:

  A CmdStanR model object.

- data:

  A list of data formatted for Stan.

- inits_list:

  Initial values for sampling.

- output_dir:

  Directory for saving output files.

- basename:

  Base name for saved output files.

- ess_check:

  Required multivariate effective sample size for convergence.

- metric:

  Metric type for HMC (`"dense_e"` or `"diag_e"`).

- refresh:

  Frequency of sampler output.

- warmup_iter:

  Number of warmup iterations.

- n_iter_chunk:

  Number of sampling iterations per chunk.

- max_iter:

  Maximum total number of iterations.

- adapt_delta:

  Target acceptance probability.

- seed:

  Random seed.

- chains:

  Number of chains.

## Value

Invisibly returns NULL. Sampling results and diagnostics are saved to
disk.

## Details

Sampling proceeds in chunks, reusing the last draw as initialization for
the next chunk. Convergence is assessed using multivariate ESS computed
via `mcmcse`.
