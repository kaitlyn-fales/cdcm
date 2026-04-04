# Sample from the CDCM posterior in checkpointed chunks

Runs posterior sampling for a compiled CDCM model in sequential chunks,
saving intermediate results to disk after each chunk and stopping when a
target multivariate effective sample size (ESS) is reached or a maximum
iteration limit is exceeded.

## Usage

``` r
dcm_sample(
  mod,
  data,
  inits_list,
  output_dir,
  basename,
  ess_check,
  refresh = 100,
  warmup_iter = 5000,
  n_iter_chunk = 1000,
  max_iter = 10^6,
  adapt_delta = 0.9,
  seed = 1234
)
```

## Arguments

- mod:

  A compiled CmdStanR model object, typically returned by
  [`compile_cdcm`](https://kaitlyn-fales.github.io/cdcm/reference/compile_cdcm.md).

- data:

  A named list of data in the format required by the Stan model,
  typically produced by
  [`get_stan_dat`](https://kaitlyn-fales.github.io/cdcm/reference/get_stan_dat.md).

- inits_list:

  A named list of initial values for model parameters, typically
  produced by
  [`get_initial_vals`](https://kaitlyn-fales.github.io/cdcm/reference/get_initial_vals.md).

- output_dir:

  Character string giving the directory where checkpointed posterior
  draws and diagnostics should be saved.

- basename:

  Character string used as the base filename for saved output files.

- ess_check:

  Numeric scalar giving the target multivariate ESS used for convergence
  assessment.

- refresh:

  Integer giving the sampling progress refresh rate passed to CmdStanR.

- warmup_iter:

  Integer giving the number of warmup iterations for the initial
  sampling phase.

- n_iter_chunk:

  Integer giving the number of post-warmup sampling iterations to run
  per chunk.

- max_iter:

  Integer giving the maximum total number of post-warmup draws allowed
  across all chunks.

- adapt_delta:

  Numeric target acceptance rate passed to CmdStanR.

- seed:

  Integer random seed for reproducibility.

## Value

Invisibly returns a list with elements:

- draws:

  A data frame containing the combined posterior draws across all
  completed sampling chunks.

- diagnostics:

  A list of NUTS diagnostics for each completed sampling chunk.

- converged:

  Logical indicating whether the ESS-based convergence criterion was
  satisfied.

- total_draws:

  Total number of post-warmup draws generated across all completed
  chunks.

If sampling is interrupted before the function completes, checkpointed
draws and diagnostics saved to disk remain available, but no final
return object is produced.

## Details

This function is designed for checkpointed **single-chain** sampling and
is especially useful for long runs on HPC systems, where jobs may need
to stop at walltime limits. For standard multi-chain sampling, use the
compiled CmdStanR model directly via `mod$sample()`.

Sampling begins with an initial warmup phase, followed by posterior
sampling in chunks. After the first chunk, the last draw is used as the
initial value for the next chunk, and the adapted step size and inverse
metric from the previous fit are reused so that sampling can continue
efficiently without repeating warmup.

After each completed chunk, posterior draws and NUTS diagnostics are
saved to disk in `output_dir` using filenames based on `basename`. This
checkpointing allows intermediate results to be recovered if sampling is
interrupted, for example due to HPC walltime limits.

Convergence is assessed using the multivariate effective sample size
(ESS), computed via `mcmcse`, on the accumulated single-chain posterior
draws. Sampling continues until the target ESS specified by `ess_check`
is reached or the maximum number of iterations is exceeded.

This function is intentionally restricted to **single-chain** sampling.
Its checkpointing logic, reuse of the last draw as initialization, reuse
of the adapted inverse metric and step size, and ESS calculation are all
based on a single continuing chain. If multiple chains are desired,
users should call the underlying CmdStanR sampler directly via
`mod$sample()` (or the compiled model object's `$sample()` method),
specifying their preferred number of chains and parallelization
settings. In that case, the checkpointed chunkwise workflow implemented
here is not used.

## Multi-chain sampling

For standard multi-chain posterior sampling, call the compiled CmdStanR
model directly rather than using `dcm_sample()`. For example:

    mod <- compile_cdcm()
    stan_dat <- get_stan_dat(data, hypothesis_idxs)
    init_list <- get_initial_vals(mod, stan_dat)

    fit <- mod$sample(
      data = stan_dat,
      init = init_list,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      seed = 1234
    )

This approach provides standard multi-chain CmdStanR sampling, but does
not include the checkpointing or chunkwise ESS-based stopping rule
implemented in `dcm_sample()`.
