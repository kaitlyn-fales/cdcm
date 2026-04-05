# cdcm: Canonical Dynamic Causal Modeling

An R package for Bayesian estimation of effective connectivity in task-based fMRI using Canonical Dynamic Causal Modeling (CDCM).

## Overview

`cdcm` provides tools for:

* Simulating fMRI time series under the CDCM framework
* Constructing Stan data for CDCM models
* Running efficient MCMC sampling via CmdStanR
* Generating initialization using Pathfinder
* Supporting reproducible effective connectivity analyses across studies

---

## Installation

```r
# Install from GitHub
devtools::install_github("kaitlyn-fales/cdcm")
```

## Requirements

* R (≥ 4.2 recommended)
* [`cmdstanr`](https://stan-dev.r-universe.dev/cmdstanr)
* CmdStan (≥ 2.35)

To install CmdStan:

```r
cmdstanr::install_cmdstan()
```

---

## Example Workflow

```r
library(cdcm)

# Define true connectivity parameters 
nu <- list(
  nu_A = c(0.4, 0.3, -0.1, 0.15),
  nu_B = c(-0.2),
  nu_C = c(0.7)
)

# Specify model structure 
idxs <- list(
  A_idxs = matrix(c(2,1,
                    1,2,
                    1,1,
                    2,2), byrow = TRUE, ncol = 2),
  B_idxs = matrix(c(1,2,1), byrow = TRUE, ncol = 3),
  C_idxs = matrix(c(1,1), byrow = TRUE, ncol = 2)
)

# Simple block design with a single stimulus
U <- matrix(rep(c(rep(0,5),rep(1,5)), 10), ncol = 1)

# Simulate data
sim <- simulate_cdcm(
  nu = nu,
  hypothesis_idxs = idxs,
  nscan = 100,
  m = 2,
  n_u = 1,
  TR = 2,
  U = U,
  SNR = 3,
  seed = 123
)

dat <- sim$simulated_data

# Compile model
mod <- compile_cdcm()

# Prepare data
stan_dat <- get_stan_dat(dat, idxs)

# Convergence criterion - Multivariate ESS (95% HPD, 5% tolerance)
minESS_check <- minESS_criterion(stan_dat, alpha = 0.05, eps = 0.05)

# Initialize via Pathfinder
init <- get_initial_vals(mod, stan_dat)

# Run sampling until multivariate ESS convergence
dir.create("Output", showWarnings = FALSE)
fit <- dcm_sample(
  mod = mod,
  data = stan_dat,
  inits_list = init,
  output_dir = "Output",
  basename = "cdcm_fit",
  ess_check = minESS_check,
)
```

---

## Models Included

* Canonical DCM (CDCM) for task-based fMRI
* Meta-analysis model for group-level inference

Additional models (e.g., multi-run meta-analysis) are planned.

---

## Notes

* CmdStan must be installed and configured before model compilation.
* Recompilation may be required when updating CmdStan:

  ```r
  compile_cdcm(force_recompile = TRUE)
  ```

---

## Citation

If you use `cdcm`, please cite the package in R with:

```r
citation("cdcm")
```


