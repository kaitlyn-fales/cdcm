# cdcm

An R package for Bayesian estimation of effective connectivity in task-based fMRI using Canonical Dynamic Causal Modeling (CDCM).

## Overview

`cdcm` provides tools for:

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

---

## Requirements

* R (≥ 4.2 recommended)
* [`cmdstanr`](https://mc-stan.org/cmdstanr/)
* CmdStan (≥ 2.35)

To install CmdStan:

```r
cmdstanr::install_cmdstan()
```

---

## Example Workflow

```r
library(cdcm)

# Compile model
mod <- compile_cdcm()

# Load preprocessed data - preprocessing functions to come in later versions
# Preprocessed data should be a list with the named elements: dat$times, dat$u, dat$y_obs
load("preprocessed_data.RData")

# Set up hypotheses - indices for where to estimate parameters in A, B, C
idxs = list(A_idxs = matrix(c(2,1,
                              1,2,
                              1,1,
                              2,2), byrow = T, ncol=2),
            B_idxs = matrix(c(2,1,2), byrow=T, ncol = 3),
            C_idxs = matrix(c(1,1), byrow=T, ncol=2))

# Prepare data
stan_dat <- get_stan_dat(dat, idxs)

# Compute number of parameters for multivariate ESS convergence
num_param <- get_num_param(stan_dat)

# Convergence check specs - 95% HPD intervals with 5% tolerance
minESS_check <- as.numeric(mcmcse::minESS(num_param, alpha = 0.05, eps = 0.05))

# Initialize via Pathfinder
init <- get_initial_vals(mod, stan_dat)

# Run sampling
dcm_sample(
  mod = mod,
  data = stan_dat,
  inits_list = init,
  output_dir = "Output",
  basename = "cdcm_fit",
  ess_check = minESS_check,
  metric = "dense_e"
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


