# Predict CDCM latent and BOLD signals from posterior draws

Predict CDCM latent and BOLD signals from posterior draws

## Usage

``` r
predict_cdcm(
  draws_df,
  hypothesis_idxs,
  nscan,
  m,
  n_u,
  TR,
  U,
  check_identifiability = FALSE,
  identifiability_tol = 1e-08
)
```

## Arguments

- draws_df:

  Posterior draws data frame. The first column is assumed to be the
  log-posterior column, the last three columns are assumed to be draw
  metadata, and the parameter columns are assumed to be ordered as:
  noise standard deviations, `nu_A`, `nu_B`, `nu_C`, `z0`, and `beta`.

- hypothesis_idxs:

  List of model index matrices containing `A_idxs`, `B_idxs`, and
  `C_idxs`.

- nscan:

  Number of scan time points.

- m:

  Number of brain regions/ROIs, denoted d in the manuscript.

- n_u:

  Number of stimulus inputs, denoted m in the manuscript.

- TR:

  Repetition time.

- U:

  Input stimulus matrix.

- check_identifiability:

  Logical. If TRUE, runs post-fit plug-in identifiability diagnostics.

- identifiability_tol:

  Numeric tolerance for identifiability diagnostics.

## Value

A list containing predicted latent states and predicted BOLD signal, and
optional identifiability check.
