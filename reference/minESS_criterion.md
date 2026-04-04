# Compute a minimum ESS threshold for CDCM sampling

Computes the minimum effective sample size (ESS) threshold corresponding
to a specified confidence level and relative precision, using the number
of model parameters implied by a CDCM Stan data list.

## Usage

``` r
minESS_criterion(data, alpha = 0.05, eps = 0.05)
```

## Arguments

- data:

  A named list of data in the format required by the Stan model,
  typically produced by
  [`get_stan_dat`](https://kaitlyn-fales.github.io/cdcm/reference/get_stan_dat.md).

- alpha:

  Significance level used in the
  [`mcmcse::minESS()`](https://rdrr.io/pkg/mcmcse/man/minESS.html)
  calculation. Defaults to `0.05`.

- eps:

  Relative precision target used in the
  [`mcmcse::minESS()`](https://rdrr.io/pkg/mcmcse/man/minESS.html)
  calculation. Defaults to `0.05`.

## Value

A numeric scalar giving the minimum ESS threshold.

## Details

This function first determines the number of model parameters implied by
the supplied Stan data list, then uses
[`mcmcse::minESS()`](https://rdrr.io/pkg/mcmcse/man/minESS.html) to
compute the minimum effective sample size required to achieve the
specified level `alpha` and relative precision `eps`.

This can be used to define the `ess_check` argument in
[`dcm_sample`](https://kaitlyn-fales.github.io/cdcm/reference/dcm_sample.md).

## Examples

``` r
if (FALSE) { # \dontrun{
stan_dat <- get_stan_dat(data, hypothesis_idxs)
ess_target <- minESS_criterion(stan_dat)
} # }
```
