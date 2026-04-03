# Path to the canonical DCM Stan program

Returns the file path to the CDCM Stan program bundled with the `cdcm`
package.

## Usage

``` r
cdcm_stan_file()
```

## Value

A character string giving the full path to the `canonical_dcm.stan`
file.

## Details

This function is primarily intended for use with
[`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html)
when compiling the model.

## Examples

``` r
if (FALSE) { # \dontrun{
stan_file <- cdcm_stan_file()
mod <- cmdstanr::cmdstan_model(stan_file)
} # }
```
