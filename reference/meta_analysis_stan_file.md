# Path to the group-level meta-analysis Stan program

Returns the file path to the group-level meta analysis model Stan
program bundled with the `cdcm` package.

## Usage

``` r
meta_analysis_stan_file()
```

## Value

A character string giving the full path to the `meta_analysis.stan`
file.

## Details

This function is primarily intended for use with
[`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html)
when compiling the model.

## Examples

``` r
if (FALSE) { # \dontrun{
stan_file <- met_analysis_stan_file()
mod <- cmdstanr::cmdstan_model(stan_file)
} # }
```
