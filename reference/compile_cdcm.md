# Compile the CDCM Stan model

Compiles the CDCM Stan program bundled with the `cdcm` package using
`cmdstanr`.

## Usage

``` r
compile_cdcm(force_recompile = FALSE, quiet = TRUE, ...)
```

## Arguments

- force_recompile:

  Logical; if `TRUE`, forces recompilation of the Stan model even if a
  compiled executable already exists.

- quiet:

  Logical; passed to
  [`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html)
  to control console output during compilation.

- ...:

  Additional arguments passed to
  [`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html).

## Value

A `CmdStanModel` object.

## Details

This function provides a convenient interface for compiling the model
without requiring users to manually manage file paths to the Stan
program.

This function requires a working CmdStan installation (version 2.35.0 or
newer) configured for use with `cmdstanr`. If CmdStan is not installed
or configured, an error will be thrown.

## Examples

``` r
if (FALSE) { # \dontrun{
mod <- compile_cdcm()

# Force recompilation if needed
mod <- compile_cdcm(force_recompile = TRUE)
} # }
```
