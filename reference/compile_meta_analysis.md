# Compile the group-level meta-analysis Stan model

Compiles the meta-analysis Stan program bundled with the `cdcm` package
using `cmdstanr`.

## Usage

``` r
compile_meta_analysis(force_recompile = FALSE, quiet = TRUE, ...)
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
