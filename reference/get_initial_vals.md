# Generate initial values for CDCM posterior sampling

Produces a list of initial values for the CDCM Stan model, optionally
using CmdStanR Pathfinder to obtain data-informed initializations.

## Usage

``` r
get_initial_vals(mod, data, pathfinder_init = NULL)
```

## Arguments

- mod:

  A compiled CmdStanR model object, typically returned by
  [`compile_cdcm`](https://kaitlyn-fales.github.io/cdcm/reference/compile_cdcm.md).

- data:

  A named list of data in the format required by the Stan model,
  typically produced by
  [`get_stan_dat`](https://kaitlyn-fales.github.io/cdcm/reference/get_stan_dat.md).

- pathfinder_init:

  Optional named list of initial values to seed the Pathfinder run. If
  `NULL`, default values are constructed internally.

## Value

A list of initial values suitable for use in Stan sampling. When
Pathfinder succeeds, this is based on the final Pathfinder draw. When
Pathfinder fails, a backup initialization list based on default values
is returned instead.

## Details

If `pathfinder_init` is not supplied, default initial values are used:
observation noise parameters are initialized at 1, neural connectivity
parameters are initialized at 0, initial neural states are initialized
at 0.1, and hemodynamic scaling parameters are initialized at 0.

The function then attempts to run CmdStanR Pathfinder and extract the
final Pathfinder draw in a format suitable for posterior sampling. If
Pathfinder fails for any reason, a backup initialization list based on
the default values is returned instead, and a message is issued.

The object returned by this function is naturally a **single-chain**
initialization list. For multi-chain sampling via `mod$sample()`, users
may either:

- run `get_initial_vals()` separately to generate one initialization
  list per chain, or

- reuse/replicate a single initialization list across chains if that is
  appropriate for the analysis.
