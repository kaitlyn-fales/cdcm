# Generate initial values using Pathfinder

Uses CmdStanR's Pathfinder algorithm to obtain initial values for MCMC
sampling. If Pathfinder fails, a fallback initialization is used.

## Usage

``` r
get_initial_vals(mod, data, pathfinder_init = NULL, num_paths = 1)
```

## Arguments

- mod:

  A CmdStanR model object.

- data:

  A list of data formatted for Stan.

- pathfinder_init:

  Optional named list of initial values for Pathfinder.

- num_paths:

  Number of Pathfinder optimization paths (default 1).

## Value

A list of initial values suitable for use in `mod$sample()`.
