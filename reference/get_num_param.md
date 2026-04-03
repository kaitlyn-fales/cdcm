# Compute number of model parameters

Calculates the total number of parameters in the canonical DCM model,
based on the dimensions of connectivity components and variance terms.
This is primarily used to determine convergence thresholds (e.g.,
required ESS).

## Usage

``` r
get_num_param(data)
```

## Arguments

- data:

  A list containing model dimensions:

  - `d_A`: Number of A parameters

  - `d_B`: Number of B parameters

  - `d_C`: Number of C parameters

  - `m`: Number of regions (ROIs)

## Value

Integer representing total number of parameters.
