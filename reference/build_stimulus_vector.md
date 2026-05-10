# Build binary block-design stimulus vector

Constructs a binary stimulus indicator vector from onset and duration
information specified in seconds.

## Usage

``` r
build_stimulus_vector(onsets, durations, nscan, TR)
```

## Arguments

- onsets:

  Numeric vector of stimulus onset times in seconds.

- durations:

  Numeric vector of stimulus durations in seconds.

- nscan:

  Number of scans/time points.

- TR:

  Repetition time in seconds.

## Value

A binary vector of length `nscan`, where 1 indicates stimulus presence
and 0 indicates absence.

## Details

Onset and duration values are converted from seconds to scan indices
using the supplied repetition time `TR`. Specifically, onset times are
mapped using `ceiling(onsets / TR)`, and durations are mapped using
`ceiling(durations / TR)`.

If a stimulus interval extends beyond `nscan`, it is truncated to the
scan range.

## Examples

``` r
# Single stimulus
build_stimulus_vector(
  onsets = c(0, 20, 40),
  durations = c(10, 10, 10),
  nscan = 30,
  TR = 2
)
#>  [1] 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0

# Multiple stimulus conditions with rest blocks
visual_stim <- build_stimulus_vector(
  onsets = c(20, 80),
  durations = c(20, 20),
  nscan = 60,
  TR = 2
)

auditory_stim <- build_stimulus_vector(
  onsets = c(50, 110),
  durations = c(20, 10),
  nscan = 60,
  TR = 2
)

task_input <- cbind(visual_stim, auditory_stim)
colnames(task_input) <- c("Visual", "Auditory")

task_input
#>       Visual Auditory
#>  [1,]      0        0
#>  [2,]      0        0
#>  [3,]      0        0
#>  [4,]      0        0
#>  [5,]      0        0
#>  [6,]      0        0
#>  [7,]      0        0
#>  [8,]      0        0
#>  [9,]      0        0
#> [10,]      1        0
#> [11,]      1        0
#> [12,]      1        0
#> [13,]      1        0
#> [14,]      1        0
#> [15,]      1        0
#> [16,]      1        0
#> [17,]      1        0
#> [18,]      1        0
#> [19,]      1        0
#> [20,]      0        0
#> [21,]      0        0
#> [22,]      0        0
#> [23,]      0        0
#> [24,]      0        0
#> [25,]      0        1
#> [26,]      0        1
#> [27,]      0        1
#> [28,]      0        1
#> [29,]      0        1
#> [30,]      0        1
#> [31,]      0        1
#> [32,]      0        1
#> [33,]      0        1
#> [34,]      0        1
#> [35,]      0        0
#> [36,]      0        0
#> [37,]      0        0
#> [38,]      0        0
#> [39,]      0        0
#> [40,]      1        0
#> [41,]      1        0
#> [42,]      1        0
#> [43,]      1        0
#> [44,]      1        0
#> [45,]      1        0
#> [46,]      1        0
#> [47,]      1        0
#> [48,]      1        0
#> [49,]      1        0
#> [50,]      0        0
#> [51,]      0        0
#> [52,]      0        0
#> [53,]      0        0
#> [54,]      0        0
#> [55,]      0        1
#> [56,]      0        1
#> [57,]      0        1
#> [58,]      0        1
#> [59,]      0        1
#> [60,]      0        0
```
