#' Build binary block-design stimulus vector
#'
#' Constructs a binary stimulus indicator vector from onset and duration
#' information specified in seconds.
#'
#' @param onsets Numeric vector of stimulus onset times in seconds.
#' @param durations Numeric vector of stimulus durations in seconds.
#' @param nscan Number of scans/time points.
#' @param TR Repetition time in seconds.
#'
#' @return A binary vector of length \code{nscan}, where 1 indicates stimulus
#'   presence and 0 indicates absence.
#'
#' @details
#' Onset and duration values are converted from seconds to scan indices using
#' the supplied repetition time \code{TR}. Specifically, onset times are mapped
#' using \code{ceiling(onsets / TR)}, and durations are mapped using
#' \code{ceiling(durations / TR)}.
#'
#' If a stimulus interval extends beyond \code{nscan}, it is truncated to the
#' scan range.
#'
#' @examples
#' # Single stimulus
#' build_stimulus_vector(
#'   onsets = c(0, 20, 40),
#'   durations = c(10, 10, 10),
#'   nscan = 30,
#'   TR = 2
#' )
#'
#' # Multiple stimulus conditions with rest blocks
#' visual_stim <- build_stimulus_vector(
#'   onsets = c(20, 80),
#'   durations = c(20, 20),
#'   nscan = 60,
#'   TR = 2
#' )
#'
#' auditory_stim <- build_stimulus_vector(
#'   onsets = c(50, 110),
#'   durations = c(20, 10),
#'   nscan = 60,
#'   TR = 2
#' )
#'
#' task_input <- cbind(visual_stim, auditory_stim)
#' colnames(task_input) <- c("Visual", "Auditory")
#'
#' task_input
#'
#' @export
build_stimulus_vector <- function(onsets, durations, nscan, TR) {

  if (!is.numeric(onsets)) {
    stop("`onsets` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(durations)) {
    stop("`durations` must be numeric.", call. = FALSE)
  }

  if (length(onsets) != length(durations)) {
    stop("`onsets` and `durations` must have the same length.", call. = FALSE)
  }

  if (!is.numeric(nscan) || length(nscan) != 1 || nscan <= 0) {
    stop("`nscan` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(TR) || length(TR) != 1 || TR <= 0) {
    stop("`TR` must be a positive number.", call. = FALSE)
  }

  onset_scans <- ceiling(onsets / TR)
  duration_scans <- ceiling(durations / TR)

  stim <- rep(0, nscan)

  for (i in seq_along(onset_scans)) {

    start_idx <- max(1, onset_scans[i])
    end_idx <- min(
      nscan,
      onset_scans[i] + duration_scans[i] - 1
    )

    stim[start_idx:end_idx] <- 1
  }

  return(stim)
}
