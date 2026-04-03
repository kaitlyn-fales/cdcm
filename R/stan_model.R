# Utility function for stan files
stan_file <- function(name) {
  system.file("stan", name, package = "cdcm")
}

#' Path to the canonical DCM Stan program
#'
#' Returns the file path to the CDCM
#' Stan program bundled with the \code{cdcm} package.
#'
#' This function is primarily intended for use with
#' \code{cmdstanr::cmdstan_model()} when compiling the model.
#'
#' @return A character string giving the full path to the
#'   \code{canonical_dcm.stan} file.
#'
#' @examples
#' \dontrun{
#' stan_file <- cdcm_stan_file()
#' mod <- cmdstanr::cmdstan_model(stan_file)
#' }
#'
#' @export
cdcm_stan_file <- function() {
  stan_file("canonical_dcm.stan")
}

#' Path to the group-level meta-analysis Stan program
#'
#' Returns the file path to the group-level meta analysis model
#' Stan program bundled with the \code{cdcm} package.
#'
#' This function is primarily intended for use with
#' \code{cmdstanr::cmdstan_model()} when compiling the model.
#'
#' @return A character string giving the full path to the
#'   \code{meta_analysis.stan} file.
#'
#' @examples
#' \dontrun{
#' stan_file <- met_analysis_stan_file()
#' mod <- cmdstanr::cmdstan_model(stan_file)
#' }
#'
#' @export
meta_analysis_stan_file <- function() {
  stan_file("meta_analysis.stan")
}

#' Compile the CDCM Stan model
#'
#' Compiles the CDCM Stan program bundled
#' with the \code{cdcm} package using \code{cmdstanr}.
#'
#' This function provides a convenient interface for compiling the model
#' without requiring users to manually manage file paths to the Stan program.
#'
#' @param force_recompile Logical; if \code{TRUE}, forces recompilation of the
#'   Stan model even if a compiled executable already exists.
#' @param quiet Logical; passed to \code{cmdstanr::cmdstan_model()} to control
#'   console output during compilation.
#' @param ... Additional arguments passed to \code{cmdstanr::cmdstan_model()}.
#'
#' @return A \code{CmdStanModel} object.
#'
#' @details
#' This function requires a working CmdStan installation (version 2.35.0 or
#' newer) configured for use with \code{cmdstanr}. If CmdStan is not installed
#' or configured, an error will be thrown.
#'
#' @examples
#' \dontrun{
#' mod <- compile_cdcm()
#'
#' # Force recompilation if needed
#' mod <- compile_cdcm(force_recompile = TRUE)
#' }
#'
#' @export
compile_cdcm <- function(force_recompile = FALSE, quiet = TRUE, ...) {
  stan_path <- cdcm_stan_file()

  if (stan_path == "") {
    stop("Could not find 'canonical_dcm.stan' in the installed cdcm package.")
  }

  if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
    stop(
      "CmdStan is not configured for cmdstanr. ",
      "Install/configure CmdStan, then try again."
    )
  }

  cmdstanr::cmdstan_model(
    stan_file = stan_path,
    force_recompile = force_recompile,
    quiet = quiet,
    ...
  )
}

#' Compile the group-level meta-analysis Stan model
#'
#' Compiles the meta-analysis Stan program bundled with the \code{cdcm}
#' package using \code{cmdstanr}.
#'
#' @param force_recompile Logical; if \code{TRUE}, forces recompilation of the
#'   Stan model even if a compiled executable already exists.
#' @param quiet Logical; passed to \code{cmdstanr::cmdstan_model()} to control
#'   console output during compilation.
#' @param ... Additional arguments passed to \code{cmdstanr::cmdstan_model()}.
#'
#' @return A \code{CmdStanModel} object.
#' @export
compile_meta_analysis <- function(force_recompile = FALSE, quiet = TRUE, ...) {
  stan_path <- meta_analysis_stan_file()

  if (stan_path == "") {
    stop("Could not find 'meta_analysis.stan' in the installed cdcm package.")
  }

  if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
    stop(
      "CmdStan is not configured for cmdstanr. ",
      "Install/configure CmdStan, then try again."
    )
  }

  cmdstanr::cmdstan_model(
    stan_file = stan_path,
    force_recompile = force_recompile,
    quiet = quiet,
    ...
  )
}
