# Utility function for stan files
stan_file <- function(name) {
  system.file("stan", name, package = "cdcm")
}
