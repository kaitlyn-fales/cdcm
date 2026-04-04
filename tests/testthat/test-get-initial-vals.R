test_that("get_initial_vals returns backup values when pathfinder fails", {
  fake_mod <- list(
    pathfinder = function(...) {
      stop("simulated pathfinder failure")
    }
  )

  toy_data <- list(
    m = 2,
    d_A = 3,
    d_B = 1,
    d_C = 1
  )

  expect_message(
    result <- get_initial_vals(fake_mod, toy_data),
    "Pathfinder failed"
  )

  expect_type(result, "list")
  expect_length(result, 1)

  expect_equal(result[[1]]$sigma, rep(1, 2))
  expect_equal(result[[1]]$nu_A, rep(0, 3))
  expect_equal(result[[1]]$nu_B, rep(0, 1))
  expect_equal(result[[1]]$nu_C, rep(0, 1))
  expect_equal(result[[1]]$z0, rep(0.1, 2))
  expect_equal(result[[1]]$beta, rep(0, 2))
})

test_that("get_initial_vals fallback ignores custom pathfinder_init and uses backup defaults", {
  fake_mod <- list(
    pathfinder = function(...) {
      stop("simulated pathfinder failure")
    }
  )

  toy_data <- list(
    m = 2,
    d_A = 2,
    d_B = 1,
    d_C = 1
  )

  custom_init <- list(
    sigma = c(9, 9),
    nu_A = c(8, 8),
    nu_B = 7,
    nu_C = 6,
    z0 = c(5, 5),
    beta = c(4, 4)
  )

  result <- get_initial_vals(fake_mod, toy_data, pathfinder_init = custom_init)

  expect_equal(result[[1]]$sigma, c(1, 1))
  expect_equal(result[[1]]$nu_A, c(0, 0))
  expect_equal(result[[1]]$z0, c(0.1, 0.1))
})

test_that("get_initial_vals errors when required data components are missing", {
  fake_mod <- list(
    pathfinder = function(...) {
      stop("simulated pathfinder failure")
    }
  )

  bad_data <- list(
    m = 2,
    d_A = 3,
    d_B = 1
  )

  expect_error(
    get_initial_vals(fake_mod, bad_data),
    "missing required component"
  )
})
