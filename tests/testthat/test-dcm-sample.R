test_that("dcm_sample errors when ess_check is missing", {
  fake_mod <- list(
    sample = function(...) stop("should not get here")
  )

  expect_error(
    dcm_sample(
      mod = fake_mod,
      data = list(),
      inits_list = list(),
      output_dir = tempdir(),
      basename = "test"
    ),
    "ess_check"
  )
})

test_that("dcm_sample errors when output_dir does not exist", {
  fake_mod <- list(
    sample = function(...) stop("should not get here")
  )

  expect_error(
    dcm_sample(
      mod = fake_mod,
      data = list(),
      inits_list = list(),
      output_dir = file.path(tempdir(), "does-not-exist"),
      basename = "test",
      ess_check = 100
    ),
    "existing directory"
  )
})
