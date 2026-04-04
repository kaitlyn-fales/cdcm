test_that("cdcm_stan_file returns an existing path", {
  path <- cdcm_stan_file()

  expect_type(path, "character")
  expect_length(path, 1)
  expect_true(file.exists(path))
})

test_that("meta_analysis_stan_file returns an existing path", {
  path <- meta_analysis_stan_file()

  expect_type(path, "character")
  expect_length(path, 1)
  expect_true(file.exists(path))
})
