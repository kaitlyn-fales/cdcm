test_that("build_stimulus_vector creates expected binary vector", {
  stim <- build_stimulus_vector(
    onsets = c(0, 20),
    durations = c(10, 10),
    nscan = 20,
    TR = 2
  )

  expected <- c(
    rep(1, 4),
    rep(0, 5),
    rep(1, 5),
    rep(0, 6)
  )

  expect_equal(stim, expected)
})

test_that("build_stimulus_vector supports multiple stimuli via cbind", {
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

  u <- cbind(visual_stim, auditory_stim)
  colnames(u) <- c("Visual", "Auditory")

  expect_equal(dim(u), c(60, 2))

  expect_equal(unname(u[1:9, ]), matrix(0, nrow = 9, ncol = 2))

  expect_equal(unname(u[10:19, "Visual"]), rep(1, 10))
  expect_equal(unname(u[10:19, "Auditory"]), rep(0, 10))

  expect_equal(unname(u[25:34, "Auditory"]), rep(1, 10))

  expect_equal(unname(u[40:49, "Visual"]), rep(1, 10))

  expect_equal(unname(u[50:54, ]), matrix(0, nrow = 5, ncol = 2))
  expect_equal(unname(u[55:59, "Auditory"]), rep(1, 5))
})

test_that("build_stimulus_vector checks input lengths", {
  expect_error(
    build_stimulus_vector(
      onsets = c(0, 20),
      durations = 10,
      nscan = 20,
      TR = 2
    ),
    "same length"
  )
})
