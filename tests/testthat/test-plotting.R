test_that("plot.single_ab_validation runs without error", {
  gold <- c("<0.25", "8", "64", ">64")
  test <- c("<0.25", "2", "16", "64")
  val <- compare_mic(gold, test)

  # expect no error (warnings are acceptable)
  expect_error(plot(val), NA)
  expect_error(plot(val, match_axes = FALSE), NA)
})

test_that("plot.multi_ab_validation runs without error", {
  gold <- c("<0.25", "8", "64", ">64", "0.5")
  test <- c("<0.25", "2", "16", "64", "1")
  ab <- c("AMK", "AMK", "CIP", "CIP", "AMK")
  val <- compare_mic(gold, test, ab = ab)

  # default should run (our multi-ab function falls back to single when faceting not requested)
  expect_error(plot(val), NA)

  # explicit faceting should also run (allow warnings)
  expect_error(plot(val, facet_wrap_ncol = 2), NA)
})
