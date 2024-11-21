context("Testing scatter.R plots")

library(notame)

test_that("mz_rt_plot uses correct data", {
  gg <- mz_rt_plot(example_set)
  expect_equal(fData(example_set), gg$data)
})

test_that("mz_rt_plot works with right objects", {
  expect_visible(mz_rt_plot(fData(example_set)))
  expect_visible(mz_rt_plot(example_set))
})

test_that("Low quality features are dropped when they should", {
  flagged <- flag_quality(example_set)
  no_low <- drop_flagged(flagged)
  gg_flagged <- mz_rt_plot(flagged, all_features = FALSE)
  gg_all <- mz_rt_plot(flagged, all_features = TRUE)

  expect_equal(gg_flagged$data, fData(no_low))
  expect_equal(gg_all$data, fData(flagged))
})
