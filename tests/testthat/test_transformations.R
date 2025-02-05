context("Testing transformations")

library(notame)

data("example_set")

test_that("Marking NAs works properly", {
  marked <- mark_nas(example_set, value = 0)
  zero_idx <- assay(example_set) == 0
  na_idx <- is.na(assay(marked))

  expect_equal(zero_idx, na_idx)
})


test_that("RF imputation works as expected", {
  marked <- mark_nas(example_set, value = 0)
  imputed <- impute_rf(marked)

  # Check that all missing values are imputed
  expect_equal(sum(is.na(assay(imputed))), 0)

  # Check that non-missing values are unchanged
  na_idx <- is.na(assay(marked))

  non_na_marked <- assay(marked)[!na_idx]
  non_na_imputed <- assay(imputed)[!na_idx]
  expect_equal(non_na_imputed, non_na_marked)
})

test_that("Flagging works as expected", {
  single_change <- matrix(1111,
    nrow = 1, ncol = 1,
    dimnames = list(
      "HILIC_neg_108_1065a2_6121",
      "Demo_1"
    )
  )
  mrg <- merge_assay(example_set, single_change)
  expect_equal(assay(mrg)[2, 1], 1111)

  block_change <- matrix(stats::runif(9),
    nrow = 3, ncol = 3,
    dimnames = list(
      rownames(assay(example_set))[4:6],
      colnames(assay(example_set))[5:7]
    )
  )
  mrg2 <- merge_assay(example_set, block_change)
  expect_true(all(assay(mrg2)[4:6, 5:7] == block_change))

  scattered_change <- matrix(stats::runif(9),
    nrow = 3, ncol = 3,
    dimnames = list(
      rownames(assay(example_set))[c(1, 5, 7)],
      colnames(assay(example_set))[c(2, 5, 8)]
    )
  )
  mrg3 <- merge_assay(example_set, scattered_change)
  expect_true(all(assay(mrg3)[c(1, 5, 7), c(2, 5, 8)] == scattered_change))

  wrong_names <- matrix(stats::runif(9), nrow = 3)
  expect_error(merge_assay(example_set, wrong_names), "Column names")

  wrong_names2 <- matrix(stats::runif(9),
    nrow = 3,
    dimnames = list(
      letters[1:3],
      colnames(assay(example_set))[5:7]
    )
  )
  expect_error(merge_assay(example_set, wrong_names2), "Row names")
})


test_that("Flagged compounds are not imputed", {
  marked <- mark_nas(example_set, 0)
  flag(marked)[c(1, 4, 6)] <- "Flagged"

  imputed <- impute_rf(marked)
  nas <- apply(assay(imputed), 1, prop_na)

  expect_true(all(nas[c(1, 4, 6)] > 0))
  expect_true(all(nas[-c(1, 4, 6)] == 0))

  # All feature parameter test
  imputed <- impute_rf(marked, all_features = TRUE)
  nas <- apply(assay(imputed), 1, prop_na)
  expect_true(all(nas == 0))
})

test_that("Inverse normalization works as expected", {
  marked <- mark_nas(example_set, value = 0)
  imputed <- impute_rf(marked)
  normalized <- inverse_normalize(imputed)

  # Ranks should be identical
  orig_ranks <- assay(imputed) %>%
    apply(1, rank) %>%
    t()
  norm_ranks <- assay(normalized) %>%
    apply(1, rank) %>%
    t()
  expect_identical(norm_ranks, orig_ranks)

  # Check that all the rows follow standard normal distribution
  # Zero means up to 0.1 accuracy
  abs_means <- abs(rowMeans(assay(normalized)))
  expect_true(all(abs_means < 0.1))
  # Unit variance
  sd_diffs <- abs(apply(assay(normalized), 1, sd) - 1)
  expect_true(all(sd_diffs < 0.1))

  # Shapiro-Wilk normality tests are non-significant
  shapiro_p <- assay(normalized) %>%
    apply(1, function(x) {
      shapiro.test(x)$p.value
    })
  expect_true(all(shapiro_p > 0.9))
})

test_that("Simple imputation works as expected", {
  marked <- mark_nas(example_set, value = 0)
  imputed <- impute_simple(marked, value = 0)

  # Check that all missing values are imputed
  expect_equal(sum(is.na(assay(imputed))), 0)

  # Check that non-missing values are unchanged
  na_idx <- is.na(assay(marked))

  non_na_marked <- assay(marked)[!na_idx]
  non_na_imputed <- assay(imputed)[!na_idx]
  expect_equal(non_na_imputed, non_na_marked)
})

test_that("Simple imputation with only one feature works", {
  marked <- example_set
  assay(marked)[2, 5] <- NA
  imputed <- impute_simple(marked, value = 0)
  #
  expect_equal(assay(imputed)[2, 5], 0)
})

test_that("PQN normalization works correctly using median of QC samples as reference", {
  data <- assay(example_set)
  # Calculate the median of QC samples
  reference <- apply(data[, example_set$QC == "QC"], 1, finite_median)
  # do the normalization
  quotient <- data / reference
  quotient_md <- apply(quotient, 2, finite_median)
  pqn_data <- t(t(data) / quotient_md)

  expect_equal(pqn_data, assay(pqn_normalization(example_set)))
})

test_that("PQN normalization works correctly using median of all samples as reference", {
  data <- assay(example_set)
  # Calculate the median of all samples
  reference <- apply(data, 1, finite_median)
  # do the normalization
  quotient <- data / reference
  quotient_md <- apply(quotient, 2, finite_median)
  pqn_data <- t(t(data) / quotient_md)

  expect_equal(pqn_data, assay(pqn_normalization(example_set, ref = "all")))
})

test_that("PQN normalization works correctly using mean of QC samples as reference", {
  data <- assay(example_set)
  # Calculate the mean of QC samples
  reference <- apply(data[, example_set$QC == "QC"], 1, finite_mean)
  # do the normalization
  quotient <- data / reference
  quotient_md <- apply(quotient, 2, finite_median)
  pqn_data <- t(t(data) / quotient_md)

  expect_equal(pqn_data, assay(pqn_normalization(example_set, method = "mean")))
})

test_that("PQN normalization works correctly using mean of all samples as reference", {
  data <- assay(example_set)
  # Calculate the mean of all samples
  reference <- apply(data, 1, finite_mean)
  # do the normalization
  quotient <- data / reference
  quotient_md <- apply(quotient, 2, finite_median)
  pqn_data <- t(t(data) / quotient_md)

  expect_equal(pqn_data, assay(pqn_normalization(example_set, ref = "all", 
                                                 method = "mean")))
})

test_that("PQN normalization works with flagged features", {
  flagged_mset <- flag_quality(example_set)
  ref_data <- assay(drop_flagged(flagged_mset))
  reference_spectrum <- apply(ref_data[, example_set$QC == "QC"], 1, finite_median)
  quotient <- ref_data / reference_spectrum
  quotient_md <- apply(quotient, 2, finite_median)

  data <- assay(flagged_mset)
  pqn_data <- t(t(data) / quotient_md)

  pqn_mset <- pqn_normalization(flagged_mset)
  expect_equal(pqn_data, assay(pqn_mset))
})

test_that("Assay control works for transformations (with drift correction)", {
  ex_set <- example_set[1:10, ]
  # New assays can be added to objects with one assay
  corrected <- correct_drift(ex_set, name = "corrected")
  # Assay not found with a single unnamed assay
  expect_error(correct_drift(ex_set, assay.type = "original"))
  # From and to can't have identical names
  names(assays(ex_set)) <- "original"
  expect_error(correct_drift(ex_set, 
    assay.type = "original", name = "original"))
  # Assay not found with multiple assays
  assay(ex_set, "alternative") <- assay(ex_set)
  expect_error(correct_drift(ex_set, 
    assay.type = "nope", name = "corrected"))
  # Only one name can be supplied to assay.type
  expect_error(correct_drift(ex_set,
    assay.type = c("original", "lel"), name = "corrected"))
  # Only one name can be supplied to transformed assay
  expect_error(correct_drift(ex_set,
    assay.type = "original", name = c("A", "B")))
  # MetaboSet works
  ms_set <- as(ex_set, "MetaboSet")
  ms_corrected <- correct_drift(ms_set)
  # assay.type throws error with MetaboSet
  expect_error(correct_drift(ms_set, assay.type = "nope", name = "corrected"))
  # # name throws error with MetaboSet
  # expect_error(correct_drift(ms_set, assay.type = "nope", name = "corrected"))
})