context("Testing merging logic")

test_that("Flag column is combined correctly when merging batches", {
    se <- toy_notame_set
    se$Batch <- sort(rep(1:3, length.out = ncol(se)))
    batches <- sapply(unique(se$Batch), \(x) se[, se$Batch == x])
    # Only one flag
    flag(batches[[1]])[1:3] <- "A"
    flag(batches[[2]])[2:3] <- "A"
    merged <- merge_notame_sets(batches, merge = "samples")
    expected_flags <- c(rep("A", 3), rep(NA, nrow(merged) - 3))
    expect_equal(flag(merged), expected_flags)
    # Different flags
    # 1: any match
    # 2: tie, so we take the first one
    # 3: most common is B
    # 4: rest are NA
    flag(batches[[1]])[1:3] <- "A"
    flag(batches[[2]])[2:3] <- "B"
    flag(batches[[3]])[3] <- "B"
    merged <- merge_notame_sets(batches, merge = "samples")
    expected_flags <- c("A", "A", "B", rep(NA, nrow(merged) - 3))
    expect_equal(flag(merged), expected_flags)
})
