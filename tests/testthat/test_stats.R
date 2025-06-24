# P-value correction ----
test_that("P-value correction works", {
  ps <- data.frame(
    x = letters,
    x_P = stats::runif(26),
    P_x = stats::runif(26),
    y_P = stats::runif(26)
  )

  adj <- .adjust_p_values(ps, flags = rep(NA, 26))

  expect_equal(adj$x_P_FDR, stats::p.adjust(ps$x_P, method = "BH"))
  expect_equal(adj[colnames(ps)], ps)
  expect_equal(ncol(adj), ncol(ps) + 2)

  adj2 <- .adjust_p_values(ps, flags = c(
    "a", "a", "a",
    rep(NA, 23)
  ))

  expect_equal(
    adj2$x_P_FDR,
    c(rep(NA, 3), stats::p.adjust(ps$x_P[4:26], method = "BH"))
  )
})

# Linear model ----
test_that("Linear model works", {
  cd <- combined_data(drop_qcs(example_set))
  lm_fit <- stats::lm(HILIC_neg_259_9623a4_4322 ~ Time,
    data = cd
  )
  smry <- summary(lm_fit)

  # Works for a simple example
  lm_res <- perform_lm(drop_qcs(example_set),
    formula_char = "Feature ~ Time"
  )

  expect_equal(lm_res$Time2_Estimate[1], smry$coefficients[2, 1])
  expect_equal(lm_res$Time2_Std_Error[1], smry$coefficients[2, 2])
  expect_equal(lm_res$Time2_t_value[1], smry$coefficients[2, 3])
  expect_equal(lm_res$Time2_P[1], smry$coefficients[2, 4])
  expect_equal(lm_res$R2[1], smry$r.squared)
  expect_equal(lm_res$Adj_R2[1], smry$adj.r.squared)

  # Works with column with only NA values
  ex_set_na <- drop_qcs(mark_nas(example_set, 0))
  assay(ex_set_na)[1:2, ] <- NA

  lm_res <- perform_lm(ex_set_na,
    formula_char = "Feature ~ Time"
  )
  expect_equal(nrow(lm_res), nrow(assay(example_set)))
  expect_equal(lm_res$Feature_ID, rownames(example_set))
  expect_true(all(is.na(lm_res[1:2, 2:ncol(lm_res)])))

  # FDR correction ignored for flagged compounds
  ex_set_na <- flag_quality(ex_set_na)
  lm_res2 <- perform_lm(ex_set_na, formula_char = "Feature ~ Group")
  flag_idx <- !is.na(flag(ex_set_na))
  expect_true(all(is.na(lm_res2$Group_P_FDR[flag_idx])))
})

test_that("Assay control works (with correlation tests)", {
  ex_set <- example_set[1:10, ]
  names(assays(ex_set)) <- c("original")
  # Don't require assay.type if object has only one assay for consistence
  lm_results <- perform_lm(drop_qcs(ex_set), 
                           formula_char = "Feature ~ Group + Time")
  
  # Assay not found
  lm_results <- expect_error(
    perform_lm(drop_qcs(ex_set), 
               formula_char = "Feature ~ Group + Time",
               assay.type = "nope"))
 
  # Assay not found with multiple assays
  assay(ex_set, "alternative") <- assay(ex_set)
  lm_results <- expect_error(
    perform_lm(drop_qcs(ex_set), 
               formula_char = "Feature ~ Group + Time",
               assay.type = "nope"))
})