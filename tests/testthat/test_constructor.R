context("Testing constructor class")

data(toy_notame_set, package = "notame")

test_that("Joining new columns to feature data preserves data", {
  # With more rows in data.frame
  new_info <- data.frame(Feature_ID = c(rownames(toy_notame_set), "additional"),
                         Feature_number = seq_len(nrow(toy_notame_set) + 1))
  with_new_info <- join_rowData(toy_notame_set, new_info)
  expect_equal(rowData(toy_notame_set), 
               rowData(with_new_info)[, 
                 which(colnames(rowData(with_new_info)) != "Feature_number")])
  
  # With fewer rows in data.frame
  new_info <- data.frame(
    Feature_ID = rownames(toy_notame_set)[1:nrow(toy_notame_set) -1],
    Feature_number = seq_len(nrow(toy_notame_set) - 1))
  with_new_info <- join_rowData(toy_notame_set, new_info)
  expect_equal(rowData(toy_notame_set)[1: nrow(toy_notame_set) - 1, ], 
               rowData(with_new_info)[1:nrow(toy_notame_set) - 1, 
                  which(colnames(rowData(with_new_info)) != "Feature_number")])

  # NAs are added to fill the difference if data.frame has fewer rows
  expect_equal(is.na(rowData(with_new_info)[nrow(toy_notame_set),
                                            "Feature_number"]),
               TRUE)
})
               
test_that("Joining new columns to pheno data preserves object", {
  # With more rows in data.frame
  new_info <- data.frame(
    Sample_ID = colnames(toy_notame_set),
    BMI = stats::runif(ncol(toy_notame_set), 22, 26))
  with_new_info <- join_colData(toy_notame_set, new_info)
  expect_equal(colData(toy_notame_set), 
               colData(with_new_info)[, 
                 which(colnames(colData(with_new_info)) != "BMI")])
  
  # With fewer rows in data.frame
  new_info <- data.frame(
    Sample_ID = colnames(toy_notame_set)[1:ncol(toy_notame_set) -1],
    BMI = seq_len(ncol(toy_notame_set) - 1))
  with_new_info <- join_colData(toy_notame_set, new_info)
  expect_equal(colData(toy_notame_set), colData(with_new_info)[, 
                 which(colnames(colData(with_new_info)) != "BMI")])

  # NAs are added to fill the difference if data.frame has fewer rows
  expect_equal(is.na(colData(with_new_info)[ncol(toy_notame_set), "BMI"]),
               TRUE)
})

