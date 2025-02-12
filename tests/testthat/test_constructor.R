context("Testing constructor class")

data(example_set, package = "notame")

test_that("Changing feature names updates all occurrences", {
  set <- as(example_set, "MetaboSet")
  rownames(set) <- paste0(letters[1:20], 1:80)
  expect_equal(featureNames(set), featureNames(assayData(set))) 
  expect_equal(featureNames(set), featureNames(featureData(set)))
  expect_equal(featureNames(set), fData(set)$Feature_ID)
})

test_that("Changing feature names only works if valid names are given", {
  set <- as(example_set, "MetaboSet")
  expect_error(featureNames(set) <- NULL)
  # Numbers are not allowed
  expect_error(featureNames(set) <- 1:80)
  # Number of new names should equal number of rows
  expect_error(featureNames(set) <- letters[1:15])
  # Duplicates are not allowed
  names <- paste0(letters[1:80])
  expect_warning(
    expect_error(featureNames(set) <- names),
    "non-unique value when setting"
  )
  # Names are not allowed to start with numbers
  names[2] <- "2a"
  expect_error(featureNames(set) <- names)
  # NAs are not allowed
  names[2] <- NA
  expect_error(featureNames(set) <- names)
})

test_that("Changing sample names updates all occurrences", {
  set <- as(example_set, "MetaboSet")
  sampleNames(set) <- tolower(sampleNames(set))
  expect_equal(sampleNames(set), sampleNames(assayData(set)))
  expect_equal(sampleNames(set), sampleNames(protocolData(set)))
  expect_equal(sampleNames(set), sampleNames(phenoData(set)))
  expect_equal(sampleNames(set), pData(set)$Sample_ID)
})

test_that("Changing sample names only works if valid names are given", {
  set <- as(example_set, "MetaboSet")
  names <- sampleNames(set)
  expect_error(sampleNames(set) <- NULL)
  # Number of new names should equal number of columns
  expect_error(sampleNames(set) <- names[1:29])
  # Duplicates are not allowed
  names[2] <- names[1]
  expect_warning(
    expect_error(sampleNames(set) <- names),
    "non-unique value"
  )
  # NAs are not allowed
  names[2] <- NA
  expect_error(sampleNames(set) <- names)
})

test_that("Joining new columns to feature data preserves data", {
  # With more rows in data.frame
  new_info <- data.frame(Feature_ID = c(rownames(example_set), "additional"),
                         Feature_number = seq_len(nrow(example_set) + 1))
  with_new_info <- join_rowData(example_set, new_info)
  expect_equal(rowData(example_set), 
               rowData(with_new_info)[, 
                 which(colnames(rowData(with_new_info)) != "Feature_number")])
  
  # With fewer rows in data.frame
  new_info <- data.frame(
    Feature_ID = rownames(example_set)[1:nrow(example_set) -1],
    Feature_number = seq_len(nrow(example_set) - 1))
  with_new_info <- join_rowData(example_set, new_info)
  expect_equal(rowData(example_set)[1: nrow(example_set) - 1, ], 
               rowData(with_new_info)[1:nrow(example_set) - 1, 
                  which(colnames(rowData(with_new_info)) != "Feature_number")])

  # NAs are added to fill the difference if data.frame has fewer rows
  expect_equal(is.na(rowData(with_new_info)[nrow(example_set),
                                            "Feature_number"]),
               TRUE)
})
               
test_that("Joining new columns to pheno data preserves object", {
  # With more rows in data.frame
  new_info <- data.frame(
    Sample_ID = colnames(example_set),
    BMI = stats::runif(ncol(example_set), 22, 26))
  with_new_info <- join_colData(example_set, new_info)
  expect_equal(colData(example_set), 
               colData(with_new_info)[, 
                 which(colnames(colData(with_new_info)) != "BMI")])
  
  # With fewer rows in data.frame
  new_info <- data.frame(
    Sample_ID = colnames(example_set)[1:ncol(example_set) -1],
    BMI = seq_len(ncol(example_set) - 1))
  with_new_info <- join_colData(example_set, new_info)
  expect_equal(colData(example_set), colData(with_new_info)[, 
                 which(colnames(colData(with_new_info)) != "BMI")])

  # NAs are added to fill the difference if data.frame has fewer rows
  expect_equal(is.na(colData(with_new_info)[ncol(example_set), "BMI"]),
               TRUE)
})

