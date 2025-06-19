context("Testing plotting and saving functions")

data(example_set, package = "notame")

# Testing save_plot helper function ----
test_that("Recursive folder creation works", {
  p <- ggplot()
  folder <- paste0(tempdir(), "\\test\\recursive")
  expect_equal(file.exists(folder), FALSE)
  file <- tempfile(tmpdir = folder, fileext = ".pdf")
  save_plot(p, file)
  expect_equal(file.exists(folder), TRUE)
  unlink(folder, recursive = TRUE)
})

test_that("Creating pdf files work", {
  p <- ggplot()
  file <- tempfile(fileext = ".pdf")
  save_plot(p, file)
  expect_equal(file.exists(file), TRUE)
  unlink(file)
})

test_that("Creating emf files work", {
  test_file_extension_helper(".emf")
})

test_that("Creating svg files work", {
  test_file_extension_helper(".svg")
})

test_that("Creating png files work", {
  test_file_extension_helper(".png")
})

test_that("Creating tiff files work", {
  test_file_extension_helper(".tiff")
})

test_that("Giving invalid file format throws error", {
  p <- ggplot()
  file <- tempfile(fileext = ".jpeg")
  expect_error(save_plot(p, file), "is not valid")
  expect_equal(file.exists(file), FALSE)
  unlink(file)
})

# Testing save functions ----
test_that("Batch plots work with and without multiple assays", {
  path <- paste0(tempdir(), "\\test\\batch_plots.pdf")
  ex_set <- example_set
  names(assays(ex_set)) <- "original"
  batch_corrected <- example_set
  assay(batch_corrected, "bcorrected") <- assay(batch_corrected)
  
  # Assay not found with one object with a single assay
  expect_error(save_batch_plots(
    orig = ex_set[1:10], file = path, assay.type1 = c("original"),
    assay.type2 = "asjf"))
  
  # Don't require assay.type if object has only one assay for consistence
  save_batch_plots(
    orig = ex_set[1, ], batch_corrected, file = path, 
    assay.type2 = "bcorrected")
  
  expect_error(save_batch_plots(
    orig = ex_set[1:10], corrected = batch_corrected[1:10],
    file = path, assay.type1 = 1))
})

