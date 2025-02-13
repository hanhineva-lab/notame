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

test_that("Subject line plots are saved without title", {
  test_plot_saving_helper(example_set, save_subject_line_plots, title = NULL, 
                          func_args = c(x = "Time", id = "Subject_ID", 
                                        color = "Group"))
})

test_that("Subject line plot naming works", {
  test_plot_saving_helper(example_set, save_subject_line_plots, 
                          title = "Metabolite_name",   
                          func_args = c(x = "Time", id = "Subject_ID"))
})

test_that("Group boxplots are saved without title", {
  test_plot_saving_helper(example_set, save_group_boxplots, title = NULL,
                          func_args = c(x= "Group", color = "Group"))
})

test_that("Group boxplot naming works", {
  test_plot_saving_helper(example_set, save_group_boxplots, 
                          title = "Feature_ID", 
                          func_args = c(x = "Group", color = "Group"))
})

test_that("Beeswarm plots are saved without title", {
  test_plot_saving_helper(example_set, save_beeswarm_plots, title = NULL,
                          func_args = c(x = "Group", color = "Group"))
})

test_that("Beeswarm plot naming works", {
  test_plot_saving_helper(example_set, save_beeswarm_plots,
                          title = "Metabolite_name",
                          func_args = c(x = "Group", color = "Group"))
})

test_that("Scatter plots are saved without title", {
  test_plot_saving_helper(example_set, save_scatter_plots, title = NULL,
                          func_args = c(x = "Injection_order", color = "Group"))
})

test_that("Scatter plot naming works", {
  test_plot_saving_helper(example_set, save_scatter_plots, title = "Feature_ID",
                          func_args = c(x = "Injection_order", color = "Group"))
})

test_that("Group lineplots are saved without title", {
  test_plot_saving_helper(example_set, save_group_lineplots, title = NULL,
                          func_args = c(x = "Time", group = "Group"))
})

test_that("Group lineplot naming works", {
  test_plot_saving_helper(example_set, save_group_lineplots,
                          title = "Feature_ID",
                          func_args = c(x = "Time", group = "Group"))
})

test_that("Batch plots work with and without multiple assays", {
  path <- paste0(tempdir(), "\\test\\batch_plots.pdf")
  ex_set <- example_set
  names(assays(ex_set)) <- c("original")

  replicates <- list(which(ex_set$QC == "QC"))
  batch_corrected <- ruvs_qc(example_set, replicates = replicates,
    name = "bcorrected")
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
      
  # MetaboSet works
  ms_set <- as(ex_set, "MetaboSet")
  ms_batch_corrected <- as(batch_corrected, "MetaboSet")
  save_batch_plots(
    orig = ms_set[1:10], corrected = ms_batch_corrected[1:10],
    file = path)

  # assay.type throws error with MetaboSet
  expect_error(save_batch_plots(
    orig = ms_set[1:10], corrected = ms_batch_corrected[1:10],
    file = path, assay.type1 = "nope"))
})

