test_plot_saving_helper <- function(func, title, func_args) {
  prefix <- paste0(tempdir(), "\\test\\")
  tmp <- example_set[1:5]
  rowData(tmp)$Metabolite_name <- c("Glucose", "Threoline", "5-AVAB", "1/2 acid", "20:0 carbon chain")
  do.call(func, c(list(drop_qcs(tmp), file_path = prefix, format = "emf", title = title), func_args))
  if (is.null(title)) {
    expect_equal(all(list.files(path = prefix) %in% paste0(rownames(tmp), ".emf")), TRUE)
  } else if (title == "Metabolite_name") {
    expect_equal(all(
      list.files(path = prefix) %in% paste0(
        gsub("[:/]", "_", rowData(tmp)$Metabolite_name),".emf")
    ), TRUE)
  } else {
    expect_equal(all(list.files(path = prefix) %in% paste0(rownames(tmp), ".emf")), TRUE)
  }
  unlink(prefix, recursive = TRUE)
}


test_file_extension_helper <- function(fileext) {
  p <- ggplot()
  file <- tempfile(fileext = fileext)
  save_plot(p, file)
  expect_equal(file.exists(file), TRUE)
  unlink(file)
}
