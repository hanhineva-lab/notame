# Check that two SummarizedExperiment objects can be combined
#
# Checks many matching criteria, basically pheno data needs to have similar 
# special columns,
# the number of samples needs to be the same and feature data need to have the
# same columns names. Throws an error if any of the criteria is not fulfilled.
#
# @param x,y SummarizedExperiment objects
.check_match <- function(x, y) {
  # Lots of checks to ensure that everything goes smoothly

  # Amount of samples must be equal
  if (nrow(colData(x)) != nrow(colData(y))) {
    warning("Unequal amount of samples")
  }
  # Resulting feature ID must be unique
  feature_id <- c(rowData(x)$Feature_ID, rowData(y)$Feature_ID)
  if (!.all_unique(feature_id)) {
    stop("Merge would result in duplicated feature ID")
  }
  common <- intersect(colnames(x), colnames(y))
  if (!identical(colData(x)[common, "Injection_order"], 
                 colData(y)[common, "Injection_order"])) {
    stop("Injection orders of common samples are not identical")
  }
  if (!identical(colData(x)$Sample_ID, colData(y)$Sample_ID)) {
    warning("Sample IDs are not identical")
    samples_x <- setdiff(colnames(x), colnames(y))
    samples_y <- setdiff(colnames(y), colnames(x))
    log_text("Merging objects with unequal amounts of samples.")
    log_text("Samples only in first object:")
    log_text(paste0(paste(samples_x, collapse = ", "), "\n"))
    log_text("Samples only in second object:")
    log_text(paste0(paste(samples_y, collapse = ", "), "\n"))
  }


  overlap_cols <- intersect(colnames(colData(x)), colnames(colData(y))) |>
    setdiff(c("Sample_ID", "Injection_order"))

  if (length(overlap_cols)) {
    for (overlap_col in overlap_cols) {
      if (!identical(colData(x)[common, overlap_col], 
                     colData(y)[common, overlap_col])) {
        stop("Columns named ", overlap_col, 
             " in pheno data have different content")
      }
    }
  }

  if (!identical(colnames(rowData(x)), colnames(rowData(y)))) {
    stop("Feature data have different column names")
  }
}

# Merge two SummarizedExperiment objects together
.merge_mode_helper <- function(x, y) {
  # Create dummy injection order if original ones differ
  common <- intersect(colnames(x), colnames(y))
  if (!identical(colData(x)[common, "Injection_order"], 
      colData(y)[common, "Injection_order"])) {
    log_text(paste0("Injection order differs between modes.",
                    "Creating dummy injection order"))
    x_modes <- unique(rowData(x)$Split)
    # Save original injection order for first mode
    if (length(x_modes) == 1) {
      colData(x)[, paste0(x_modes, "_Injection_order")] <- x$Injection_order
    }
    dummy_injection <- as.numeric(-seq_along(colData(x)$Sample_ID))
    names(dummy_injection) <- x$Sample_ID
    colData(x)$Injection_order <- dummy_injection
    # Save original injection order in other modes
    colData(y)[, paste0(unique(rowData(y)$Split), "_Injection_order")] <-
      y$Injection_order
    # Update dummy injection
    y_in_x <- y$Sample_ID %in% x$Sample_ID
    new_io <- seq(from = min(dummy_injection) - 1, length.out = sum(!y_in_x))
    names(new_io) <- y$Sample_ID[!y_in_x]
    dummy_injection <- append(dummy_injection, new_io)

    colData(y)$Injection_order <- dummy_injection[match(colData(y)$Sample_ID,
                                                        names(dummy_injection))]
    log_text("Dummy injection order (row numbers) created")
  }
  # Check that the match is ok
  .check_match(x, y)

  merged_coldata <- dplyr::full_join(as.data.frame(colData(x)), 
                                     as.data.frame(colData(y)),
                                     by = intersect(colnames(colData(x)),
                                                    colnames(colData(y)))) |>
    S4Vectors::DataFrame()
  rownames(merged_coldata) <- merged_coldata$Sample_ID
  merged_rowdata <- rbind(rowData(x), rowData(y)) |>
    S4Vectors::DataFrame()
  if (identical(colnames(assay(x)), colnames(assay(y)))) {
    merged_assay <- rbind(assay(x), assay(y))
  } else {
    merged_assay <- dplyr::bind_rows(as.data.frame(assay(x)),
                                     as.data.frame(assay(y))) |> 
      as.matrix()
    rownames(merged_assay) <- rownames(merged_rowdata)
  }

  merged_object <- SummarizedExperiment(assays = merged_assay, 
                                        colData = merged_coldata,
                                        rowData = merged_rowdata)

  merged_object
}

# Convert objects in ... to a list
# OR if a list is given in the first place, preserve that list
.to_list <- function(...) {
  # Combine the objects to a list
  objects <- list(...)
  # If a list is given in the first place, it should move to top level
  if (length(objects) == 1) {
    if (is(objects[[1]], "list")) {
      objects <- objects[[1]]
    }
  }

  objects
}

#' Merge SummarizedExperiment objects together
#'
#' Merges two or more SummarizedExperiment objects together. Can be used to 
#' merge analytical modes or batches.
#'
#' @param ... \code{
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' objects or a list of objects
#' @param merge what to merge? features is used for combining analytical modes,
#' samples is used for batches
#' @param assay.type character, assay to be used in case of multiple assays. 
#' The same assay needs to be present in all objects to be merged, and the 
#' resultant object contains this single assay.
#'
#' @return A merged SummarizedExperiment object.
#'
#' @details When merging samples, sample IDs that begin with "QC" or "Ref" are 
#' combined so that they have running numbers on them. This means that if both 
#' batches have samples called "QC_1", this will not result in an error,
#' but the sample IDs will be adjusted so that they are unique.
#'
#' @examples
#' # Merge analytical modes
#' data(hilic_neg_sample, hilic_pos_sample, rp_neg_sample, rp_pos_sample)
#' merged <- merge_notame_sets(
#'   hilic_neg_sample, hilic_pos_sample,
#'   rp_neg_sample, rp_pos_sample
#' )
#' # Merge batches
#' batch1 <- toy_notame_set[, toy_notame_set$Batch == 1]
#' batch2 <- toy_notame_set[, toy_notame_set$Batch == 2]
#' merged <- merge_notame_sets(batch1, batch2, merge = "samples")
#'
#' @export
merge_notame_sets <- function(..., merge = c("features", "samples"),
                              assay.type = NULL) {
  merge <- match.arg(merge)
  # Combine the objects to a list
  objects <- .to_list(...)
  # Class check
  if (!all(vapply(objects, class, character(1)) == "SummarizedExperiment")) {
    stop("The arguments should only contain SummarizedExperiment objects")
  }
  # Check assay.type and prepare objects for merge
  from_list <- lapply(objects, .get_from_name, assay.type)
  if (!length(unique(from_list)) == 1) {
    stop("The same assay must be present in all objects or alternatively,
          use objects with a single assay")
  } else {
    from <- unlist(unique(from_list))
  }
  objects <- lapply(objects, function(object) {
    assays(object) <- assays(object)[from]
    object
  })

  # Choose merging function
  if (merge == "features") {
    merge_fun <- .merge_mode_helper
  } else {
    merge_fun <- .merge_batch_helper
  }
  # Merge objects together one by one
  merged <- NULL
  for (object in objects) {
    if (is.null(merged)) {
      merged <- object
    } else {
      merged <- merge_fun(merged, object)
    }
  }
  merged
}

.rowdata_batch_helper <- function(fx, fy) {
  non_identical_cols <- !identical(colnames(fx), colnames(fy))
  if (non_identical_cols) {
    only_x_cols <- setdiff(colnames(fx), colnames(fy))
    only_x <- fx[c("Feature_ID", only_x_cols)]
    fx[only_x_cols] <- NULL
    only_y_cols <- setdiff(colnames(fy), colnames(fx))
    only_y <- fy[c("Feature_ID", only_y_cols)]
    fy[only_y_cols] <- NULL
  }

  # Combine common features: all NAs in fx are replaced by a value from fy
  common_features <- intersect(fx$Feature_ID, fy$Feature_ID)
  for (cf in common_features) {
    na_idx <- is.na(fx[cf, ])
    fx[cf, na_idx] <- fy[cf, na_idx]
  }
  new_features <- setdiff(fy$Feature_ID, fx$Feature_ID)

  merged_rowdata <- rbind(fx, fy[new_features, ])

  if (non_identical_cols) {
    merged_rowdata <-dplyr::left_join(merged_rowdata, only_x, 
                                      by = "Feature_ID") |>
      dplyr::left_join(only_y, by = "Feature_ID")
    rownames(merged_rowdata) <- merged_rowdata$Feature_ID
  }
  merged_rowdata
}


.merge_batch_helper <- function(x, y) {
  merged_coldata <- rbind(colData(x), colData(y))
  if (anyDuplicated(merged_coldata$Sample_ID)) {
    log_text(paste0("Found duplicated sample IDs when merging, ",
                    "renaming QC and Ref samples"))
    qc_idx <- grepl("QC", merged_coldata$Sample_ID)
    merged_coldata$Sample_ID[qc_idx] <- 
      paste0("QC_", seq_len(sum(grepl("QC", merged_coldata$Sample_ID))))
    ref_idx <- grepl("Ref", merged_coldata$Sample_ID)
    merged_coldata$Sample_ID[ref_idx] <- 
      paste0("Ref_", seq_len(sum(grepl("Ref", merged_coldata$Sample_ID))))
  }

  rownames(merged_coldata) <- merged_coldata$Sample_ID
  merged_coldata <- merged_coldata |>
    S4Vectors::DataFrame()

  if (identical(rownames(assay(x)), rownames(assay(y)))) {
    merged_assay <- cbind(assay(x), assay(y))
    colnames(merged_assay) <- rownames(merged_coldata)
  } else {
    merged_assay <- dplyr::bind_rows(as.data.frame(t(assay(x))),
                                     as.data.frame(t(assay(y)))) |> 
     t()
    colnames(merged_assay) <- rownames(merged_coldata)
  }

  merged_rowdata <- .rowdata_batch_helper(rowData(x), rowData(y)) |>
    S4Vectors::DataFrame()

  merged_object <- SummarizedExperiment(assays = merged_assay, 
                                        colData = merged_coldata,
                                        rowData = merged_rowdata)

  merged_object
}