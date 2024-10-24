library(dplyr)
devtools::load_all()
# Create a toy dataset for use in examples

set.seed(38)
n_features <- 20
modes <- c("HILIC_neg", "HILIC_pos", "RP_neg", "RP_pos")
feature_data_modes <- list()
for (mode in modes) {
  feature_data <- data.frame(
    Split = mode,
    Alignment = seq_len(n_features),
    Average_Mz = stats::runif(n_features, 100, 500),
    Average_Rt_min = stats::runif(n_features, 0.5, 8),
    Column = strsplit(mode, "_")[[1]][1], 
    Ion_mode = strsplit(mode, "_")[[1]][2],
    stringsAsFactors = FALSE)
  
  # Create Feature ID
  round_mz <- as.numeric(feature_data$Average_Mz) %>%
    round(digits = 4) %>%
    as.character() %>%
    gsub("[.]", "_", .)
  round_rt <- as.numeric(feature_data$Average_Rt_min) %>%
    round(digits = 4) %>%
    as.character() %>%
    gsub("[.]", "_", .)
  feature_data$Feature_ID <- paste0(mode, "_", round_mz, "a", round_rt)
  feature_data <- dplyr::select(feature_data, "Feature_ID", everything())

  rownames(feature_data) <- feature_data$Feature_ID
  feature_data_modes[[mode]] <- feature_data
}

# Pheno data
n_samples <- 50
# Make batch as if paired samples were measured in different batches
qc_idx <- seq(1, n_samples / 2, length.out = 5) %>% round()
subject_ids <- as.character(seq_len(n_samples / 2))
group <- sample(LETTERS[1:2], n_samples / 2, replace = TRUE)
time <- as.character(rep(c(1, 2), each = n_samples / 2))
batch <- rep(c(1, 2), each = n_samples / 2)
group[qc_idx] <- "QC"
subject_ids[qc_idx] <- "QC"
time[c(qc_idx, qc_idx + n_samples / 2)] <- "QC"
qc <- ifelse(group == "QC", "QC", "Sample")
HILIC_neg_Datafile <- paste0("HILIC_neg_", seq_len(n_samples))
HILIC_pos_Datafile <- paste0("HILIC_pos_", seq_len(n_samples))
RP_neg_Datafile <- paste0("RP_neg_", seq_len(n_samples))
RP_pos_Datafile <- paste0("RP_pos_", seq_len(n_samples))
pheno_data <- data.frame(Injection_order = seq_len(n_samples),
                         Sample_ID = paste0("Demo_", seq_len(n_samples)),
                         Subject_ID = subject_ids, Group = factor(group), 
                         QC = factor(qc), Time = factor(time), Batch = batch,
                         HILIC_neg_Datafile, HILIC_pos_Datafile,
                         RP_neg_Datafile, RP_pos_Datafile)

rownames(pheno_data) <- pheno_data$Sample_ID


# Assay data

# Random means for each feature per batch/time point, higher overall abundance 
# for batch 2
set.seed(2024)
means_b1 <- stats::runif(n_features * length(modes), 3000, 33000)
set.seed(2024)
means_b2 <- stats::runif(n_features * length(modes), 4000, 44000)

# Normally distributed data around the mean
assay_data_b1 <- t(sapply(means_b1, function(x) {
  rnorm(n_samples/2, x, 0.3 * x)
}))

## Less variance for second batch with higher overall abundances
assay_data_b2 <- t(sapply(means_b2, function(x) {
  rnorm(n_samples/2, x, 0.2 * x)
}))

# Get feature means of biological samples (pooling)
qc_boolean <- pheno_data$QC[seq_len(n_samples / 2)] == "QC"
means_qc_b1 <- apply(assay_data_b1[, qc_boolean], FUN = mean, MARGIN = 1)
means_qc_b2 <- apply(assay_data_b2[, qc_boolean], FUN = mean, MARGIN = 1)

# Less variance in QC samples
assay_data_b1[, qc_boolean] <- t(sapply(means_qc_b1, function(x) {
  rnorm(sum(qc_boolean), x, 0.1 * x)
}))
assay_data_b2[, qc_boolean] <- t(sapply(means_qc_b2, function(x) {
  rnorm(sum(qc_boolean), x, 0.1 * x)
}))

# Add drift effect to the data
coefs <- stats::runif(n_features * length(modes), 0.4, 0.9) * 
  sample(c(-1, 1), n_features * length(modes), replace = TRUE)

## Add drift for batch 1
for (i in seq_len(nrow(assay_data_b1))) {
  set.seed(2024)
  if (rnorm(1) > 0) {
    assay_data_b1[i, ] <- assay_data_b1[i, ] + means_qc_b1[i] * coefs[i] * log(pheno_data$Injection_order[seq_len(n_samples/2)]) * 0.1
  } else {
    assay_data_b1[i, ] <- assay_data_b1[i, ] + means_qc_b1[i] * coefs[i] * 0.1 * pheno_data$Injection_order[seq_len(n_samples/2)]
  }
}

## Add drift for batch 2
for (i in seq_len(nrow(assay_data_b2))) {
  set.seed(2024)
  if (rnorm(1) > 0) {
    assay_data_b2[i, ] <- assay_data_b2[i, ] + means_qc_b2[i] * coefs[i] * log(pheno_data$Injection_order[seq_len(n_samples/2)]) * 0.1
  } else {
    assay_data_b2[i, ] <- assay_data_b2[i, ] + means_qc_b2[i] * coefs[i] * 0.1 * pheno_data$Injection_order[seq_len(n_samples/2)]
  }
}

## Add random noise for batch 1
for (i in seq_len(nrow(assay_data_b1))) {
  set.seed(2024)
  assay_data_b1[i, ] <- assay_data_b1[i, ] + rnorm(n_samples / 2, 0, 0.02 * means_qc_b1[i])
}

## Add random noise for batch 2
for (i in seq_len(nrow(assay_data_b2))) {
  set.seed(2024)
  assay_data_b2[i, ] <- assay_data_b2[i, ] + rnorm(n_samples / 2, 0, 0.02 * means_qc_b2[i])
}
assay_data <- cbind(assay_data_b1, assay_data_b2)

## Convert to positive values (artifact from the logic above)
assay_data <- abs(assay_data)

# Set random indexes to zero (for missing values)
n_missing <- 500
row_zeros <- sample(seq_len(nrow(assay_data)), n_missing, replace = TRUE)
col_zeros <- sample(seq_len(ncol(assay_data)), n_missing, replace = TRUE)
for (i in seq_len(n_missing)) {
  assay_data[row_zeros[i], col_zeros[i]] <- 0
}  

# Set dimension names
rownames(assay_data) <- unlist(lapply(feature_data_modes, rownames))
colnames(assay_data) <- rownames(pheno_data)

feature_data <- bind_rows(feature_data_modes)
  
# Construct objects, with all modes and separately
example_set <- construct_metabosets(
  exprs = assay_data, pheno_data = pheno_data, feature_data = feature_data,
  group_col = "Group", time_col = "Time", subject_col = "Subject_ID",
  split_data = FALSE)

hilic_neg_sample <- construct_metabosets(
  exprs = assay_data[feature_data$Split == "HILIC_neg", ],
  pheno_data = pheno_data, 
  feature_data = feature_data[feature_data$Split == "HILIC_neg", ],
  group_col = "Group", time_col = "Time", subject_col = "Subject_ID",
  split_data = FALSE)

hilic_pos_sample <- construct_metabosets(
  exprs = assay_data[feature_data$Split == "HILIC_pos", ],
  pheno_data = pheno_data, 
  feature_data = feature_data[feature_data$Split == "HILIC_pos", ],
  group_col = "Group", time_col = "Time", subject_col = "Subject_ID", 
  split_data = FALSE)
  
rp_neg_sample <- construct_metabosets(
  exprs = assay_data[feature_data$Split == "RP_neg", ],
  pheno_data = pheno_data, 
  feature_data = feature_data[feature_data$Split == "RP_neg", ],
  group_col = "Group", time_col = "Time", subject_col = "Subject_ID",
  split_data = FALSE)
  
rp_pos_sample <- construct_metabosets(
  exprs = assay_data[feature_data$Split == "RP_pos", ],
  pheno_data = pheno_data, 
  feature_data = feature_data[feature_data$Split == "RP_pos", ],
  group_col = "Group", time_col = "Time", subject_col = "Subject_ID",
  split_data = FALSE)

usethis::use_data(example_set, overwrite = TRUE)
usethis::use_data(hilic_neg_sample, overwrite = TRUE)
usethis::use_data(hilic_pos_sample, overwrite = TRUE)
usethis::use_data(rp_neg_sample, overwrite = TRUE)
usethis::use_data(rp_pos_sample, overwrite = TRUE)

