# Load necessary libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(dplyr)
library(permute)
library(data.table)
library(fields)
library(pROC)

# Set working directory
script_dir <- "./modeling/ppp"
setwd(script_dir)

# Load custom functions
source("calc_mark_prob.R")
source("calc_dist_optimized.R")

# Parse command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]
tissue <- args[2]
intensity_type <- args[3]
cluster_num <- as.numeric(args[4])
radius <- as.numeric(args[5])
hr <- as.numeric(args[6])
resample_percent <- as.numeric(args[7])
random_seed <- as.numeric(args[8])

# Define TMC based on tissue type
TMC <- ifelse(tissue %in% c("SI", "LI"), "Stanford", "Florida")

# Define data directory

# Create matrices for radius and hr
r_matrix <- matrix(rep(radius, cluster_num^2), nrow = cluster_num)
hr_matrix <- matrix(rep(hr, cluster_num^2), nrow = cluster_num)

# Initialize lists for deviance calculations
deviance_train_list <- c()
deviance_val_list <- c()
deviance_test_list <- c()

# Load model files
load(file.path(data_dir, TMC, tissue, sprintf("fmla_%d_100-%d_%d_%s_across3_self_quad_d_no_between_dummy.Rda", 
                                              cluster_num, radius, hr, intensity_type)))
load(file.path(data_dir, TMC, tissue, sprintf("family_%d_100-%d_%d_%s_across3_self_quad_d_no_between_dummy.Rda", 
                                              cluster_num, radius, hr, intensity_type)))
load(file.path(data_dir, TMC, tissue, sprintf("coef_%d_100-%d_%d_%s_across3_self_quad_d_no_between_dummy.Rda", 
                                              cluster_num, radius, hr, intensity_type)))

# Extract coefficients for interaction terms
coef_names <- names(coef)
G <- lapply(0:(cluster_num-1), function(i) {
  unlist(lapply(seq(100, radius, 100), function(r) {
    unlist(lapply(0:(cluster_num-1), function(j) {
      interaction_name <- if (i <= j) {
        sprintf("X%d x X%d x %d", i, j, r)
      } else {
        sprintf("X%d x X%d x %d", j, i, r)
      }
      coef[grep(interaction_name, coef_names)]
    }))
  }))
})

# Compute intensity coefficients B
B <- c(coef[1], sapply(2:cluster_num, function(i) coef[i] + coef[1]))

# Define possible cell types
cell_types <- 0:4

# Load simulated pattern data
simulated_pattern_path <- file.path(data_dir, TMC, tissue, "random_resample", 
  sprintf("random_simulated_pattern_%d_500_%d_%s_across3_resample5_percentage_%d_seed_%d_self_quad_d_no_between_dummy.Rda", 
          cluster_num, hr, intensity_type, resample_percent, random_seed))
load(simulated_pattern_path)

# Convert cell pattern to data frame
cell_pattern_df <- data.frame(marked_simulated_pattern[[1]])
cell_pattern_window <<- marked_simulated_pattern[[1]]$window

# Predict cell types and intensities
marks <- c()
intensities <- matrix(nrow = 0, ncol = length(B))

for (point_idx in 1:nrow(cell_pattern_df)) {
  point_current <- cell_pattern_df[point_idx, ]
  counts_current <- t(apply(point_current[, 1:2], 1, calc_dist_optimized, cell_pattern_df, cell_pattern_window))
  mark_current <- apply(counts_current, 1, calc_mark_prob, B, G, cell_types)
  marks <- c(marks, mark_current[[1]][[1]])
  intensities <- rbind(intensities, mark_current[[1]][[2]])
}

# Save predicted cell types and intensities
pred_cell_pattern_df <- data.frame(x = cell_pattern_df$x, y = cell_pattern_df$y, marks = marks)
pred_intensities <- data.frame(intensities)

output_base_path <- file.path(data_dir, TMC, tissue, "random_resample")
filename_pred_cell <- file.path(output_base_path, sprintf("pred_cell_type_%d_500_%d_%s_across3_resample5_percentage_%d_seed_%d_self_quad_d_no_between_dummy.Rda",
                                                          cluster_num, hr, intensity_type, resample_percent, random_seed))
filename_pred_intensity <- file.path(output_base_path, sprintf("pred_intensities_%d_500_%d_%s_across3_resample5_percentage_%d_seed_%d_self_quad_d_no_between_dummy.Rda",
                                                               cluster_num, hr, intensity_type, resample_percent, random_seed))

save(pred_cell_pattern_df, file = filename_pred_cell)
save(pred_intensities, file = filename_pred_intensity)

# Calculate AUC metrics
true_labels <- cell_pattern_df$marks
preds <- pred_intensities

# Macro-averaged AUC
total_auc <- sapply(levels(true_labels), function(level) {
  auc(roc(response = as.integer(true_labels == level), predictor = preds[, level]))
})
macro_avg_auc <- mean(total_auc)

# Weighted Macro-Averaged AUC
original_freq <- prop.table(table(true_labels))
weighted_macro_avg_auc <- sum(sapply(levels(true_labels), function(level) {
  auc(roc(response = as.integer(true_labels == level), predictor = preds[, level])) * original_freq[level]
}))

# Micro-averaged AUC
all_preds <- unlist(preds)
all_true <- unlist(lapply(levels(true_labels), function(level) as.integer(true_labels == level)))
micro_avg_auc <- auc(roc(response = all_true, predictor = all_preds))[1]

# Print and save AUC results
print(sprintf("Macro-Averaged AUC: %f", macro_avg_auc))
print(sprintf("Weighted Macro-Averaged AUC: %f", weighted_macro_avg_auc))
print(sprintf("Micro-Averaged AUC: %f", micro_avg_auc))

aucroc <- list(macro_avg_auc, weighted_macro_avg_auc, micro_avg_auc, total_auc)

filename_auc <- file.path(output_base_path, sprintf("AUCROC_%d_500_%d_%s_across3_resample5_percentage_%d_seed_%d_self_quad_d_no_between_dummy.Rda",
                                                    cluster_num, hr, intensity_type, resample_percent, random_seed))
save(aucroc, file = filename_auc)

print(aucroc)
