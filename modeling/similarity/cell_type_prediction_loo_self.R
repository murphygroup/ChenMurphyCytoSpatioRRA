# Load required libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(spatstat.geom)
library(dplyr)
library(proxy)
library(fields)
library(pROC)
library(data.table)

# Set random seed for reproducibility
set.seed(1)

# Get command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments
train_tissue <- args[2]
test_tissue <- args[3]
intensity_type <- args[4]
cluster_num <- as.numeric(args[5])
radius <- as.numeric(args[6])
hr <- as.numeric(args[7])
img_idx1 <- as.numeric(args[8])

# Determine Tissue Mapping Center (TMC) for train and test tissues
TMC1 <- if (train_tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'
TMC2 <- if (test_tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Load model parameters (formula, family, coefficients) for the trained model
load(file.path(data_dir, TMC1, train_tissue, paste0('fmla_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_looidx_', img_idx1, '_across3_self_quad_d_no_between_dummy.Rda')))
load(file.path(data_dir, TMC1, train_tissue, paste0('family_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_looidx_', img_idx1, '_across3_self_quad_d_no_between_dummy.Rda')))
load(file.path(data_dir, TMC1, train_tissue, paste0('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_looidx_', img_idx1, '_across3_self_quad_d_no_between_dummy.Rda')))

# Extract model coefficients and interactions
G <- list()
coef_names <- names(coef)

for (i in 0:(cluster_num - 1)) {
  coef_current <- c()
  for (r in seq(100, radius, 100)) {
    for (j in 0:(cluster_num - 1)) {
      mark1 <- as.character(i)
      mark2 <- as.character(j)
      interaction_name <- if (i <= j) {
        paste0('X', mark1, 'x', 'X', mark2, 'x', r)
      } else {
        paste0('X', mark2, 'x', 'X', mark1, 'x', r)
      }
      coef_current <- c(coef_current, coef[grep(interaction_name, coef_names, value = FALSE)])
    }
  }
  G[[i + 1]] <- coef_current
}

# Compute intensity terms
B <- c()
for (i in 1:cluster_num) {
  B_current <- if (i == 1) coef[i] else coef[i] + coef[1]
  B <- c(B, B_current)
}

# Define function to calculate cell type probabilities
calc_mark_prob <- function(counts) {
  m_prob <- sapply(seq_along(types), function(m) {
    m_int <- as.integer(m)
    interact_term <- sum(G[[m_int + 1]] * counts)
    intensity_term <- B[m_int + 1]
    exp(intensity_term + interact_term)
  })
  
  m_prob <- m_prob / sum(m_prob)
  final_type <- sample(types, 1, prob = m_prob)
  return(list(final_type, m_prob))
}

# Iterate through test images
for (img_idx2 in 1:30) {
  
  quad_file <- file.path(data_dir, TMC2, test_tissue, paste0('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx2, '_self_quad_d_no_between_dummy.Rda'))
  
  if (file.exists(quad_file)) {
    if (train_tissue != test_tissue || img_idx1 == img_idx2) {
      
      aucroc_file <- file.path(data_dir, TMC1, train_tissue, paste0('AUCROC_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_quad_d_no_between_dummy.Rda'))
      pred_file <- file.path(data_dir, TMC1, train_tissue, paste0('pred_intensities_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_quad_d_no_between_dummy.Rda'))
      
      if (!file.exists(aucroc_file) && !file.exists(pred_file)) {
        print(img_idx2)
        
        # Load quadrature scheme and cell list
        load(quad_file)
        load(file.path(data_dir, TMC2, test_tissue, paste0('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx2, '.Rda')))
        
        Y <- Quad_all_all$moadf$.mpl.Y[which(Quad_all_all$moadf$.mpl.Y != 0)]
        
        cell_pattern <- cell_list[[1]]
        cell_pattern_df <- data.frame(cell_pattern)
        cell_pattern_window <<- cell_pattern$window
        
        # Predict cell types and intensities
        marks <- c()
        intensities <- c()
        for (point_current_idx in 1:nrow(cell_pattern_df)) {
          point_current <- cell_pattern_df[point_current_idx, ]
          counts_current <- t(apply(point_current[, 1:2], 1, calc_dist_optimized, cell_pattern_df))
          mark_current <- apply(counts_current, 1, calc_mark_prob)
          marks <- c(marks, mark_current[[1]][[1]])
          intensities <- rbind(intensities, mark_current[[1]][[2]])
        }
        
        # Save predicted cell types
        pred_cell_type <- data.frame(x = cell_pattern_df$x, y = cell_pattern_df$y, marks = marks)
        save(pred_cell_type, file = file.path(data_dir, TMC1, train_tissue, paste0('pred_cell_type_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_quad_d_no_between_dummy.Rda')))
        
        # Save predicted intensities
        pred_intensities <- data.frame(intensities)
        save(pred_intensities, file = pred_file)
        
        # Calculate AUC-ROC metrics
        true <- cell_pattern_df$marks
        preds <- pred_intensities
        
        num_classes <- length(levels(true))
        total_auc <- sapply(seq_len(num_classes), function(i) {
          true_binary <- as.numeric(true == levels(true)[i])
          auc(roc(response = true_binary, predictor = preds[, i]))
        })
        
        macro_avg_auc <- mean(total_auc)
        print(paste("Macro-Averaged AUC:", macro_avg_auc))
        
        original_freq <- table(cell_pattern_df$marks) / length(cell_pattern_df$marks)
        weighted_macro_avg_auc <- sum(sapply(seq_len(num_classes), function(i) {
          true_binary <- as.numeric(true == levels(true)[i])
          auc(roc(response = true_binary, predictor = preds[, i])) * original_freq[i]
        }))
        print(paste("Weighted Macro-Averaged AUC:", weighted_macro_avg_auc))
        
        all_preds <- unlist(preds)
        all_true <- unlist(lapply(seq_len(num_classes), function(i) as.numeric(true == levels(true)[i])))
        micro_avg_auc <- auc(roc(response = all_true, predictor = all_preds))[1]
        print(paste("Micro-Averaged AUC:", micro_avg_auc))
        
        aucroc <- list(macro_avg_auc, weighted_macro_avg_auc, micro_avg_auc, total_auc)
        save(aucroc, file = aucroc_file)
      }
    }
  }
}
