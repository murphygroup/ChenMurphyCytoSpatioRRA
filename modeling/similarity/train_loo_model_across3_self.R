# Load required libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(dplyr)
library(permute)
library(data.table)

# Get command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments
train_tissue <- args[2]
intensity_type <- args[3]
num <- as.numeric(args[4])
r <- as.numeric(args[5])
hr <- as.numeric(args[6])

# Determine Tissue Mapping Center (TMC) based on tissue type
TMC <- if (train_tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Function to get the number of cells in a dataset
get_cell_num <- function(individual_data) {
  SUBSET <- individual_data$.mpl.SUBSET
  Z <- (individual_data$.mpl.Y != 0)
  sum(Z & SUBSET)
}

# Function to ensure all cell types are included in the dataset
fill_missing_marks <- function(pattern, marks, example) {
  pattern <- rbindlist(list(pattern, example), fill = TRUE)[-nrow(pattern), ]
  for (mark in marks) {
    if (!(mark %in% pattern$marks)) {
      pattern <- rbindlist(list(pattern, example), fill = TRUE)
      pattern[nrow(pattern), 1:4] <- 0
      pattern[nrow(pattern), 7:ncol(pattern)] <- 0
      pattern[nrow(pattern), 5] <- mark
    }
  }
  return(pattern)
}

# Function to remove unused cell types
remove_redundant_marks <- function(pattern, training_marks) {
  for (mark in unique(pattern$marks)) {
    if (!(mark %in% training_marks)) {
      pattern <- pattern[-which(pattern$marks == mark), ]
      for (j in ncol(pattern):7) {
        if (grepl(paste0('X', mark), colnames(pattern)[j], fixed = TRUE)) {
          pattern <- pattern[, -j, with = FALSE]
        }
      }
    }
  }
  return(pattern)
}

# Function to generate a formula string for modeling
get_formula <- function(interactions) {
  formula_string <- ".mpl.Y ~ marks"
  for (interaction in interactions) {
    formula_string <- paste0(formula_string, " + ", interaction)
  }
  return(as.formula(formula_string))
}

# Function to split IDs into folds for cross-validation
get_split_ID <- function(all_ID, fold) {
  region_length <- round(length(all_ID) / fold)
  region_seq <- seq(1, length(all_ID), region_length)
  region_seq <- c(region_seq, tail(all_ID, n = 1) + 1)
  return(region_seq)
}

# Define matrices for interactions
r_matrix <- matrix(rep(r, num * num), nrow = num)
hr_matrix <- matrix(rep(hr, num * num), nrow = num)

# Initialize lists for deviance tracking
deviance_train_list <- c()
deviance_val_list <- c()
deviance_test_list <- c()

# Load the dataset
quad_file <- file.path(data_dir, TMC, train_tissue, paste0('quad_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda'))
load(quad_file)

# Get unique image IDs and create train-test splits
all_img_ID <- unique(Quad_all_all$moadf$pattern_ID)
all_split <- get_split_ID(all_img_ID, length(all_img_ID))

# Loop through outer folds
for (outer_idx in 1:(length(all_split) - 1)) {
  print(outer_idx)
  
  # Define train-test split
  pattern_test <- Quad_all_all$moadf[Quad_all_all$moadf$pattern_ID >= all_split[outer_idx] & Quad_all_all$moadf$pattern_ID < all_split[outer_idx + 1], ]
  pattern_train <- Quad_all_all$moadf[Quad_all_all$moadf$pattern_ID < all_split[outer_idx] | Quad_all_all$moadf$pattern_ID >= all_split[outer_idx + 1], ]
  
  coef_file <- file.path(data_dir, TMC, train_tissue, paste0('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda'))
  
  # Train model if coefficients do not already exist
  if (!file.exists(coef_file)) {
    print('Training...')
    
    # Prepare training dataset
    Quad_training <- Quad_all_all
    rm(Quad_all_all)
    gc()
    
    Quad_training$moadf <- pattern_train
    formula_interaction <- colnames(pattern_train)[7:(ncol(pattern_train) - 1)]
    formula_interaction <- formula_interaction[-which(formula_interaction == "pattern_ID")]
    formula <- get_formula(formula_interaction)
    
    # Train model
    model_train <- mppm.fit.ppp(Data = Quad_training, formula, interaction = MultiStraussHard(iradii = r_matrix, hradii = hr_matrix))
    print(model_train$Fit$FIT$coefficients)
    
    # Save model parameters
    save(model_train$Fit$fmla, file = file.path(data_dir, TMC, train_tissue, paste0('fmla_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda')))
    save(model_train$Fit$FIT$family, file = file.path(data_dir, TMC, train_tissue, paste0('family_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda')))
    save(model_train$Fit$FIT$coefficients, file = coef_file)
    
    # Cleanup
    rm(model_train)
    gc()
  }
  
  # Compute number of cells in train and test datasets
  cell_num_file <- file.path(data_dir, TMC, train_tissue, paste0('cell_num_loo_', num, '_', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda'))
  
  if (!file.exists(cell_num_file)) {
    cell_num_train <- get_cell_num(pattern_train)
    cell_num_test <- get_cell_num(pattern_test)
    cell_num <- list(cell_num_train, cell_num_test)
    save(cell_num, file = cell_num_file)
  }
}
