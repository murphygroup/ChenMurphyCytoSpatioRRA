# Load required libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(dplyr)
library(permute)
library(data.table)


# Set script directory
script_dir <- './modeling/ppp/'
setwd(script_dir)

# Load required functions
source('mppm.fit.ppp.R')
source('glm.prep.R')
source('glm.fit.prep.R')
source('evaluate.ppp.R')
source('evalPairPotential.ppp.R')
source('get_devi.R')

# Function to construct the formula for model fitting
get_formula <- function(interactions) {
  formula_string <- ".mpl.Y ~ marks"
  for (interaction in interactions) {
    formula_string <- paste(formula_string, interaction, sep = " + ")
  }
  return(as.formula(formula_string))
}

# Function to compute max log pseudo-likelihood (PL)
get_maxlogpl <- function(individual_data, coef, fmla, family) {
  individual_data <<- individual_data
  .mpl.W <- individual_data$.mpl.W
  caseweight <- individual_data$caseweight
  gcontrol <- list()
  
  ctrl <- do.call(glm.control, resolve.defaults(gcontrol, list(maxit = 50)))
  devi_list <- glm.prep(fmla, family = quasi(link = "log", variance = "mu"), 
                        weights = .mpl.W * caseweight, data = individual_data, 
                        start = coef, subset = (individual_data$.mpl.SUBSET == "TRUE"), 
                        control = ctrl)
  
  W <- with(individual_data, .mpl.W * caseweight)
  SUBSET <- individual_data$.mpl.SUBSET
  Z <- (individual_data$.mpl.Y != 0)
  
  # Filter non-NA deviance values
  devi_list_all <- devi_list[SUBSET][!is.na(devi_list[SUBSET])]
  devi_list_real <- devi_list[Z & SUBSET][!is.na(devi_list[Z & SUBSET])]
  devi_list_dummy <- devi_list[!Z & SUBSET][!is.na(devi_list[!Z & SUBSET])]
  
  # Compute cell-type-specific deviances
  devi_list_cell_type_all_list <- c()
  devi_list_cell_type_real_list <- c()
  devi_list_cell_type_dummy_list <- c()
  
  for (c in 0:4) {
    devi_list_cell_type_all <- devi_list_all[individual_data$marks == c]
    devi_list_cell_type_real <- devi_list_real[individual_data$marks == c & Z]
    devi_list_cell_type_dummy <- devi_list_dummy[individual_data$marks == c & !Z]
    
    devi_list_cell_type_all <- devi_list_cell_type_all[!is.na(devi_list_cell_type_all)]
    devi_list_cell_type_real <- devi_list_cell_type_real[!is.na(devi_list_cell_type_real)]
    devi_list_cell_type_dummy <- devi_list_cell_type_dummy[!is.na(devi_list_cell_type_dummy)]
    
    devi_list_cell_type_all_list <- c(devi_list_cell_type_all_list, sum(devi_list_cell_type_all) / length(devi_list_cell_type_all))
    devi_list_cell_type_real_list <- c(devi_list_cell_type_real_list, sum(devi_list_cell_type_real) / length(devi_list_cell_type_real))
    devi_list_cell_type_dummy_list <- c(devi_list_cell_type_dummy_list, sum(devi_list_cell_type_dummy) / length(devi_list_cell_type_dummy))
  }
  
  rm(devi_list)
  gc()
  
  return(list(
    'devi_all_avg' = sum(devi_list_all) / length(devi_list_all),
    'devi_real_avg' = sum(devi_list_real) / length(devi_list_real),
    'devi_dummy_avg' = sum(devi_list_dummy) / length(devi_list_dummy),
    'devi_all_cell_type_avg' = devi_list_cell_type_all_list,
    'devi_real_cell_type_avg' = devi_list_cell_type_real_list,
    'devi_dummy_cell_type_avg' = devi_list_cell_type_dummy_list
  ))
}

# Read command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]  # Updated to take data directory as an argument
tissue <- args[2]
intensity_type <- args[3]
num <- as.numeric(args[4])
r <- args[5]
hr <- as.numeric(args[6])

# Determine Tissue Mapping Center (TMC)
TMC <- if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Define quad file path
quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))

# Load quad data
load(quad_file)

# Define interaction matrices
r_matrix <- matrix(rep(500, num * num), nrow = num)
hr_matrix <- matrix(rep(hr, num * num), nrow = num)

# Train the model and save coefficients
print('Training model...')
formula_interaction <- colnames(Quad_all_all$moadf)[7:(ncol(Quad_all_all$moadf) - 1)]
formula_interaction <- formula_interaction[formula_interaction != "pattern_ID"]
formula <- get_formula(formula_interaction)

# Fit the model using MultiStraussHard interaction
model_train <- mppm.fit.ppp(Data = Quad_all_all, formula, interaction = MultiStraussHard(iradii = r_matrix, hradii = hr_matrix))

# Extract model parameters
print(model_train$Fit$FIT$coefficients)
fmla <- model_train$Fit$fmla
family <- model_train$Fit$FIT$family
coef <- model_train$Fit$FIT$coefficients

# Clean up memory
rm(model_train)
gc()

# Save model parameters
save(fmla, file = file.path(data_dir, TMC, tissue, paste0('fmla_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))
save(family, file = file.path(data_dir, TMC, tissue, paste0('family_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))
save(coef, file = file.path(data_dir, TMC, tissue, paste0('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))

# Compute deviance and save results
print('Computing deviance...')
deviance_ori <- get_devi(Quad_all_all$moadf, coef, fmla, family)
save(deviance_ori, file = file.path(data_dir, TMC, tissue, paste0('devi_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))
print(deviance_ori)
