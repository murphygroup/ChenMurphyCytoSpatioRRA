# Load required libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(dplyr)
library(permute)
library(data.table)

# Set script directory
script_dir = './modeling/ppp/'
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

# Read command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]
tissue <- args[2]
intensity_type <- args[3]
num <- as.numeric(args[4])
r <- as.numeric(args[5])
hr <- as.numeric(args[6])

# Define data directory

# Determine Tissue Mapping Center (TMC)
TMC <- if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Define quad file path
quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))

# Load quad data
load(quad_file)

# Define interaction matrices
r_matrix <- matrix(rep(500, num * num), nrow = num)
hr_matrix <- matrix(rep(hr, num * num), nrow = num)

# Train the model and save coefficients if needed
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
save(fmla, file = file.path(data_dir, TMC, tissue, paste0('fmla_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))
save(family, file = file.path(data_dir, TMC, tissue, paste0('family_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))
save(coef, file = file.path(data_dir, TMC, tissue, paste0('coef_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))

# Compute deviance and save results
print('Computing deviance...')
deviance_ori <- get_devi(Quad_all_all$moadf, coef, fmla, family)
save(deviance_ori, file = file.path(data_dir, TMC, tissue, paste0('devi_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))
print(deviance_ori)
