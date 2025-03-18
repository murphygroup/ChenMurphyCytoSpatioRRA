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

# Function to construct formula for model fitting
get_formula = function(interactions) {
  formula_string = ".mpl.Y ~ marks"
  for (interaction in interactions) {
    formula_string = paste(formula_string, interaction, sep = " + ")
  }
  return(as.formula(formula_string))
}

# Read command-line arguments
args = commandArgs(TRUE)
data_dir = args[1]
tissue = args[2]
intensity_type = args[3]
num = as.numeric(args[4])
r = as.numeric(args[5])
hr = as.numeric(args[6])
shuf_num = as.numeric(args[7])

# Determine Tissue Mapping Center (TMC) based on tissue type
TMC = if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Check if model coefficients file exists, otherwise train a new model
coef_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('coef_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda'))

if (!file.exists(coef_file)) {
  
  # Load shuffled data file for model training
  shuf_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda'))
  load(shuf_file)
  
  print('training...')
  
  # Extract relevant data for model training
  Quad_all_all$moadf = Quad_all_all$moadf[Quad_all_all$moadf$pattern_ID == 1,]
  formula_interaction = colnames(Quad_all_all$moadf)[7:(ncol(Quad_all_all$moadf) - 1)]
  
  print(formula_interaction)
  
  # Remove "pattern_ID" from formula interactions
  formula_interaction = formula_interaction[formula_interaction != "pattern_ID"]
  formula = get_formula(formula_interaction)
  
  # Train the model
  model_train = mppm.fit.ppp(Data = Quad_all_all, formula, interaction = MultiStrauss(radii = r_matrix))
  
  # Extract and save model parameters
  print(model_train$Fit$FIT$coefficients)
  fmla = model_train$Fit$fmla
  family = model_train$Fit$FIT$family
  coef = model_train$Fit$FIT$coefficients
  
  # Clean up memory
  rm(model_train)
  gc()
  
  # Save model components
  save(fmla, file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('fmla_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda')))
  save(family, file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('family_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda')))
  save(coef, file = coef_file)
  
} else {
  # Load existing model parameters
  load(file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('fmla_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda')))
  load(file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('family_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda')))
  load(coef_file)
  print(coef)
}

# Compute deviance for shuffled model against original data
print('shuf on ori...')

# Load original data file
ori_file = file.path(data_dir, TMC, tissue, paste0('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))
load(ori_file)

# Calculate and save deviance
deviance_shuf_on_ori = get_devi(Quad_all_all$moadf, coef, fmla, family)
save(deviance_shuf_on_ori, file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('devi_shuf_on_ori_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda')))
print(deviance_shuf_on_ori)
