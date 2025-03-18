library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(dplyr)
library(permute)
library(data.table)
library(gtools)

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


# Determine Tissue Mapping Center (TMC)
TMC = if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Load model parameters
load(file.path(data_dir, TMC, tissue, 'random_shuf', paste0('fmla_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda')))
load(file.path(data_dir, TMC, tissue, 'random_shuf', paste0('family_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda')))
load(file.path(data_dir, TMC, tissue, 'random_shuf', paste0('coef_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda')))

print('fitting..')

# Set seed for reproducibility
set.seed(shuf_num)

# Get list of shuffled data files
ori_shuf_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda'))
shuf_file_list = Sys.glob(file.path(data_dir, TMC, tissue, 'random_shuf', paste0('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_*_self_quad_d_no_between_dummy.Rda')))
shuf_file_list = mixedsort(shuf_file_list)[1:100]  # Sort and take the first 100 shuffled files
print(shuf_file_list)

# Randomly select a different shuffle file
shuf_file = sample(shuf_file_list, 1)

# Extract shuffle number from filename
path_split = strsplit(shuf_file, split = '_')[[1]]
other_shuf_num = path_split[length(path_split) - 5]

# Ensure `shuf_num` and `other_shuf_num` are not the same
if (shuf_num == other_shuf_num) {
  file.remove(file.path(data_dir, TMC, tissue, 'random_shuf', paste0('devi_shuf_on_another_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', shuf_num, '_on_', other_shuf_num, '_self_quad_d_no_between_dummy.Rda')))
  shuf_file = sample(shuf_file_list, 1)  # Resample another shuffle file
  path_split = strsplit(shuf_file, split = '_')[[1]]
  other_shuf_num = path_split[length(path_split) - 5]
}

# Compute deviance if not already computed
deviance_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('devi_shuf_on_another_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', shuf_num, '_on_', other_shuf_num, '_self_quad_d_no_between_dummy.Rda'))

if (!file.exists(deviance_file)) {
  print(shuf_file)
  load(shuf_file)  # Load the shuffled dataset
  training = Quad_all_all$moadf
  deviance_shuf = get_devi(training, coef, fmla, family)  # Compute deviance
  
  # Save deviance results
  save(deviance_shuf, file = deviance_file)
  print(deviance_shuf)
} else {
  # Load existing deviance file and save it under a new name
  load(file.path(data_dir, TMC, tissue, 'random_shuf', paste0('devi_shuf_on_other_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', shuf_num, '_on_', other_shuf_num, '_self_quad_d_no_between_dummy.Rda')))
  save(deviance_shuf, file = deviance_file)
  print('done')
}
