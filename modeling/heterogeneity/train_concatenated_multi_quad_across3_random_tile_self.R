library(spatstat)
library(spatstat.utils)
library(spatstat.data)
# library(ggplot2)
library(dplyr)
library(permute)
library(data.table)
# script_dir = '/home/hrchen/Documents/Research/hubmap/ppm/script'
script_dir = '/home/haoranch/projects/HuBMAP/ppm/script'
# script_dir = '/Users/hrchen/Documents/Research/ppm'

setwd(script_dir)
source('mppm.fit.ppp.R')
source('glm.prep.R')
source('glm.fit.prep.R')
source('evaluate.ppp.R')
source('evalPairPotential.ppp.R')
source('get_devi.R')
get_formula = function(interactions){
  formula_string = ".mpl.Y ~ marks"
  for (interaction in interactions) {
    formula_string = paste(formula_string, interaction, sep = " + ")
  }
  return(as.formula(formula_string))
}


args = commandArgs(TRUE)
tissue = args[1]
intensity_type = args[2]
num = as.numeric(args[3])
#r = as.numeric(args[4])
r = args[4]
hr = as.numeric(args[5])
img_idx = as.numeric(args[6])
tile_side = args[7]
tile_num = args[8]
# tissue = 'LI'
# intensity_type = 'mean'
# num = 5
# r = 150
# shuf_num = 1

# data_dir = '/data2/PPM/HUBMAP_DATA_new/Stanford'
data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'

if (tissue == 'LI' | tissue == 'SI'){
  TMC = 'Stanford'
} else{
  TMC = 'Florida'
}
# data_dir = '/Users/hrchen/Documents/Research/ppm'
quad_file = file.path(data_dir, TMC, tissue, 'random_tile', paste('quad_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = ''))

load(quad_file)

cell_pattern_file = file.path(data_dir, TMC, tissue, 'random_tile', paste('cell_list_', num, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '.Rda', sep = ''))

load(cell_pattern_file)

cell_pattern_window = cell_list[[1]]$window
rect = c(cell_pattern_window$xrange[1]+500, cell_pattern_window$yrange[1]+500, cell_pattern_window$xrange[2]-500, cell_pattern_window$yrange[2]-500)
training = Quad_all_all$moadf
training = training[which(training$x >= rect[1] & training$x <= rect[3] & training$y >= rect[2] & training$y <= rect[4]),]
training = training[complete.cases(training),]

Quad_all_all$moadf = training

#print(training[1:5,])
r_matrix = matrix(rep(500, num*num), nrow = num)
hr_matrix = matrix(rep(hr, num*num), nrow = num)

if (!file.exists(file.path(data_dir, TMC, tissue, 'random_tile', paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))){
#if (T){
  print('training...')
  formula_interaction = colnames(Quad_all_all$moadf)[7:(ncol(Quad_all_all$moadf)-1)]
  formula_interaction = formula_interaction[-which(formula_interaction == "pattern_ID")]
  formula = get_formula(formula_interaction)
  model_train = mppm.fit.ppp(Data=Quad_all_all, formula, interaction=MultiStraussHard(iradii=r_matrix, hradii=hr_matrix)) 
  print(model_train$Fit$FIT$coefficients)
  fmla = model_train$Fit$fmla
  family = model_train$Fit$FIT$family
  coef <- model_train$Fit$FIT$coefficients
  rm(model_train)
  gc()
  save(fmla, file = file.path(data_dir, TMC, tissue, 'random_tile', paste('fmla_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
  save(family, file = file.path(data_dir, TMC, tissue, 'random_tile', paste('family_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
  save(coef, file = file.path(data_dir, TMC, tissue, 'random_tile', paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
}
print('fitting...')
load(file = file.path(data_dir, TMC, tissue, 'random_tile', paste('fmla_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, 'random_tile', paste('family_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, 'random_tile', paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
deviance_ori = get_devi(training, coef, fmla, family)
save(deviance_ori, file =file.path(data_dir, TMC, tissue, 'random_tile', paste('devi_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
print(deviance_ori)          
