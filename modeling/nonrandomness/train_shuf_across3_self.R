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
r = as.numeric(args[4])
hr = as.numeric(args[5])
shuf_num = as.numeric(args[6])

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
#if (shuf_num == 0){
#ori_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3', '.Rda', sep = ''))
#} else{
#}
#load(ori_file)
#Quad_all_all$moadf = data_shuffled 
#rm(data_shuffled)
#gc()
#if (T){
if (!file.exists(file.path(data_dir, TMC, tissue, 'random_shuf', paste('coef_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))){
  shuf_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = ''))
#}
#load(ori_file)
load(shuf_file)
  print('training...')
  Quad_all_all$moadf = Quad_all_all$moadf[Quad_all_all$moadf$pattern_ID == 1,]
  formula_interaction = colnames(Quad_all_all$moadf)[7:(ncol(Quad_all_all$moadf)-1)]
  print(formula_interaction)
  formula_interaction = formula_interaction[-which(formula_interaction == "pattern_ID")]
  formula = get_formula(formula_interaction)
  model_train = mppm.fit.ppp(Data=Quad_all_all, formula, interaction=MultiStrauss(radii=r_matrix))
  print(model_train$Fit$FIT$coefficients)
  fmla = model_train$Fit$fmla
  family = model_train$Fit$FIT$family
  coef <- model_train$Fit$FIT$coefficients
  rm(model_train)
  gc()
  save(fmla, file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('fmla_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
  save(family, file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('family_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
  save(coef, file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('coef_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
} else{
  load(file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('fmla_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
  load(file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('family_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
  load(file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('coef_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
  print(coef)
}

#if (!file.exists(file =file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))){
#deviance_shuf = get_devi(Quad_all_all$moadf, coef, fmla, family)
#save(deviance_shuf, file =file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
#print(deviance_shuf)          
#}

#if (!file.exists(file =file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_ori_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))){
if(T){
print('shuf on ori...')
ori_file = file.path(data_dir, TMC, tissue, paste('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = ''))
load(ori_file)
deviance_shuf_on_ori = get_devi(Quad_all_all$moadf, coef, fmla, family)
save(deviance_shuf_on_ori, file =file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_ori_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
print(deviance_shuf_on_ori)    
}
