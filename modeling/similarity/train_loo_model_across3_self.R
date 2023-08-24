library(spatstat)
library(spatstat.utils)
library(spatstat.data)
# library(ggplot2)
library(dplyr)
library(permute)
library(data.table)
# script_dir = '/home/hrchen/Documents/Research/hubmap/ppm/script'
script_dir = '/home/haoranch/projects/HuBMAP/ppm/script'
setwd(script_dir)
source('mppm.fit.ppp.R')
source('glm.prep.R')
source('glm.fit.prep.R')
source('get_devi.R')

get_cell_num = function(individual_data){
  SUBSET <- individual_data$.mpl.SUBSET
  Z <- (individual_data$.mpl.Y != 0)
  cell_num_data = dim(individual_data[Z & SUBSET,])[1]
  return(cell_num_data)
}



fill_missing_marks = function(pattern, marks, example){
  pattern = rbindlist(list(pattern, example), fill = T)
  pattern = pattern[-nrow(pattern), ]
  for (mark in marks) {
    if (length(which(pattern$marks == mark)) == 0){
      # pattern[nrow(pattern) + 1,] = pattern[nrow(pattern),]
      pattern = rbindlist(list(pattern, example), fill = T)
      pattern[nrow(pattern), 1:4] = 0
      pattern[nrow(pattern), 7:ncol(pattern)] = 0
      pattern[nrow(pattern), 5] = mark
    }
  }
  return(pattern)
}

remove_redundant_marks = function(pattern, training_marks){
  for (mark in unique(pattern$marks)) {
    if (!(mark %in% training_marks)) {
      pattern = pattern[-which(pattern$marks == mark),]
      for (j in ncol(pattern):7){
        if (grepl(paste('X', mark, sep = ''), colnames(pattern)[j], fixed = TRUE)){
          pattern = pattern[, -j]
        }
      }
    }
  }
  return(pattern)
}

get_formula = function(interactions){
  formula_string = ".mpl.Y ~ marks"
  for (interaction in interactions) {
    formula_string = paste(formula_string, interaction, sep = " + ")
  }
  return(as.formula(formula_string))
}

get_split_ID = function(all_ID, fold){
  region_length = round(length(all_ID) / fold)
  region_seq = c()
  for (i in seq(1, length(all_ID), region_length)) {
    region_seq = c(region_seq, all_ID[i])
  }
  region_seq = c(region_seq, tail(all_ID, n=1)+1)
  return(region_seq)
}

args = commandArgs(TRUE)
train_tissue = args[1]
intensity_type = args[2]
num = as.numeric(args[3])
r = as.numeric(args[4])
hr = as.numeric(args[5])
# train_tissue = 'LI'
# intensity_type = 'mean'
# num = 5
# r = 50

if (train_tissue == 'SI' | train_tissue == 'LI') {
  TMC = 'Stanford'
} else{
  TMC = 'Florida'
}

data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'




r_matrix = matrix(rep(r, num*num), nrow = num)
hr_matrix = matrix(rep(hr, num*num), nrow = num)

deviance_train_list = c()
deviance_val_list = c()
deviance_test_list = c()

load(file.path(data_dir, TMC, train_tissue, paste('quad_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
all_img_ID = unique(Quad_all_all$moadf$pattern_ID)
all_split = get_split_ID(all_img_ID, length(all_img_ID))

for (outer_idx in 1:(length(all_split)-1)) {
  print(outer_idx)
#for (outer_idx in 2:2){
  pattern_test = Quad_all_all$moadf[which(Quad_all_all$moadf$pattern_ID >= all_split[outer_idx] & Quad_all_all$moadf$pattern_ID < all_split[outer_idx+1]), ]
  pattern_train = Quad_all_all$moadf[which(Quad_all_all$moadf$pattern_ID < all_split[outer_idx] | Quad_all_all$moadf$pattern_ID >= all_split[outer_idx+1]), ]
  training_ID = unique(pattern_train$pattern_ID)
  training_split = get_split_ID(training_ID, length(training_ID))
print(file.path(data_dir, TMC, train_tissue, paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3.Rda', sep = '')))  
if (!file.exists(file.path(data_dir, TMC, train_tissue, paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))){
    #if (T){
    print('training...')
    Quad_training = Quad_all_all
    rm(Quad_all_all)
    gc()  
    Quad_training$moadf = pattern_train
    formula_interaction = colnames(pattern_train)[7:(ncol(pattern_train)-1)]
    formula_interaction = formula_interaction[-which(formula_interaction == "pattern_ID")]
    formula = get_formula(formula_interaction)
    
    model_train = mppm.fit.ppp(Data=Quad_training, formula, interaction=MultiStraussHard(iradii=r_matrix, hradii=hr_matrix))
    print(model_train$Fit$FIT$coefficients)
    rm(Quad_training)
    gc()  
    fmla = model_train$Fit$fmla
    family = model_train$Fit$FIT$family
    coef <- model_train$Fit$FIT$coefficients
    rm(model_train)
    gc()
    # save(model_train, file = file.path(data_dir, TMC, train_tissue, paste('model_', num, '_', r, '_', hr, '_', intensity_type, '_', outer_idx, '.Rda', sep = '')))
    save(fmla, file = file.path(data_dir, TMC, train_tissue, paste('fmla_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
    save(family, file = file.path(data_dir, TMC, train_tissue, paste('family_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
    save(coef, file = file.path(data_dir, TMC, train_tissue, paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
  }
  #if (file.exists(file.path(data_dir, TMC, train_tissue, paste('devi_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = ''))) == F){
#    if (T){
#    print('fitting...')
#    load(file = file.path(data_dir, TMC, train_tissue, paste('fmla_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
#    load(file = file.path(data_dir, TMC, train_tissue, paste('family_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
#    load(file = file.path(data_dir, TMC, train_tissue, paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
    
    #deviance_train = get_devi(pattern_train, coef, fmla, family)
    
    #deviance_test = get_devi(pattern_test, coef, fmla, family)      
    # pattern_test2[is.na(pattern_test2)] = 0
    # pattern_test2 = remove_redundant_marks(pattern_test2, unique(pattern_train$marks))
    # pattern_test2 = fill_missing_marks(pattern_test2, unique(pattern_train$marks), pattern_train_example)
    # pattern_test2[is.na(pattern_test2)] = 0
    # deviance_test2 = get_devi(pattern_test2)
    #results = c(deviance_train, deviance_test)
    #print(results)
    
    #save(results, file = file.path(data_dir, TMC, train_tissue, paste('devi_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
    load(file.path(data_dir, TMC, train_tissue, paste('quad_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
  }
  if (file.exists(file.path(data_dir, TMC, train_tissue, paste('cell_num_loo_', num, '_', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = ''))) == F){
    #if (T){
    cell_num_train = get_cell_num(pattern_train)
    cell_num_test = get_cell_num(pattern_test)
    cell_num = list(cell_num_train, cell_num_test)
    save(cell_num, file = file.path(data_dir, TMC, train_tissue, paste('cell_num_loo_', num, '_100-', r, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
  }
}




