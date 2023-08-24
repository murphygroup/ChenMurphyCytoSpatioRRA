library(spatstat)
library(spatstat.utils)
library(spatstat.data)
# library(ggplot2)
library(dplyr)
library(permute)
library(data.table)
library(biglm)
library(speedglm)
# script_dir = '/home/hrchen/Documents/Research/hubmap/ppm/script'
script_dir = '/home/haoranch/projects/HuBMAP/ppm/script'
# script_dir = '/Users/hrchen/Documents/Research/ppm'

setwd(script_dir)
#source('mppm.fit.ppp.biglm.R')
#source('mppm.fit.ppp.speedglm.R')

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

get_maxlogpl = function(individual_data, coef, fmla, family){

  individual_data <<- individual_data
  .mpl.W = individual_data$.mpl.W
  caseweight = individual_data$caseweight
  gcontrol = list()

  ctrl <- do.call(glm.control, resolve.defaults(gcontrol,
                                                list(maxit = 50)))
  devi_list <- glm.prep(fmla, family = quasi(link = "log", variance = "mu"), weights = .mpl.W * caseweight, data = individual_data, start = coef, subset = (individual_data$.mpl.SUBSET == "TRUE"), control = ctrl)
  W <- with(individual_data, .mpl.W * caseweight)
  SUBSET <- individual_data$.mpl.SUBSET
  Z <- (individual_data$.mpl.Y != 0)
  # print(sum(devi_list)/length(devi_list))
  # print(length(devi_list))
  devi_list_all = devi_list[SUBSET]
  devi_list_real = devi_list[Z & SUBSET]
  devi_list_dummy = devi_list[!Z & SUBSET]

  devi_list_all = devi_list_all[!is.na(devi_list_all)]
  devi_list_real = devi_list_real[!is.na(devi_list_real)]
  devi_list_dummy = devi_list_dummy[!is.na(devi_list_dummy)]

  devi_list_cell_type_all_list = c()
  devi_list_cell_type_real_list = c()
  devi_list_cell_type_dummy_list = c()

  for (c in 0:4){
    devi_list_cell_type_all = devi_list_all[individual_data$marks == c]
    devi_list_cell_type_real = devi_list_all[individual_data$marks == c & Z]
    devi_list_cell_type_dummy = devi_list_all[individual_data$marks == c & !Z]

    devi_list_cell_type_all = devi_list_cell_type_all[!is.na(devi_list_cell_type_all)]
    devi_list_cell_type_real = devi_list_cell_type_real[!is.na(devi_list_cell_type_real)]
    devi_list_cell_type_dummy = devi_list_cell_type_dummy[!is.na(devi_list_cell_type_dummy)]

    devi_list_cell_type_all_list = c(devi_list_cell_type_all_list, sum(devi_list_cell_type_all) / length(devi_list_cell_type_all))
    devi_list_cell_type_real_list = c(devi_list_cell_type_real_list, sum(devi_list_cell_type_real) / length(devi_list_cell_type_real))
    devi_list_cell_type_dummy_list = c(devi_list_cell_type_dummy_list, sum(devi_list_cell_type_dummy) / length(devi_list_cell_type_dummy))
  }
  rm(devi_list)
  gc()

  return(list('devi_all_avg' = sum(devi_list_all) / length(devi_list_all),
              'devi_real_avg' = sum(devi_list_real) / length(devi_list_real),
              'devi_dummy_avg' = sum(devi_list_dummy) / length(devi_list_dummy),
              'devi_all_cell_type_avg' = devi_list_cell_type_all_list,
              'devi_real_cell_type_avg' = devi_list_cell_type_real_list,
              'devi_dummy_cell_type_avg' = devi_list_cell_type_dummy_list))
}


args = commandArgs(TRUE)
tissue = args[1]
intensity_type = args[2]
num = as.numeric(args[3])
#r = as.numeric(args[4])
r = args[4]
hr = as.numeric(args[5])

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
quad_file = file.path(data_dir, TMC, tissue, paste('quad_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = ''))

load(quad_file)
r_matrix = matrix(rep(500, num*num), nrow = num)
hr_matrix = matrix(rep(hr, num*num), nrow = num)

if (!file.exists(file.path(data_dir, TMC, tissue, paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))){
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
  save(fmla, file = file.path(data_dir, TMC, tissue, paste('fmla_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
  save(family, file = file.path(data_dir, TMC, tissue, paste('family_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
  save(coef, file = file.path(data_dir, TMC, tissue, paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
}
print('fitting...')
load(file = file.path(data_dir, TMC, tissue, paste('fmla_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, paste('family_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, paste('coef_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
deviance_ori = get_devi(Quad_all_all$moadf, coef, fmla, family)
save(deviance_ori, file =file.path(data_dir, TMC, tissue, paste('devi_', num, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
print(deviance_ori)          
