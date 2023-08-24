library(spatstat)
library(spatstat.utils)
library(spatstat.data)
# library(ggplot2)
library(dplyr)
library(permute)
library(data.table)
library(gtools)
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

# data_dir = '/data/PPM/HUBMAP_DATA_new/Stanford'
data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'

if (tissue == 'LI' | tissue == 'SI'){
  TMC = 'Stanford'
} else{
  TMC = 'Florida'
}



load(file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('fmla_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('family_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('coef_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))

print('fitting..')
set.seed(shuf_num)
#for (other_shuf_num in c(1:10)){
ori_shuf_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = ''))
shuf_file_list = Sys.glob(file.path(data_dir, TMC, tissue, 'random_shuf', paste('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_*_self_dummy_grid_eps_20.Rda', sep = '')))
shuf_file_list = mixedsort(shuf_file_list)[1:100]
print(shuf_file_list)
shuf_file = sample(shuf_file_list, 1) 

path_split = strsplit(shuf_file, split = '_')[[1]]
other_shuf_num = path_split[length(path_split)-5]

if (shuf_num == other_shuf_num){
file.remove(file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_another_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', shuf_num, '_on_', other_shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
shuf_file = sample(shuf_file_list, 1)
path_split = strsplit(shuf_file, split = '_')[[1]]
other_shuf_num = path_split[length(path_split)-5]
}

#  if (ori_shuf_file != shuf_file){
if(T){
    if (!file.exists(file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_another_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', shuf_num, '_on_', other_shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))){
      print(shuf_file)
      load(shuf_file)
      training = Quad_all_all$moadf
      deviance_shuf = get_devi(training, coef, fmla, family)
      save(deviance_shuf, file =file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_another_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', shuf_num, '_on_', other_shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
      print(deviance_shuf)
    } else{
load(file =file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_other_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', shuf_num, '_on_', other_shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))

save(deviance_shuf, file =file.path(data_dir, TMC, tissue, 'random_shuf', paste('devi_shuf_on_another_shuf_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', shuf_num, '_on_', other_shuf_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
print('done')
}
  }

