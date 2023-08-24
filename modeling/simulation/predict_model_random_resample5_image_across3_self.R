library(spatstat)
library(spatstat.utils)
library(spatstat.data)
# library(ggplot2)
library(dplyr)
library(permute)
library(data.table)
library(fields)
library(pROC)

# script_dir = '/home/hrchen/Documents/Research/hubmap/ppm/script'
script_dir = '/home/haoranch/projects/HuBMAP/ppm/script'
setwd(script_dir)

source('calc_mark_prob.R')
source('calc_dist_optimized.R')



args = commandArgs(TRUE)
tissue = args[1]
intensity_type = args[2]
cluster_num = as.numeric(args[3])
radius = as.numeric(args[4])
hr = as.numeric(args[5])
img_idx = as.numeric(args[6])
# tissue = 'LN'
# intensity_type = 'total'
# cluster_num = 5
# radius = 500
# hr = 1
# img_idx = 15


if (tissue == 'SI' | tissue == 'LI') {
  TMC = 'Stanford'
} else{
  TMC = 'Florida'
}


data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
# data_dir = '/data/PPM/HUBMAP_DATA_new'




r_matrix = matrix(rep(radius, cluster_num*cluster_num), nrow = cluster_num)
hr_matrix = matrix(rep(hr, cluster_num*cluster_num), nrow = cluster_num)

deviance_train_list = c()
deviance_val_list = c()
deviance_test_list = c()

load(file = file.path(data_dir, TMC, tissue, paste('fmla_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, paste('family_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, paste('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_dummy_grid_eps_20.Rda', sep = '')))



# cell_pattern = cell_list[[1]]
# cell_pattern_df = data.frame(cell_pattern)
# cell_pattern_window <<- cell_pattern$window


G = list()
coef_names = names(coef)
for (i in 0:(cluster_num-1)){
  coef_current = c()
  for (r in seq(100, radius, 100)){
    for (j in 0:(cluster_num-1)){
      mark1 = as.character(i)
      mark2 = as.character(j)
      if (i <= j){
        interaction_name = paste('X', mark1, 'x', 'X', mark2, 'x', r, sep = '')
      } else{
        interaction_name = paste('X', mark2, 'x', 'X', mark1, 'x', r, sep = '')
      }
      coef_current = c(coef_current, coef[grep(interaction_name, coef_names, value=F)])
      # print(coef_current)
      
    }
  }
  G[[i+1]] = coef_current
}

B = c()
for (i in 1:cluster_num){
  if (i == 1){
    B_current = coef[i]
  } else{
    B_current = coef[i] + coef[1]
  }
  B = c(B, B_current)
}


cell_types = c(0, 1, 2, 3, 4)

for (resample_percent in c(0, 50, 100, 150, 200, 250, 300, 350, 400)){
  if (file.exists(file.path(data_dir, TMC, tissue, 'random_resample', paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))){
    # if(!file.exists(file.path(data_dir, TMC, tissue, 'random_resample', paste('AUCROC_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))){
    if (T){
      if(!file.exists(file.path(data_dir, TMC, tissue, 'random_resample', paste('pred_intensities_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))){
        load(file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))
        cell_pattern_df = data.frame(marked_simulated_pattern[[1]])
        cell_pattern_window <<- marked_simulated_pattern[[1]]$window
        
        marks = c()
        intensities = c()
        for (point_current_idx in 1:nrow(cell_pattern_df)){
          # print(point_current_idx)
          point_current = cell_pattern_df[point_current_idx, ]
          counts_current = t(apply(point_current[,1:2], 1, calc_dist_optimized, cell_pattern_df, cell_pattern_window))
          mark_current = apply(counts_current, 1, calc_mark_prob, B, G, cell_types)
          marks = c(marks, mark_current[[1]][[1]])
          intensities = rbind(intensities, mark_current[[1]][[2]])
        }
        
        pred_cell_pattern_df = data.frame(x = cell_pattern_df$x, y = cell_pattern_df$y, marks = marks)
        pred_intensities = data.frame(intensities)
        
        filename1 = file.path(data_dir, TMC, tissue, 'random_resample', paste('pred_cell_type_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))
        save(pred_cell_pattern_df, file = filename1)
        
        filename2 = file.path(data_dir, TMC, tissue, 'random_resample', paste('pred_intensities_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))
        save(pred_intensities, file = filename2)
      } 
      else{
        filename2 = file.path(data_dir, TMC, tissue, 'random_resample', paste('pred_intensities_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))
        load(filename2)
      } 

      load(file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))
      cell_pattern_df = data.frame(marked_simulated_pattern[[1]])
      # Calculate macro-averaged AUC
      true = cell_pattern_df$marks
      preds = pred_intensities
      
      total_auc <- c()
      num_classes <- length(levels(true))
      for (i in seq_len(num_classes)) {
        true_binary <- ifelse(true == levels(true)[i], 1, 0)
        auc <- auc(roc(response = true_binary, predictor = preds[,i]))
        total_auc <- c(total_auc, auc)
      }
      macro_avg_auc <- mean(total_auc)
      print(paste("Macro-Averaged AUC: ", macro_avg_auc))
      
      original_freq = table(cell_pattern_df$marks) / length(cell_pattern_df$marks)
      weighted_macro_avg_auc <- 0
      num_classes <- length(levels(true))
      for (i in seq_len(num_classes)) {
        true_binary <- ifelse(true == levels(true)[i], 1, 0)
        auc <- auc(roc(response = true_binary, predictor = preds[,i])) * original_freq[i]
        # print(auc)
        weighted_macro_avg_auc <- weighted_macro_avg_auc + auc
      }
      print(paste("Weighted Macro-Averaged AUC: ", weighted_macro_avg_auc))
      
      # Calculate micro-averaged AUC
      all_preds <- c()
      all_true <- c()
      for (i in seq_len(num_classes)) {
        true_binary <- ifelse(true == levels(true)[i], 1, 0)
        all_preds <- c(all_preds, preds[,i])
        all_true <- c(all_true, true_binary)
      }
      micro_avg_auc <- auc(roc(response = all_true, predictor = all_preds))[1]
      print(paste("Micro-Averaged AUC: ", micro_avg_auc))
      
      aucroc = list(macro_avg_auc, weighted_macro_avg_auc, micro_avg_auc, total_auc)
      filename3 = file.path(data_dir, TMC, tissue, 'random_resample', paste('AUCROC_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))
      save(aucroc, file = filename3)
      
      print(aucroc) 
      
      # load(file.path(data_dir, TMC, tissue, 'random_resample', paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))
      # pattern_test = Quad_all_all$moadf
      # print(pattern_test[1:5,])
      # results = get_devi(pattern_test, coef, fmla, family)
      # print(resample_percent)
      # print(results)
      # save(results, file = file.path(data_dir, TMC, tissue, 'random_resample', paste('devi_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_', tissue, '_', tissue, '_across3_', img_idx, '_resample5_', resample_percent,'_self_dummy_grid_eps_20.Rda', sep = '')))
    }
  }
}


