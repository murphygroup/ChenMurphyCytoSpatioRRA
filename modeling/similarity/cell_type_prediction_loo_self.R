library(spatstat)
library(spatstat.utils)
library(spatstat.data)
# library(spatstat.random)
library(spatstat.geom)
library(dplyr)
library(proxy)
library(fields)
library(pROC)

set.seed(1)
script_dir = '/home/haoranch/projects/HuBMAP/ppm/script'
# script_dir = '/home/hrchen/Documents/Research/hubmap/ppm/script'

# setwd(script_dir)
# source('rmh.default_new.R')
# source('rmhEngine_new.R')


get_random_sampled_marks = function(types, num, prob){
  sampled_marks = sample(types, num, TRUE, prob)
  return(sampled_marks)
}

library(data.table)

calc_dist_optimized = function(new_points, existing_points){
  existing_points$marks = as.integer(existing_points$marks)
  existing_points <- as.data.table(existing_points) # convert to data.table
  names(existing_points) <- c("x", "y", "marks") # change according to your column names
  
  new_points = data.frame(t(new_points))
  x = new_points[1, 1]
  y = new_points[1, 2]
  xmin <- max(0, x-radius)
  ymin <- max(0, y-radius)
  xmax <- min(cell_pattern_window$xrange[2], x+radius)
  ymax <- min(cell_pattern_window$yrange[2], y+radius)
  
  # efficient subsetting using data.table
  existing_points_within_range <- existing_points[x > xmin & x < xmax & y > ymin & y < ymax]
  
  eu_dist <- t(rdist(new_points, existing_points_within_range[,1:2]))
  
  # Add eu_dist as a column to existing_points_within_range
  existing_points_within_range[, eu_dist := eu_dist]
  
  existing_points_within_range <- existing_points_within_range[eu_dist < 500 & eu_dist > 0, .(marks, eu_dist)]
  
  mark_count_row = as.integer(existing_points_within_range[, marks])
  mark_count_col = existing_points_within_range[, eu_dist] %/% 100 + 1
  mark_count = table(mark_count_row, mark_count_col)
  mark_count_matrix = matrix(0, nrow = 5, ncol = 5)
  # for (k in 1:nrow(mark_count)) {
  #   for (l in 1:ncol(mark_count))
  #   mark_count_matrix[as.numeric(rownames(mark_count)[k]),as.numeric(colnames(mark_count)[l])] <- mark_count[k,l]
  # }
  mark_count_matrix[as.matrix(expand.grid(as.numeric(rownames(mark_count)), as.numeric(colnames(mark_count))))] <- as.vector(mark_count)
  # print(mark_count_matrix)
  
  return(mark_count_matrix)
}



calc_dist = function(new_points, existing_points){
  new_points = data.frame(t(new_points))
  x = new_points[1, 1]
  y = new_points[1, 2]
  xmin = max(0, x-radius)
  ymin = max(0, y-radius)
  xmax = min(cell_pattern_window$xrange[2], x+radius)
  ymax = min(cell_pattern_window$yrange[2], y+radius)
  existing_points_within_range = existing_points[which(existing_points[,1] > xmin & existing_points[,1] < xmax & existing_points[,2] > ymin & existing_points[,2] < ymax),]
  eu_dist = rdist(new_points, existing_points_within_range[,1:2])
  existing_points_within_range = cbind(existing_points_within_range, c(eu_dist))
  existing_points_within_range = existing_points_within_range[which(existing_points_within_range[,4]<500 & existing_points_within_range[,4] > 0),3:4]
  mark_count = matrix(0, nrow = 5, ncol = 5)
  for (i in 1:dim(existing_points_within_range)[1]){
    mark_count_row = as.integer(existing_points_within_range[i, 1])
    mark_count_col = existing_points_within_range[i, 2] %/% 100 + 1
    # print(existing_points_within_range[i, 2])
    mark_count[mark_count_row, mark_count_col]  = mark_count[mark_count_row, mark_count_col] + 1
  }
  return(mark_count)
}


args = commandArgs(TRUE)
train_tissue = args[1]
test_tissue = args[2]
intensity_type = args[3]
cluster_num = as.numeric(args[4])
radius = as.numeric(args[5])
hr = as.numeric(args[6])
img_idx1 = as.numeric(args[7])
# train_tissue = 'LI'
# intensity_type = 'mean'
# num = 5
# r = 50

if (train_tissue == 'SI' | train_tissue == 'LI') {
  TMC1 = 'Stanford'
} else{
  TMC1 = 'Florida'
}

if (test_tissue == 'SI' | test_tissue == 'LI') {
  TMC2 = 'Stanford'
} else{
  TMC2 = 'Florida'
}


data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
# data_dir = '/data/PPM/HUBMAP_DATA_new'


load(file = file.path(data_dir, TMC1, train_tissue, paste('fmla_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_looidx_', img_idx1, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC1, train_tissue, paste('family_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_looidx_', img_idx1, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC1, train_tissue, paste('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_looidx_', img_idx1, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))





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

# B[2] = -11 

# types = c('0', '1', '2', '3', '4')
types = c(0, 1, 2, 3, 4)


calc_mark_prob = function(counts){
  m_prob = c()
  for (m in types){
    m_int = as.integer(m)+1
    G_type = G[[m_int]]
    
    # counts_type = counts[-seq(m_int, m_int+20, 5)]
    counts_type = counts
    interact_term = sum(G_type * counts_type)
    intensity_term = B[m_int]
    
    m_prob = c(m_prob, unname(exp(intensity_term + interact_term)))
    # print(intensity_term + interact_term)
    # m_prob = c(m_prob, unname(interact_term))
    # print(m_int)
    # print(G_type)
    # print(counts_type)
    # print(c(exp(interact_term), exp(intensity_term)))
    
  }
  # print(m_prob)
  
  m_prob = m_prob / sum(m_prob)
  # print(m_prob)
  # print(types)
  final_type = sample(types, 1, prob = m_prob)
  # return(final_type)
  # final_type = which.max(m_prob)-1
  # final_type = which.min(abs(m_prob - Y[point_current_idx]))-1
  return(list(final_type, m_prob))
}


for (img_idx2 in 1:30){
  if (file.exists(file.path(data_dir, TMC2, test_tissue, paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx2, '_self_dummy_grid_eps_20.Rda', sep = '')))){
    if (train_tissue != test_tissue | img_idx1 == img_idx2){
      if (!file.exists(file.path(data_dir, TMC1, train_tissue, paste('AUCROC_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_dummy_grid_eps_20.Rda', sep = '')))){
        # if (!file.exists(file.path(data_dir, TMC1, train_tissue, paste('pred_intensities_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_dummy_grid_eps_20.Rda', sep = '')))){
        if (T){
          print(img_idx2)
          load(file = file.path(data_dir, TMC2, test_tissue, paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx2, '_self_dummy_grid_eps_20.Rda', sep = '')))
          load(file.path(data_dir, TMC2, test_tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx2, '.Rda', sep = '')))
          
          Y = Quad_all_all$moadf$.mpl.Y[which(Quad_all_all$moadf$.mpl.Y != 0)]
          
          cell_pattern = cell_list[[1]]
          cell_pattern_df = data.frame(cell_pattern)
          cell_pattern_window <<- cell_pattern$window
          
          
          marks = c()
          intensities = c()
          for (point_current_idx in 1:nrow(cell_pattern_df)){
            # print(point_current_idx)
            point_current = cell_pattern_df[point_current_idx, ]
            counts_current = t(apply(point_current[,1:2], 1, calc_dist_optimized, cell_pattern_df))
            mark_current = apply(counts_current, 1, calc_mark_prob)
            marks = c(marks, mark_current[[1]][[1]])
            intensities = rbind(intensities, mark_current[[1]][[2]])
          }
          
          pred_cell_type = data.frame(x = cell_pattern_df$x, y = cell_pattern_df$y, marks = marks)
          
          print(dim(pred_cell_type))
          filename1 = file.path(data_dir, TMC1, train_tissue, paste('pred_cell_type_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_dummy_grid_eps_20.Rda', sep = ''))
          save(pred_cell_type, file = filename1)
          
          pred_intensities = data.frame(intensities)
          filename2 = file.path(data_dir, TMC1, train_tissue, paste('pred_intensities_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_dummy_grid_eps_20.Rda', sep = ''))
          save(pred_intensities, file = filename2)
        }
        else{
          load(file.path(data_dir, TMC2, test_tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx2, '.Rda', sep = '')))
          cell_pattern = cell_list[[1]]
          cell_pattern_df = data.frame(cell_pattern)
          
          filename2 = file.path(data_dir, TMC1, train_tissue, paste('pred_intensities_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_dummy_grid_eps_20.Rda', sep = ''))
          # pred_intensities = read.csv(filename2)
          load(filename2)
          #print(pred_intensities)
          # save(pred_intensities, file = filename2) 
        }
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
        filename2 = file.path(data_dir, TMC1, train_tissue, paste('AUCROC_loo_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', train_tissue, '_', img_idx1, '_', test_tissue, '_', img_idx2, '_self_dummy_grid_eps_20.Rda', sep = ''))
        save(aucroc, file = filename2)
      }
    }
  }
}

