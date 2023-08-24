library(spatstat)
library(spatstat.utils)
library(spatstat.data)
# library(spatstat.random)
library(spatstat.geom)
library(dplyr)
library(proxy)

script_dir = '/home/haoranch/projects/HuBMAP/ppm/script'
# script_dir = '/home/hrchen/Documents/Research/hubmap/ppm/script'

# setwd(script_dir)
# source('rmh.default_new.R')
# source('rmhEngine_new.R')


get_random_sampled_marks = function(types, num, prob){
  sampled_marks = sample(types, num, TRUE, prob)
  return(sampled_marks)
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
  eu_dist = dist(new_points, existing_points_within_range[,1:2])
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

# 
args = commandArgs(TRUE)
tissue = args[1]
intensity_type = args[2]
cluster_num = as.numeric(args[3])
radius = as.numeric(args[4])
hr = as.numeric(args[5])
img_idx = as.numeric(args[6])
resample_percent = as.numeric(args[7])
# 
# tissue = 'SI'
# intensity_type = 'total'
# cluster_num = 5
# radius = 500
# hr = 1
# img_idx = 11
# resample_percent = 200

if (tissue == 'SI' | tissue == 'LI') {
  TMC = 'Stanford'
} else{
  TMC = 'Florida'
}
data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
# data_dir = '/data2/PPM/HUBMAP_DATA_new'

set.seed(3)

load(file = file.path(data_dir, TMC, tissue, paste('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_dummy_grid_eps_20.Rda', sep = '')))
load(file = file.path(data_dir, TMC, tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx, '.Rda', sep = '')))


# load('/data2/PPM/HUBMAP_DATA_new/Stanford/SI/coef_5_100-500_1_total_across3.Rda')
# load('/data2/PPM/HUBMAP_DATA_new/Stanford/SI/family_5_100_1_total_across3.Rda')
# load('/data2/PPM/HUBMAP_DATA_new/Stanford/SI/fmla_5_100_1_total_across3.Rda')
# load('/data2/PPM/HUBMAP_DATA_new/Stanford/SI/cell_list_5_total_across3_1.Rda')


cell_pattern = cell_list[[1]]
cell_pattern_df = data.frame(cell_pattern)
# mark_types = unique(cell_pattern$marks)
# mark_prob = 

cell_pattern_window <<- cell_pattern$window


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
    B_current = exp(coef[i])
  } else{
    B_current = exp(coef[i] + coef[1])
  }
  B = c(B, B_current)
}

types = c('0', '1', '2', '3', '4')


calc_mark_prob = function(counts){
  m_prob = c()
  for (m in types){
    m_int = as.integer(m)+1
    G_type = G[[m_int]]
    
    # counts_type = counts[-seq(m_int, m_int+20, 5)]
    counts_type = counts
    interact_term = exp(sum(G_type * counts_type))
    intensity_term = B[m_int]
    m_prob = c(m_prob, unname(intensity_term * interact_term))
    # m_prob = c(m_prob, unname(interact_term))
    # print(m_int)
    # print(G_type)
    # print(counts_type)
    # print(intensity_term)
    # print(interact_term)
  }
  # print(m_prob)
  m_prob = m_prob / sum(m_prob)
  # print(m_prob)
  return(factor(sample(types, 1, prob = m_prob), levels = types))
}

# counting = t(apply(cell_pattern_df[1:5000, 1:2], 1, calc_dist, cell_pattern_df))
# calc_dist(cell_pattern_df[1,1:2], cell_pattern_df)
# calc_mark_prob(counting[1,])

# test = apply(counting, 1, calc_mark_prob)

# resampling parameters and initialization

#model_mark_prob = B / sum(B)

model_mark_prob = (dplyr::count(cell_pattern_df, marks)$n) / dim(cell_pattern_df)[1]
# cell_pattern_split_list = split(cell_pattern_df, sample(1:batch_num, nrow(cell_pattern_df), replace=T))
pattern_exist = cell_pattern_df
pattern_exist[, 3] = as.factor(sample(types, nrow(cell_pattern_df), replace=T, model_mark_prob))
cell_pattern_num = dim(cell_pattern_df)[1]
resample_num = ceiling(cell_pattern_num * resample_percent / 100)
print(resample_num)
if (resample_num != 0){
  for (i in 1:resample_num){
    # print(i)
    point_current_idx = sample(1:cell_pattern_num, 1)
    # print(c(i, point_current_idx))
    point_current =  pattern_exist[point_current_idx, ]
    counts_current = t(apply(point_current[,1:2], 1, calc_dist, pattern_exist))
    mark_prob_current = apply(counts_current, 1, calc_mark_prob)
    point_current[, 3] = mark_prob_current
    pattern_exist[point_current_idx, ] = point_current
  }
}


marked_simulated_pattern = list()
marked_simulated_pattern[[1]] = ppp(c(pattern_exist[,1]), c(pattern_exist[,2]), window = cell_pattern_window, marks = factor(pattern_exist[,3], levels = types))

filename1 = file.path(data_dir, TMC, tissue, 'resample', paste('marked_simulated_pattern_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))

save(marked_simulated_pattern, file = filename1)


