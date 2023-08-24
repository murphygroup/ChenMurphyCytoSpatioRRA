library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)

# data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
# args = commandArgs(TRUE)
# tissue = args[1]
# intensity_type = args[2]
# n = as.numeric(args[3])
# r = as.numeric(args[4])
# img_idx = args[6]
# cluster_num = 5

get_tile_pattern = function(original_pattern, rect){
  original_coords_x = original_pattern$x
  original_coords_y = original_pattern$y
  original_marks = original_pattern$marks
  original_pattern_df = data.frame(x = original_coords_x, y = original_coords_y, marks = original_marks)
  original_pattern_tile_df = original_pattern_df[which(original_pattern_df$x >= rect[1] & original_pattern_df$x <= rect[3] & original_pattern_df$y >= rect[2] & original_pattern_df$y <= rect[4]),]
  original_pattern_tile = ppp(c(original_pattern_tile_df[,1]), c(original_pattern_tile_df[,2]), window = owin(c(rect[1], rect[3]), c(rect[2], rect[4])), marks = factor(original_pattern_tile_df[,3]))
  return(original_pattern_tile)
}

tiling = function(tissue, tile_side, total_tile_num){

  tile_cell_num_df <- data.frame(matrix(ncol = 5, nrow = 0))
  
  colnames(tile_cell_num_df) <- c('img_idx', 'tile_idx', 'cell_num', 'cell_type_num', 'keep')
  data_dir = '/data/PPM/HUBMAP_DATA_new'
  if (tissue == 'LI' | tissue == 'SI'){
    TMC = 'Stanford'
  } else{
    TMC = 'Florida'
  }
  dir.create(file.path(data_dir, TMC, tissue, 'random_tile'), showWarnings = FALSE)
  cell_list_files = Sys.glob(file.path(data_dir, TMC, tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3_*.Rda', sep = '')))
  # print(cell_list_files)
  cell_list_file_num = length(cell_list_files)
  for (img_idx in 1:cell_list_file_num){
    # for (img_idx in 8:8){
    cell_list_file = cell_list_files[img_idx]
    file_name = basename(cell_list_file)
    load(cell_list_file)
    cell_pattern = cell_list[[1]]
    
    
    xmin_image = 0
    ymin_image = 0
    xmax_image = cell_pattern$window$xrange[2]
    ymax_image = cell_pattern$window$yrange[2]
    
    xmin_tile_range = xmin_image + 500 + tile_side / 2
    ymin_tile_range = ymin_image + 500 + tile_side / 2
    xmax_tile_range = xmax_image - 500 - tile_side / 2
    ymax_tile_range = ymax_image - 500 - tile_side / 2
    sampled_tile_num = 0
    all_sampled_tile_num = 0
    while (sampled_tile_num < total_tile_num & all_sampled_tile_num < 1000){
      all_sampled_tile_num = all_sampled_tile_num + 1
      xcenter_tile = sample(xmin_tile_range:xmax_tile_range, 1)
      ycenter_tile = sample(ymin_tile_range:ymax_tile_range, 1)
      
      xmin_tile = xcenter_tile - 500 - tile_side / 2
      xmax_tile = xcenter_tile + 500 + tile_side / 2
      ymin_tile = ycenter_tile - 500 - tile_side / 2
      ymax_tile = ycenter_tile + 500 + tile_side / 2    
      
      tile_rect = c(xmin_tile, ymin_tile, xmax_tile, ymax_tile)
      tile_real_rect = c(xmin_tile + 500, ymin_tile + 500, xmax_tile - 500, ymax_tile - 500)
      
      cell_pattern_tile = get_tile_pattern(cell_pattern, tile_rect)
      cell_real_pattern_tile = get_tile_pattern(cell_pattern, tile_real_rect)
      
      cell_type_num = length(sort(unique(cell_real_pattern_tile$marks)))
      
      cell_list = list()
      cell_list[['tile']] = cell_pattern_tile 
      print(cell_pattern_tile)
      # print(split(cell_list[[1]]))
      
      if (cell_pattern_tile$n > cell_pattern$n / total_tile_num / 5){
        if (cell_type_num == cluster_num){
          sampled_tile_num = sampled_tile_num + 1
          filename1 = file.path(data_dir, TMC, tissue, 'random_tile', paste(substr(file_name, 1, nchar(file_name)-4), '_tile_', tile_side, '_', sampled_tile_num, '.Rda', sep = ''))
          save(cell_list, file=filename1)
        }
    }
  # print(test2 / cell_list_file_num / t / t)
  # return(tile_cell_num_df)
    }
  }
}
cluster_num = 5
intensity_type = "total"
set.seed(1)
tissue_list = c('SPLEEN','THYMUS','LI','SI')
# tissue_list = c('LN')
total_tile_num = 5
cell_pattern_tile_num = c()
  for (tissue in tissue_list) {
    for (tile_side in c(5000,2500)){
    # for (tissue in c('SPLEEN')) {
      tiling(tissue, tile_side, total_tile_num)
    }
}



