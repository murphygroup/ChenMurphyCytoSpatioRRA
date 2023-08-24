library(spatstat)
library(spatstat.utils)
library(spatstat.data)


print('Reading cell distributions...')

get_cells = function (tissue, intensity_type, cluster_num) {
  if (tissue == 'LI' | tissue == 'SI'){
    TMC = 'Stanford'
  } else{
    TMC = 'Florida'
  }
  args = commandArgs(TRUE)
  data_dir = args[1]
  cell_type_files = Sys.glob(file.path(data_dir, TMC, tissue, '*', paste('*_cluster_', as.character(cluster_num), '_label_', intensity_type, '_across3.csv', sep = '')))
  print(cell_type_files)
  cell_list_all = list()
  cell_pp_idx = 1
  for (i in 1:length(cell_type_files)) {
    #print(i)
    cell_type_file = cell_type_files[i]
    cell_type = read.csv(cell_type_file)
    cell_type_labels = sort(unique(cell_type[, 2]))
    cell_type_num = length(cell_type_labels)
    cell_type_label_num = c()
    for (c in cell_type_labels){
      cell_type_label_num = c(cell_type_label_num, length(which(cell_type[, 2] == c)))
    }

    if (cell_type_num == cluster_num){
      cell_list = list()
      region_name = paste(strsplit(basename(cell_type_file), split = '_')[[1]][1], '_stitched_expressions.ome.tiff', sep = '')
      print(cell_type_files[i])
      cell_center_file = file.path(dirname(as.character(cell_type_file)), paste(region_name, '-cell_centers.csv', sep = ''))
      img_shape_file = file.path(dirname(as.character(cell_type_file)), paste(region_name, 'img_shape.txt', sep = '-'))
      img_name = paste(basename(dirname(as.character(cell_type_file))), region_name, sep = '/')
      cell_center = read.csv(cell_center_file)
      cell_info = merge(cell_center, cell_type)
      img_shape = c(read.delim(img_shape_file, header = F, sep = '\n')[[1]])

      cell_pp = ppp(c(cell_info[,2]), c(cell_info[,3]), window = owin(c(0, img_shape[1]), c(0, img_shape[2])), marks = factor(cell_type[,2]))
      cell_list[[img_name]] = cell_pp
      cell_list_all[[img_name]] = cell_pp
      print(cell_pp_idx)
      print(cell_list[[img_name]])
      print(cell_type_label_num)

      filename1 = file.path(data_dir, TMC, tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', cell_pp_idx ,'.Rda', sep = ''))
#       save(cell_list, file=filename1)
      cell_pp_idx = cell_pp_idx + 1
    }
  }
       filename1 = file.path(data_dir, TMC, tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3.Rda', sep = ''))
#        save(cell_list_all, file=filename1)

}

for (cluster_num in 5:5){
  for (tissue in c('LI','SI','THYMUS','LN','SPLEEN')) {
    for (intensity_type in c('total')) {
      get_cells(tissue, intensity_type, cluster_num)
    }
  }
}

