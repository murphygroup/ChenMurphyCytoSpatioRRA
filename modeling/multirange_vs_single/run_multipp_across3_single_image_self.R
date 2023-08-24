args = commandArgs(TRUE)
data_dir = argv[1]
for (radius in c(100,200,300,400,500)){
  for (tissue in c('LN','SPLEEN','LI','SI','THYMUS')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }
          for (img_idx in 1:30){
            if (file.exists(file.path(data_dir, TMC, tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx, '.Rda', sep = '')))){
             if (file.exists(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'_self_dummy_grid_eps_20.Rda', sep = ''))) == F){
print(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'_self_dummy_grid_eps_20.Rda', sep = '')))

             system(paste('Rscript ./modeling/multirange_vs_single/multipp_across3_single_image_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))
              }
            }
          }
        }
      }
    }
  }
}
