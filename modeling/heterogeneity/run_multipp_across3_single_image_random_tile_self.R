args = commandArgs(TRUE)
data_dir = argv[1]
for (tile_num in 1:5){

for (radius in c(100,200,300,400,500)){
  for (tissue in c('LI','SI','LN','THYMUS','SPLEEN')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }
          for (img_idx in 1:30){
            for (tile_side in c(5000,2500)){
             #print(file.path(data_dir, TMC, tissue, 'random_tile', paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '.Rda', sep = '')))
              if (file.exists(file.path(data_dir, TMC, tissue, 'random_tile', paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '.Rda', sep = '')))){
              if (file.exists(file.path(data_dir, TMC, tissue, 'random_tile', paste('quad_', cluster_num, '_100-500_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = ''))) == F){
            system(paste('Rscript ./modeling/heterogeneity/multipp_across3_single_image_random_tile_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', tile_side, ' ', tile_num, ' &', sep = ''))
}
}
              }
            }
          } 
        }
      }
    }
  }
}
