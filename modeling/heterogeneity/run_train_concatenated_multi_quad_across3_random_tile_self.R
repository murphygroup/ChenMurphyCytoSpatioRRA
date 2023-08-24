args = commandArgs(TRUE)
data_dir = argv[1]
for (radius in c(500)){
  for (tissue in c('SI','LI','THYMUS','SPLEEN','LN')) {

    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }
          for (img_idx in 2:2){
            for (tile_side in c(2500)){
            for (tile_num in 4:4){
		if (file.exists(file.path(data_dir, TMC, tissue, 'random_tile', paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))){
if(T){

                print(file.path(data_dir, TMC, tissue, 'random_tile', paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))

                system(paste('Rscript ./modeling/heterogeneity/train_concatenated_multi_quad_across3_random_tile_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', tile_side, ' ', tile_num, ' &', sep = ''))
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
