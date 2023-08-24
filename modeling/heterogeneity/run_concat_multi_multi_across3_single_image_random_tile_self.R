args = commandArgs(TRUE)
data_dir = argv[1]
for (radius in c(500)){
  for (tissue in c('LI','SI','THYMUS','LN','SPLEEN')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }
          for (img_idx in 1:30){
            for (tile_num in 1:5){
            for (tile_side in c(2500,5000)){
		if (file.exists(file.path(data_dir, TMC, tissue, 'random_tile', paste('quad_', cluster_num, '_', '100', '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))){
                if (!file.exists(file.path(data_dir, TMC, tissue, 'random_tile', paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))){
                print(file.path(data_dir, TMC, tissue, 'random_tile', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_tile_', tile_side, '_', tile_num, '_self_dummy_grid_eps_20.Rda', sep = '')))
                #print(paste('srun -p model1,model2,model3,pool1,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'_self.out', sep = '')), ' -n 1 -c 4 --mem 46G Rscript multipp_across3_single_image.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))

                system(paste('Rscript ./modeling/heterogeneity/concat_multi_multi_quad_single_image_across3_random_tile_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', tile_side, ' ', tile_num, ' &', sep = ''))
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
