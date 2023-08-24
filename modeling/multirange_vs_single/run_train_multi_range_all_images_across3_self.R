
args = commandArgs(TRUE)
data_dir = argv[1]

  for (tissue in c('LN','SPLEEN','THYMUS','LI','SI')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        if (T){
          radius = 500

        if (tissue == 'SI' | tissue == 'LI') {
          TMC = 'Stanford'
        } else{
          TMC = 'Florida'
        }

         hr = 1
         if (file.exists(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))) {
if(T){	 
           print(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
           system(paste('Rscript ./modeling/multirange_vs_single/train_multi_range_all_images_across3_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' &', sep = ''))
	    }
          }
        }
      }
    }
  }
