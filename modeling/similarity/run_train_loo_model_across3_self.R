
args = commandArgs(TRUE)
data_dir = args[1]

  for (tissue in c('LI','SI','LN','SPLEEN','THYMUS')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
	for (radius in c(500)){

        if (tissue == 'SI' | tissue == 'LI') {
          TMC = 'Stanford'
        } else{
          TMC = 'Florida'
        }

        hr = 1
        #if(T){
            if (file.exists(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))){
             #if (!file.exists(file.path(data_dir, TMC, tissue, paste('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_looidx_', outer_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))){ 
                print(file.path(data_dir, TMC, tissue, paste('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
                system(paste('Rscript ./modeling/similarity/train_loo_model_across3_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' &', sep = ''))

              #}
            #}
          }
        }
      }
    }
  }

