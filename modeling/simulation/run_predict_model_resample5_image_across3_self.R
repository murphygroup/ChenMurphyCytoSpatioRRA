
args = commandArgs(TRUE)
data_dir = argv[1]
tissue_list = c('LN','SPLEEN','THYMUS','LI','SI')
#tissue_list = c('LN')
  for (tissue in tissue_list) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
	for (radius in c(500)){
        for (img_idx in 1:30){
        if (tissue == 'SI' | tissue == 'LI') {
          TMC = 'Stanford'
        } else{
          TMC = 'Florida'
        }

        hr = 1
            if (file.exists(file.path(data_dir, TMC, tissue, paste('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_dummy_grid_eps_20.Rda', sep = '')))){
            if (!file.exists(file.path(data_dir, TMC, tissue, 'resample', paste('AUCROC_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_', tissue, '_', tissue, '_across3_', img_idx, '_resample5_400_self_dummy_grid_eps_20.Rda', sep = '')))){
#if(T){
system(paste('Rscript ./modeling/simulation/predict_model_resample5_image_across3_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))
               }
           } 
          }
    }
  }
}
}
