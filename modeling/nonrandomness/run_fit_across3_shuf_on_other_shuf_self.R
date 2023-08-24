args = commandArgs(TRUE)
data_dir = args[1]

for (shuffle_num in 1:2500) {
  for (tissue in c('LN','SPLEEN','THYMUS','SI','LI')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
	for (radius in c(100)){
          if (tissue == 'LI' | tissue == 'SI'){
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }

         hr = 1
         if (file.exists(file.path(data_dir, TMC, tissue, 'random_shuf', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuffle_num, '_self_dummy_grid_eps_20.Rda', sep = '')))) {
           if (file.exists(file.path(data_dir, TMC, tissue, 'random_shuf', paste('coef_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuffle_num, '_self_dummy_grid_eps_20.Rda', sep = '')))) {
           #if (T){
           print(file.path(data_dir, TMC, tissue, paste('fit_shuf_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_', shuffle_num, '.out', sep = '')))
             system(paste('Rscript ./modeling/nonrandomness/fit_across3_shuf_on_other_shuf_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', shuffle_num, ' &', sep = ''))
	    }
          }
        }
      }
    }
  }
}
