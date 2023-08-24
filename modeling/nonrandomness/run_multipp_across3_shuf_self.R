args = commandArgs(TRUE)
data_dir = args[1]
for (radius in c(100)) {
  for (tissue in c('LN','SI','THYMUS','SPLEEN')){
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (shuf_num in 1:2500){
          hr = 1
          if (tissue == 'LI' | tissue == 'SI'){
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          } 

        dir.create(file.path(data_dir, TMC, tissue, 'random_shuf'), showWarnings = FALSE)
        print(file.path(data_dir, TMC, tissue, 'random_shuf'))
	if (file.exists(file.path(data_dir, TMC, tissue, 'random_shuf', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_dummy_grid_eps_20.Rda', sep = ''))) == F){

          system(paste('Rscript ./modeling/nonrandomness/multipp_across3_shuf_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', shuf_num, sep = ''))
        }
      }
    }
  }
}
}
