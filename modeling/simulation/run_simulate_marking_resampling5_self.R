args = commandArgs(TRUE)
data_dir = argv[1]
for (radius in c(500)){
  for (tissue in c('LN','SPLEEN','THYMUS','LI','SI')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }
              dir.create(file.path(data_dir, TMC, tissue, 'resample'), showWarnings = FALSE)

          for (img_idx in 1:30){
           if (file.exists(file.path(data_dir, TMC, tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx, '.Rda', sep = '')))){

	      for (resample_percent in c(0, 50, 100,150,200,250,300,350,400)){
                 output_path = file.path(data_dir, TMC, tissue, 'resample', paste('marked_simulated_pattern_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))
            if (!file.exists(output_path)){
	      print(output_path)
              #print(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'.Rda', sep = '')))
              #print(paste('srun -p model1,model2,model3,pool1,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'.out', sep = '')), ' -n 1 -c 4 --mem 46G Rscript multipp_across3_single_image.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))

           system(paste('Rscript ./modeling/simulation/simulate_marking_resampling5_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', resample_percent, ' &', sep = ''))
}              
		}

}
            }
          } 
        }
      }
    }
  }
