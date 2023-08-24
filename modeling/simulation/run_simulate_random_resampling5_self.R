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
          for (img_idx in c(1:30)){
           if (file.exists(file.path(data_dir, TMC, tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx, '.Rda', sep = '')))){
              dir.create(file.path(data_dir, TMC, tissue, 'random_resample'), showWarnings = FALSE)

	      for (resample_percent in c(0,100,200,300,400)){
                 output_path = file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))
            #if (!file.exists(output_path)){
if(T){
	      print(output_path)
              #print(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'.Rda', sep = '')))
              #print(paste('srun -p model1,model2,model3,pool1,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'.out', sep = '')), ' -n 1 -c 4 --mem 46G Rscript multipp_across3_single_image.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))

          system(paste('Rscript ./modeling/simulation/simulate_random_resampling5_self_same_pattern.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', resample_percent, ' 10000 &', sep = ''))
}              
		}

}
            }
          } 
        }
      }
    }
  
}
