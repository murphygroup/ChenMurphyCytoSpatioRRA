data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
for (radius in c(500)){
  for (tissue in c('SI','LI','LN','SPLEEN','THYMUS')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }
           if (file.exists(file.path(data_dir, TMC, tissue, paste('cell_list_', cluster_num, '_', intensity_type, '_across3.Rda', sep = '')))){

	      for (resample_percent in c(0,100,200,300,400)){
            
#if (T){
#file_path = file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_resample5_percentage_', resample_percent, '_seed_', random_seed, '_self_quad_d_no_between_dummy.Rda', sep = ''))

# file_path = file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_resample5_percentage_', resample_percent, '_seed_', random_seed, '_self_quad_d_no_between_dummy.Rda', sep = ''))
              for (cell_num in c(10000,15000,20000)){

file_path = file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_resample5_', resample_percent, '_cell_num_', cell_num, '_self_global_same_pattern_quad_d_no_between_dummy.Rda', sep = ''))
            #if (file.exists(file_path)){
#file_path = file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_resample5_global_', resample_percent, '_', random_seed, '.out', sep = ''))

#if (T){
if (!file.exists(file_path) || (Sys.time() - file.info(file_path)$mtime) > as.difftime(60, units = "days")) {


print(file_path)
              #print(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'.Rda', sep = '')))
              #print(paste('srun -p model1,model2,model3,pool1,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'.out', sep = '')), ' -n 1 -c 4 --mem 46G Rscript multipp_across3_single_image.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))

           system(paste('srun -p model1,model2,pool1,model3,model4,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_resample5_global_', resample_percent, '_same_pattern_', cell_num,'.out', sep = '')), ' -n 1 -c 4 --mem 20G Rscript simulate_random_resampling5_global_self_same_pattern.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', resample_percent, ' ', cell_num, ' &', sep = ''))
}              
		}
}
}
            }
          } 
       }
    }
}
