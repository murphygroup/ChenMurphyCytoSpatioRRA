data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
for (img_idx in 1:30){
for (radius in c(100,200,300,400,500)){
  for (tissue in c('LN','SPLEEN','THYMUS','LI','SI')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }
              #for (resample_percent in seq(0, batch_num, 2)){
              for (resample_percent in c(0,50,100,150,200,250,300,350,400)){

if (file.exists(file.path(data_dir, TMC, tissue, 'resample', paste('marked_simulated_pattern_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))){
              
              if (file.exists(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_100-500_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))) == F){
if (file.exists(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius ,'_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))) == F){
print(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '.out', sep = '')))
#if (T){
              #if (file.exists(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', total_resample5_num, '_', resample_percent, '.out', sep = ''))) == F){
	      #print(paste('srun -p pool1,model1,model2,model3,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', total_resample5_num, '_', resample_percent, '.out', sep = '')), ' -n 1 -c 4 --mem 46G Rscript multipp_across3_single_image_resample.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', total_resample5_num, ' ', resample_percent, ' &', sep = ''))	
              system(paste('srun -p pool1,model1,model2,model3,pool3-bigmem,gpu -t 12:00:00 -o ', file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self.out', sep = '')), ' -n 1 -c 8 --mem 46G Rscript multipp_across3_single_image_resample5_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', resample_percent, ' &', sep = ''))              
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

