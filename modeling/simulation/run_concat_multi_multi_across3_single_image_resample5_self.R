data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
for (radius in c(500)){
  for (tissue in c('SPLEEN','LN','THYMUS','LI','SI')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }
          for (img_idx in 1:30){
            for (resample_percent in c(0,50,100,150,200,250,300,350,400)){
                #print(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))

if (file.exists(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))){
                #print(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))
if (!file.exists(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))){
#if (T) {
                print(file.path(data_dir, TMC, tissue, 'resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))
                #print(paste('srun -p model1,model2,model3,pool1,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'_self.out', sep = '')), ' -n 1 -c 4 --mem 46G Rscript multipp_across3_single_image.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))

                system(paste('srun -p short1,model1,pool1,gpu,model2,model3,model4 -o ', file.path(data_dir, TMC, tissue, 'resample', paste('concat_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self.out', sep = '')), ' -n 1 -c 16 --mem 120G Rscript concat_multi_multi_quad_single_image_across3_resample5_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', resample_percent, ' &', sep = ''))
}
               }
              }
            }
          } 
        }
      }
    }
  }

