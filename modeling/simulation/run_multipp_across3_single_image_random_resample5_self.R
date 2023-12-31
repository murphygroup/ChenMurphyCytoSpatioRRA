args = commandArgs(TRUE)
data_dir = argv[1]
for (img_idx in 1:30){
for (resample_percent in c(0, 50, 100, 150, 200,250,300,350,400)){

for (radius in c(100,200,300,400,500)){
  for (tissue in c('SPLEEN','LN','THYMUS','LI','SI')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        for (hr in c(1)) {
          if (tissue == 'LI' | tissue == 'SI') {
            TMC = 'Stanford'
          } else{
            TMC = 'Florida'
          }

if (file.exists(file.path(data_dir, TMC, tissue, 'random_resample', paste('random_simulated_pattern_', cluster_num, '_500_1_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))){
              if (file.exists(file.path(data_dir, TMC, tissue, 'random_resample', paste('quad_', cluster_num, '_100-500_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))) == F){
              if (file.exists(file.path(data_dir, TMC, tissue, 'random_resample', paste('quad_', cluster_num, '_', radius ,'_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = ''))) == F){

print(file.path(data_dir, TMC, tissue, 'random_resample', paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_', resample_percent, '_self_dummy_grid_eps_20.Rda', sep = '')))
             system(paste('Rscript multipp_across3_single_image_random_resample5_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' ', resample_percent, ' &', sep = ''))


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
