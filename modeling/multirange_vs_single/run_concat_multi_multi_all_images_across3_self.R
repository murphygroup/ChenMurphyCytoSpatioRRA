args = commandArgs(TRUE)
data_dir = argv[1]
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
          #if (T){
          if (file.exists(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = ''))) == F){
              print(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))
              #print(paste('srun -p model1,model2,model3,pool1,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx ,'.out', sep = '')), ' -n 1 -c 4 --mem 46G Rscript multipp_across3_single_image.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))

            system(paste('Rscript ./modeling/multirange_vs_single/concat_multi_multi_all_images_across3_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' &', sep = ''))
              

           
          } 
        }
      }
    }
  }
}
