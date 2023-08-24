
args = commandArgs(TRUE)
data_dir = argv[1]



for (cluster_num in c(5)) {
  for (radius in c(100,200,300,400,500)) {  
    for (tissue in c('SPLEEN','LN','THYMUS','LI','SI')) {
      for (intensity_type in c('total')) {
        if (tissue == 'SI' | tissue == 'LI') {
          TMC = 'Stanford'
        } else{
          TMC = 'Florida'
        }
        for (hr in c(1)) {
          if (file.exists(file.path(data_dir, TMC, tissue, paste('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_50.Rda', sep = ''))) == F){
          #if (T){
            if (T){
              print(file.path(data_dir, TMC, tissue, paste('combine_', cluster_num, '_', radius, '_', hr, '_',  intensity_type, '_across3.Rda', sep = '')))
              #system(paste('srun -p model1,model2,model3,model4,pool3-bigmem,gpu,pool1,interactive -o ', file.path(data_dir, TMC, tissue, paste('combine_', cluster_num, '_', radius, '_', hr, '_',  intensity_type, '_across3_self.out', sep = '')), ' -n 1 -c 4 --mem 46G Rscript combine_quad_across3_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' &', sep = '')) 
          system(paste('Rscript ./modeling/multirange_vs_single/train_models_nested.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' &', sep = ''))
            }
          }
        }
      }
    }
  }
}
