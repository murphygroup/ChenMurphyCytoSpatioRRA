
args = commandArgs(TRUE)
data_dir = args[1]
tissue_list1 <- tissue_list2 <- c('LI','SI','LN','SPLEEN','THYMUS')
#tissue_list1 = c('SPLEEN')
#tissue_list2 = c('LI')
  for (tissue1 in tissue_list1) {
    for (tissue2 in tissue_list2){
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
	for (radius in c(500)){

        if (tissue1 == 'SI' | tissue1 == 'LI') {
          TMC = 'Stanford'
        } else{
          TMC = 'Florida'
        }

        hr = 1
        for (img_idx in 1:30){
            if (file.exists(file.path(data_dir, TMC, tissue1, paste('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_looidx_', img_idx, '_across3_self_dummy_grid_eps_20.Rda', sep = '')))){
if(T){
                
system(paste('Rscript ./modeling/similarity/cell_type_prediction_loo_self.R ', tissue1, ' ', tissue2, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))
              
} 
            }
          }
        }
      }
    }
  }
}
