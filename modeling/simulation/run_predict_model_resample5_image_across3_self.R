
data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
tissue_list = c('LN','SPLEEN','THYMUS','LI','SI')
#tissue_list = c('LN')
  for (tissue in tissue_list) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
	for (radius in c(500)){
        for (img_idx in 1:30){
	#if (T){
# for (radius in c(20)) {
#   for (tissue in c('SPLEEN')) {
#     for (cluster_num in 5:5) {
#       for (intensity_type in c('mean')) {
          #print(paste('Running', tissue, cluster_num, intensity_type, radius, sep = ' '))
          #print(paste('srun -p model3,model4 -o ', file.path(data_dir, tissue, paste('shuf_', cluster_num, '_', radius, '_', intensity_type, '_', shuffle_num, '_self.out', sep = '')), ' -n 1 -c 1 --mem 21G Rscript get_likelihood.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', shuffle_num, ' &', sep = ''))
          # print(file.path(data_dir, tissue, paste('quad_', cluster_num, '_', radius, '_', intensity_type, '_self_dummy_grid_eps_20.Rda', sep = '')))
         #if (TRUE){
        if (tissue == 'SI' | tissue == 'LI') {
          TMC = 'Stanford'
        } else{
          TMC = 'Florida'
        }

        hr = 1
            if (file.exists(file.path(data_dir, TMC, tissue, paste('coef_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_dummy_grid_eps_20.Rda', sep = '')))){
            if (!file.exists(file.path(data_dir, TMC, tissue, 'resample', paste('AUCROC_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_', tissue, '_', tissue, '_across3_', img_idx, '_resample5_400_self_dummy_grid_eps_20.Rda', sep = '')))){
#if(T){
print(paste('srun -p short1,interactive,pool1,pool3-bigmem,model1,model2,model3 -o ', file.path(data_dir, TMC, tissue, 'resample', paste('predict_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_self.out', sep = '')), ' -n 1 -c 4 --mem 21G Rscript predict_model_resampled5_image_across3.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = ''))
system(paste('srun -p pool1,model1,model2,model3,model4,pool3-bigmem -t 72:00:00 -o ', file.path(data_dir, TMC, tissue, 'resample', paste('predict_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_resample5_self.out', sep = '')), ' -n 1 -c 4 --mem 21G Rscript predict_model_resample5_image_across3_self.R ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &', sep = '')) 
               }
           } 
          }
    }
  }
}
}
