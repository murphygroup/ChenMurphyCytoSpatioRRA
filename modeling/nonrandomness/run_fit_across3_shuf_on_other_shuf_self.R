# Set data directory
args = commandArgs(TRUE)
data_dir = args[1]

for (tissue in c('LN', 'SPLEEN', 'THYMUS', 'SI', 'LI')) {
  for (cluster_num in 5:5) {
    for (intensity_type in c('total')) {
      for (radius in c(100)) {
        
        # Determine Tissue Mapping Center (TMC)
        TMC = if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'
        
        # Set shuffle range dynamically
        shuf_range = if (tissue %in% c('LI', 'SI')) 1:2500 else 1:150
        
        for (shuffle_num in shuf_range) {
          hr = 1
          
          # Define required file paths
          quad_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuffle_num, '_self_quad_d_no_between_dummy.Rda'))
          coef_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('coef_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuffle_num, '_self_quad_d_no_between_dummy.Rda'))
          output_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('model_shuf_on_other_shuf_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', shuffle_num, '_self.out'))
          
          # Check if required quad and coefficient files exist
          if (file.exists(quad_file) && file.exists(coef_file)) {
            print(file.path(data_dir, TMC, tissue, paste0('fit_shuf_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_', shuffle_num, '.out')))
            
            # Submit job using system call
            system(paste0(
              'srun -p pool1,model1,model2,model3,model4,pool1,pool3-bigmem,gpu,interactive -o ',
              output_file,
              ' -n 1 -c 4 --mem 46G Rscript ./modeling/nonrandomness/fit_across3_shuf_on_other_shuf_self.R ',
              data_dir, ' ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', shuffle_num, ' &'
            ))
          }
        }
      }
    }
  }
}
