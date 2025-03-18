# Set data directory
args = commandArgs(TRUE)
data_dir = args[1]

# Define parameter ranges
tissues <- c('LN', 'SPLEEN', 'THYMUS', 'LI', 'SI')
cluster_nums <- c(5)
intensity_types <- c('total')
radii <- c(100, 200, 300, 400, 500)
hr_values <- c(1)

# Iterate over all parameter combinations
for (tissue in tissues) {
  for (cluster_num in cluster_nums) {
    for (intensity_type in intensity_types) {
      for (radius in radii) {
        
        # Determine Tissue Mapping Center (TMC)
        TMC <- if (tissue %in% c('SI', 'LI')) 'Stanford' else 'Florida'
        
        for (hr in hr_values) {
          # Define the quad file path to check existence
          quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))
          
          # Proceed if the quad file exists
          if (file.exists(quad_file)) {
            
            print(quad_file)
            
            system(paste0(
              'srun -p model1,model2,model3,model4,pool1,pool3-bigmem,gpu -o ',
              file.path(data_dir, TMC, tissue, paste0('model_single_range_all_images_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_self.out')),
              ' -n 1 -c 4 --mem 85G Rscript ./modeling/heterogeneity/train_single_range_all_images_across3_self.R ',
              data_dir, ' ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' &'
            ))
          }
        }
      }
    }
  }
}
