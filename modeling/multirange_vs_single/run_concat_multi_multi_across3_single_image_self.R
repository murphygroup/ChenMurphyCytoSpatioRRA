# Set data directory
args = commandArgs(TRUE)
data_dir = args[1]

# Define parameter ranges
radii <- c(500)
tissues <- c('LN', 'SPLEEN', 'THYMUS', 'LI', 'SI')
cluster_nums <- c(5)
intensity_types <- c('total')
hr_values <- c(1)
img_indices <- 1:30  # Range of image indices

# Iterate over all parameter combinations
for (radius in radii) {
  for (tissue in tissues) {
    for (cluster_num in cluster_nums) {
      for (intensity_type in intensity_types) {
        for (hr in hr_values) {
          
          # Determine Tissue Mapping Center (TMC)
          TMC <- if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'
          
          for (img_idx in img_indices) {
            # Define the original quad file path to check existence
            original_quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_quad_d_no_between_dummy.Rda'))
            
            # Proceed only if the original file exists
            if (file.exists(original_quad_file)) {
              
              # Define new quad file path with "100-" radius notation
              new_quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_quad_d_no_between_dummy.Rda'))
              
              # Check if new file is missing or older than 7 days
              if (!file.exists(new_quad_file) || (Sys.time() - file.info(new_quad_file)$mtime) > as.difftime(7, units = "days")) {
                
                # Print file path for logging
                print(original_quad_file)
                
                system(paste0(
                  'srun -p short1,model1,model2,model3,model4,pool1,pool3-bigmem -t 2:00:00 -o ',
                  file.path(data_dir, TMC, tissue, paste0('concat_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self.out')),
                  ' -n 1 -c 4 --mem 46G Rscript ./modeling/heterogeneity/concat_multi_multi_quad_single_image_across3_self.R ',
                  data_dir, ' ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &'
                ))
              }
            }
          }
        }
      }
    }
  }
}
