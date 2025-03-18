# Set data directory
args = commandArgs(TRUE)
data_dir = args[1]

# Define parameter ranges
radii <- c(500)
tissues <- c('SPLEEN', 'LN', 'THYMUS', 'LI', 'SI')
cluster_nums <- c(5)
intensity_types <- c('total')
hr_values <- c(1)

# Iterate over all parameter combinations
for (radius in radii) {
  for (tissue in tissues) {
    for (cluster_num in cluster_nums) {
      for (intensity_type in intensity_types) {
        for (hr in hr_values) {
          
          # Determine Tissue Mapping Center (TMC)
          TMC <- if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'
          
          # Define the quad file path to check existence and modification time
          quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', cluster_num, '_100-', radius, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))
          
          # Check if the file is missing or older than 6 days
          if (!file.exists(quad_file) || (Sys.time() - file.info(quad_file)$mtime) > as.difftime(6, units = "days")) {
            
            # Print file path for logging
            print(file.path(data_dir, TMC, tissue, paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda')))
            
            system(paste0(
              'srun -p model3,model4,short1,pool1,gpu,interactive -o ',
              file.path(data_dir, TMC, tissue, paste0('concat_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.out')),
              ' -n 1 -c 4 --mem 80G Rscript ./modeling/heterogeneity/concat_multi_multi_all_images_across3_self.R ',
              data_dir, ' ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' &'
            ))
          } 
        }
      }
    }
  }
}
