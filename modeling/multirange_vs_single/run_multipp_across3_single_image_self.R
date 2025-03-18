# Read command-line arguments
args = commandArgs(TRUE)
data_dir = args[1]

# Define parameter ranges
radii <- c(100, 200, 300, 400, 500)
tissues <- c('SI', 'LN', 'LI', 'THYMUS', 'SPLEEN')
cluster_num <- 5
intensity_types <- c('total')
hr_values <- c(1)

# Iterate over all parameter combinations
for (radius in radii) {
  for (tissue in tissues) {
    for (intensity_type in intensity_types) {
      for (hr in hr_values) {
        
        # Determine Tissue Mapping Center (TMC)
        TMC <- if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'
        
        # Process 30 image indices
        for (img_idx in 1:30) {
          # Define file path for checking existence
          cell_list_file <- file.path(data_dir, TMC, tissue, paste0('cell_list_', cluster_num, '_', intensity_type, '_across3_', img_idx, '.Rda'))
          
          if (file.exists(cell_list_file)) {
            # Define output file paths
            quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_quad_d_no_between_dummy.Rda'))
            output_file <- file.path(data_dir, TMC, tissue, paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_quad_d_no_between_dummy.out'))
            
            # Submit SLURM job
            system(paste0(
              'srun -p pool1,model1,model2,model3,model4,pool3-bigmem -t 12:00:00 -o ',
              output_file,
              ' -n 1 -c 4 --mem 120G Rscript ./modeling/heterogeneity/multipp_across3_single_image_self_quad_d_no_between_dummy.R ',
              data_dir, ' ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', img_idx, ' &'
            ))
          }
        }
      }
    }
  }
}
