# Read command-line arguments
args = commandArgs(TRUE)
data_dir = args[1]

# Define parameter ranges
cluster_nums <- c(5)
radii <- c(100, 200, 300, 400, 500)
tissues <- c('SPLEEN', 'LN', 'THYMUS', 'LI', 'SI')
intensity_types <- c('total')
hr_values <- c(1)

# Iterate over all parameter combinations
for (cluster_num in cluster_nums) {
  for (radius in radii) {  
    for (tissue in tissues) {
      for (intensity_type in intensity_types) {
        
        # Determine Tissue Mapping Center (TMC)
        TMC <- if (tissue %in% c('SI', 'LI')) 'Stanford' else 'Florida'
        
        for (hr in hr_values) {
          # Define file path to check existence and modification time
          file_path <- file.path(data_dir, TMC, tissue, paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))
          
          # Check if file is missing or older than 20 days
          if (!file.exists(file_path) || (Sys.time() - file.info(file_path)$mtime) > as.difftime(20, units = "days")) {
            
            print(file.path(data_dir, TMC, tissue, paste0('combine_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3.Rda')))
            
            # Submit SLURM job
            system(paste0(
              'srun -p pool1,model1,model2,model3,model4,pool3-bigmem,gpu,pool1,interactive -o ',
              file.path(data_dir, TMC, tissue, paste0('combine_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_self.out')),
              ' -n 1 -c 4 --mem 46G Rscript ./modeling/heterogeneity/combine_quad_across3_self.R ',
              data_dir, ' ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' &'
            ))
          }
        }
      }
    }
  }
}
