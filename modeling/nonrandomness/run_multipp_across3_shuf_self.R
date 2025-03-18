args = commandArgs(TRUE)
data_dir = args[1]

for (radius in c(100)) {
  for (tissue in c('LN', 'SI', 'THYMUS', 'SPLEEN', 'LI')) {
    for (cluster_num in 5:5) {
      for (intensity_type in c('total')) {
        
        # Set shuffling range based on tissue type
        shuf_range = if (tissue %in% c('LI', 'SI')) 1:2500 else 1:150
        
        for (shuf_num in shuf_range) {
          hr = 1
          TMC = if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'
          
          # Create directory if not exists
          dir.create(file.path(data_dir, TMC, tissue, 'random_shuf'), showWarnings = FALSE)
          
          # Define file path
          file_path = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda'))
          
          # Check if file is missing or older than 30 days, then run system command
          if (!file.exists(file_path) || (Sys.time() - file.info(file_path)$mtime) > as.difftime(30, units = "days")) {
            print(file_path)
            
            system(paste0(
              'srun -p model3,model1,model2,pool3-bigmem,pool1,gpu -t 12:00:00 -o ',
              file.path(data_dir, TMC, tissue, 'random_shuf', paste0('quad_', cluster_num, '_', radius, '_', intensity_type, '_across3_shuf_', shuf_num, '_self.out')),
              ' -n 1 -c 4 --mem 64G Rscript ./modeling/nonrandomness/multipp_across3_shuf_self.R ',
              data_dir, ' ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', shuf_num, ' &'
            ))
          }
        }
      }
    }
  }
}
