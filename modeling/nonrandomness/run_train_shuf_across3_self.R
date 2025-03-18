args = commandArgs(TRUE)
data_dir = args[1]

for (tissue in c('SI', 'LI', 'THYMUS', 'LN', 'SPLEEN')) {
  for (cluster_num in 5:5) {
    for (intensity_type in c('total')) {
      for (radius in c(100)) {
        
        # Set Tissue Mapping Center (TMC) based on tissue type
        TMC = if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'
        
        # Set shuffle range dynamically based on tissue type
        shuf_range = if (tissue %in% c('LI', 'SI')) 1:2500 else 1:150
        
        for (shuffle_num in shuf_range) {
          hr = 1
          
          # Define file paths
          quad_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('quad_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuffle_num, '_self_quad_d_no_between_dummy.Rda'))
          devi_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('devi_shuf_on_ori_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuffle_num, '_self_quad_d_no_between_dummy.Rda'))
          output_file = file.path(data_dir, TMC, tissue, 'random_shuf', paste0('model_shuf_', cluster_num, '_', radius, '_', hr, '_', intensity_type, '_across3_shuf_', shuffle_num, '_self_quad_d_no_between_dummy.out'))
          
          # Check if the required quad file exists
          if (file.exists(quad_file)) {
            
            # Check if devi_file does not exist or is older than 7 days
            if (!file.exists(devi_file) || (Sys.time() - file.info(devi_file)$mtime) > as.difftime(7, units = "days")) {
              print(quad_file)
              
              # Submit job using system call
              system(paste0(
                'srun -p model1,model2,model3,model4,pool1,pool3-bigmem,short1,gpu,interactive -o ', 
                output_file,
                ' -n 1 -c 8 --mem 46G Rscript ./modeling/nonrandomness/train_shuf_across3_self_quad_d_no_between_dummy.R ',
                data_dir, ' ', tissue, ' ', intensity_type, ' ', cluster_num, ' ', radius, ' ', hr, ' ', shuffle_num, ' &'
              ))
            }
          }
        }
      }
    }
  }
}
