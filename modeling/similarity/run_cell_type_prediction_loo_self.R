args = commandArgs(TRUE)
data_dir = args[1]

# Configuration
tissue_list1 <- c("LI", "SI", "LN", "SPLEEN", "THYMUS")
tissue_list2 <- c("LI", "SI", "LN", "SPLEEN", "THYMUS")
cluster_nums <- 5:5
intensity_types <- c("total")
radii <- c(500)
hr <- 1
img_idx_range <- 1:30

for (tissue1 in tissue_list1) {
  for (tissue2 in tissue_list2) {
    for (cluster_num in cluster_nums) {
      for (intensity_type in intensity_types) {
        for (radius in radii) {
          
          # Determine TMC based on tissue type
          TMC <- ifelse(tissue1 %in% c("SI", "LI"), "Stanford", "Florida")
          
          for (img_idx in img_idx_range) {
            
            # File paths
            coef_file <- file.path(data_dir, TMC, tissue1, 
                                   paste("coef_", cluster_num, "_100-", radius, "_", hr, "_", 
                                         intensity_type, "_looidx_", img_idx, 
                                         "_across3_self_quad_d_no_between_dummy.Rda", sep = ""))
            
            devi_file <- file.path(data_dir, TMC, tissue1, 
                                   paste("devi_loo_", cluster_num, "_100-", radius, "_", hr, "_", 
                                         intensity_type, "_", tissue1, "_", img_idx, "_", tissue2, 
                                         "_across3_self_quad_d_no_between_dummy.Rda", sep = ""))
            
            predict_log_file <- file.path(data_dir, TMC, tissue1, 
                                          paste("predict_", cluster_num, "_100-", radius, "_", hr, "_", 
                                                intensity_type, "_across3_", tissue1, "_", img_idx, "_", 
                                                tissue2, "_self.out", sep = ""))
            
            # Check if coef file exists before running the prediction
            if (file.exists(coef_file)) {
              print(devi_file)
              
              system(paste(
                "srun -p pool1 -t 72:00:00 -o", predict_log_file, 
                "-n 1 -c 4 --mem 20G Rscript ./similarity/cell_type_prediction_loo_self.R", 
                data_dir, ' ', tissue1, tissue2, intensity_type, cluster_num, radius, hr, img_idx, "&", sep = " "
              ))
            }
          }
        }
      }
    }
  }
}
