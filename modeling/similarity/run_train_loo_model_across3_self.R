args = commandArgs(TRUE)
data_dir = args[1]

# Configuration
tissues <- c("SPLEEN")
cluster_nums <- 5:5
intensity_types <- c("total")
radii <- c(500)
hr <- 1

for (tissue in tissues) {
  for (cluster_num in cluster_nums) {
    for (intensity_type in intensity_types) {
      for (radius in radii) {
        
        # Determine TMC based on tissue type
        TMC <- ifelse(tissue %in% c("SI", "LI"), "Stanford", "Florida")
        
        # File paths
        quad_file <- file.path(data_dir, TMC, tissue, 
                               paste("quad_", cluster_num, "_100-", radius, "_", hr, "_", 
                                     intensity_type, "_across3_self_quad_d_no_between_dummy.Rda", sep = ""))
        
        coef_file <- file.path(data_dir, TMC, tissue, 
                               paste("coef_", cluster_num, "_100-", radius, "_", hr, "_", 
                                     intensity_type, "_across3_self_quad_d_no_between_dummy.Rda", sep = ""))
        
        log_file <- file.path(data_dir, TMC, tissue, 
                              paste("model_loo_", cluster_num, "_100-", radius, "_", hr, "_", 
                                    intensity_type, "_across3_self.out", sep = ""))
        
        # Check if the quad file exists
        if (file.exists(quad_file)) {
          print(coef_file)
          
          # Run the model training script with appropriate resources
          system(paste(
            "srun -p gpu,model3,model4,pool3-bigmem -t 72:00:00 -o", log_file, 
            "-n 1 -c 8 --mem 185G Rscript ./similarity/train_loo_model_across3_self_quad_d_no_between_dummy.R", 
            data_dir, ' ', tissue, intensity_type, cluster_num, radius, hr, "&", sep = " "
          ))
        }
      }
    }
  }
}
