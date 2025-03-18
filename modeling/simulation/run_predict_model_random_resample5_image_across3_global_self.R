# Set data directory
args = commandArgs(TRUE)
data_dir = args[1]

# Configuration
random_seeds <- 1:5
radii <- c(500)
tissue_list <- c("SPLEEN", "LN", "LI", "THYMUS", "SI")
cluster_nums <- 5:5
intensity_types <- c("total")
hr_values <- c(1)
resample_percents <- c(100, 200, 300, 400, 500, 600, 700)

idx <- 0

for (random_seed in random_seeds) {
  for (radius in radii) {
    for (tissue in tissue_list) {
      for (cluster_num in cluster_nums) {
        for (intensity_type in intensity_types) {
          for (hr in hr_values) {
            
            # Determine TMC based on tissue type
            TMC <- ifelse(tissue %in% c("SI", "LI"), "Stanford", "Florida")
            
            for (resample_percent in resample_percents) {
              
              # Define file paths
              simulated_pattern_path <- file.path(data_dir, TMC, tissue, "random_resample",
                                                  paste0("random_simulated_pattern_", cluster_num, "_500_", hr, "_", intensity_type, 
                                                         "_across3_resample5_percentage_", resample_percent, "_seed_", random_seed, 
                                                         "_self_quad_d_no_between_dummy.Rda"))
              
              pred_cell_type_path <- file.path(data_dir, TMC, tissue, "random_resample",
                                               paste0("pred_cell_type_", cluster_num, "_500_", hr, "_", intensity_type, 
                                                      "_across3_resample5_percentage_", resample_percent, "_seed_", random_seed, 
                                                      "_self_quad_d_no_between_dummy.Rda"))
              
              # Execute only if simulated pattern exists but pred_cell_type does not
              if (file.exists(simulated_pattern_path) && !file.exists(pred_cell_type_path)) {
                
                print(paste("Processing:", pred_cell_type_path))
                
                # Run the first job
                quad_out_path <- file.path(data_dir, TMC, tissue, "random_resample",
                                           paste0("quad_", cluster_num, "_", radius, "_", hr, "_", intensity_type, 
                                                  "_across3_resample5_global_", resample_percent, "_", random_seed, ".out"))
                
                system(paste("srun -p model1,model2,pool1,model3,model4,pool3-bigmem -t 12:00:00 -o", 
                             quad_out_path, 
                             "-n 1 -c 4 --mem 46G Rscript ./simulation/multipp_across3_single_image_random_resample_global_self_quad_d_no_between_dummy.R", 
                             tissue, intensity_type, cluster_num, radius, hr, resample_percent, random_seed, "&"))
                
                # Run the second job
                pred_out_path <- file.path(data_dir, TMC, tissue, "random_resample",
                                           paste0("pred_", cluster_num, "_", radius, "_", hr, "_", intensity_type, 
                                                  "_across3_resample5_global_", resample_percent, "_", random_seed, ".out"))
                
                system(paste("srun -p model1,model2,pool1,model3,model4,pool3-bigmem -t 12:00:00 -o", 
                             pred_out_path, 
                             "-n 1 -c 4 --mem 46G Rscript ./simulation/predict_model_random_resample5_image_across3_global_self.R", 
                             data_dir, ' ', tissue, intensity_type, cluster_num, radius, hr, resample_percent, random_seed, "&"))
                
                idx <- idx + 1
              }
            }
          }
        }
      }
    }
  }
}

print(idx)
