data_dir <- "/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new"

# Configuration
tissue_list <- c("LI", "SI", "THYMUS", "LN", "SPLEEN")
cluster_nums <- 5:5
intensity_types <- c("total")
radii <- c(500)
hr_values <- c(1)
img_idx_range <- 1:30
resample_percents <- c(0, 100, 200, 300, 400)

for (radius in radii) {
  for (tissue in tissue_list) {
    for (cluster_num in cluster_nums) {
      for (intensity_type in intensity_types) {
        for (hr in hr_values) {
          
          # Determine TMC based on tissue type
          TMC <- ifelse(tissue %in% c("SI", "LI"), "Stanford", "Florida")
          
          for (img_idx in img_idx_range) {
            cell_list_path <- file.path(data_dir, TMC, tissue, paste("cell_list_", cluster_num, "_", intensity_type, "_across3_", img_idx, ".Rda", sep = ""))
            
            # Check if cell list file exists before proceeding
            if (file.exists(cell_list_path)) {
              
              # Create 'random_resample' directory if it doesn't exist
              dir.create(file.path(data_dir, TMC, tissue, "random_resample"), showWarnings = FALSE)
              
              for (resample_percent in resample_percents) {
                output_path <- file.path(data_dir, TMC, tissue, "random_resample",
                                         paste("random_simulated_pattern_", cluster_num, "_", radius, "_", hr, "_", intensity_type,
                                               "_across3_", img_idx, "_resample5_", resample_percent, "_self_quad_d_no_between_dummy.Rda", sep = ""))
                
                # Always print output path (kept from original script)
                print(output_path)
                
                system(paste(
                  "srun -p short1,model1,model2,pool1,pool3-bigmem,interactive,gpu -t 2:00:00 -o",
                  file.path(data_dir, TMC, tissue, "random_resample",
                            paste("random_simulated_pattern_", cluster_num, "_", radius, "_", hr, "_", intensity_type,
                                  "_across3_", img_idx, "_resample5_", resample_percent, "_self_same_pattern.out", sep = "")),
                  "-n 1 -c 4 --mem 22G Rscript simulate_random_resampling5_self_same_pattern.R",
                  tissue, intensity_type, cluster_num, radius, hr, img_idx, resample_percent, "10000 &",
                  sep = " "
                ))
              }
            }
          } 
        }
      }
    }
  }
}
