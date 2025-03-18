data_dir <- "/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new"

# Configuration
tissue_list <- c("LN", "SPLEEN", "THYMUS", "LI", "SI")
cluster_nums <- 5:5
intensity_types <- c("total")
radii <- c(500)
hr_values <- c(1)
img_idx_range <- 1:1
resample_percents <- c(0, 50, 100, 150, 200, 250, 300, 350, 400)

for (radius in radii) {
  for (tissue in tissue_list) {
    for (cluster_num in cluster_nums) {
      for (intensity_type in intensity_types) {
        for (hr in hr_values) {
          
          # Determine TMC based on tissue type
          TMC <- ifelse(tissue %in% c("SI", "LI"), "Stanford", "Florida")
          
          # Create 'resample' directory if it doesn't exist
          dir.create(file.path(data_dir, TMC, tissue, "resample"), showWarnings = FALSE)
          
          for (img_idx in img_idx_range) {
            cell_list_path <- file.path(data_dir, TMC, tissue, paste("cell_list_", cluster_num, "_", intensity_type, "_across3_", img_idx, ".Rda", sep = ""))
            
            # Check if cell list file exists before proceeding
            if (file.exists(cell_list_path)) {
              
              for (resample_percent in resample_percents) {
                output_path <- file.path(data_dir, TMC, tissue, "resample",
                                         paste("marked_simulated_pattern_", cluster_num, "_", radius, "_", hr, "_", intensity_type,
                                               "_across3_", img_idx, "_resample5_", resample_percent, "_self_quad_d_no_between_dummy.Rda", sep = ""))
                
                # Run only if output file doesn't exist
                if (!file.exists(output_path)) {
                  print(output_path)
                  
                  system(paste(
                    "srun -p gpu,model1,model2,model3,pool1,pool3-bigmem,model4 -t 72:00:00 -o",
                    file.path(data_dir, TMC, tissue, "resample",
                              paste("marked_simulated_pattern_", cluster_num, "_", radius, "_", hr, "_", intensity_type,
                                    "_across3_", img_idx, "_resample5_", resample_percent, "_self.out", sep = "")),
                    "-n 1 -c 4 --mem 20G Rscript simulate_marking_resampling5_self.R",
                    tissue, intensity_type, cluster_num, radius, hr, img_idx, resample_percent, "&",
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
}
