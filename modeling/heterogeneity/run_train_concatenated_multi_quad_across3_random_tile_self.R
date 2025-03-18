args = commandArgs(TRUE)
data_dir = args[1]
# Configuration
tissue_list <- c("SI", "LI", "THYMUS", "SPLEEN", "LN")
cluster_nums <- 5:5
intensity_types <- c("total")
radii <- c(500)
hr_values <- c(1)
img_idx_range <- 1:30
tile_sides <- c(2500, 5000)
tile_nums <- 1:5

for (radius in radii) {
  for (tissue in tissue_list) {
    for (cluster_num in cluster_nums) {
      for (intensity_type in intensity_types) {
        for (hr in hr_values) {
          
          # Determine TMC based on tissue type
          TMC <- ifelse(tissue %in% c("SI", "LI"), "Stanford", "Florida")
          
          for (img_idx in img_idx_range) {
            for (tile_side in tile_sides) {
              for (tile_num in tile_nums) {
                
                # Define file paths
                quad_file <- file.path(data_dir, TMC, tissue, "random_tile",
                                       paste("quad_", cluster_num, "_100-", radius, "_", hr, "_", intensity_type, 
                                             "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, 
                                             "_self_quad_d_no_between_dummy.Rda", sep = ""))
                
                devi_file <- file.path(data_dir, TMC, tissue, "random_tile",
                                       paste("devi_", cluster_num, "_100-", radius, "_", hr, "_", intensity_type, 
                                             "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, 
                                             "_self_quad_d_no_between_dummy.Rda", sep = ""))
                
                output_file <- file.path(data_dir, TMC, tissue, "random_tile",
                                         paste("train_", cluster_num, "_", radius, "_", hr, "_", intensity_type, 
                                               "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, 
                                               "_self.out", sep = ""))
                
                # Check if the quad file exists and devi file does not exist
                if (file.exists(quad_file) && !file.exists(devi_file)) {
                  
                  print(quad_file)
                  
                  system(paste(
                    "srun -p gpu,short1,model1,model2,pool1,pool3-bigmem -o", 
                    output_file, 
                    "-n 1 -c 8 --mem 46G Rscript ./heterogeneity/train_concatenated_multi_quad_across3_random_tile_self.R", 
                    data_dir, ' ', tissue, intensity_type, cluster_num, radius, hr, img_idx, tile_side, tile_num, "&", 
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
