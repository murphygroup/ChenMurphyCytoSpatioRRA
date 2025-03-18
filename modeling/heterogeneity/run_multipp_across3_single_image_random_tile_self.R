args = commandArgs(TRUE)
data_dir = args[1]
# Configuration
tissue_list <- c("LI", "SI", "LN", "THYMUS", "SPLEEN")
cluster_nums <- 5:5
intensity_types <- c("total")
radii <- c(100, 200, 300, 400, 500)
hr_values <- c(1)
img_idx_range <- 1:30
tile_sides <- c(2500, 5000)
tile_nums <- 1:5

for (tile_num in tile_nums) {
  for (radius in radii) {
    for (tissue in tissue_list) {
      for (cluster_num in cluster_nums) {
        for (intensity_type in intensity_types) {
          for (hr in hr_values) {
            
            # Determine TMC based on tissue type
            TMC <- ifelse(tissue %in% c("SI", "LI"), "Stanford", "Florida")
            
            for (img_idx in img_idx_range) {
              for (tile_side in tile_sides) {
                
                # Define file paths
                cell_list_file <- file.path(data_dir, TMC, tissue, "random_tile", 
                                            paste("cell_list_", cluster_num, "_", intensity_type, "_across3_", 
                                                  img_idx, "_tile_", tile_side, "_", tile_num, ".Rda", sep = ""))
                
                quad_file <- file.path(data_dir, TMC, tissue, "random_tile", 
                                       paste("quad_", cluster_num, "_", radius, "_", hr, "_", intensity_type, 
                                             "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, 
                                             "_self_quad_d_no_between_dummy.Rda", sep = ""))
                
                output_file <- file.path(data_dir, TMC, tissue, "random_tile", 
                                         paste("quad_", cluster_num, "_", radius, "_", hr, "_", intensity_type, 
                                               "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, 
                                               "_self_quad_d_no_between_dummy.out", sep = ""))
                
                # Check if cell list file exists
                if (file.exists(cell_list_file)) {
                  
                  # Check if quad file exists or is older than 60 days
                  if (!file.exists(quad_file) || (Sys.time() - file.info(quad_file)$mtime) > as.difftime(60, units = "days")) {
                    
                    print(quad_file)
                    
                    system(paste(
                      "srun -p short1,pool1,model1,model2,model3,pool3-bigmem,interactive -t 2:00:00 -o", 
                      output_file, 
                      "-n 1 -c 4 --mem 64G Rscript ./heterogeneity/multipp_across3_single_image_random_tile_self_no_between_dummy.R", 
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
}
