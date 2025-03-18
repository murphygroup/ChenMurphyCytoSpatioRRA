# Define the data directory
args = commandArgs(TRUE)
data_dir = args[1]

# Define parameters
radius_values <- c(500)
tissues <- c("SI", "LI", "LN", "SPLEEN", "THYMUS")
cluster_num <- 5
intensity_type <- "total"
hr <- 1
resample_percents <- c(0, 100, 200, 300, 400, 500, 600, 700)
random_seeds <- 1:1

# Iterate over parameter combinations
for (radius in radius_values) {
  for (tissue in tissues) {
    
    # Assign TMC based on tissue type
    TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")
    
    # Check if the cell list file exists
    cell_list_file <- file.path(data_dir, TMC, tissue, sprintf("cell_list_%d_%s_across3.Rda", cluster_num, intensity_type))
    
    if (file.exists(cell_list_file)) {
      
      for (resample_percent in resample_percents) {
        for (random_seed in random_seeds) {
          
          # Define the output file path
          output_file <- file.path(
            data_dir, TMC, tissue, "random_resample",
            sprintf(
              "random_simulated_pattern_%d_%d_%d_%s_across3_resample5_percentage_%d_seed_%d_self_quad_d_no_between_dummy.Rda",
              cluster_num, radius, hr, intensity_type, resample_percent, random_seed
            )
          )
          
          # Check if the file exists and its age
          if (!file.exists(output_file) || (Sys.time() - file.info(output_file)$mtime) > as.difftime(60, units = "days")) {
            
            # Print the file path for debugging
            print(output_file)
            
            # Submit the job
            system(paste(
              "srun -p model1,model2,pool1,model3,model4,pool3-bigmem -t 72:00:00 -o",
              file.path(
                data_dir, TMC, tissue, "random_resample",
                sprintf(
                  "random_simulated_pattern_%d_%d_%d_%s_across3_resample5_global_%d_%d_same_pattern.out",
                  cluster_num, radius, hr, intensity_type, resample_percent, random_seed
                )
              ),
              "-n 1 -c 4 --mem 20G Rscript ./simulation/simulate_random_resampling5_global_self_same_pattern.R",
              data_dir, ' ', tissue, intensity_type, cluster_num, radius, hr, resample_percent, random_seed, "&"
            ))
          }
        }
      }
    }
  }
}
