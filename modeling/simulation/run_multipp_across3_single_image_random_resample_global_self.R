# Define data directory
args = commandArgs(TRUE)
data_dir = args[1]

# Initialize index counter
idx = 0

# Loop through random seeds
for (random_seed in 1:5) {
  # Loop through radius values
  for (radius in c(100, 200, 300, 400, 500)) {
    # Loop through tissue types
    for (tissue in c("SPLEEN", "LN", "LI", "THYMUS", "SI")) {
      for (cluster_num in 5:5) {
        for (intensity_type in c("total")) {
          for (hr in c(1)) {
            
            # Assign TMC based on tissue type
            TMC = if (tissue == "LI" | tissue == "SI") "Stanford" else "Florida"
            
            # Loop through resample percentages
            for (resample_percent in c(0, 100, 200, 300, 400, 500, 600, 700)) {
              
              # Define input file path
              file_path = file.path(
                data_dir, TMC, tissue, "random_resample",
                paste(
                  "random_simulated_pattern_", cluster_num, "_500_", hr, "_",
                  intensity_type, "_across3_resample5_percentage_", resample_percent, 
                  "_seed_", random_seed, "_self_quad_d_no_between_dummy.Rda",
                  sep = ""
                )
              )
              
              # Check if the file exists
              if (file.exists(file_path)) {
                
                # Define output file path
                file_path = file.path(
                  data_dir, TMC, tissue, "random_resample",
                  paste(
                    "quad_", cluster_num, "_", radius, "_", hr, "_", intensity_type, 
                    "_across3_resample5_percentage_", resample_percent, 
                    "_seed_", random_seed, "_self_quad_d_no_between_dummy.Rda",
                    sep = ""
                  )
                )
                
                # Check if output file does not exist or is older than 45 days
                if (!file.exists(file_path) || (Sys.time() - file.info(file_path)$mtime) > as.difftime(45, units = "days")) {
                  
                  # Submit job to SLURM
                  system(
                    paste(
                      "srun -p model1,model2,pool1,model3,model4,pool3-bigmem -t 12:00:00 -o",
                      file.path(
                        data_dir, TMC, tissue, "random_resample",
                        paste(
                          "quad_", cluster_num, "_", radius, "_", hr, "_", intensity_type, 
                          "_across3_resample5_global_", resample_percent, "_", random_seed, ".out",
                          sep = ""
                        )
                      ),
                      "-n 1 -c 4 --mem 64G Rscript ./modeling/simulation/multipp_across3_single_image_random_resample_global_self_quad_d_no_between_dummy.R",
                      data_dir, tissue, intensity_type, cluster_num, radius, hr, resample_percent, random_seed, "&",
                      sep = " "
                    )
                  )
                }
              }
            }
          }
        }
      }
    }
  }
}
