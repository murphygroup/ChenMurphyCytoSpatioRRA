# Define data directory
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments

# Parse command-line arguments
args <- commandArgs(TRUE)
tissue <- args[2]
intensity_type <- args[3]
n <- as.numeric(args[4])
r <- as.numeric(args[5])
hr <- as.numeric(args[6])
img_idx <- as.numeric(args[7])
tile_side <- args[8]
tile_num <- args[9]

# Determine TMC location
TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")

# Initialize concatenated_quad variable
for (radii in seq(100, r, 100)) {
  
  # Load data for current radius
  quad_file <- file.path(
    data_dir, TMC, tissue, "random_tile",
    paste0("quad_", n, "_", radii, "_", hr, "_", intensity_type, 
           "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, "_self_quad_d_no_between_dummy.Rda")
  )
  
  load(quad_file)
  current_quad <- Quad_all$moadf
  
  # Remove "is_real" column if present
  current_quad <- current_quad[, !colnames(current_quad) %in% "is_real"]
  
  print(current_quad[1:5, ])  # Print first 5 rows for verification
  
  # Free memory
  rm(Quad_all)
  gc()
  
  # Rename columns with radius suffix
  colnames(current_quad)[7:(ncol(current_quad) - 2)] <- 
    paste0(colnames(current_quad)[7:(ncol(current_quad) - 2)], "x", radii)
  
  # Initialize or update concatenated_quad
  if (radii == 100) {
    last_current_quad <- current_quad[, 7:(ncol(current_quad) - 2)]
    concatenated_quad <- current_quad[, 1:(ncol(current_quad) - 2)]
    pattern_ID <- current_quad[, (ncol(current_quad) - 1)]
    caseweight <- current_quad[, ncol(current_quad)]
  } else {
    current_quad <- current_quad[, 7:(ncol(current_quad) - 2)]
    concatenated_quad <- cbind(concatenated_quad, current_quad - last_current_quad)
    last_current_quad <- current_quad
  }
}

# Free memory
rm(current_quad)
gc()

# Append pattern_ID and caseweight to concatenated_quad
concatenated_quad <- cbind(concatenated_quad, pattern_ID, caseweight)

# Load the final Quad_all file
final_quad_file <- file.path(
  data_dir, TMC, tissue, "random_tile",
  paste0("quad_", n, "_", r, "_", hr, "_", intensity_type, 
         "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, "_self_quad_d_no_between_dummy.Rda")
)

load(final_quad_file)
Quad_all_all <- Quad_all
Quad_all_all$moadf <- concatenated_quad

# Print first 5 rows for verification
print(concatenated_quad[1:5, ])

# Save updated Quad_all_all object
output_file <- file.path(
  data_dir, TMC, tissue, "random_tile",
  paste0("quad_", n, "_100-", r, "_", hr, "_", intensity_type, 
         "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, "_self_quad_d_no_between_dummy.Rda")
)
save(Quad_all_all, file = output_file)
