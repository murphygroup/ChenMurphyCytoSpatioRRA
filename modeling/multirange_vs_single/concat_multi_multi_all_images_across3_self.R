# Set the data directory from command-line argument
args <- commandArgs(TRUE)
data_dir <- args[1]  # Updated to take data directory as an argument
tissue <- args[2]
intensity_type <- args[3]
n <- as.numeric(args[4])
r <- as.numeric(args[5])
hr <- as.numeric(args[6])

# Determine Tissue Mapping Center (TMC)
TMC <- if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Initialize variables for concatenated quadrature data
concatenated_quad <- NULL
last_current_quad <- NULL
pattern_ID <- NULL
caseweight <- NULL

# Iterate over radius values from 100 to r, in steps of 100
for (radii in seq(100, r, 100)) {
  
  # Load quadrature data for the current radius
  quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', n, '_', radii, '_1_total_across3_self_quad_d_no_between_dummy.Rda'))
  load(quad_file)
  
  # Extract relevant data from loaded quadrature object
  current_quad <- Quad_all_all$moadf
  print(colnames(current_quad))  # Print column names for debugging
  
  # Free memory by removing large object
  rm(Quad_all_all)
  gc()
  
  # Append radius identifier to column names (excluding last two columns)
  for (c in 7:(ncol(current_quad) - 2)) {
    colnames(current_quad)[c] <- paste0(colnames(current_quad)[c], "x", radii)
  }
  
  # Initialize concatenated dataset at first iteration
  if (radii == 100) {
    last_current_quad <- current_quad[, 7:(ncol(current_quad) - 2)]
    concatenated_quad <- current_quad[, 1:(ncol(current_quad) - 2)]
    pattern_ID <- current_quad[, ncol(current_quad) - 1]
    caseweight <- current_quad[, ncol(current_quad)]
  } else {
    # Compute difference from the previous radius dataset and concatenate
    current_quad <- current_quad[, 7:(ncol(current_quad) - 2)]
    concatenated_quad <- cbind(concatenated_quad, current_quad - last_current_quad)
    last_current_quad <- current_quad
  }
}

# Free memory after processing
rm(current_quad)
gc()

# Append pattern_ID and caseweight columns back to the concatenated dataset
concatenated_quad <- cbind(concatenated_quad, pattern_ID, caseweight)

# Load the final quadrature object and replace its moadf data
final_quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', n, '_', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))
load(final_quad_file)
Quad_all_all$moadf <- concatenated_quad

# Save the updated quadrature object
output_file <- file.path(data_dir, TMC, tissue, paste0('quad_', n, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))
save(Quad_all_all, file = output_file)

# Print column names of the final dataset for verification
print(colnames(concatenated_quad))
