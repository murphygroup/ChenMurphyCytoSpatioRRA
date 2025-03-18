# Load required libraries
library(spatstat)

# Read command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]  # Updated to take data directory as an argument
tissue <- args[2]
intensity_type <- args[3]
num <- as.numeric(args[4])
r <- as.numeric(args[5])
hr <- as.numeric(args[6])

# Determine Tissue Mapping Center (TMC)
TMC <- if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Initialize container for combined quadrature data
Quad_all_all <- NULL

# Loop over image indices from 1 to 30
for (i in 1:30) {
  print(i)  # Print current image index
  
  # Construct file path for the quadrature data
  quad_image_file <- file.path(data_dir, TMC, tissue, paste0('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', i, '_self_quad_d_no_between_dummy.Rda'))
  
  # Check if the file exists before loading
  if (file.exists(quad_image_file)) {
    load(file = quad_image_file)  # Load quadrature data
    Quad_all$moadf$pattern_ID <- i  # Assign pattern ID
    
    # Remove "is_real" column if it exists
    if ("is_real" %in% colnames(Quad_all$moadf)) {
      Quad_all$moadf$is_real <- NULL
    }
    
    # If first iteration, initialize Quad_all_all
    if (is.null(Quad_all_all)) {
      Quad_all_all <- Quad_all
    } else {
      # Append new data to the existing quadrature object
      Quad_all_all$moadf <- rbind(Quad_all_all$moadf, Quad_all$moadf)
      Quad_all_all$Info$rownames <- c(Quad_all_all$Info$rownames, as.character(i))
      
      rownames(Quad_all$Inter$interaction) <- as.character(i)
      Quad_all_all$Inter$interaction <- rbind.hyperframe(Quad_all_all$Inter$interaction, Quad_all$Inter$interaction)
      
      Quad_all_all$npat <- Quad_all_all$npat + 1
      rownames(Quad_all$data) <- as.character(i)
      Quad_all_all$data <- rbind.hyperframe(Quad_all_all$data, Quad_all$data)
      
      levels(Quad_all_all$data$id) <- c(levels(Quad_all_all$data$id), as.character(i))
      Quad_all_all$data$id[i] <- i
      
      Quad_all_all$Y[[i]] <- Quad_all$Y[[1]]
      levels(Quad_all_all$datadf$id) <- c(levels(Quad_all_all$datadf$id), as.character(i))
      Quad_all_all$datadf <- rbind(Quad_all_all$datadf, i)
    }
    
    # Free memory after each iteration
    rm(Quad_all)
    gc()
  }
}

# Save the final combined quadrature data
quad_save_file <- file.path(data_dir, TMC, tissue, paste0('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_quad_d_no_between_dummy.Rda'))
save(Quad_all_all, file = quad_save_file)
