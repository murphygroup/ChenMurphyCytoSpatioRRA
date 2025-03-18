# Load required libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)

# Function to extract a tile from the original pattern
get_tile_pattern <- function(original_pattern, rect) {
  original_pattern_df <- data.frame(
    x = original_pattern$x,
    y = original_pattern$y,
    marks = original_pattern$marks
  )
  
  # Extract tile based on given coordinates
  tile_df <- original_pattern_df[
    original_pattern_df$x >= rect[1] & original_pattern_df$x <= rect[3] & 
      original_pattern_df$y >= rect[2] & original_pattern_df$y <= rect[4], 
  ]
  
  # Create a point pattern object for the tile
  tile_pattern <- ppp(
    tile_df$x, tile_df$y, 
    window = owin(c(rect[1], rect[3]), c(rect[2], rect[4])),
    marks = factor(tile_df$marks)
  )
  
  return(tile_pattern)
}

# Function to generate tiles from microscopy images
tiling <- function(tissue, tile_side, total_tile_num) {
  args <- commandArgs(TRUE)
  data_dir <- args[1]  # Set data directory from arguments
  
  # Determine the corresponding TMC
  TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")
  
  # Create directory for storing tile data if it doesn't exist
  dir.create(file.path(data_dir, TMC, tissue, "random_tile"), showWarnings = FALSE)
  
  # Get list of cell list files
  cell_list_files <- Sys.glob(
    file.path(data_dir, TMC, tissue, paste0("cell_list_", cluster_num, "_", intensity_type, "_across3_*.Rda"))
  )
  
  for (img_idx in seq_along(cell_list_files)) {
    cell_list_file <- cell_list_files[img_idx]
    file_name <- basename(cell_list_file)
    load(cell_list_file)
    
    cell_pattern <- cell_list[[1]]
    
    # Image boundaries
    xmin_image <- 0
    ymin_image <- 0
    xmax_image <- cell_pattern$window$xrange[2]
    ymax_image <- cell_pattern$window$yrange[2]
    
    # Define valid tiling range within the image
    xmin_tile_range <- xmin_image + 500 + tile_side / 2
    ymin_tile_range <- ymin_image + 500 + tile_side / 2
    xmax_tile_range <- xmax_image - 500 - tile_side / 2
    ymax_tile_range <- ymax_image - 500 - tile_side / 2
    
    sampled_tile_num <- 0
    all_sampled_tile_num <- 0
    
    while (sampled_tile_num < total_tile_num & all_sampled_tile_num < 1000) {
      all_sampled_tile_num <- all_sampled_tile_num + 1
      
      # Randomly select tile center
      xcenter_tile <- sample(xmin_tile_range:xmax_tile_range, 1)
      ycenter_tile <- sample(ymin_tile_range:ymax_tile_range, 1)
      
      # Define tile boundaries
      tile_rect <- c(
        xcenter_tile - 500 - tile_side / 2,
        ycenter_tile - 500 - tile_side / 2,
        xcenter_tile + 500 + tile_side / 2,
        ycenter_tile + 500 + tile_side / 2
      )
      
      tile_real_rect <- c(
        tile_rect[1] + 500, tile_rect[2] + 500,
        tile_rect[3] - 500, tile_rect[4] - 500
      )
      
      # Extract tile patterns
      cell_pattern_tile <- get_tile_pattern(cell_pattern, tile_rect)
      cell_real_pattern_tile <- get_tile_pattern(cell_pattern, tile_real_rect)
      
      # Count unique cell types in the tile
      cell_type_num <- length(sort(unique(cell_real_pattern_tile$marks)))
      
      # Save valid tiles
      if (cell_pattern_tile$n > cell_pattern$n / total_tile_num / 5 & cell_type_num == cluster_num) {
        sampled_tile_num <- sampled_tile_num + 1
        filename <- file.path(
          data_dir, TMC, tissue, "random_tile",
          paste0(substr(file_name, 1, nchar(file_name) - 4), "_tile_", tile_side, "_", sampled_tile_num, ".Rda")
        )
        save(list(tile = cell_pattern_tile), file = filename)
      }
    }
  }
}

# Set parameters
cluster_num <- 5
intensity_type <- "total"
set.seed(1)

# List of tissues to process
tissue_list <- c("SPLEEN", "THYMUS", "LN", "SI", "LI")

# Perform tiling for each tissue and tile size
for (tissue in tissue_list) {
  for (tile_side in c(5000, 2500)) {
    tiling(tissue, tile_side, total_tile_num = 5)
  }
}
