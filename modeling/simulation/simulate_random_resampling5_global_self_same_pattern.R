# Load necessary libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(spatstat.geom)
library(dplyr)
library(proxy)

# Define script directory
script_dir <- "./modeling/ppp"

# Function to randomly sample marks based on given probabilities
get_random_sampled_marks <- function(types, num, prob) {
  return(sample(types, num, TRUE, prob))
}

# Function to calculate distance between new points and existing points
calc_dist <- function(new_points, existing_points) {
  new_points <- data.frame(t(new_points))
  x <- new_points[1, 1]
  y <- new_points[1, 2]
  
  xmin <- max(0, x - radius)
  ymin <- max(0, y - radius)
  xmax <- min(cell_pattern_random_window$xrange[2], x + radius)
  ymax <- min(cell_pattern_random_window$yrange[2], y + radius)
  
  existing_points_within_range <- existing_points[
    which(existing_points[,1] > xmin & existing_points[,1] < xmax & 
            existing_points[,2] > ymin & existing_points[,2] < ymax), ]
  
  eu_dist <- dist(new_points, existing_points_within_range[,1:2])
  existing_points_within_range <- cbind(existing_points_within_range, eu_dist)
  existing_points_within_range <- existing_points_within_range[
    which(existing_points_within_range[,4] < 500 & existing_points_within_range[,4] > 0), 3:4]
  
  mark_count <- matrix(0, nrow = 5, ncol = 5)
  
  for (i in 1:nrow(existing_points_within_range)) {
    mark_count_row <- as.integer(existing_points_within_range[i, 1])
    mark_count_col <- existing_points_within_range[i, 2] %/% 100 + 1
    mark_count[mark_count_row, mark_count_col] <- mark_count[mark_count_row, mark_count_col] + 1
  }
  
  return(mark_count)
}

# Function to efficiently calculate distances using data.table
calc_dist_optimized <- function(new_points, existing_points) {
  existing_points$marks <- as.integer(existing_points$marks)
  existing_points <- as.data.table(existing_points)
  names(existing_points) <- c("x", "y", "marks")
  
  new_points <- data.frame(t(new_points))
  x <- new_points[1, 1]
  y <- new_points[1, 2]
  
  xmin <- max(0, x - radius)
  ymin <- max(0, y - radius)
  xmax <- min(cell_pattern_window$xrange[2], x + radius)
  ymax <- min(cell_pattern_window$yrange[2], y + radius)
  
  existing_points_within_range <- existing_points[x > xmin & x < xmax & y > ymin & y < ymax]
  eu_dist <- t(rdist(new_points, existing_points_within_range[,1:2]))
  
  existing_points_within_range[, eu_dist := eu_dist]
  existing_points_within_range <- existing_points_within_range[eu_dist < 500 & eu_dist > 0, .(marks, eu_dist)]
  
  mark_count <- table(as.integer(existing_points_within_range[, marks]), existing_points_within_range[, eu_dist] %/% 100 + 1)
  mark_count_matrix <- matrix(0, nrow = 5, ncol = 5)
  
  mark_count_matrix[as.matrix(expand.grid(as.numeric(rownames(mark_count)), as.numeric(colnames(mark_count))))] <- as.vector(mark_count)
  
  return(mark_count_matrix)
}

# Read command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]
tissue <- args[2]
intensity_type <- args[3]
cluster_num <- as.numeric(args[4])
radius <- as.numeric(args[5])
hr <- as.numeric(args[6])
resample_percent <- as.numeric(args[7])
cell_num <- as.numeric(args[8])

# Assign TMC based on tissue type
TMC <- ifelse(tissue %in% c("SI", "LI"), "Stanford", "Florida")

# Define data directory

# Set random seed
set.seed(3)

# Load coefficients
coef_file <- file.path(data_dir, TMC, tissue, sprintf("coef_%d_100-%d_%d_%s_across3_self_quad_d_no_between_dummy.Rda", cluster_num, radius, hr, intensity_type))
load(coef_file)
print(coef)

# Construct interaction matrix G
G <- list()
coef_names <- names(coef)

for (i in 0:(cluster_num - 1)) {
  coef_current <- c()
  for (r in seq(100, radius, 100)) {
    for (j in 0:(cluster_num - 1)) {
      mark1 <- as.character(i)
      mark2 <- as.character(j)
      
      interaction_name <- if (i <= j) {
        sprintf("X%sxX%sx%d", mark1, mark2, r)
      } else {
        sprintf("X%sxX%sx%d", mark2, mark1, r)
      }
      
      coef_current <- c(coef_current, coef[grep(interaction_name, coef_names, value = FALSE)])
    }
  }
  G[[i + 1]] <- coef_current
}

# Construct B matrix
B <- exp(coef)
B[-1] <- exp(coef[-1] + coef[1])

# Define possible types
types <- as.character(0:4)

# Function to calculate mark probabilities
calc_mark_prob <- function(counts) {
  m_prob <- sapply(seq_along(types), function(m_int) {
    G_type <- G[[m_int]]
    interact_term <- exp(sum(G_type * counts))
    intensity_term <- B[m_int]
    intensity_term * interact_term
  })
  
  m_prob <- m_prob / sum(m_prob)
  print(m_prob)
  
  return(factor(sample(types, 1, prob = m_prob), levels = types))
}

# Load cell list for multiple images
cell_pattern_all_images <- NULL
size_all_images <- 0

for (i in 1:20) {
  file_dir <- file.path(data_dir, TMC, tissue, sprintf("cell_list_%d_%s_across3_%d.Rda", cluster_num, intensity_type, i))
  if (file.exists(file_dir)) {
    load(file_dir)
    if (is.null(cell_pattern_all_images)) {
      cell_pattern_all_images <- data.frame(cell_list[[1]])
      size_all_images <- cell_list[[1]]$window$xrange[2] * cell_list[[1]]$window$yrange[2]
    } else {
      cell_pattern_all_images <- rbind(cell_pattern_all_images, data.frame(cell_list[[1]]))
      size_all_images <- size_all_images + cell_list[[1]]$window$xrange[2] * cell_list[[1]]$window$yrange[2]
    }
  }
}

# Compute model mark probabilities
model_mark_prob <- dplyr::count(cell_pattern_all_images, marks)$n / nrow(cell_pattern_all_images)

# Generate random cell pattern
img_size <- 2000
intensity <- cell_num / (img_size^2)
file_path <- file.path(data_dir, TMC, tissue, "random_resample", sprintf("random_same_pattern_cell_num_%d.Rda", cell_num))

cell_pattern_random <- rpoispp(intensity, win = c(c(0, img_size), c(0, img_size)))
save(cell_pattern_random, file = file_path)

# Convert pattern to dataframe
cell_pattern_random_df <- data.frame(cell_pattern_random)
cell_pattern_random_window <<- cell_pattern_random$window

# Apply resampling
pattern_exist <- cell_pattern_random_df
pattern_exist[, 3] <- as.factor(sample(types, nrow(cell_pattern_random_df), replace = TRUE, model_mark_prob))

# Perform resampling based on resample_percent
resample_num <- ceiling(nrow(cell_pattern_random_df) * resample_percent / 100)
print(resample_num)

if (resample_num > 0) {
  for (i in 1:resample_num) {
    point_idx <- sample(nrow(cell_pattern_random_df), 1)
    counts_current <- t(apply(pattern_exist[point_idx, 1:2, drop = FALSE], 1, calc_dist, pattern_exist))
    pattern_exist[point_idx, 3] <- apply(counts_current, 1, calc_mark_prob)
  }
}

# Save results
save(pattern_exist, file = file_path)
write.csv(pattern_exist, gsub(".Rda", ".csv", file_path), row.names = FALSE)
