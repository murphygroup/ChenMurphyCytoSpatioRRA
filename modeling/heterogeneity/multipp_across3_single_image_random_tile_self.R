# Load necessary libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(ggplot2)
library(dplyr)
library(permute)
library(data.table)

# Set script directory
script_dir <- "./modeling/ppp/"
setwd(script_dir)

# Load required source files
source('mppm.ppp.R')
source('mppm.quad.ppp.R')
source('bt.frame.ppp.R')
source('quadscheme.ppp.R')
source('default.dummy.ppp.R')
source('default.n.tiling.ppp.R')
source('mpl.engine.ppp.R')
source('mpl.prepare.ppp.R')
source('evalInteraction.ppp.R')
source('quadBlockSizes.ppp.R')
source('evalInterEngine.ppp.R')
source('evaluate.ppp.R')
source('evalPairPotential.ppp.R')

# Function to calculate the average number of particles across tomograms
get_avg_particle_num <- function(tomograms) {
  particle_num_total <- sum(sapply(tomograms, function(t) t$n))
  return(floor(particle_num_total / length(tomograms)))
}

print("Generating quadrature scheme...")

# Parse command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]  # Dynamically passed instead of hardcoded
tissue <- args[2]
intensity_type <- args[3]
n <- as.numeric(args[4])
r <- as.numeric(args[5])
hr <- as.numeric(args[6])
img_idx <- args[7]
tile_side <- args[8]
tile_num <- args[9]

# Determine TMC location
TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")

print(paste("Processing:", tissue, intensity_type, n, r, sep = " "))

# Load cell list data
cell_list_file <- file.path(data_dir, TMC, tissue, "random_tile", 
                            paste0("cell_list_", n, "_", intensity_type, "_across3_", 
                                   img_idx, "_tile_", tile_side, "_", tile_num, ".Rda"))
load(cell_list_file)

P <- cell_list
print(paste("Number of images:", length(P)))

# Compute average particle count
avg_particle_num <- get_avg_particle_num(P)

# Prepare interaction list
interaction <- list()
for (i in seq_along(P)) {
  num <- length(unique(P[[i]]$marks))
  ir_matrix <- matrix(rep(r, num * num), nrow = num)
  hr_matrix <- matrix(rep(hr, num * num), nrow = num)
  interaction[[i]] <- MultiStraussHard(iradii = ir_matrix, hradii = hr_matrix)
}

# Construct hyperframe
H <- hyperframe(Y = P)
I <- hyperframe(Interaction = interaction)

# Generate quadrature scheme
Quad_all <- mppm.quad.ppp(Y ~ marks, data = H, interaction = I, 
                          average_number = avg_particle_num, run_quad = FALSE)

# Save output
output_file <- file.path(data_dir, TMC, tissue, "random_tile", 
                         paste0("quad_", n, "_", r, "_", hr, "_", intensity_type, 
                                "_across3_", img_idx, "_tile_", tile_side, "_", 
                                tile_num, "_self_quad_d_no_between_dummy.Rda"))
save(Quad_all, file = output_file)

# Cleanup memory
gc()

# Print warnings (if any)
print(warnings())
