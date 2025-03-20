# Load necessary libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(ggplot2)
library(dplyr)
library(permute)
library(data.table)

# Set working directory
script_dir <- "./modeling/ppp"
setwd(script_dir)

# Load required functions
source_files <- c(
  "mppm.ppp.R", "mppm.quad.ppp.R", "bt.frame.ppp.R", "quadscheme.ppp.R",
  "default.dummy.ppp.R", "default.n.tiling.ppp.R", "mpl.engine.ppp.R",
  "mpl.prepare.ppp.R", "evalInteraction.ppp.R", "quadBlockSizes.ppp.R",
  "evalInterEngine.ppp.R", "evaluate.ppp.R", "evalPairPotential.ppp.R"
)
lapply(source_files, source)

# Function to compute the average number of particles
get_avg_particle_num <- function(tomograms) {
  total_particles <- sum(sapply(tomograms, function(t) t$n))
  avg_particles <- floor(total_particles / length(tomograms))
  return(avg_particles)
}

print("Generating quadrature scheme...")

# Define data directory

# Parse command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]
tissue <- args[2]
intensity_type <- args[3]
n <- as.numeric(args[4])
r <- as.numeric(args[5])
hr <- as.numeric(args[6])
resample_percent <- as.numeric(args[7])
random_seed <- as.numeric(args[8])

# Assign TMC based on tissue type
TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")

# Print parameters
print(sprintf("%s %s %d %d", tissue, intensity_type, n, r))

# Load the simulated pattern data
simulated_pattern_file <- file.path(
  data_dir, TMC, tissue, "random_resample",
  sprintf(
    "random_simulated_pattern_%d_500_%d_%s_across3_resample5_percentage_%d_seed_%d_self_quad_d_no_between_dummy.Rda",
    n, hr, intensity_type, resample_percent, random_seed
  )
)
load(simulated_pattern_file)

# Store the loaded pattern
P <- marked_simulated_pattern

# Print the number of images
print(sprintf("Image num = %d", length(P)))

# Compute the average number of particles
avg_particle_num <- get_avg_particle_num(P)

# Define interaction models
interaction <- lapply(P, function(p) {
  num <- length(unique(p$marks))
  ir_matrix <- matrix(rep(r, num * num), nrow = num)
  hr_matrix <- matrix(rep(hr, num * num), nrow = num)
  MultiStraussHard(iradii = ir_matrix, hradii = hr_matrix)
})

# Create hyperframes for model fitting
H <- hyperframe(Y = P)
I <- hyperframe(Interaction = interaction)

# Perform quadrature scheme fitting
Quad_all <- mppm.quad.ppp(
  Y ~ marks, data = H, interaction = I,
  average_number = avg_particle_num, run_quad = FALSE
)

# Define the output file path
quad_output_file <- file.path(
  data_dir, TMC, tissue, "random_resample",
  sprintf(
    "quad_%d_%d_%d_%s_across3_resample5_percentage_%d_seed_%d_self_quad_d_no_between_dummy.Rda",
    n, r, hr, intensity_type, resample_percent, random_seed
  )
)

# Filter and clean quadrature points
Quad_all$moadf <- Quad_all$moadf %>%
  group_by(pattern_ID) %>%
  mutate(is_real = (.mpl.Y != 0)) %>%
  group_by(x, y, .add = TRUE) %>%
  filter(any(is_real)) %>%
  ungroup()

# Save the processed quadrature scheme
save(Quad_all, file = quad_output_file)

# Run garbage collection to free up memory
gc()

# Print any warnings
print(warnings())
