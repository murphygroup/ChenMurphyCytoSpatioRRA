# Load required libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(ggplot2)
library(dplyr)
library(permute)
library(data.table)

# Set script directory
script_dir <- './modeling/ppp/'
setwd(script_dir)

# Load required functions
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

# Function to compute the average number of particles per tomogram
get_avg_particle_num <- function(tomograms) {
  particle_num_total <- 0
  tomogram_num <- length(tomograms)
  for (i in 1:tomogram_num) {
    particle_num_total <- particle_num_total + tomograms[[i]]$n
  }
  return(floor(particle_num_total / tomogram_num))
}

# Print status message
print('Generating quadrature scheme...')

# Read command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]  # Updated to take data directory as an argument
tissue <- args[2]
intensity_type <- args[3]
n <- as.numeric(args[4])
r <- as.numeric(args[5])
hr <- as.numeric(args[6])
img_idx <- args[7]

# Determine Tissue Mapping Center (TMC)
TMC <- if (tissue %in% c('LI', 'SI')) 'Stanford' else 'Florida'

# Print input parameters
print(paste(tissue, intensity_type, n, r, sep = ' '))

# Load cell list data for the specified tissue and image index
cell_list_file <- file.path(data_dir, TMC, tissue, paste0('cell_list_', n, '_', intensity_type, '_across3_', img_idx, '.Rda'))
load(cell_list_file)
P <- cell_list

# Print number of images loaded
print(paste('Image count:', length(P)))

# Compute average particle number
avg_particle_num <- get_avg_particle_num(P)

# Initialize interaction list
interaction <- list()
for (i in 1:length(P)) {
  num <- length(unique(P[[i]]$marks))
  ir_matrix <- matrix(rep(r, num * num), nrow = num)
  hr_matrix <- matrix(rep(hr, num * num), nrow = num)
  interaction[[i]] <- MultiStraussHard(iradii = ir_matrix, hradii = hr_matrix)
}

# Create hyperframes for input data and interactions
H <- hyperframe(Y = P)
I <- hyperframe(Interaction = interaction)

# Compute quadrature scheme
Quad_all <- mppm.quad.ppp(Y ~ marks, data = H, interaction = I, average_number = avg_particle_num, run_quad = FALSE)

# Save output quadrature data
quad_file <- file.path(data_dir, TMC, tissue, paste0('quad_', n, '_', r, '_', hr, '_', intensity_type, '_across3_', img_idx, '_self_dummy_d_eps_1000.Rda'))
save(Quad_all, file = quad_file)

# Perform garbage collection
gc()

# Print any warnings generated during execution
print(warnings())
