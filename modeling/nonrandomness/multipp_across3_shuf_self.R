library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(ggplot2)
library(dplyr)
library(permute)
library(data.table)

# Set script directory
script_dir = './modeling/ppp/'
setwd(script_dir)

# Load required function scripts
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

# Function to compute the average number of particles across tomograms
get_avg_particle_num = function(tomograms) {
  particle_num_total = 0
  tomogram_num = length(tomograms)
  
  for (i in 1:tomogram_num) {
    particle_num_total = particle_num_total + tomograms[[i]]$n
  }
  
  return(floor(particle_num_total / tomogram_num))
}

print('Generating quadrature scheme...')

# Set data directory
data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'

# Read command line arguments
args = commandArgs(TRUE)
data_dir = args[1]
tissue = args[2]
intensity_type = args[3]
n = as.numeric(args[4])
r = as.numeric(args[5])
hr = as.numeric(args[6])
shuf_num = as.numeric(args[7])

# Determine Tissue Mapping Center (TMC) based on tissue type
if (tissue == 'LI' | tissue == 'SI') {
  TMC = 'Stanford'
} else {
  TMC = 'Florida'
}

print(paste(tissue, intensity_type, n, r, sep = ' '))

# Load cell list data
load(file.path(data_dir, TMC, tissue, paste('cell_list_', n, '_', intensity_type, '_across3.Rda', sep = '')))
P = cell_list_all

# Shuffle marks within each pattern using a fixed seed
for (p in 1:length(P)) {
  set.seed(shuf_num)
  P[[p]]$marks = sample(P[[p]]$marks, replace = TRUE)
}

print(paste('image num = ', length(P), sep = ''))

# Compute average number of particles per tomogram
avg_particle_num = get_avg_particle_num(P)

# Initialize interaction list and validity check
interaction = list()
valid = TRUE

# Construct interaction model for each pattern
for (i in 1:length(P)) {
  num = length(unique(P[[i]]$marks))
  print(c(num, levels(P[[i]]$marks)))
  
  r_matrix = matrix(rep(r, num * num), nrow = num)
  hr_matrix = matrix(rep(hr, num * num), nrow = num)
  interaction[[i]] = MultiStraussHard(iradii = r_matrix, hradii = hr_matrix)
  
  # Check if the number of unique marks matches the expected number
  if (num != n) {
    valid = FALSE
  }
}

# Proceed if valid, otherwise print a warning
if (valid) {
  H = hyperframe(Y = P)
  I = hyperframe(Interaction = interaction)
  print(H)
  
  # Generate quadrature scheme for point process modeling
  Quad_all_all = mppm.quad.ppp(Y ~ marks, data = H, interaction = I, average_number = avg_particle_num, run_quad = FALSE)
  
  # Save the quadrature scheme
  filename = file.path(data_dir, TMC, tissue, 'random_shuf', paste('quad_', n, '_', r, '_', hr, '_', intensity_type, '_across3_shuf_', shuf_num, '_self_quad_d_no_between_dummy.Rda', sep = ''))
  save(Quad_all_all, file = filename)
  
  # Run garbage collection
  gc()
} else {
  print('Missing cell type after shuffling')
}

# Print any warnings that occurred
print(warnings())
