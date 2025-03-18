# Load necessary libraries
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(dplyr)
library(permute)
library(data.table)

# Set script directory dynamically
script_dir <- "./modeling/ppp/"
setwd(script_dir)

# Load required source files
source("mppm.fit.ppp.R")
source("glm.prep.R")
source("glm.fit.prep.R")
source("evaluate.ppp.R")
source("evalPairPotential.ppp.R")
source("get_devi.R")

# Function to construct formula dynamically
get_formula <- function(interactions) {
  formula_string <- ".mpl.Y ~ marks"
  for (interaction in interactions) {
    formula_string <- paste(formula_string, interaction, sep = " + ")
  }
  return(as.formula(formula_string))
}

# Parse command-line arguments
args <- commandArgs(TRUE)
data_dir <- args[1]  # Dynamically passed instead of hardcoded
tissue <- args[2]
intensity_type <- args[3]
num <- as.numeric(args[4])
r <- args[5]  # Keeping as character if intended
hr <- as.numeric(args[6])
img_idx <- as.numeric(args[7])
tile_side <- args[8]
tile_num <- args[9]

# Determine TMC location
TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")

# Load quadrature scheme
quad_file <- file.path(data_dir, TMC, tissue, "random_tile", 
                       paste0("quad_", num, "_100-", r, "_", hr, "_", intensity_type, 
                              "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, "_self_quad_d_no_between_dummy.Rda"))
load(quad_file)

# Load cell pattern
cell_pattern_file <- file.path(data_dir, TMC, tissue, "random_tile", 
                               paste0("cell_list_", num, "_", intensity_type, "_across3_", 
                                      img_idx, "_tile_", tile_side, "_", tile_num, ".Rda"))
load(cell_pattern_file)

# Define region of interest (excluding 500-pixel boundary)
cell_pattern_window <- cell_list[[1]]$window
rect <- c(cell_pattern_window$xrange[1] + 500, 
          cell_pattern_window$yrange[1] + 500, 
          cell_pattern_window$xrange[2] - 500, 
          cell_pattern_window$yrange[2] - 500)

# Filter training data within the region
training <- Quad_all_all$moadf
training <- training[training$x >= rect[1] & training$x <= rect[3] & 
                       training$y >= rect[2] & training$y <= rect[4],]
training <- training[complete.cases(training),]

# Update training data
Quad_all_all$moadf <- training

# Define interaction matrices
r_matrix <- matrix(rep(500, num * num), nrow = num)
hr_matrix <- matrix(rep(hr, num * num), nrow = num)

# Train model
print("Training model...")
formula_interaction <- colnames(training)[7:(ncol(training) - 1)]
formula_interaction <- formula_interaction[formula_interaction != "pattern_ID"]
formula <- get_formula(formula_interaction)

model_train <- mppm.fit.ppp(Data = Quad_all_all, formula, 
                            interaction = MultiStraussHard(iradii = r_matrix, hradii = hr_matrix))

# Extract model components
print(model_train$Fit$FIT$coefficients)
fmla <- model_train$Fit$fmla
family <- model_train$Fit$FIT$family
coef <- model_train$Fit$FIT$coefficients

# Clean up memory
rm(model_train)
gc()

# Save trained model components
save(fmla, file = file.path(data_dir, TMC, tissue, "random_tile", 
                            paste0("fmla_", num, "_100-", r, "_", hr, "_", intensity_type, 
                                   "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, "_self_quad_d_no_between_dummy.Rda")))
save(family, file = file.path(data_dir, TMC, tissue, "random_tile", 
                              paste0("family_", num, "_100-", r, "_", hr, "_", intensity_type, 
                                     "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, "_self_quad_d_no_between_dummy.Rda")))
save(coef, file = file.path(data_dir, TMC, tissue, "random_tile", 
                            paste0("coef_", num, "_100-", r, "_", hr, "_", intensity_type, 
                                   "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, "_self_quad_d_no_between_dummy.Rda")))

# Fit model and calculate deviance
print("Fitting model...")
deviance_ori <- get_devi(training, coef, fmla, family)
save(deviance_ori, file = file.path(data_dir, TMC, tissue, "random_tile", 
                                    paste0("devi_", num, "_100-", r, "_", hr, "_", intensity_type, 
                                           "_across3_", img_idx, "_tile_", tile_side, "_", tile_num, "_self_quad_d_no_between_dummy.Rda")))

# Output deviance
print(deviance_ori)
