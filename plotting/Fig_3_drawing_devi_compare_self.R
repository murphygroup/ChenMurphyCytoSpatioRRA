# Load required libraries
library(ggplot2)
library(ggbreak)
library(patchwork)
library(ggraph)
library(grid)

# Define tissue list and data directory
tissue_list <- c("SPLEEN", "LN", "THYMUS", "SI", "LI")
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments

# Initialize empty plot variable
p_total <- NULL

# Loop through each tissue type
for (tissue in tissue_list) {
  
  # Assign corresponding TMC based on tissue type
  TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")
  
  # Load deviance data
  file_path <- file.path(data_dir, TMC, tissue, "devi_5_100-500_1_total_across3_self_quad_d_no_between_dummy.Rda")
  load(file_path)
  
  # Initialize vectors
  devi_real <- deviance_ori$devi_real_avg
  devi_dummy <- deviance_ori$devi_dummy_avg
  devi_all <- deviance_ori$devi_all_avg
  
  # Iterate over different radius values
  for (i in seq(100, 500, 100)) {
    file_path <- file.path(data_dir, TMC, tissue, paste0("devi_5_", i, "_1_total_across3_self_quad_d_no_between_dummy.Rda"))
    load(file_path)
    
    devi_real <- c(devi_real, deviance_ori$devi_real_avg)
    devi_dummy <- c(devi_dummy, deviance_ori$devi_dummy_avg)
    devi_all <- c(devi_all, deviance_ori$devi_all_avg)
  }
  
  # Prepare data for visualization
  radii_range <- c("100-500", seq(100, 500, 100))
  devi_data <- data.frame(
    radii = radii_range,
    all = devi_all,
    real = devi_real,
    dummy = devi_dummy
  )
  
  devi_data_long <- reshape2::melt(devi_data, id.vars = "radii", variable.name = "Cell", value.name = "devi")
  devi_data_long$radii <- factor(devi_data_long$radii, levels = unique(devi_data_long$radii))
  
  # Generate plot for current tissue
  p <- ggplot(devi_data_long, aes(x = radii, y = devi, group = Cell, colour = Cell)) +
    geom_line() +
    geom_point() +
    facet_wrap(~Cell, scales = "free_y", ncol = 1) +
    scale_x_discrete("Strauss radii", guide = guide_axis(angle = 45)) +
    ylab("Average deviance per cell") +
    ggtitle(tissue) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black"),
      strip.text = element_blank(),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.key.size = unit(1, "cm")
    )
  
  # Combine plots
  p_total <- if (is.null(p_total)) p else p_total + p
}

# Layout adjustment
p_total <- p_total + 
  plot_layout(ncol = 3, guides = "collect") + 
  guide_area() + 
  guides(colour = guide_legend(title = NULL))

# Remove redundant y-axis labels
for (i in c(2, 3, 5)) {
  p_total[[i]] <- p_total[[i]] + theme(axis.title.y = element_blank())
}

# Save the final plot
output_path <- "./fig/Fig_3_all_tissue_devi_100-500_rerun.png"
png(filename = output_path, width = 4500, height = 3000, res = 500)
print(p_total)

# Add text labels using grid
grid.text("Strauss radii", x = 0.5, y = 0.025, gp = gpar(fontsize = 15))
grid.text("Average deviance per cell", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 15))

dev.off()
