# Load required libraries
library(ggplot2)
library(patchwork)
library(grid)
library(reshape2)

# Define tissue list and data directory
tissue_list <- c("SPLEEN", "LN", "THYMUS", "SI", "LI")
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments

# Initialize empty plot variable
p_total <- NULL

# Initialize dataframe to save supplementary data
supplementary_data <- data.frame()

# Loop through each tissue type
for (tissue in tissue_list) {
  
  # Assign corresponding TMC based on tissue type
  TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")
  
  # Load initial deviance data with CI
  file_path <- file.path(data_dir, TMC, tissue, "devi_5_100-500_1_total_across3_self_quad_d_no_between_dummy_SE_CI.Rda")
  load(file_path)
  
  # Initialize vectors with initial radius values
  devi_real_avg <- c(deviance_ori$devi_real_avg)
  devi_dummy_avg <- c(deviance_ori$devi_dummy_avg)
  devi_all_avg <- c(deviance_ori$devi_all_avg)
  
  # CI upper and lower initialization
  devi_real_upper <- c(deviance_ori$devi_real_ci["upper_CI"])
  devi_real_lower <- c(deviance_ori$devi_real_ci["lower_CI"])
  devi_dummy_upper <- c(deviance_ori$devi_dummy_ci["upper_CI"])
  devi_dummy_lower <- c(deviance_ori$devi_dummy_ci["lower_CI"])
  devi_all_upper <- c(deviance_ori$devi_all_ci["upper_CI"])
  devi_all_lower <- c(deviance_ori$devi_all_ci["lower_CI"])
  
  # Iterate over radius values (100,200,...,500)
  for (i in seq(100, 500, 100)) {
    file_path <- file.path(data_dir, TMC, tissue, paste0("devi_5_", i, "_1_total_across3_self_quad_d_no_between_dummy_SE_CI.Rda"))
    load(file_path)
    
    devi_real_avg <- c(devi_real_avg, deviance_ori$devi_real_avg)
    devi_dummy_avg <- c(devi_dummy_avg, deviance_ori$devi_dummy_avg)
    devi_all_avg <- c(devi_all_avg, deviance_ori$devi_all_avg)
    
    devi_real_upper <- c(devi_real_upper, deviance_ori$devi_real_ci["upper_CI"])
    devi_real_lower <- c(devi_real_lower, deviance_ori$devi_real_ci["lower_CI"])
    
    devi_dummy_upper <- c(devi_dummy_upper, deviance_ori$devi_dummy_ci["upper_CI"])
    devi_dummy_lower <- c(devi_dummy_lower, deviance_ori$devi_dummy_ci["lower_CI"])
    
    devi_all_upper <- c(devi_all_upper, deviance_ori$devi_all_ci["upper_CI"])
    devi_all_lower <- c(devi_all_lower, deviance_ori$devi_all_ci["lower_CI"])
  }
  
  # Prepare data for plotting with CI
  radii_range <- c("100-500", seq(100, 500, 100))
  devi_data <- data.frame(
    Tissue = tissue,
    Radii = radii_range,
    all = devi_all_avg,
    real = devi_real_avg,
    dummy = devi_dummy_avg,
    all_upper = devi_all_upper,
    all_lower = devi_all_lower,
    real_upper = devi_real_upper,
    real_lower = devi_real_lower,
    dummy_upper = devi_dummy_upper,
    dummy_lower = devi_dummy_lower
  )
  
  # Append to supplementary data
  supplementary_data <- rbind(supplementary_data, devi_data)
  
  # Convert data to long format
  devi_data_long <- melt(devi_data[, c("Radii", "all", "real", "dummy")], id.vars = "Radii",
                         variable.name = "Cell", value.name = "devi")
  
  devi_data_long$upper_CI <- melt(devi_data[, c("Radii", "all_upper", "real_upper", "dummy_upper")], id.vars = "Radii")$value
  devi_data_long$lower_CI <- melt(devi_data[, c("Radii", "all_lower", "real_lower", "dummy_lower")], id.vars = "Radii")$value
  
  devi_data_long$Radii <- factor(devi_data_long$Radii, levels = unique(devi_data_long$Radii))
  
  # Plot with smaller points and narrower black error bars
  p <- ggplot(devi_data_long, aes(x = Radii, y = devi, colour = Cell)) +
    geom_point(size=1.5) +
    geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.15, colour = "black") +
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
      legend.key.size = unit(1, "cm"),
      panel.grid = element_blank(),
      panel.background = element_blank()
    )
  
  # Combine plots
  p_total <- if (is.null(p_total)) p else p_total + p
}

# Save supplementary data
write.csv(supplementary_data, file="Fig_3_supplementary_data_table.csv", row.names = FALSE)

# Adjust plot layout
p_total <- p_total + plot_layout(ncol = 3, guides = "collect") + guide_area() + guides(colour = guide_legend(title = NULL))

# Save the final plot
output_path <- "Fig_3_all_tissue_devi_100-500_CI.png"
png(filename = output_path, width = 4500, height = 3000, res = 500)
print(p_total)

# Add text labels
grid.text("Strauss radii", x = 0.5, y = 0.025, gp = gpar(fontsize = 15))
grid.text("Average deviance per cell", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 15))

dev.off()
