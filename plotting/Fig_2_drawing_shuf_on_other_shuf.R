# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define tissue list and data directory
tissue_list <- c("SPLEEN", "THYMUS", "LN", "SI", "LI")
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments

# Initialize data frames
devi_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(devi_df) <- c("devi", "cat", "tissue")

devi_hline_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(devi_hline_df) <- c("devi", "cat", "tissue")

# Track index positions
total_idx <- 1
total_hline_idx <- 1

# Process each tissue
for (tissue in tissue_list) {
  TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")
  
  # Process shuffled models on other shuffled patterns
  devi_shuf_on_another_shuf_file_list <- Sys.glob(file.path(
    data_dir, TMC, tissue, "random_shuf",
    "devi_shuf_on_another_shuf_5_100_1_total_across3_*_on_*"
  ))
  
  for (file in devi_shuf_on_another_shuf_file_list) {
    load(file)
    devi_df[total_idx, ] <- c(deviance_shuf[[1]], "cell-type-shuffled models on\nother cell type-shuffled patterns", tissue)
    total_idx <- total_idx + 1
  }
  
  # Process shuffled models on original patterns
  devi_shuf_on_ori_file_list <- Sys.glob(file.path(
    data_dir, TMC, tissue, "random_shuf",
    "devi_shuf_on_ori_5_100_1_total_across3_shuf_*"
  ))
  
  for (file in devi_shuf_on_ori_file_list) {
    load(file)
    devi_df[total_idx, ] <- c(deviance_shuf_on_ori[[1]], "cell type-shuffled models\non original patterns", tissue)
    total_idx <- total_idx + 1
  }
}

# Convert data to numeric and factor for proper visualization
devi_df$devi <- as.numeric(devi_df$devi)
devi_df$tissue <- factor(devi_df$tissue, levels = tissue_list)

# Apply log transformation
devi_df$devi <- log(devi_df$devi)
devi_hline_df$devi <- log(as.numeric(devi_hline_df$devi))

# Compute boxplot stats for y-axis limits
sts <- unlist(lapply(tissue_list, function(tissue_ylim) {
  devi_df_ylim <- devi_df[which(devi_df$tissue == tissue_ylim), 1]
  boxplot.stats(devi_df_ylim)$stats
}))

# Create the plot
p <- ggplot(data = devi_df, aes(x = tissue, y = devi, color = cat)) +
  geom_boxplot(position = "identity", outlier.shape = NA) +
  coord_cartesian(ylim = c(min(sts), max(sts))) +
  xlab("Tissue") +
  ylab("log(Average deviance per cell)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 0.3, colour = "black", linetype = 1),
    legend.title = element_blank(),
    legend.position = "top",
    text = element_text(size = 14)
  )

# Save the plot
ggsave(filename = "Fig_2_shuf_on_other_shuf_no_dummy.png",
       plot = p, dpi = 500, width = 8, height = 5)
