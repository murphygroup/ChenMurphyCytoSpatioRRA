# Load required libraries
library(ggplot2)
library(ggbreak)
library(patchwork)
library(ggraph)
library(ggforce)
library(dplyr)

# Define parameters
radii <- 500
aucroc_idx <- 2
tissue_list <- c("SPLEEN", "THYMUS", "LN", "SI", "LI")
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments

# Initialize total plot variable
p_total <- NULL

# Loop through each training tissue
for (tissue1 in tissue_list) {
  
  # Assign corresponding TMC for training tissue
  TMC1 <- ifelse(tissue1 %in% c("LI", "SI"), "Stanford", "Florida")
  
  # Initialize dataframes
  AUCROC_dataframe <- data.frame(AUCROC=NA, training_tissue=NA, test_tissue=NA, training_idx=NA, test_idx=NA, tissue_type=NA)
  AUCROC_idx <- 1
  
  # Loop through each test tissue
  for (tissue2 in tissue_list) {
    
    # Assign corresponding TMC for test tissue
    TMC2 <- ifelse(tissue2 %in% c("LI", "SI"), "Stanford", "Florida")
    
    # Same tissue: Leave-one-out (LOO) AUCROC
    if (tissue1 == tissue2) {
      for (img_idx1 in 1:30) {
        AUCROC1_dir <- file.path(data_dir, TMC1, tissue1, paste0("AUCROC_loo_5_100-", radii, "_1_total_across3_", tissue1, "_", img_idx1, "_", tissue1, "_", img_idx1, "_self_quad_d_no_between_dummy.Rda"))
        
        if (file.exists(AUCROC1_dir)) {
          load(AUCROC1_dir)
          AUCROC_dataframe[AUCROC_idx, ] <- c(aucroc[[aucroc_idx]], tissue1, tissue1, img_idx1, img_idx1, "held out images\nof the same tissue")
          AUCROC_idx <- AUCROC_idx + 1
        }
      }
    }
    
    # Different tissues: Cross-tissue AUCROC
    if (tissue1 != tissue2) {
      for (img_idx1 in 1:30) {
        for (img_idx2 in 1:30) {
          AUCROC2_dir <- file.path(data_dir, TMC1, tissue1, paste0("AUCROC_loo_5_100-", radii, "_1_total_across3_", tissue1, "_", img_idx1, "_", tissue2, "_", img_idx2, "_self_quad_d_no_between_dummy.Rda"))
          
          if (file.exists(AUCROC2_dir)) {
            load(AUCROC2_dir)
            AUCROC_dataframe[AUCROC_idx, ] <- c(aucroc[[aucroc_idx]], tissue1, tissue2, img_idx1, img_idx2, "images of other tissues")
            AUCROC_idx <- AUCROC_idx + 1
          }
        }
      }
    }
  }
  
  # Convert AUCROC column to numeric
  AUCROC_dataframe$AUCROC <- as.numeric(AUCROC_dataframe$AUCROC)
  
  # Generate plot for current tissue
  p <- ggplot(AUCROC_dataframe, aes(x = reorder(test_tissue, -AUCROC, FUN = mean), y = AUCROC, fill = tissue_type)) +
    geom_violin() +
    stat_summary(fun = mean, geom = "point", shape = 4, size = 0.8, color = "black") +
    xlab('') +
    ylab('Weighted Macro AUCROC') +
    theme_minimal(base_size = 11) +
    ggtitle(tissue1) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    coord_cartesian(ylim = c(0.48, 0.87))
  
  # Combine plots
  p_total <- if (is.null(p_total)) p else p_total + p
}

# Adjust layout
p_total <- p_total + 
  plot_layout(ncol = 3, guides = "collect") + 
  guide_area()

# Remove redundant y-axis labels
for (i in c(2, 3, 5)) {
  p_total[[i]] <- p_total[[i]] + theme(axis.title.y = element_blank())
}

# Save plot
output_path <- "Fig_4_AUCROC_tissue_vs_other_tissues_AUCROC_no_dummy.png"
ggsave(filename = output_path, plot = p_total, dpi = 500, width = 6, height = 5)
