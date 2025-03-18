# Load required libraries
library(ggplot2)
library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(patchwork)
library(grid)
library(gridExtra)

# Define parameters
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments
tissue_list <- c("SPLEEN", "LN", "THYMUS", "SI", "LI")
resample_percent_list <- c(0, 100, 200, 300, 400, 500, 600, 700)

# Initialize data frames
mark_simulateion_df <- data.frame(matrix(ncol=5, nrow=0))
colnames(mark_simulateion_df) <- c("tissue", "seed", "devi", "resample_percent", "type")
mark_simulateion_df_idx <- 1

mark_ori_df <- data.frame(matrix(ncol=3, nrow=0))
colnames(mark_ori_df) <- c("tissue", "seed", "devi")
mark_ori_df_idx <- 1

# Loop over tissues
for (tissue in tissue_list) {
  for (seed in 1:30) {
    for (resample_percent in resample_percent_list) {
      TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")
      
      file_dir <- file.path(data_dir, TMC, tissue, "random_resample", 
                            paste("AUCROC_5_500_1_total_across3_resample5_percentage_", resample_percent, 
                                  "_seed_", seed, "_self_quad_d_no_between_dummy.Rda", sep=""))
      
      if (file.exists(file_dir)) {
        load(file_dir)
        mark_simulateion_df[mark_simulateion_df_idx, ] <- c(tissue, seed, aucroc[[2]], resample_percent, "simulated")
        mark_simulateion_df_idx <- mark_simulateion_df_idx + 1
      }
    }
  }
}

# Convert data types
mark_simulateion_df$devi <- as.numeric(mark_simulateion_df$devi)
mark_simulateion_df$resample_percent <- factor(mark_simulateion_df$resample_percent, levels=resample_percent_list)

# Generate plots
p_total <- NULL
for (tissue in tissue_list) {
  mark_simulateion_tissue_df <- mark_simulateion_df[which(mark_simulateion_df$tissue == tissue), ]
  
  p <- ggplot(data=mark_simulateion_tissue_df, aes(x=resample_percent, y=devi, group=seed)) +
    geom_line() + geom_point(size=0.5) +
    scale_x_discrete("Resampling Fraction", labels=as.character(0:7)) +
    ylab("Weighted Macro AUCROC") +
    ggtitle(tissue) +
    theme(text=element_text(size=12),
          panel.background=element_blank(),
          axis.line=element_line(colour="black"),
          strip.background=element_blank(),
          strip.text.x=element_blank(),
          plot.title=element_text(hjust=0.5)) +
    ylim(0.495, 0.65)
  
  if (is.null(p_total)) {
    p_total <- p
  } else {
    p_total <- p_total + p
  }
}

p_total <- p_total + plot_layout(ncol=3, guides="collect") + guide_area()
p_total[[2]] <- p_total[[2]] + theme(axis.title.y=element_blank())
p_total[[3]] <- p_total[[3]] + theme(axis.title.y=element_blank())
p_total[[5]] <- p_total[[5]] + theme(axis.title.y=element_blank())

ggsave(filename=file.path("./fig", "Fig_S1_AUCROC_simulation_random_resample5_new_datasets_new.png"), 
       plot=p_total, dpi=500, width=8, height=8)
