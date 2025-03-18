# Load required libraries
library(igraph)
library(ggplot2)
library(ggplotify)
library(gridExtra)

# Define parameters
tissue_list <- c("SPLEEN", "THYMUS", "LN", "SI", "LI")
cluster_num <- 5
args <- commandArgs(TRUE)
data_dir <- args[1]  # Set data directory from arguments

# Initialize variables
coef_intensity_all <- NULL
cell_num_tissues <- c()

# Set up PNG output
png(file.path("./fig/Fig_6_cell_type_net_across3_quad_d_no_between_dummy.png"),
    width = 3000, height = 2000, res = 300)

par(mfrow = c(2, 3), mar = c(1, 1, 2, 0))

# Process each tissue
for (tissue in tissue_list) {
  
  TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")
  
  # Load coefficient file
  coef_file <- file.path(data_dir, TMC, tissue, "coef_5_100-500_1_total_across3_self_quad_d_no_between_dummy.Rda")
  load(coef_file)
  
  # Process intensity coefficients
  coef_intensity_tissue <- coef[1:5]
  B <- coef_intensity_tissue
  B[-1] <- B[-1] + coef_intensity_tissue[1]
  coef_intensity_tissue <- B
  
  coef_intensity_all <- if (is.null(coef_intensity_all)) coef_intensity_tissue else rbind(coef_intensity_all, coef_intensity_tissue)
  
  # Load cell list file
  cell_list_file <- file.path(data_dir, TMC, tissue, "cell_list_5_total_across3.Rda")
  load(cell_list_file)
  
  # Compute total cell count
  cell_num_images <- sum(sapply(cell_list_all, function(x) x$n))
  cell_num_tissues <- c(cell_num_tissues, cell_num_images)
}

# Normalize intensity coefficients
coef_intensity_all <- exp(coef_intensity_all)
sum_vector <- rowSums(coef_intensity_all)
tissue_normalized_matrix <- sweep(coef_intensity_all, 1, sum_vector, "/")
sum_vector_celltype <- colSums(tissue_normalized_matrix)
coef_intensity_all_normalized <- sweep(tissue_normalized_matrix, 2, sum_vector_celltype, "/")

# Process each tissue again
for (tissue in tissue_list) {
  TMC <- ifelse(tissue %in% c("LI", "SI"), "Stanford", "Florida")
  tissue_idx <- which(tissue_list == tissue)
  
  # Load coefficient file
  coef_file <- file.path(data_dir, TMC, tissue, "coef_5_100-500_1_total_across3_self_quad_d_no_between_dummy.Rda")
  load(coef_file)
  
  coef_interaction <- coef[6:length(coef)]
  cell_type_list <- c("cytotoxic T cell", "other cells", "lymphocyte of B lineage", "proliferating T cell", "CD4-positive T cell")
  
  cell_type_links <- data.frame(from=character(), to=character(), weight=numeric(), sign=character(), tissue=character())
  cell_type_self_interaction_links <- data.frame(from=character(), to=character(), weight=numeric(), sign=character(), tissue=character())
  cell_type_nodes <- data.frame(cell_type=character(), influence=numeric(), tissue=character())
  
  # Process cell-cell interactions
  for (i in 0:(length(cell_type_list)-2)) {
    for (j in (i+1):(length(cell_type_list)-1)) {
      for (r in seq(100, 500, 100)) {
        interaction_name <- paste0("InteractionmarkX", i, "xX", j, "x", r)
        coef_idx <- grep(interaction_name, names(coef_interaction))
        current_coef <- coef_interaction[coef_idx]
        
        sign_color <- ifelse(sign(current_coef) == 1, "darkblue", "darkred")
        cell_type_links <- rbind(cell_type_links, data.frame(from=cell_type_list[i+1], to=cell_type_list[j+1], weight=abs(current_coef), sign=sign_color, tissue=tissue))
      }
    }
  }
  
  # Process self-interactions
  for (i in 0:(length(cell_type_list)-1)) {
    interaction_name <- paste0("X", i)
    current_influence <- 0
    for (j in seq_along(coef_interaction)) {
      if (grepl(interaction_name, names(coef_interaction)[j], fixed=TRUE) &&
          grepl(paste0("InteractionmarkX", i, "xX", i), names(coef_interaction)[j], fixed=TRUE)) {
        
        sign_color <- ifelse(sign(coef_interaction[j]) == 1, "darkblue", "darkred")
        cell_type_self_interaction_links <- rbind(cell_type_self_interaction_links, data.frame(from=cell_type_list[i+1], to=paste0(cell_type_list[i+1], "_2"), weight=abs(coef_interaction[j]), sign=sign_color, tissue=tissue))
        current_influence <- current_influence + coef_interaction[j]
      }
    }
    cell_type_nodes <- rbind(cell_type_nodes, data.frame(cell_type=cell_type_list[i+1], influence=current_influence, tissue=tissue))
  }
  
  # Normalize node and edge weights
  node_weight <- as.numeric(cell_type_nodes$influence)
  shifted_node_weight <- exp(log10(node_weight) - min(log10(node_weight)) + 1)
  cell_type_nodes$influence <- shifted_node_weight
  
  edge_weight <- as.numeric(cell_type_links$weight)
  shifted_edge_weight <- exp(log10(edge_weight) - min(log10(edge_weight)) + 1)
  cell_type_links$weight <- shifted_edge_weight / 0.3e2
  
  # Generate network plot
  net <- graph_from_data_frame(d=cell_type_links, vertices=cell_type_nodes, directed=FALSE)
  color <- c("green4", "magenta", "cyan", "red1", "yellow2")
  V(net)$color <- color
  V(net)$size <- shifted_node_weight * 10
  E(net)$width <- shifted_edge_weight / 0.3e2
  E(net)$edge.color <- E(net)$sign
  
  plot(net, edge.color=E(net)$sign, vertex.label=NA)
  title(main=tissue, cex.main=2.5)
}

dev.off()

# Save self-interaction network
png(file.path("./fig/Fig_S4_cell_type_self_interaction_net_across3_quad_d_no_between_dummy.png"),
    width = 3000, height = 2500, res = 300)

par(mfrow = c(2, 3), mar = c(1, 1, 2, 0))

for (tissue in tissue_list) {
  cell_type_nodes <- subset(cell_type_nodes, tissue == tissue)
  cell_type_links <- subset(cell_type_self_interaction_links, tissue == tissue)
  
  net <- graph_from_data_frame(d=cell_type_links, vertices=cell_type_nodes, directed=FALSE)
  
  color <- rep(c("green4", "magenta", "cyan", "red1", "yellow2"), 2)
  V(net)$color <- color
  V(net)$size <- shifted_node_weight * 10
  E(net)$width <- shifted_edge_weight / 0.3e2
  E(net)$edge.color <- E(net)$sign
  
  plot(net, edge.color=E(net)$sign, vertex.label=NA)
  title(main=tissue, cex.main=2.5)
}

dev.off()
