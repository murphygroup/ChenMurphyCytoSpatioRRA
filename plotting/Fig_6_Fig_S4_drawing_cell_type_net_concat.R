library(igraph)
library(ggplot2)
library(ggplotify)
library(gridExtra)

args <- commandArgs(TRUE)
data_dir <- args[1]  


tissue_list = c('SPLEEN','THYMUS', 'LN', 'SI', 'LI')
plot_list <- list()
cluster_num = 5
png(paste('Fig_6_cell_type_net_across3_quad_d_no_between_dummy.png', sep = ''), width = 3000, height = 2000, res = 300)
par(mfrow = c(2, 3), mar = c(1, 1, 2, 0))
cell_num_tissues = c()
for (tissue in tissue_list){
  if (tissue == 'LI' | tissue == 'SI'){
    TMC <- 'Stanford'
  } else{
    TMC <- 'Florida'
  }
  load(paste(data_dir,'/', TMC, '/', tissue, '/coef_5_100-500_1_total_across3_self_quad_d_no_between_dummy.Rda', sep = ''))
  coef_intensity_tissue <- coef[1:5]
  B <- c()
  for (i in 1:cluster_num){
    if (i == 1){
      B_current <- coef_intensity_tissue[i]
    } else{
      B_current <- coef_intensity_tissue[i] + coef_intensity_tissue[1]
    }
    B <- c(B, B_current)
  }
  coef_intensity_tissue <- B
  if (tissue == tissue_list[1]){
    coef_intensity_all = coef_intensity_tissue
  } else{
    coef_intensity_all = rbind(coef_intensity_all, coef_intensity_tissue)
  }
  
  load(paste(data_dir,'/', TMC, '/', tissue, '/cell_list_5_total_across3.Rda', sep = ''))
  cell_num_images = 0
  for (image_idx in 1:length(cell_list_all)){
    cell_num_images = cell_num_images + cell_list_all[[image_idx]]$n
  }
  cell_num_tissues = c(cell_num_tissues, cell_num_images)
}

coef_intensity_all = exp(coef_intensity_all)
sum_vector <- rowSums(coef_intensity_all)
tissue_normalized_matrix <- sweep(coef_intensity_all, 1, sum_vector, "/")
sum_vector_celltype <- colSums(tissue_normalized_matrix)
coef_intensity_all_normalized <- sweep(tissue_normalized_matrix, 2, sum_vector_celltype, "/")
for (tissue in tissue_list){
  tissue_idx = which(tissue_list == tissue)
  if (tissue == 'LI' | tissue == 'SI'){
    TMC = 'Stanford'
  } else{
    TMC = 'Florida'
  }
  
  
  load(paste(data_dir,'/', TMC, '/', tissue, '/coef_5_100-500_1_total_across3_self_quad_d_no_between_dummy.Rda', sep = ''))
  coef_interaction = coef[6:length(coef)]
  cell_type_num = 5
  cell_type_list = c('cytotoxic T cell', 'other cells', 'lymphocyte of B lineage', 'proliferating T cell', 'CD4-positive T cell')
  
  for (i in 0:(cell_type_num-1)){
    for (j in (i):(cell_type_num-1)){
      for (r in seq(100, 500, 100)){
        current_interaction_name = paste('InteractionmarkX', i, 'xX', j, 'x', r, sep = '')
        coef_interaction_idx = grep(current_interaction_name, names(coef_interaction))
        current_coef = coef_interaction[coef_interaction_idx]
        coef_interaction[coef_interaction_idx] = current_coef
      }
    }
  }

  cell_type_links = data.frame(matrix(ncol = 5))
  colnames(cell_type_links) = c('from', 'to', 'weight', 'sign', 'tissue')
  cell_type_self_interaction_links = data.frame(matrix(ncol = 5))
  colnames(cell_type_self_interaction_links) = c('from', 'to', 'weight', 'sign', 'tissue')
  cell_type_nodes = data.frame(matrix(ncol = 3, nrow = cell_type_num))
  colnames(cell_type_nodes) = c('cell_type', 'influence', 'tissue')
  
  interaction_names = names(coef_interaction)
  idx = 1
  for (i in 0:(cell_type_num-2)){
    for (j in (i+1):(cell_type_num-1)){
      for (r in seq(100, 500, 100)){
        current_interaction_name = paste('InteractionmarkX', i, 'xX', j, 'x', r, sep = '')
        current_coef = coef_interaction[grep(current_interaction_name,names(coef_interaction))]

        if (sign(current_coef) == 1){
          current_sign = 'darkblue'
        } else{
          current_sign = 'darkred'
        }
        cell_type_links[idx,] = c(cell_type_list[i+1], cell_type_list[j+1], abs(current_coef), current_sign, tissue)
        idx = idx + 1
      }
    }
  }
  
  idx = 1
  for (i in 0:(cell_type_num-1)){
    current_interaction_name = paste('X', i, sep = '')
    current_influence = 0
    for (j in 1:length(coef_interaction)){
      if (grepl(current_interaction_name, interaction_names[j], fixed=TRUE) 
          & grepl(paste('InteractionmarkX', i, 'xX', i, sep = ''), interaction_names[j], fixed=TRUE)){
        if (sign(coef_interaction[j]) == 1){
          current_sign = 'darkblue'
        } else{
          current_sign = 'darkred'
        }
        cell_type_self_interaction_links[idx,] = c(cell_type_list[i+1], paste(cell_type_list[i+1], '2', sep='_'), abs(coef_interaction[j]), current_sign, tissue)
        current_influence = current_influence + coef_interaction[j]
        idx = idx + 1
      }
    }
    cell_type_nodes[(i+1),] = c(cell_type_list[i+1], current_influence, tissue)
  }
  
  cell_type_self_interaction_nodes = cell_type_nodes
  
  x = c(0, 0, 0, 0, 0)
  y = c(4, 3, 2, 1, 0)
  cell_type_self_interaction_nodes = cbind(cell_type_self_interaction_nodes, x, y)
  cell_type_self_interaction_nodes_copy = cell_type_self_interaction_nodes
  cell_type_self_interaction_nodes_copy$cell_type = paste(cell_type_self_interaction_nodes_copy$cell_type,2,sep='_')
  cell_type_self_interaction_nodes_copy$x = cell_type_self_interaction_nodes_copy$x + 1
  cell_type_self_interaction_nodes = rbind(cell_type_self_interaction_nodes, cell_type_self_interaction_nodes_copy)
  
  
  x = c(0, -0.951, 0.951, -0.588, 0.588)
  y = c(1, 0.309, 0.309, -0.809, -0.809)
  cell_type_nodes = cbind(cell_type_nodes, x, y)
  if (tissue == tissue_list[1]){
    cell_type_nodes_all_tissues = cell_type_nodes
    cell_type_links_all_tissues = cell_type_links
    cell_type_self_interaction_nodes_all_tissues = cell_type_self_interaction_nodes
    cell_type_self_interaction_links_all_tissues = cell_type_self_interaction_links
  } else{
    cell_type_nodes_all_tissues = rbind(cell_type_nodes_all_tissues, cell_type_nodes)
    cell_type_links_all_tissues = rbind(cell_type_links_all_tissues, cell_type_links)
    cell_type_self_interaction_nodes_all_tissues = rbind(cell_type_self_interaction_nodes_all_tissues, cell_type_self_interaction_nodes)
    cell_type_self_interaction_links_all_tissues = rbind(cell_type_self_interaction_links_all_tissues, cell_type_self_interaction_links)
  }
}


node_weight = as.numeric(cell_type_nodes_all_tissues$influence)
positive_node_weight <- node_weight - min(node_weight)
scaled_node_weight <- log10(node_weight)
min_value = min(scaled_node_weight)
shifted_node_weight <- exp(scaled_node_weight - min_value + 1)

cell_type_nodes_all_tissues$influence = shifted_node_weight
  
edge_weight = as.numeric(cell_type_links_all_tissues$weight)
scaled_edge_weight <- log10(edge_weight)
min_value = min(scaled_edge_weight)

shifted_edge_weight <- exp(scaled_edge_weight - min_value + 1)


cell_type_links_all_tissues$weight = shifted_edge_weight

for (tissue in tissue_list){
  cell_type_nodes = cell_type_nodes_all_tissues[which(cell_type_nodes_all_tissues$tissue == tissue),]
  cell_type_links = cell_type_links_all_tissues[which(cell_type_links_all_tissues$tissue == tissue),]
  net <- graph_from_data_frame(d=cell_type_links, vertices=cell_type_nodes, directed=F) 
  color = c('green4', 'magenta', 'cyan', 'red1', 'yellow2')
  V(net)$color = color
  node_weight = as.numeric(V(net)$influence)
  V(net)$size <- node_weight * 10
  
  edge_weight = as.numeric(E(net)$weight) 
  E(net)$width <- edge_weight / 0.3e2
  E(net)$edge.color <- E(net)$sign
  plot(net, edge.color = E(net)$sign, vertex.label=NA)
  title(main = tissue,  cex.main=2.5)

}
dev.off()



png(paste('Fig_S4_cell_type_self_interaction_net_across3_quad_d_no_between_dummy.png', sep = ''), width = 3000, height = 2500, res = 300)
par(mfrow = c(2, 3), mar = c(1, 1, 2, 0))


node_weight = as.numeric(cell_type_self_interaction_nodes_all_tissues$influence)
edge_weight = as.numeric(cell_type_self_interaction_links_all_tissues$weight)
scaled_edge_weight <- log10(edge_weight)
min_value = min(scaled_edge_weight)
shifted_edge_weight <- exp(scaled_edge_weight - min_value + 1)
cell_type_self_interaction_links_all_tissues$weight <- shifted_edge_weight / 0.3e2

for (tissue in tissue_list){
  cell_type_nodes = cell_type_self_interaction_nodes_all_tissues[which(cell_type_self_interaction_nodes_all_tissues$tissue == tissue),]
  cell_type_links = cell_type_self_interaction_links_all_tissues[which(cell_type_self_interaction_links_all_tissues$tissue == tissue),]
  net <- graph_from_data_frame(d=cell_type_links, vertices=cell_type_nodes, directed=F) 
  color = c('green4', 'magenta', 'cyan', 'red1', 'yellow2', 'green4', 'magenta', 'cyan', 'red1', 'yellow2')
  V(net)$color = color
  node_weight = as.numeric(V(net)$influence)
  scaled_node_weight <- log10(node_weight)
  min_value = min(scaled_node_weight)
  shifted_node_weight <- exp(scaled_node_weight - min_value + 1)
  V(net)$size <- shifted_node_weight  * 10
  edge_weight = as.numeric(E(net)$weight)
  E(net)$width <- edge_weight
  E(net)$edge.color <- E(net)$sign
  plot(net, edge.color = E(net)$sign, vertex.label=NA)
  title(main = tissue,  cex.main=2.5)

}
dev.off()



