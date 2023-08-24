import scanpy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

data_dir = sys.argv[1]
# get info from cellar data
script_dir = os.path.dirname(os.path.abspath(__file__))
cellar_info = scanpy.read_h5ad(f'{script_dir}/CODEX_Florida_19-003-lymph-node-R2_cellar.h5ad')
cellar_channel_intensity = cellar_info.X
cell_labels = cellar_info.obs
channel_names_cellar = cellar_info.var_names

# extract cell type marker intensity profile from Cellar
cell_annotations = cell_labels.annotations.unique()
cellar_marker_profile = {}
for annotation in cell_annotations:
	if annotation != '':
		cell_indices = cell_labels.index[cell_labels['annotations'] == annotation].tolist()
		cell_indices = np.array(cell_indices).astype(int)
		annotation_intensity = pd.DataFrame(cellar_channel_intensity).loc[cell_indices, :]
		annotation_intensity.columns = channel_names_cellar
		annotation_intensity.index = range(annotation_intensity.shape[0])
		cellar_marker_profile[annotation] = annotation_intensity
	
# collect cell type marker intensity profiler from
nuclei_intensity = pd.read_csv(f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/reg1_stitched_expressions.ome.tiff-nuclei_channel_total.csv')
nuclei_boundary_intensity = pd.read_csv(f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/reg1_stitched_expressions.ome.tiff-nucleus_boundaries_channel_total.csv')
cell_intensity = pd.read_csv(f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/reg1_stitched_expressions.ome.tiff-cell_channel_total.csv')
cell_boundary_intensity = pd.read_csv(f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/reg1_stitched_expressions.ome.tiff-cell_boundaries_channel_total.csv')

cell_total_intensity = nuclei_intensity.loc[:, nuclei_intensity.columns != 'ID'] + \
nuclei_boundary_intensity.loc[:, nuclei_boundary_intensity.columns != 'ID'] + \
cell_intensity.loc[:, cell_intensity.columns != 'ID'] + \
cell_boundary_intensity.loc[:, cell_boundary_intensity.columns != 'ID']

channel_names_CODEX = cell_total_intensity.columns

common_channels = ['CD11c', 'CD21', 'CD4', 'CD8', 'Ki67']
common_channel_num = len(common_channels)
# check marker intensity distribution
cellar_channel_intensity_df = pd.DataFrame(cellar_channel_intensity, columns=channel_names_cellar)
cellar_channel_intensity_common_channels = cellar_channel_intensity_df.loc[:, cellar_channel_intensity_df.columns.isin(common_channels)]
cellar_channel_intensity_common_channels = cellar_channel_intensity_common_channels.reindex(sorted(cellar_channel_intensity_common_channels.columns), axis=1)


for cluster_num in [39]:
	cell_clustering = pd.read_csv(f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/reg1_stitched_expressions.ome.tiff_individual_tissue_cluster_' + str(cluster_num) + '_label_total_5_tissues_individual_z_score.csv')
	cell_clusters = cell_clustering.Cluster.unique()
	

	ppm_marker_profile = {}
	for cluster in cell_clusters:
		if cluster != '':
			cell_indices = cell_clustering.index[cell_clustering['Cluster'] == cluster].tolist()
			# for cell_index in cell_indices:
			cell_indices = np.array(cell_indices).astype(int)
			cluster_intensity = pd.DataFrame(cell_total_intensity).loc[cell_indices, :]
			cluster_intensity.columns = channel_names_CODEX
			cluster_intensity.index = range(cluster_intensity.shape[0])
			ppm_marker_profile[cluster] = cluster_intensity
			

	

	cluster_dir = f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/' + str(cluster_num) + '_clusters_test' + '/'
	if not os.path.exists(cluster_dir):
		os.makedirs(cluster_dir)
	
	for key in cellar_marker_profile.keys():
		print(key)

		current_cell_type_profile = cellar_marker_profile[key]
		print(current_cell_type_profile.shape)
		current_cell_type_profile_common_channels = current_cell_type_profile.loc[:, current_cell_type_profile.columns.isin(common_channels)]
		current_cell_type_profile_common_channels = current_cell_type_profile_common_channels.reindex(sorted(current_cell_type_profile_common_channels.columns), axis=1)
		marker_means = current_cell_type_profile_common_channels.mean(axis=0).sort_index()
		sns.boxplot(data=current_cell_type_profile_common_channels, showfliers=False)
		plt.xlabel('Markers')
		plt.ylabel('Cell total intensity')
		plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
		plt.tight_layout()
		plt.annotate('Num of Cells\n= ' + str(current_cell_type_profile_common_channels.shape[0]), xy=(0.02, 0.9), xycoords='axes fraction')
		plt.savefig(cluster_dir + '/' + key + '_common_' + str(common_channel_num) + '_channels.png', dpi=500)
		plt.close()
		pd.DataFrame.to_csv(current_cell_type_profile_common_channels, cluster_dir + '/' + key + '_common_' + str(common_channel_num) + '_channels.csv', index=False)
	
	
	
	
	for key in sorted(ppm_marker_profile.keys()):
		print(key)
		current_cell_type_profile = ppm_marker_profile[key]
		print(current_cell_type_profile.shape)
		current_cell_type_profile_common_channels = current_cell_type_profile.loc[:, current_cell_type_profile.columns.isin(common_channels)]
		marker_means = current_cell_type_profile_common_channels.mean(axis=0)
		sns.boxplot(data = current_cell_type_profile_common_channels, showfliers=False)
		plt.xlabel('Markers')
		plt.ylabel('Cell total intensity')
		plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
		plt.tight_layout()
		plt.annotate('Num of Cells\n= ' + str(current_cell_type_profile_common_channels.shape[0]), xy=(0.02, 0.9), xycoords='axes fraction')
		plt.savefig(cluster_dir + '/' + str(cluster_num) + '_cluster_' + str(key) + '_common_' + str(common_channel_num) + '_channels.png', dpi=500)
		plt.close()
		pd.DataFrame.to_csv(current_cell_type_profile_common_channels, cluster_dir + '/' + str(cluster_num) + '_cluster_' + str(key) + '_common_' + str(common_channel_num) + '_channels.csv', index=False)



	
