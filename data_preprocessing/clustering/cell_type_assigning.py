import sys

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import bz2
import pickle
import os
from os.path import join
from sklearn.metrics.pairwise import rbf_kernel, linear_kernel, laplacian_kernel, chi2_kernel, sigmoid_kernel, polynomial_kernel
from sklearn.metrics.pairwise import euclidean_distances
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler as SS

def convert_profile_to_pairwise(profile):
	channel_num = len(profile.columns)
	pairwise_profile = pd.DataFrame()
	for i in range(channel_num):
		for j in range(channel_num):
			if i != j:
				current_pairwise_profile = profile.loc[:, profile.columns[i]] / profile.loc[:, profile.columns[j]]
				current_pairwise_name = profile.columns[i] + '/' + profile.columns[j]
				pairwise_profile[current_pairwise_name] = current_pairwise_profile
	return pairwise_profile

def gaussian_kernel(x, y):
	k_h = np.exp((-(np.linalg.norm(x-y) ** 2)) / 10000)
	return k_h

def kernel_test_statstic(profile1, profile2):
	kernel_stat = 0
	for i in range(profile1.shape[0]):
		for j in range(profile2.shape[0]):
			# print(i, j, gaussian_kernel(profile1[i], profile2[j]))
			current_gaussian = gaussian_kernel(profile1[i], profile2[j])
			if current_gaussian == np.nan:
				print(current_gaussian)
			kernel_stat += current_gaussian
	return kernel_stat

def calculate_kernel_stat(profile1, profile2):
	kernel_stat2 = np.sum(rbf_kernel(profile2.values, profile1.values, gamma=1))
	kernel_stat1 = np.sum(rbf_kernel(profile1.values, gamma=1))
	kernel_stat3 = np.sum(rbf_kernel(profile2.values, gamma=1))
	profile1_cell_num = profile1.values.shape[0]
	profile2_cell_num = profile2.values.shape[0]
	kernel_stat = (1 / (profile1_cell_num ** 2)) * kernel_stat1 - (2 / (profile1_cell_num * profile2_cell_num)) * kernel_stat2 + (1 / (profile2_cell_num ** 2)) * kernel_stat3
	return kernel_stat

data_dir = sys.argv[1]
c = 39
for cluster_num in range(c, c+1):
	cluster_dir = f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/' + str(cluster_num) + '_clusters_test' + '/'
	cell_type_list = ['cytotoxic T cell', 'endothelial cell', 'lymphocyte of B lineage', 'proliferating T cell', 'CD4-positive T cell']
	cell_type_list_vis = ['cytotoxic T cell', 'other cells', 'B cell', 'proliferating T cell', 'CD4-positive T cell']

	for cell_type in cell_type_list:
		cellar_profile = pd.read_csv(join(cluster_dir, cell_type + '_common_5_channels.csv'))
		cellar_profile_cell_num = cellar_profile.values.shape[0]
		cellar_profile = cellar_profile.div(cellar_profile.sum(axis=1), axis=0)
		cellar_profile.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
		exist_cluster_list = []
		original_kernel_stat_list = []
		for i in range(cluster_num):
			SPRM_profile_file = join(cluster_dir, str(cluster_num) + '_cluster_' + str(i) + '_common_5_channels.csv')
			if os.path.exists(SPRM_profile_file):
				SPRM_profile = pd.read_csv(SPRM_profile_file)
				SPRM_profile = SPRM_profile.div(SPRM_profile.sum(axis=1), axis=0)
				SPRM_profile.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
				SPRM_profile_cell_num = SPRM_profile.values.shape[0]
				original_kernel_stat = calculate_kernel_stat(cellar_profile, SPRM_profile)
				original_kernel_stat_list.append(original_kernel_stat)
				print(i, original_kernel_stat, SPRM_profile_cell_num)
				exist_cluster_list.append(i)
		if cell_type == cell_type_list[0]:
			cell_type_kernel_stat = np.array(original_kernel_stat_list)
		else:
			cell_type_kernel_stat = np.vstack((cell_type_kernel_stat, original_kernel_stat_list))
	
	cell_type_kernel_stat_dataframe = pd.DataFrame(cell_type_kernel_stat)
	cell_type_kernel_stat_dataframe.columns = exist_cluster_list
	cell_type_kernel_stat_dataframe.index = cell_type_list_vis
	
	ax = sns.heatmap(cell_type_kernel_stat_dataframe, cmap='magma', vmin=0, vmax=1)
	plt.xticks(np.arange(0.5, len(cell_type_kernel_stat_dataframe.columns) + 0.5, 1),
	           cell_type_kernel_stat_dataframe.columns, size=9, rotation=90)
	# plt.yticks(np.arange(0.5, len(cell_type_kernel_stat_dataframe.index) + 0.5, 1),
	#            cell_type_kernel_stat_dataframe.index, size=11)
	
	cbar = ax.collections[0].colorbar
	cbar.ax.tick_params(labelsize=13)
	plt.xticks(size=9)
	plt.yticks(size=11)
	# plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
	plt.xlabel('cluster number', size=13)
	plt.ylabel('cell type', size=13)
	plt.tight_layout()
	# plt.show()
	plt.savefig(join(cluster_dir, 'gaussian_kernel_' + str(c) + '.png'), dpi=500)
	plt.clf()
	plt.close()



for cluster_num in range(c, c+1):
	cluster_dir = f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/' + str(cluster_num) + '_clusters_test' + '/'
	cell_type_kernel_stat_pairwise = np.zeros((cluster_num, cluster_num))
	for i in range(cluster_num-1):
		for j in range(i+1, cluster_num):
			SPRM_profile_file1 = join(cluster_dir, str(cluster_num) + '_cluster_' + str(i) + '_common_5_channels.csv')
			SPRM_profile_file2 = join(cluster_dir, str(cluster_num) + '_cluster_' + str(j) + '_common_5_channels.csv')
			if os.path.exists(SPRM_profile_file1) and os.path.exists(SPRM_profile_file2):
				SPRM_profile1 = pd.read_csv(SPRM_profile_file1)
				SPRM_profile2 = pd.read_csv(SPRM_profile_file2)
				# SPRM_profile = np.log(SPRM_profile)
				SPRM_profile1 = SPRM_profile1.div(SPRM_profile1.sum(axis=1), axis=0)
				SPRM_profile2 = SPRM_profile2.div(SPRM_profile2.sum(axis=1), axis=0)
				SPRM_profile1.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
				SPRM_profile2.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
				SPRM_profile1_cell_num = SPRM_profile1.values.shape[0]
				SPRM_profile2_cell_num = SPRM_profile2.values.shape[0]
				# original_kernel_stat = np.mean(euclidean_distances(cellar_profile, SPRM_profile))
				original_kernel_stat = calculate_kernel_stat(SPRM_profile1, SPRM_profile2)
				# if cellar_profile_cell_num < SPRM_profile_cell_num:
				# 	SPRM_profile_sampled = SPRM_profile.sample(n=cellar_profile_cell_num, random_state=i)
				# 	original_kernel_stat = euclidean_distances(cellar_profile, SPRM_profile_sampled)
				#
				# else:
				# 	cellar_profile_sampled = cellar_profile.sample(n=SPRM_profile_cell_num, random_state=i)
				# 	original_kernel_stat = euclidean_distances(cellar_profile_sampled, SPRM_profile)
				cell_type_kernel_stat_pairwise[i, j] = original_kernel_stat
				print(i, original_kernel_stat, SPRM_profile_cell_num)
	
	for i in range(cluster_num-1):
		for j in range(i+1, cluster_num):
			cell_type_kernel_stat_pairwise[j, i] = cell_type_kernel_stat_pairwise[i, j]
	
	
	cell_type_kernel_stat_pairwise_dataframe = pd.DataFrame(cell_type_kernel_stat_pairwise)
	zero_labels = cell_type_kernel_stat_pairwise_dataframe.index[cell_type_kernel_stat_pairwise_dataframe.eq(0).all(axis=1)].to_list()
	
	cell_type_kernel_stat_pairwise_dataframe = cell_type_kernel_stat_pairwise_dataframe.loc[:, (cell_type_kernel_stat_pairwise_dataframe != 0).any(axis=0)]
	cell_type_kernel_stat_pairwise_dataframe = cell_type_kernel_stat_pairwise_dataframe.loc[(cell_type_kernel_stat_pairwise_dataframe != 0).any(axis=1)]
	
	cell_type_kernel_stat_pairwise = cell_type_kernel_stat_pairwise_dataframe.values
	
	
	ax = sns.heatmap(cell_type_kernel_stat_pairwise_dataframe, cmap='magma', mask = np.triu(cell_type_kernel_stat_pairwise), vmin=0, vmax=1)

	plt.xticks(np.arange(0.5, len(cell_type_kernel_stat_pairwise_dataframe.columns) + 0.5, 1),
	           cell_type_kernel_stat_pairwise_dataframe.columns, size=11, rotation=90)
	
	plt.yticks(np.arange(0.5, len(cell_type_kernel_stat_dataframe.columns) + 0.5, 1),
	           cell_type_kernel_stat_pairwise_dataframe.columns, size=11)
	

		
	cbar = ax.collections[0].colorbar
	cbar.ax.tick_params(labelsize=13)
	
	plt.xticks(size=11)
	plt.yticks(rotation=0)
	
	# plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
	plt.xlabel('cluster number', size=13)
	plt.ylabel('cluster number', size=13)
	plt.tight_layout()
	# plt.show()
	plt.savefig(join(cluster_dir, 'gaussian_kernel_pairwise_' + str(c) + '.png'), dpi=500)
	plt.clf()
	plt.close()



	cell_type_kernel_stat_together = np.vstack((cell_type_kernel_stat, cell_type_kernel_stat_pairwise)).T
	model = KMeans(n_clusters=5, random_state=1).fit(cell_type_kernel_stat_together)
	
	# align cluster num from KMeans to real num
	labels = model.labels_.astype(int)
	label_new_assignment = np.zeros((len(labels)), dtype=int)
	cell_type_new_assignment = cell_type_kernel_stat_dataframe.idxmin(axis='columns')
	for l in range(len(cell_type_new_assignment)):
		current_assignment = l
		current_cluster = cell_type_new_assignment.values[l]
		previous_cluster = labels[cell_type_kernel_stat_dataframe.columns.to_list().index(current_cluster)]
		label_new_assignment[labels == previous_cluster] = current_assignment
	

	
	cluster_num_index = np.linspace(0, c-1, c).astype(int)
	cluster_num_index = [elem for elem in cluster_num_index if elem not in zero_labels]
	cluster_num_index_labels = np.vstack((cluster_num_index, label_new_assignment)).T
	np.save(join(cluster_dir, 'total_cell_intensity_5_tissues_individual_z_score_' + str(c) + '_clusters_to_5_clusters_labels.npy'), cluster_num_index_labels)


	# cluster_num_index_labels = np.load(f'{data_dir}/total_cell_intensity_5_tissues_individual_z_score_' + str(c) + '_clusters_to_5_clusters_labels.npy')
	cell_type_cell_num = np.zeros(5)
	for i in cluster_num_index:
		cluster_dir = f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/' + str(cluster_num) + '_clusters_test' + '/'
		SPRM_profile_cluster = pd.read_csv(join(cluster_dir, str(cluster_num) + '_cluster_' + str(i) + '_common_5_channels.csv'))
		cell_num = len(SPRM_profile_cluster)
		row_index = [index for index, content in enumerate(cluster_num_index_labels) if content[0] == i]
		cell_type_cell_num[cluster_num_index_labels[row_index, 1]] += cell_num
	np.savetxt(join(cluster_dir, 'cell_type_num.txt'), cell_type_cell_num.astype(int), fmt='%i')
	print(cell_type_cell_num)
	