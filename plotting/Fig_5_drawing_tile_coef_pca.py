import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import seaborn as sns
from scipy.stats import iqr
import sys

if __name__ == '__main__':
	data_dir = sys.argv[1]
	tissue_list = ['SPLEEN', 'THYMUS', 'LN', 'SI', 'LI']

	
	fig, axes = plt.subplots(2, 3, figsize=(14, 10))
	plot_pos = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1]]
	tile_side_list = ['2500x2500 tiles',
	                  '5000x5000 tiles',
	                  'original images']
	tile_side_original_list = ['2500', '5000', '1']
	hetero_df = pd.DataFrame(columns = ["Heterogeneity", "Tissue", "Tile type"])
	for tile_side_idx in range(len(tile_side_list)):
		tile_side_original = tile_side_original_list[tile_side_idx]
		tile_side = tile_side_list[tile_side_idx]
		for tissue in tissue_list:
			if tissue == 'LI' or tissue == 'SI':
				TMC = 'Stanford'
			else:
				TMC = 'Florida'
			if tile_side == 'original images':
				coef_df = pd.read_csv(Path(data_dir, TMC, tissue, 'coef_tissue.csv'))
			else:
				coef_df = pd.read_csv(Path(data_dir, TMC, tissue, 'random_tile', 'coef_tile_' + tile_side_original + '.csv'))


			coef_df['tissue'] = tissue
			coef_df['side'] = tile_side
			if tissue == tissue_list[0] and tile_side_original == '2500':
				coef = coef_df
			else:
				coef = pd.concat([coef, coef_df])
			
	coef = coef.fillna(0)

	scaler = StandardScaler().fit(coef.values[:, :-2])
	coef_ss = scaler.transform(coef.values[:, :-2])
	pca = PCA(n_components=2).fit(coef_ss)
	coef_pca = pca.transform(coef_ss)
	coef['PC1'] = coef_pca[:, 0]
	coef['PC2'] = coef_pca[:, 1]
	
	for tile_side in tile_side_list:
		for tissue in tissue_list:
			if tissue == 'LI' or tissue == 'SI':
				TMC = 'Stanford'
			else:
				TMC = 'Florida'

			coef_tissue = coef[(coef.tissue == tissue) & (coef.side==tile_side)].values[:, :-4]
			coef_tissue_ori_df = pd.read_csv(Path(data_dir, TMC, tissue, 'random_tile', 'coef_5_100-500_1_total_across3.csv'))
			coef_tissue_ori = coef_tissue_ori_df.values.reshape(1, -1)
			
			dist_matrix = euclidean_distances(coef_tissue, coef_tissue_ori)
			hetero = np.median(dist_matrix)
			hetero_df.loc[len(hetero_df.index)] = [hetero, tissue, tile_side]


	plot_idx = -1
	for tissue in tissue_list:
		if tissue == 'LI' or tissue == 'SI':
			TMC = 'Stanford'
		else:
			TMC = 'Florida'
		plot_idx += 1
		coef_tissue = coef[coef.tissue == tissue]
		sns.set_palette("bright")
		
		if tissue == tissue_list[0]:
			ax1 = sns.scatterplot(ax=axes[plot_idx // 3, plot_idx % 3], data=coef_tissue, x='PC1', y='PC2', hue='side', legend='auto', alpha=0.7)
			ax1.legend(loc='upper left')
			plt.setp(ax1.get_legend().get_texts(), fontsize='14')  # for legend text
			ax1.tick_params(axis='both', labelsize=14)
		
		else:
			ax1 = sns.scatterplot(ax=axes[plot_idx // 3, plot_idx % 3], data=coef_tissue, x='PC1', y='PC2', hue='side', legend=False, alpha=0.7)
			ax1.tick_params(axis='both', labelsize=14)

		coef_tissue_ori_df = pd.read_csv(Path(data_dir, TMC, tissue, 'random_tile', 'coef_5_100-500_1_total_across3.csv'))
		coef_tissue_ori = coef_tissue_ori_df.values.reshape(1, -1)
		coef_tissue_ori_ss = scaler.transform(coef_tissue_ori)
		coef_tissue_ori_pca = pca.transform(coef_tissue_ori_ss)
		coef_tissue_ori_pca_df = pd.DataFrame(columns=['PC1', 'PC2', 'type'])
		coef_tissue_ori_pca_df.loc[len(coef_tissue_ori_pca_df.index)] = [coef_tissue_ori_pca[0, 0], coef_tissue_ori_pca[0, 1], 'all image']
		if tissue == tissue_list[0]:
			ax2 = sns.scatterplot(ax=axes[plot_idx // 3, plot_idx % 3], data=coef_tissue_ori_pca_df, x='PC1', y='PC2', legend='full', marker='X', color='r', alpha=0.7, label='all original images')
			ax2.legend(loc='upper left')
			plt.setp(ax2.get_legend().get_texts(), fontsize='14')  # for legend text
			ax2.tick_params(axis='both', labelsize=14)
		
		
		else:
			ax2 = sns.scatterplot(ax=axes[plot_idx // 3, plot_idx % 3], data=coef_tissue_ori_pca_df, x='PC1', y='PC2', legend=False, marker='X', color='r', alpha=0.7, label='all original images')
			ax2.tick_params(axis='both', labelsize=14)

		axes[plot_idx // 3, plot_idx % 3].set_title(str(tissue), size=16)
		axes[plot_idx // 3, plot_idx % 3].text(-31, 31, "(" + chr(plot_idx + 65) + ")", size=16)
		
		axes[plot_idx // 3, plot_idx % 3].set_xlim(-30, 30)
		axes[plot_idx // 3, plot_idx % 3].set_ylim(-30, 30)
		if plot_idx == 2 or plot_idx == 3 or plot_idx == 4:
			axes[plot_idx // 3, plot_idx % 3].set_xlabel('PC1', fontsize=16)
		else:
			axes[plot_idx // 3, plot_idx % 3].set_xlabel("")

		if plot_idx == 0 or plot_idx == 3:
			axes[plot_idx // 3, plot_idx % 3].set_ylabel('PC2', fontsize=16)
		else:
			axes[plot_idx // 3, plot_idx % 3].set_ylabel("")
	
	plot_idx += 1
	
	hetero_df1 = hetero_df[hetero_df['Tile type'] == 'original images']
	hetero_df2 = hetero_df[hetero_df['Tile type'] == '5000x5000 tiles']
	hetero_df3 = hetero_df[hetero_df['Tile type'] == '2500x2500 tiles']
	
	hetero_df1['Tile type'] = 'original\nimages'
	hetero_df2['Tile type'] = '5000x5000\ntiles'
	hetero_df3['Tile type'] = '2500x2500\ntiles'
	
	
	hetero_df = pd.concat([hetero_df1, hetero_df2, hetero_df3])
	sns.set_palette("bright")
	
	# Create the line plot with 'tab10' color palette
	ax = sns.lineplot(ax=axes[plot_idx // 3, plot_idx % 3], data=hetero_df, x='Tile type', y='Heterogeneity',
	                  hue="Tissue", palette="tab10", legend='auto')
	legend = ax.get_legend()
	legend.set_title(None)  # Remove the legend title
	plt.setp(ax.get_legend().get_texts(), fontsize='14')  # for legend text

	axes[plot_idx // 3, plot_idx % 3].text(-0.12, 58.5, "(" + chr(plot_idx + 65) + ")", size=16)
	
	axes[plot_idx // 3, plot_idx % 3].set_xlabel('')
	axes[plot_idx // 3, plot_idx % 3].set_ylabel('Heterogeneity', fontsize=14)
	plt.tick_params(axis='both', which='major', labelsize=14)
	
	plt.savefig('Fig_5_tile_all_self_quad_d_no_between_dummy.png', dpi=500)
	plt.clf()
	plt.close()

