import sys
import numpy as np
import pandas as pd
import glob
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from os.path import join
def get_common_markers(files):
	common_markers = []
	for file in files:
		if file == files[0]:
			common_markers = pd.read_csv(file).columns.str.lower()
		else:
			common_markers = np.intersect1d(common_markers, pd.read_csv(file).columns.str.lower())
	common_markers = common_markers[common_markers != 'ID']
	common_markers = np.insert(common_markers, 0, 'ID', axis=0)
	return common_markers

def get_intensity_matrix(files):
	common_markers = ['ID', 'CD11c', 'CD21', 'CD4', 'CD8', 'Ki67']
	print(common_markers)
	intensity_image_pieces = []
	cell_num = []
	for file in files:
		print(file)
		intensity = pd.read_csv(file)[common_markers].values
		intensity_image_pieces.append(intensity)
		cell_num.append(intensity.shape[0])
	intensity_image_stacked = np.vstack(intensity_image_pieces)
	cell_id = intensity_image_stacked[:, 0].astype(int)
	return intensity_image_stacked[:, 1:], cell_id, cell_num
		

if __name__ == '__main__':
	data_dir = sys.argv[1]
	for intensity_type in ['total']:
		all_feature_matrix_z_5_tissues_pieces = []
		cell_file_list_5_tissues = []
		cell_id_5_tissues = []
		cell_num_5_tissues = []
		tissue_type = ['LN', 'THYMUS', 'SPLEEN', 'LI', 'SI']
		for tissue in tissue_type:
			cytoplasm_intensity_file_list = sorted(glob.glob(join(data_dir, '**', tissue, '**', '*cell_channel_' + intensity_type + '.csv'), recursive=True))
			print(tissue)
			print(len(cytoplasm_intensity_file_list))
			# print(cytoplasm_intensity_file_list)
			cytoplasm_feature_matrix, cell_id, cell_num = get_intensity_matrix(cytoplasm_intensity_file_list)
			cell_file_list_5_tissues = cell_file_list_5_tissues + cytoplasm_intensity_file_list
			cell_id_5_tissues = cell_id_5_tissues + cell_id.tolist()
			cell_num_5_tissues = cell_num_5_tissues + cell_num
			nuclear_intensity_file_list = sorted(glob.glob(join(data_dir, '**', tissue, '**', '*nuclei_channel_' + intensity_type + '.csv'), recursive=True))
			nuclear_feature_matrix, _, _ = get_intensity_matrix(nuclear_intensity_file_list)
			
			nuclear_membrane_intensity_file_list = sorted(glob.glob(join(data_dir, '**',  tissue, '**', '*nucleus_boundaries_channel_' + intensity_type + '.csv'), recursive=True))
			nuclear_membrane_feature_matrix, _, _ = get_intensity_matrix(nuclear_membrane_intensity_file_list)
			
			cell_membrane_intensity_file_list = sorted(glob.glob(join(data_dir, '**', tissue, '**', '*cell_boundaries_channel_' + intensity_type + '.csv'), recursive=True))
			cell_membrane_feature_matrix, _, _ = get_intensity_matrix(cell_membrane_intensity_file_list)
			
			all_feature_matrix = cell_membrane_feature_matrix + cytoplasm_feature_matrix + nuclear_membrane_feature_matrix + nuclear_feature_matrix
			
			ss = StandardScaler()
			all_feature_matrix_z = ss.fit_transform(all_feature_matrix)
			all_feature_matrix_z_5_tissues_pieces.append(all_feature_matrix_z)
		all_feature_matrix_z = np.vstack(all_feature_matrix_z_5_tissues_pieces)
		cc = 5
		for c in range(cc, cc+46):
			print(c)
			model = KMeans(n_clusters=c, random_state=100).fit(all_feature_matrix_z)
			labels = model.labels_.astype(int)
			idx = 0
			for i in range(len(cell_file_list_5_tissues)):
				save_file_name = cell_file_list_5_tissues[i].split('-')[0] + '_individual_tissue_cluster_' + str(c) + '_label_' + intensity_type + '_5_tissues_individual_z_score.csv'
				cell_id_image = cell_id_5_tissues[idx:(idx+cell_num_5_tissues[i])]
				cell_label_image = labels[idx:(idx+cell_num_5_tissues[i])]
				cell_type_total_num = []
				for a in range(0, c):
					cell_type_total_num.append(sum(cell_label_image == a))
				save_matrix = pd.DataFrame(data=np.vstack((cell_id_image, cell_label_image)).T, columns=['ID', 'Cluster'])
				# save_matrix.to_csv(save_file_name, index=False)
				idx += cell_num_5_tissues[i]
				print(save_file_name)
				print(cell_type_total_num)