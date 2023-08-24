import sys

import numpy as np
import glob
import pandas as pd
import os
from os.path import join
import sys
if __name__ == '__main__':
	data_dir = sys.argv[1]
	cluster_num = '39'
	cluster_dir = f'{data_dir}/Florida/LN/800f1703d81373c58ca5ca0b76e52d79/' + str(cluster_num) + '_clusters_test' + '/'

	label_matching_checklist = np.load(f'{cluster_dir}/total_cell_intensity_5_tissues_individual_z_score_{cluster_num}_clusters_to_5_clusters_labels.npy')
	for tissue in ['LN','SPLEEN','THYMUS','LI','SI']:
		label_dir_list = sorted(glob.glob(f'{data_dir}/**/{tissue}/**/reg*_stitched_expressions.ome.tiff_individual_tissue_cluster_{cluster_num}_label_total_5_tissues_individual_z_score.csv', recursive=True))
		for label_dir in label_dir_list:
			region_idx = os.path.basename(label_dir).split('_')[0][-1]
			new_csv_dir = join(os.path.dirname(label_dir), f'reg{region_idx}_stitched_expressions.ome.tiff_individual_tissue_cluster_5_label_total_across3.csv')
			original_labels = pd.read_csv(label_dir)
			for c in range(label_matching_checklist.shape[0]):
				original_labels['Cluster'].replace([label_matching_checklist[c, 0]], -label_matching_checklist[c, 1], inplace=True)
			original_labels['Cluster'] = -original_labels['Cluster']
			original_labels['Cluster'] = original_labels['Cluster'].where(original_labels['Cluster'] >= 0, 1)
			print(new_csv_dir)
			pd.DataFrame.to_csv(original_labels, new_csv_dir, index=False)