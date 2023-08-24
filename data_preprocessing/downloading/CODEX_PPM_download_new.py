import numpy as np
import os
from os.path import join
import sys
from skimage.io import imread

def data_download(id, tissue, TMC):
	img_dir = join(data_dir, TMC, tissue, id)
	if not os.path.exists(img_dir):
		os.makedirs(img_dir)
	BASE_URL = 'https://g-d00e7b.09193a.5898.dn.glob.us/'
	
	def download_files(prefix, img_dir, id):
		file_names = [
			'cell_centers.csv',
			'img_binary.png',
			'nuclei_channel_total.csv',
			'cell_channel_total.csv',
			'cell_boundaries_channel_total.csv',
			'nucleus_boundaries_channel_total.csv'
		]
		
		for file_name in file_names:
			if file_name != 'img_binary.png':
				output_file = join(img_dir, f'reg1_stitched_expressions.ome.tiff-{file_name}')
				download_url = f'{BASE_URL}{id}/sprm_outputs/{prefix}-{file_name}?download=1'
			else:
				output_file = join(img_dir, f'reg1_stitched_expressions.ome.tiff_{file_name}')
				download_url = f'{BASE_URL}{id}/sprm_outputs/{prefix}_{file_name}?download=1'
			os.system(f'wget -O {output_file} {download_url}')
	
	if TMC == 'Florida':
		download_files('reg1_stitched_expressions.ome.tiff', img_dir, id)
		
		# Check if the binary image is empty and if so, re-download with a different prefix
		if os.stat(join(img_dir, 'reg1_stitched_expressions.ome.tiff_img_binary.png')).st_size == 0:
			download_files('reg001_expr.ome.tiff', img_dir, id)
	else:
		download_files('reg001_expr.ome.tiff', img_dir, id)
	
	image = imread(join(img_dir,'reg1_stitched_expressions.ome.tiff_img_binary.png'))
	np.savetxt(join(img_dir,'reg1_stitched_expressions.ome.tiff-img_shape.txt'), image.shape, fmt='%i')
	
if __name__ == '__main__':
	script_dir = os.path.dirname(os.path.abspath(__file__))
	tissue_list = ['LN', 'LI', 'SI', 'SPLEEN', 'THYMUS']
	data_dir = sys.argv[1]
	for tissue in tissue_list:
		if tissue == 'LI' or tissue == 'SI':
			TMC = 'Stanford'
		else:
			TMC = 'Florida'
		checklist = np.loadtxt(join(script_dir, 'checklist_CODEX_PPM_' + tissue + '_new.txt'), dtype=str)
		img_num = checklist.shape[0]
		for idx in range(img_num):
			data_download(checklist[idx], tissue, TMC)

