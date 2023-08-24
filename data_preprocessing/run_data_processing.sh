#!/bin/bash

python3 ./data_preprocessing/downloading/CODEX_PPM_download_new.py $1
python3 ./data_preprocessing/clustering/clustering_5_tissues_individual_z_score.py $1
python3 ./data_preprocessing/clustering/cellar_annotation_preprocessing.py $1
python3 ./data_preprocessing/clustering/cell_type_assigning.py $1
python3 ./data_preprocessing/clustering/clusters_combining.py $1
