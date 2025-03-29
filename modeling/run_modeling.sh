#!/bin/bash

Rscript ./modeling/nonrandomness/read_cells_single_image_new_datasets.R $1
Rscript ./modeling/nonrandomness/run_multipp_across3_shuf_self.R $1
Rscript ./modeling/nonrandomness/run_train_shuf_across3_self.R $1
Rscript ./modeling/nonrandomness/run_fit_across3_shuf_on_other_shuf_self.R $1

Rscript ./modeling/multirange_vs_single/run_multipp_across3_single_image_self.R $1
Rscript ./modeling/multirange_vs_single/run_combine_quad_across3_self.R $1
Rscript ./modeling/multirange_vs_single/run_concat_multi_multi_all_images_across3_self.R $1
Rscript ./modeling/multirange_vs_single/run_concat_multi_multi_across3_single_image_self.R $1
Rscript ./modeling/multirange_vs_single/run_train_multi_range_all_images_across3_self.R $1
Rscript ./modeling/multirange_vs_single/run_train_single_range_all_images_across3_self.R $1

Rscript ./modeling/similarity/run_train_loo_model_across3_self.R $1
Rscript ./modeling/similarity/run_cell_type_prediction_loo_self.R $1

Rscript ./modeling/heterogeneity/image_random_tiling.R $1
Rscript ./modeling/heterogeneity/run_multipp_across3_single_image_random_tile_self.R $1
Rscript ./modeling/heterogeneity/run_concat_multi_multi_across3_single_image_random_tile_self.R $1
Rscript ./modeling/heterogeneity/run_train_concatenated_multi_quad_across3_random_tile_self.R $1

Rscript ./modeling/simulation/run_simulate_random_resampling5_self_same_pattern_global.R $1
Rscript ./modeling/simulation/run_multipp_across3_single_image_random_resample_global_self.R $1
Rscript ./modeling/simulation/run_predict_model_random_resample5_image_across3_global_self.R $1




