#!/bin/bash

bash ./data_preprocessing/run_data_processing.sh  $1
bash ./modeling/run_modeling.sh $1
