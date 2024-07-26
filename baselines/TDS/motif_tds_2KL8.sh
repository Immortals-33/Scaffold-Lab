#!/bin/bash

source activate se3_tds

python ./experiments/inference_particle_filter.py \
       inference.motif_scaffolding.test_name=2KL8
