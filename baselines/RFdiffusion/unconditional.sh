#!/bin/bash
# Unconditional generation for RFdiffusion
# This would perform unconditioanl generation in range (100, 1100, 100), 100 protein backbones for each length
# You can adapt this script to generate arbitrary numbers and lengths of protein backbones
# Replace $RFdiffusion_dir of your root directory of RFdiffusion

for i in {1..10}
do

    j=$((i * 100))
    output_prefix="unconditional/length_$j/length_$j"
    $RFdiffusion_dir/scripts/run_inference.py \
    "inference.output_prefix=$output_prefix" \
    "contigmap.contigs=[$j-$j]" \
    "inference.num_designs=100"
done
