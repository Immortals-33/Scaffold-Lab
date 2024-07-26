#!/bin/bash

elements=("1PRW" "1BCF" "5TPN" "5IUS" "3IXT" "5YUI" "1QJG" "1YCR" \
         "7MRX_60" "7MRX_85" "7MRX_128" "4JHW" "4ZYP" "5WN9" \
         "6VW1" "5TRV_short" "5TRV_medium" "5TRV_long" "6E6R_short" \
         "6E6R_medium" "6E6R_long" "6EXZ_short" "6EXZ_medium" "6EXZ_long")

for element in "${elements[@]}"; do
    sed "s/2KL8/$element/g" motif_tds_2KL8.sh > "motif_tds_$element.sh"
    sh "motif_tds_$element.sh"
done
