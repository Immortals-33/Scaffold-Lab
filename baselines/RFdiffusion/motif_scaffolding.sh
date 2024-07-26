#!/bin/bash
# This would run the 25-cases motif-scaffolding benchmark for RFdiffusion
# The contig are randomly sampled based on the contig provided by the RFdiffusion paper,
  # i.e. The length of scaffolds of ome cases would vary with different samples
# We generate 100 samples for each case
# Replace $RFdiffusion_dir of your root directory of RFdiffusion

# Reference: Supplementary Methods Table 9 of https://www.nature.com/articles/s41586-023-06415-8

# 1PRW
echo "Start designing 1PRW..."

output_prefix="motif_scaffolding/1PRW/1PRW"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/1PRW.pdb \
"contigmap.contigs=[5-20/A16-35/10-25/A52-71/5-20]" \
inference.num_designs=100
echo "1PRW is done!"
    
# 1BCF    
echo "Start designing 1BCF..."

output_prefix="motif_scaffolding/1BCF/1BCF"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/1BCF.pdb \
"contigmap.contigs=[8-15/A92-99/16-30/A123-130/16-30/A47-54/16-30/A18-25/8-15]" \
inference.num_designs=100

echo "1BCF is done!"

# 5TPN
echo "Start designing 5TPN..."

output_prefix="motif_scaffolding/5TPN/5TPN"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/5TPN.pdb \
"contigmap.contigs=[10-40/A163-181/10-40]" \
inference.num_designs=100

echo "5TPN is done!"

# 5IUS
echo "Start designing 5IUS..."

output_prefix="motif_scaffolding/5IUS/5IUS"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/5IUS.pdb \
"contigmap.contigs=[0-30/A119-140/15-40/A63-82/0-30]" \
inference.num_designs=100

echo "5IUS is done!"

# 3IXT
echo "Start designing 3IXT..."
output_prefix="motif_scaffolding/3IXT/3IXT"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/3IXT.pdb \
"contigmap.contigs=[10-40/P254-277/10-40]" \
inference.num_designs=100

echo "3IXT is done!"

# 5YUI
echo "Start designing 5YUI..."

output_prefix="motif_scaffolding/5YUI/5YUI"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/5YUI.pdb \
"contigmap.contigs=[5-30/A93-97/5-20/A118-120/10-35/A198-200/10-30]" \
inference.num_designs=100

echo "5YUI is done!"
    
# 1QJG (Note that the motif map of single residue motif should be add the same number like "A38-38", otherwise it'll cause error
echo "Start designing 1QJG..."

output_prefix="motif_scaffolding/1QJG/1QJG"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/1QJG.pdb \
"contigmap.contigs=[10-20/A38-38/15-30/A14-14/15-30/A99-99/10-20]" \
inference.num_designs=100

echo "1QJG is done!"

# 1YCR
echo "Start designing 1YCR..."

output_prefix="motif_scaffolding/1YCR/1YCR"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/1YCR.pdb \
"contigmap.contigs=[10-40/B19-27/10-40]" \
inference.num_designs=100

echo "1YCR is done!"

# 2KL8
echo "Start designing 2KL8..."

output_prefix="motif_scaffolding/2KL8/2KL8"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/2KL8.pdb \
"contigmap.contigs=[A1-7/20/A28-79]" \
inference.num_designs=100

echo "2KL8 is done!"

# 4JHW
echo "Start designing 4JHW..."

output_prefix="motif_scaffolding/4JHW/4JHW"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/4JHW.pdb \
"contigmap.contigs=[10-25/F196-212/15-30/F63-69/10-25]" \
inference.num_designs=100

echo "4JHW is done!"

# 4ZYP (Note that we specify the total length to be 30-50!)
echo "Start designing 4ZYP..."
output_prefix="motif_scaffolding/4ZYP/4ZYP"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/4ZYP.pdb \
"contigmap.contigs=[10-40/A422-436/10-40]" \
"contigmap.length=30-50" \
inference.num_designs=100
echo "1PRW is done!"
    
# 5WN9 (Note that we specify the total length to be 35-50!)
echo "Start designing 5WN9..."

output_prefix="motif_scaffolding/5WN9/5WN9"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/5WN9.pdb \
"contigmap.contigs=[10-40/A170-189/10-40]" \
"contigmap.length=35-50" \
inference.num_designs=100

echo "5WN9 is done!"

# 6VW1 (This is not included in the benchmark of Scaffold-Lab, though we left it here for a reference)
echo "Start designing 6VW1..."

output_prefix="motif_scaffolding/6VW1/6VW1"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/6VW1.pdb \
"contigmap.contigs=[E400-510/0 20-30/A24-42/4-10/A64-82/0-5]" \
"contigmap.length=62-83" \
inference.num_designs=100

echo "6VW1 is done!"

# 7MRX_60
echo "Start designing 7MRX_60..."

output_prefix="motif_scaffolding/7MRX_60/7MRX_60"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/7MRX.pdb \
"contigmap.contigs=[0-38/B25-46/0-38]" \
"contigmap.length=60-60" \
inference.num_designs=100
echo "7MRX_60 is done!"
    
# 7MRX_85
echo "Start designing 7MRX_85..."

output_prefix="motif_scaffolding/7MRX_85/7MRX_85"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/7MRX.pdb \
"contigmap.contigs=[0-63/B25-46/0-63]" \
"contigmap.length=85-85" \
inference.num_designs=100

echo "7MRX_85 is done!"

# 7MRX_128
echo "Start designing 7MRX_128..."

output_prefix="motif_scaffolding/7MRX_128/7MRX_128"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/7MRX.pdb \
"contigmap.contigs=[0-122/B25-46/0-122]" \
"contigmap.length=128-128" \
inference.num_designs=100

echo "7MRX_128 is done!"

# 5TRV_short
echo "Start designing 5TRV..."

output_prefix="motif_scaffolding/5TRV_short/5TRV_short"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/5TRV.pdb \
"contigmap.contigs=[0-35/A45-65/0-35]" \
"contigmap.length=56-56" \
inference.num_designs=100

echo "5TRV_short is done!"

# 5TRV_medium
echo "Start designing 5TRV_medium..."

output_prefix="motif_scaffolding/5TRV_medium/5TRV_medium"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/5TRV.pdb \
"contigmap.contigs=[0-65/A45-65/0-65]" \
"contigmap.length=86-86" \
inference.num_designs=100

echo "5TRV_short is done!"

# 5TRV_long
echo "Start designing 5TRV_long..."

output_prefix="motif_scaffolding/5TRV_long/5TRV_long"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/5TRV.pdb \
"contigmap.contigs=[0-95/A45-65/0-95]" \
"contigmap.length=116-116" \
inference.num_designs=100

echo "5TRV_long is done!"

# 6E6R_short
echo "Start designing 6E6R_short..."
output_prefix="motif_scaffolding/6E6R_short/6E6R_short"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/6E6R.pdb \
"contigmap.contigs=[0-35/A23-35/0-35]" \
"contigmap.length=48-48" \
inference.num_designs=100
echo "6E6R_short is done!"
    
# 6E6R_medium
echo "Start designing 6E6R_medium..."

output_prefix="motif_scaffolding/6E6R_medium/6E6R_medium"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/6E6R.pdb \
"contigmap.contigs=[0-65/A23-35/0-65]" \
"contigmap.length=78-78" \
inference.num_designs=100

echo "6E6R_medium is done!"

# 6E6R_long
echo "Start designing 6E6R_long..."

output_prefix="motif_scaffolding/6E6R_long/6E6R_long"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/6E6R.pdb \
"contigmap.contigs=[0-95/A23-35/0-95]" \
"contigmap.length=108-108" \
inference.num_designs=100

echo "6E6R_long is done!"

# 6EXZ_short
echo "Start designing 6EXZ_short..."

output_prefix="motif_scaffolding/6EXZ_short/6EXZ_short"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/6EXZ.pdb \
"contigmap.contigs=[0-35/A557-571/0-35]" \
"contigmap.length=50-50" \
inference.num_designs=100

echo "6EXZ_short is done!"

# 6EXZ_medium
echo "Start designing 6EXZ_medium..."

output_prefix="motif_scaffolding/6EXZ_medium/6EXZ_medium"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/6EXZ.pdb \
"contigmap.contigs=[0-65/A557-571/0-65]" \
"contigmap.length=80-80" \
inference.num_designs=100

echo "6EXZ_medium is done!"

# 6EXZ_long
echo "Start designing 6EXZ_long..."

output_prefix="motif_scaffolding/6EXZ_long/6EXZ_long"
$RFdiffusion_dir/scripts/run_inference.py \
inference.output_prefix=$output_prefix \
inference.input_pdb=../input_pdbs/6EXZ.pdb \
"contigmap.contigs=[0-95/A557-571/0-95]" \
"contigmap.length=110-110" \
inference.num_designs=100

echo "6EXZ_long is done!"