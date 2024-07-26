***

## Genie

### Environment Set-up

```bash
# Download project and weights
git clone https://github.com/aqlaboratory/genie.git
cd genie
pip install -e . # Install Genie as a package

# Create a virtual environment and set-up
conda env create -n genie python==3.9 # The Python version here is not tightly restricted
source activate genie
pip install -e . # Install Genie as a package
```

### Unconditional Generation

```bash
source activate genie
# Make sure you're now under the root directory of Genie
# We adapted the original sample script of Genie to "adapted_sample.py", making it capable of sampling dynamically by length. 
cp adapted_sample.py genie/ # Copy the adapted script into the directory where original sampling script was stored 

# Change `min_length`, `max_length` and `interval` inside "adapted_sample.py" to sample the lengths you want
# By default the command below will sample between length in range(50, 300, 50) with 100 samples for each length
python genie/adapted_sample.py \
       -n scope_l_256 \
       -r weights \
       --batch_size 5 \
       --num_batches 20 # The number of generation per length is (batch_size * num_batches). Following this rule, you can adapt these parameters based on your own hardware resources.
       
# After finishing of the run, the coordinate files will be stored under $Genie_root_directory/opt/
# We now use "convert_coords_to_pdb.py" to convert coordinate files into PDB files and rename them
python convert_coords_to_pdb.py
```



> Reference and acknowledgements: [Genie original repository](https://github.com/aqlaboratory/genie)