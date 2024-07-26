***

## FrameDiff

### Environment Set-up

```bash
# Download project and weights
git clone https://github.com/jasonkyuyim/se3_diffusion.git
cd se3_diffusion

# Create a virtual environment and set-up
conda env create -f se3.yml
source activate se3
pip install -e . # Install FrameDiff as a package
```

### Unconditional Generation

```bash
source activate se3
# Make sure you're now under the root directory of FrameDiff

# Change length parameters using command-line Hydra
# By default the command below will sample between length in range(50, 500, 50) with 100 samples for each length
python experiments/inference_se3_diffusion.py \
       inference.samples.samples_per_length=100 \
       inference.samples.min_length=50 \
       inference.samples.max_length=500 \
       inference.samples.length_step=50
       
# After finishing of the run, gather PDBs inside inference_outputs for further analysis
```



> Reference and acknowledgements: [FrameDiff original repository](https://github.com/jasonkyuyim/se3_diffusion)