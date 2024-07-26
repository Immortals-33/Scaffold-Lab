***

## FrameFlow

### Environment Set-up

```bash
# Download project and weights
git clone https://github.com/microsoft/protein-frame-flow.git
cd frame-flow
wget https://zenodo.org/records/12776473/files/weights.tar.gz?download=1
tar -xzvf weights.tar.gz

# Create a virtual environment and set-up
conda env create -f fm.yml
source activate fm
pip install torch-scatter -f https://data.pyg.org/whl/torch-2.0.0+cu117.html
pip install -e . # Install FrameFlow as a package
```

### Unconditional Generation

```bash
source activate se3
# Make sure you're now under the root directory of FrameFlow

# Change length parameters using command-line Hydra
# By default the command below will sample length in [50, 100] with 100 samples for each length
python -W ignore experiments/inference_se3_flows.py \
       inference.samples.samples_per_length=100 \
       inference.samples.length_subset=[50, 100]
       
# After finishing of the run, gather PDBs inside inference_outputs for further analysis
```



> Reference and acknowledgements: [FrameFlow original repository](https://github.com/microsoft/protein-frame-flow)