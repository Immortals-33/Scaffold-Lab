***

## RFdiffusion

### Environment Set-up

```bash
# Download project and weights
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion
mkdir models && cd models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
cd - # Back to root directory

# Set-up environment
conda env create -f env/SE3nv.yml

conda activate SE3nv
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../.. # change into the root directory of the repository
pip install -e . # install the rfdiffusion module from the root of the repository
```

### Unconditional Generation

```bash
source activate SE3nv
sh unconditional.sh
```

### Motif-scaffolding

```bash
# Generate 
source activate SE3nv
sh motif_scaffolding.sh
```

> We recommend users to go through the example scripts before running experiments. When running motif-scaffolding benchmark, it is strongly recommended to be familiar with the [grammar of contig](https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#motif-scaffolding).
>
> Reference and acknowledgements: [RFdiffusion original repository](https://github.com/RosettaCommons/RFdiffusion)
