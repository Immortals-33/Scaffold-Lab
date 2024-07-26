***

## TDS

### Environment Set-up

```bash
# Download project and weights
git clone https://github.com/blt2114/twisted_diffusion_sampler.git
cd twisted_diffusion_sampler/protein_exp # Root directory of protein experiments

# Create a virtual environment and set-up
conda env create -f se3_tds.yml
source activate se3_tds
pip install -e . # Install TDS as a package
pip install numba
```

### Motif-scaffolding

```bash
source activate se3_tds
# Make sure you're now under the root directory of protein_exp inside TDS

# First run 2KL8 for as an example
sh motif_tds_2KL8.sh
# Then run the benchmark in batch
       
# Since TDS has its own contig specifications and PDB naming format, so we recommend users to gather and concatenate the 'sc_results.csv' inside each case for further analysis,
# which is mostly compatible with the evaluation outputs of Scaffold-Lab
```



> Reference and acknowledgements: [TDS original repository](https://github.com/blt2114/twisted_diffusion_sampler)