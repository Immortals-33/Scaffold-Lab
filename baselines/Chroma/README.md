***

## Chroma

### Environment Set-up

```bash
# Download project and weights
git clone https://github.com/generatebio/chroma.git
cd chroma

# Set-up environment
conda env create -n chroma python==3.9 # The Python version here is not tightly restricted

source activate chroma
pip install -e chroma # Locally installed, recommended
pip install generate-chroma # Alternatively, directly install Chroma as a pip package

# You need to get an API key from offical Chroma before starting benchmarking
```

### Unconditional Generation

```bash
source activate chroma
# Make sure you're now under the root directory of Chroma
# By default, this would generate protein backbones with length (50 + range(100, 1100, 100)) and 100 samples per length
python unconditional.py
```

### Motif-scaffolding

*Chroma* does not provide an explicit example for motif-scaffolding. We utilized its `SubStructureConditioner` to perform the motif-scaffolding task to test its performance. As the protein system of *Chroma* is different from others (also because some disconnected positions in some PDB files) we did not use a general indexing system to perform motif-scaffolding tasks directly. Instead, we curated the benchmark scripts by manually indexing the motif positions for each protein.  

The overall framework of motif-scaffolding in Chroma is stored as a Python class `MotifScaffolding()` inside `motif_scaffolding.py`. Inside the `chroma_benchmark` folder are 24 python files for 24 cases (**6VW1** is a protein complex which could only be handled by a subset of methods, therefore excluded from our benchmarking procedure), where the contigs are manually curated to suite the length range of the original RFdiffusion benchmark. 

By default, we generate length-variable proteins for each case. (i.e. The length of protein may vary from different samples within a pre-specified range within certain cases.)

```bash
source activate chroma
# Take 2KL8 as an example
python benchmark/2KL8.py -n 100 -p "2KL8" # Here -n ~ the number of samples you want, -p ~ the name of protein as prefix
```

> We recommend users to go through the example scripts before running experiments. When running motif-scaffolding benchmark, it is strongly recommended to be familiar with the [grammar of contig](https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#motif-scaffolding).
>
> Reference and acknowledgements: [Chroma original repository](https://github.com/generatebio/chroma)