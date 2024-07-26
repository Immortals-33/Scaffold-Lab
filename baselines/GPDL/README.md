***

## GPDL

### Environment Set-up

```bash
# Download project and weights
git clone https://github.com/sirius777coder/GPDL.git
cd GPDL

# Create a virtual environment and set-up
conda create -n gpdl python=3.8
conda activate gpdl
conda install pip
conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
conda install -c conda-forge biotite
pip install "fair-esm[esmfold]"
pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'

# The full generation procedure of GPDL includes:
 # - Initial generation by GPDL-Inpainting
 # - Fixed-backbone design by ESM-IF1
 # - Backbone refined by GPDL-Hallucination
# Here the second step (ESM-IF1 sequence design) needs a separate virtual environment different from which of ESM.
conda create -n esm_if python=3.9
conda activate esm_if
conda install pytorch cudatoolkit=11.3 -c pytorch
conda install pyg -c pyg -c conda-forge
conda install pip
pip install biotite
pip install git+https://github.com/facebookresearch/esm.git
```

### Motif-scaffolding

```bash
# Make sure you're now under the root directory of GPDL

# First run 2KL8 for as an example
cd gpdl_sample
sh sub_2KL8.sh
```



> NOTE: This guidelines is WIP. We're working on updating the full generation pipeline.
>
> Reference and acknowledgements: [GPDL original repository](https://github.com/sirius777coder/GPDL)