# Scaffold-Lab: A Unified Framework for Evaluating Protein Backbone Generation Methods

***

Official implementation for [**_Scaffold-Lab: Critical Evaluation and Ranking of Protein Backbone Generation Methods in A Unified Framework_**](https://www.biorxiv.org/content/10.1101/2024.02.10.579743v2). 

> Note: This repository is currently **under development**. The complete code and data regards to the manuscript will be released soon. However, please feel free to raise any question or discussion either via [issue](https://github.com/Immortals-33/Scaffold-Lab/issues) or email <h2knight@sjtu.edu.cn>!



## Description

**_Scaffold-Lab_** is the first unified framework for evaluating different protein backbone generation methods.  

We present the benchmark for both unconditional generation and conditional generation in terms of *designability*, *diversity*, *novelty*, *efficiency* and *structural properties*. Currently evaluated methods are listed below:

### Unconditional Generation

- *RFdiffusion*: [Paper](https://www.nature.com/articles/s41586-023-06415-8) | [Code](https://github.com/RosettaCommons/RFdiffusion)
- *Chroma*: [Paper](https://www.nature.com/articles/s41586-023-06728-8) | [Code](https://github.com/generatebio/chroma)
- *FrameDiff*: [Paper](https://openreview.net/forum?id=m8OUBymxwv) | [Code](https://github.com/jasonkyuyim/se3_diffusion)
- *FrameFlow*: [Paper](https://arxiv.org/abs/2310.05297) | [Code](https://github.com/microsoft/protein-frame-flow)
- *Genie*: [Paper](https://arxiv.org/abs/2301.12485) | [Code](https://github.com/aqlaboratory/genie)

### Conditional Generation

- *GPDL*: [Paper](https://www.biorxiv.org/content/10.1101/2023.10.26.564121v1) | [Code](https://github.com/sirius777coder/GPDL)
- *TDS*: [Paper](https://arxiv.org/abs/2306.17775) | [Code](https://github.com/blt2114/twisted_diffusion_sampler)
- *RFdiffusion*: [Paper](https://www.nature.com/articles/s41586-023-06415-8) | [Code](https://github.com/RosettaCommons/RFdiffusion)
- *Chroma*: [Paper](https://www.nature.com/articles/s41586-023-06728-8) | [Code](https://github.com/generatebio/chroma)

***



## Table of Contents

* [Description](https://github.com/Immortals-33/Scaffold-Lab?tab=readme-ov-file#description)
* [Installation](https://github.com/Immortals-33/Scaffold-Lab?tab=readme-ov-file#installation)
* [Outline](https://github.com/Immortals-33/Scaffold-Lab?tab=readme-ov-file#outline)
* [Usage](https://github.com/Immortals-33/Scaffold-Lab?tab=readme-ov-file#usage)
  * [Unconditional Generation](https://github.com/Immortals-33/Scaffold-Lab?tab=readme-ov-file#unconditional-generation-1)
  * [Conditional Generation (Motif-scaffolding)](https://github.com/Immortals-33/Scaffold-Lab?tab=readme-ov-file#conditional-generation-motif-scaffolding)

***

## Installation

To quickly set up an environment, just simply run:

```bash
conda create -f scaffold-lab.yml
source activate scaffold-lab
# You may also need to install some dependencies manually in certain cases
pip install hydra-core --upgrade
pip install hydra-joblib-launcher --upgrade
pip install ml-collections GPUtil hjson h5py
```



***

## Outline

Here is a guide about how you can go through this repository. We aim to provide an easy-to-use evaluation pipeline as well as maximize the utility of individual scripts. Let's go through the structure of this repository as a start:

* `scaffold_lab`: This is the main directory to run different evaluations described in our paper.

* `analysis`: Scripts for calculating several metrics, including *diversity*, *novelty* and *structural properties*.  
* `baselines`: In order to generate protein backbones directly inside this repository, you may find the code of different methods baselines for [unconditional generation](https://github.com/Immortals-33/Scaffold-Lab#unconditional-generation) and [conditional generation](https://github.com/Immortals-33/Scaffold-Lab#unconditional-generation) then clone their repository under this content. it is highly recommended to run inference for different baselines **inside their own virtual environment** for potential conflicts of environmental dependencies.
  * Inside the `experiment` folder we provide scripts for performing motif-scaffolding experiments by *Chroma* using its `SubstrctureConditioner`. Refer the script for detailed information if you want.
* `config`: We place different configuration settings of [Hydra](https://github.com/facebookresearch/hydra) here to organize for evaluations. **Hydra** is a hierarchical configuration framework to help users systematize different experimental settings. Though it might be confusing when you first get in touch with it, it is a powerful tool to help you perform experiments efficiently with different combinations of parameters, for example, _the number of sequences to generate_. We recommend readers to [Docs](https://hydra.cc/docs/intro/) for advanced usage.

***



## Usage
### Unconditional Generation

Let's start by running a simple evaluation here: 

```bash
python scaffold_lab/unconditional/refolding.py 
```

 This performs a simple refolding analysis for the proteins we put inside `demo/unconditional/`.

***

### Conditional Generation (Motif-scaffolding)

To run a minimal version on motif-scaffolding task, simply run:

```bash
python scaffold_lab/motif_scaffolding/motif_refolding.py
```

This performs a evaluation on `demo/motif_scaffolding/2KL8/` where the outputs would be saved under `outputs/2KL8/`.
