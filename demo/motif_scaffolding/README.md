***

## Provide Detailed Information for Motif-scaffolding Task

When running refolding on motif-scaffolding task under this demo, we use `motif_info.csv` to place information, whose format follows some certain s. The following is a brief guideline:

* `pdb_name`: The PDB ID of native protein. It is used to extract the native motif of PDBs under `reference_pdbs/`. 

* `sample_num`: We recommend rename the proteins to be evaluated as `${pdb_name}_${sample_number}.pdb`. e.g. In this demo, `2KL8_2.pdb`, `2KL8_33.pdb`. This is particularly important when sample length-variable scaffolds for a single case, for different samples correspond to different motif regions.

* `contig`: The information on motifs and scaffolds of each designed protein. The contig grammar is adapted from [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#motif-scaffolding) but with a little bit modification. To be short:

  * Consider each motif and each scaffold is a **part**. Different parts are split by `/`.
  * The parts don't starts with letters correspond to **scaffold**. Different from the grammar in _RFdiffusion_, we use two same numbers and concatenate them using `-`. e.g. `20-20` represents a 20-residue scaffold.
  * The parts starts with letters correspond to **motif**, where the letters correspond to the chain harboring the native motif. e.g. In this demo, `A1-7` represents the motif as residue 1-7 of chain A in _2KL8_.

  Together is a complete contig to indicate the positions of motifs and scaffolds in both the designed and native proteins. e.g. `A1-7/20-20/A28-79` means the sample first place residue 1-7 of chain A in _2KL8_, then connected by a 20-residue scaffold and finally with another motifs of residue 28-79 of chain A in _2KL8_. 

* `motif_indices`: A str-like list that represents the indexes of motif part. This is used to extract the designed motif and calculated against the native-one. Note: This is a **1-based** index system.

* `mask`: Similar to `motif-indices`, this is a `[N,1]` boolean list that represent the motif information, where `False` means scaffold and `True` means motif. 

Altogether, the main script runs as the following logic:

* Read an input folder with files named like `${pdb_name}_${sample_num}.pdb`;
* Read the csv file to iterate each protein and extract necessary information;
* Running refolding and calculate metrics using these information.



> This is an initial version and we're working on simplify the grammar to make it more user-friendly. Bug fixes and PR are also welcome!