import os
import tree
import time
import numpy as np
import hydra
import torch
import subprocess
import re
import logging
import pandas as pd
import sys
import shutil
import GPUtil
from pathlib import Path
from typing import *
from omegaconf import DictConfig, OmegaConf

import esm
from biotite.sequence.io import fasta

import rootutils
path = rootutils.find_root(search_from='./', indicator=[".git", "setup.cfg"])
rootutils.set_root(
    path=path, # path to the root directory
    project_root_env_var=True, # set the PROJECT_ROOT environment variable to root directory
    dotenv=True, # load environment variables from .env if exists in root directory
    pythonpath=True, # add root directory to the PYTHONPATH (helps with imports)
    cwd=True, # change current working directory to the root directory (helps with filepaths)
)

from analysis import utils as au
from data import structure_utils as su


class Refolder:
    
    def __init__(
        self,
        conf:DictConfig,
        conf_overrides: Dict=None
        ):
        
        self._log = logging.getLogger(__name__)
        
        OmegaConf.set_struct(conf, False)
        
        self._conf = conf
        self._infer_conf = conf.inference
        self._sample_conf = self._infer_conf.samples
        
        self._rng = np.random.default_rng(self._infer_conf.seed)
        
        # Set-up accelerator
        if torch.cuda.is_available():
            if self._infer_conf.gpu_id is None:
                available_gpus = ''.join(
                    [str(x) for x in GPUtil.getAvailable(
                        order="memory", limit = 8)]
                )
                self.device = f'cuda:{available_gpus[0]}'
            else:
                self.device = f'cuda:{self._infer_conf.gpu_id}'
        else:
            self.device = 'cpu'
        self._log.info(f'Using device: {self.device}')
        
        
        # Set-up directories
        #self._weights_path = self._infer_conf.weights_path
        output_dir = self._infer_conf.output_dir

        self._output_dir = output_dir
        os.makedirs(self._output_dir, exist_ok=True)
        self._pmpnn_dir = self._infer_conf.pmpnn_dir
        self._sample_dir = self._infer_conf.backbone_pdb_dir
        self._CA_only = self._infer_conf.CA_only
        
        # Configs for motif-scaffolding
        self._motif_csv = self._infer_conf.motif_csv_dir
        self._input_pdbs_dir = self._infer_conf.input_pdbs_dir
        
        #Save config
        config_folder = os.path.basename(Path(self._output_dir))
        config_path = os.path.join(self._output_dir, f"{config_folder}.yaml")
        with open(config_path, 'w') as f:
            OmegaConf.save(config=self._conf, f=f)
        self._log.info(f'Saving self-consistency config to {config_path}')
        
        # Load models and experiment
        if 'cuda' in self.device:
            self._folding_model = esm.pretrained.esmfold_v1().eval()
        elif self.device == 'cpu': # ESMFold is not supported for half-precision model when running on CPU
            self._folding_model = esm.pretrained.esmfold_v1().float().eval()
        self._folding_model = self._folding_model.to(self.device)
        

    def get_csv_data(
        self,
        csv_info: pd.DataFrame,
        pdb_name: str,
        sample_num: Union[str, int]
    ):
        csv_info = pd.read_csv(csv_info)
        csv_info['sample_num'] = csv_info['sample_num'].astype(int)
        sample_item = csv_info[(csv_info['pdb_name'] == pdb_name) & (csv_info['sample_num'] == int(sample_num))]
        if not sample_item.empty:
            return(
                sample_item['contig'].iloc[0],
                sample_item['mask'].iloc[0],
                sample_item['motif_indices'].iloc[0],
            )
        
    def motif_indices_to_contig(self, motif_indices: str):
        if motif_indices.startswith('[') and motif_indices.endswith(']'):
            motif_indices = motif_indices.strip('[]').split(', ')
            try:
                motif_indices = [int(index) for index in motif_indices]
            except ValueError as e:
                raise ValueError(f"Error converting motif_indices_str to list of integers: {e}")

            sorted_indices = sorted(motif_indices)
            contig = ""
            range_start = None

            for i, index in enumerate(sorted_indices):
                if range_start is None:
                    range_start = index
                if i == len(sorted_indices) - 1 or sorted_indices[i + 1] != index + 1:
                    if contig:
                        contig += "/"
                    if range_start == index:
                        contig += f"A{range_start}"
                    else:
                        contig += f"A{range_start}-{index}"
                    range_start = None
            return contig
        else:
            raise ValueError(f"Invalid input for motif_indices_to_contig: {motif_indices}")

    def motif_indices_to_fixed_positions(self, motif_indices, chain='A'):
        # Converts motif indices to the fixed positions string format
        motif_indices = motif_indices.strip('[]').split(', ')
        motif_indices = sorted([int(index) for index in motif_indices])
        fixed_positions = ' '.join(str(idx) for idx in motif_indices)
        return f"{fixed_positions}"
    
    def run_sampling(self):
        # Run ProteinMPNN

        for pdb_file in os.listdir(self._sample_dir):
            if ".pdb" in pdb_file:
                backbone_name = os.path.splitext(pdb_file)[0]
                sample_num = backbone_name.split("_")[-1]
                parts = backbone_name.split('_')
                backbone_name = parts[0] if len(parts) == 2 else '_'.join(parts[:-1])                
                contig, mask, motif_indices = self.get_csv_data(self._motif_csv, backbone_name, sample_num)
                
                # Deal with contig
                if '6VW1' not in pdb_file:
                    reference_contig = '/'.join(re.findall(r'[A-Za-z]+\d+-\d+', contig)) 
                design_contig = self.motif_indices_to_contig(motif_indices)
                
                # Handle complex case for PDB 6VW1
                if backbone_name == '6VW1':
                    reference_contig = "A24-42/A64-82"
                    parts_6VW1 = design_contig.split("/")
                    design_contig = '/'.join(parts_6VW1[:-1])
                    chain_B = parts_6VW1[-1]
                    start, end = map(int, chain_B[1:].split("-"))
                    chain_B_indices = list(range(start, end + 1))
                
                if '_' in backbone_name: # Handle length-variable design for different PDB cases
                    reference_pdb = os.path.join(self._input_pdbs_dir, f'{backbone_name.split("_")[0]}.pdb')
                else:
                    reference_pdb = os.path.join(self._input_pdbs_dir, f'{backbone_name}.pdb')
                design_pdb = os.path.join(self._sample_dir, pdb_file)
                
                # Extract motif and calculate motif-RMSD
                reference_motif = au.motif_extract(reference_contig, reference_pdb, atom_part="CA")
                design_motif = au.motif_extract(design_contig, design_pdb, atom_part="CA")
                rms = au.rmsd(reference_motif, design_motif)
                
                # Save outputs
                basename_dir = os.path.basename(os.path.normpath(self._sample_dir))
                backbone_dir = os.path.join(self._output_dir, basename_dir, f'{backbone_name}_{sample_num}')
                os.makedirs(backbone_dir, exist_ok=True)
                self._log.info(f'Running self-consistency on {backbone_name}')
                shutil.copy2(os.path.join(self._sample_dir, pdb_file), backbone_dir)
                print(f'copied {pdb_file} to {backbone_dir}')
                
                #seperate_pdb_folder = os.path.join(backbone_dir, backbone_name)
                pdb_path = os.path.join(backbone_dir, pdb_file)
                sc_output_dir = os.path.join(backbone_dir, 'self_consistency')
                os.makedirs(sc_output_dir, exist_ok=True)
                shutil.copy(pdb_path, os.path.join(
                    sc_output_dir, os.path.basename(pdb_path)))
                
                if backbone_name == '6VW1':
                    _ = self.run_self_consistency(
                    sc_output_dir,
                    pdb_path,
                    motif_mask=np.array(eval(mask)),
                    motif_indices=motif_indices,
                    rms=rms,
                    complex_motif=chain_B_indices
                )
                else:
                    _ = self.run_self_consistency(
                        sc_output_dir,
                        pdb_path,
                        motif_mask=np.array(eval(mask)),
                        motif_indices=motif_indices,
                        rms=rms
                    )
                self._log.info(f'Done sample: {pdb_path}')
    
    def run_self_consistency(
            self,
            decoy_pdb_dir: str,
            reference_pdb_path: str,
            motif_mask: Optional[np.ndarray]=None,
            motif_indices: Optional[List]=None,
            rms: Optional[float]=None,
            complex_motif: Optional[List]=None
            ):
        """Run self-consistency on design proteins against reference protein.
        
        Args:
            decoy_pdb_dir: directory where designed protein files are stored.
            reference_pdb_path: path to reference protein file
            motif_mask: Optional mask of which residues are the motif.

        Returns:
            Writes ProteinMPNN outputs to decoy_pdb_dir/seqs
            Writes ESMFold outputs to decoy_pdb_dir/esmf
            Writes results in decoy_pdb_dir/sc_results.csv
        """

        # Run ProteinMPNN
        output_path = os.path.join(decoy_pdb_dir, "parsed_pdbs.jsonl")
        process = subprocess.Popen([
            'python',
            f'{self._pmpnn_dir}/helper_scripts/parse_multiple_chains.py',
            f'--input_path={decoy_pdb_dir}',
            f'--output_path={output_path}',
        ])
        _ = process.wait()
        num_tries = 0
        ret = -1
        pmpnn_args = [
            'python',
            f'{self._pmpnn_dir}/protein_mpnn_run.py',
            '--out_folder',
            decoy_pdb_dir,
            '--jsonl_path',
            output_path,
            '--num_seq_per_target',
            str(self._sample_conf.seq_per_sample),
            '--sampling_temp',
            '0.1',
            '--seed',
            '33',
            '--batch_size',
            '10',
        ]
        if self._infer_conf.gpu_id is not None:
            pmpnn_args.append('--device')
            pmpnn_args.append(str(self._infer_conf.gpu_id))
        if self._CA_only == True:
            pmpnn_args.append('--ca_only')
        
        # Fix desired motifs    
        if motif_indices is not None:
            fixed_positions = self.motif_indices_to_fixed_positions(motif_indices)
            chains_to_design = "A"
            # This is particularlly for 6VW1
            if complex_motif is not None:
                motif_indices = motif_indices.strip('[]').split(', ')
                motif_indices = sorted([int(index) for index in motif_indices])
                motif_indices = [element for element in motif_indices if element not in complex_motif]
                complex_motif = " ".join(map(str, complex_motif)) # List2str
                fixed_positions = " ".join(map(str, motif_indices)) # List2str
                # fixed_positions = self.motif_indices_to_fixed_positions(motif_indices)
                #complex_motif = self.motif_indices_to_fixed_positions(complex_motif)
                print(motif_indices)
                print(fixed_positions)
                fixed_positions = fixed_positions + ", " + complex_motif
                print(fixed_positions)
                chains_to_design = "A B"
            path_for_fixed_positions = os.path.join(decoy_pdb_dir, "fixed_pdbs.jsonl")

            
            subprocess.call([
                'python',
                os.path.join(self._pmpnn_dir, 'helper_scripts/make_fixed_positions_dict.py'),
                '--input_path', output_path,
                '--output_path', path_for_fixed_positions,
                '--chain_list', chains_to_design,
                '--position_list', fixed_positions
            ])
            
            pmpnn_args.extend([
                '--chain_id_jsonl', os.path.join(decoy_pdb_dir, "assigned_pdbs.jsonl"),
                '--fixed_positions_jsonl', path_for_fixed_positions
            ])
            
        while ret < 0:
            try:
                process = subprocess.Popen(
                    pmpnn_args,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT
                )
                ret = process.wait()
            except Exception as e:
                num_tries += 1
                self._log.info(f'Failed ProteinMPNN. Attempt {num_tries}/5')
                torch.cuda.empty_cache()
                if num_tries > 4:
                    raise e
        mpnn_fasta_path = os.path.join(
            decoy_pdb_dir,
            'seqs',
            os.path.basename(reference_pdb_path).replace('.pdb', '.fa')
        )

        # Run ESMFold on each ProteinMPNN sequence and calculate metrics.
        mpnn_results = {
            'tm_score': [],
            'sample_path': [],
            'header': [],
            'sequence': [],
            'rmsd': [],
            'pae': [],
            'ptm': [],
            'plddt': [],
            'length': [],
            'backbone_motif_rmsd': []
        }
        if motif_mask is not None:
            # Only calculate motif RMSD if mask is specified.
            mpnn_results['motif_rmsd'] = []
        esmf_dir = os.path.join(decoy_pdb_dir, 'esmf')
        os.makedirs(esmf_dir, exist_ok=True)
        fasta_seqs = fasta.FastaFile.read(mpnn_fasta_path)
        sample_feats = su.parse_pdb_feats('sample', reference_pdb_path)
        
        for i, (header, string) in enumerate(fasta_seqs.items()):

            # Run ESMFold
            self._log.info(f'Running ESMfold......')
            esmf_sample_path = os.path.join(esmf_dir, f'sample_{i}.pdb')
            _, full_output = self.run_folding(string, esmf_sample_path)
            esmf_feats = su.parse_pdb_feats('folded_sample', esmf_sample_path)
            sample_seq = su.aatype_to_seq(sample_feats['aatype'])

            # Calculate scTM of ESMFold outputs with reference protein
            _, tm_score = su.calc_tm_score(
                sample_feats['bb_positions'], esmf_feats['bb_positions'],
                sample_seq, sample_seq)
            rmsd = su.calc_aligned_rmsd(
                sample_feats['bb_positions'], esmf_feats['bb_positions'])
            pae = torch.mean(full_output['predicted_aligned_error']).item()
            ptm = full_output['ptm'].item()
            plddt = full_output['mean_plddt'].item()
            if motif_mask is not None:
                sample_motif = sample_feats['bb_positions'][motif_mask]
                of_motif = esmf_feats['bb_positions'][motif_mask]
                motif_rmsd = su.calc_aligned_rmsd(
                    sample_motif, of_motif)
                mpnn_results['motif_rmsd'].append(f'{motif_rmsd:.3f}')
            if rms is not None:
                mpnn_results['backbone_motif_rmsd'].append(f'{rms:.3f}')
            mpnn_results['rmsd'].append(f'{rmsd:.3f}')
            mpnn_results['tm_score'].append(f'{tm_score:.3f}')
            mpnn_results['sample_path'].append(esmf_sample_path)
            mpnn_results['header'].append(header)
            mpnn_results['sequence'].append(string)
            mpnn_results['pae'].append(f'{pae:.3f}')
            mpnn_results['ptm'].append(f'{ptm:.3f}')
            mpnn_results['plddt'].append(f'{plddt:.3f}')
            mpnn_results['length'].append(len(string))
            #mpnn_results['mpnn_score'].append(f'{mpnn_score:.3f}')

        # Save results to CSV
        csv_path = os.path.join(decoy_pdb_dir, 'sc_results.csv')
        mpnn_results = pd.DataFrame(mpnn_results)
        mpnn_results.to_csv(csv_path)

    def run_folding(self, sequence, save_path):
        """
        Run ESMFold on sequence.
        TBD: Add options for OmegaFold and AlphaFold2.
        """
        with torch.no_grad():
            output = self._folding_model.infer(sequence)
            output_dict = {key: value.cpu() for key, value in output.items()}
            output = self._folding_model.output_to_pdb(output)
        with open(save_path, "w") as f:
            f.write(output[0])
        return output, output_dict  
    
@hydra.main(version_base=None, config_path="../../config", config_name="motif_scaffolding.yaml")
def run(conf: DictConfig) -> None:
    
    print('Starting refolding for motif-scaffolding task......')
    start_time = time.time()
    refolder = Refolder(conf)
    refolder.run_sampling()
    elapsed_time = time.time() - start_time
    print(f"Finished in {elapsed_time:.2f}s. Voila!")
    
if __name__ == '__main__':
    run()
