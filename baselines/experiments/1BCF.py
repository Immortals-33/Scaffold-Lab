import os 
import logging
import typing as T
import csv
import pandas as pd
import argparse
import random
from typing import *

import torch

from chroma import Chroma, Protein, conditioners
from benchmark.scripts import motif_scaffolding


def create_parser():
    parser = argparse.ArgumentParser(description="Motif-scaffolding with Chroma")
    parser.add_argument(
        "-n",
        "--num",
        type=int,
        help="Numbers of proteins to sample"
    )
    parser.add_argument(
        "-p",
        "--pdb",
        type=str,
        help="Name of PDB"
    )
    return parser

def system_setup(
    seed: int = 33
):
    design = motif_scaffolding.MotifScaffolding('1BCF')
    
    config_info, total_length, components, sampled_contig = design.parse_contig(
        contig="8-15/A92-99/16-30/A123-130/16-30/A47-54/16-30/A18-25/8-15"
    )
    
    motif_idx = design.get_motif_indices(sampled_contig)
    
    index_tuple = design.adjust_motif(
        "(self.C==1)",
        "[start - 1]",
        "[end - 1] + 1"
    )
    
    new_protein = design.build_scaffold(
        contig=config_info,
        total_length=total_length,
        components=components,
        index=index_tuple 
    )
    
    motif_seq, mask_aa = design.create_sequence_mask(
        protein=new_protein,
        total_length=total_length,
        index=motif_idx
    )
    
    print(motif_idx, total_length)
    print(motif_seq)
    print(mask_aa)
    
    final_protein = design.scaffold_sample(
        motif_seq, 
        mask_aa,
        new_protein,
        seed
    )
    
    
    
    return final_protein, sampled_contig, motif_idx, total_length
    
def sample_motif_scaffolding(
    name: str = None,
    num_samples: int = 100
):
    os.makedirs(f'../motif_scaffolding/{name}', exist_ok=True)
    result_file = os.path.join(f'../motif_scaffolding/{name}', 'backbone_results.csv')
    
    results_df = pd.DataFrame(columns=[
        'name',
        'number',
        'contig',
        'motif_indices',
        'ELBO',
        'seq',
        'total_length'
    ])
    
    for sample_number in range(1, num_samples + 1):
        filename = f'{name}_{sample_number}.pdb'
        sample_file = os.path.join(os.path.join(f'../motif_scaffolding/{name}', filename))
        if not os.path.exists(sample_file):
            protein, contig, motif, length = system_setup(random.randint(1, 10000))
            protein.to_PDB(os.path.join(f'../motif_scaffolding/{name}', filename))
            sequence = protein.sequence()
            elbo = f'{chroma.score(protein)['elbo'].score:.3f}'
            
            # Write results
            results_df = results_df._append({
                'name': name,
                'number': sample_number,
                'contig': contig,
                'motif_indices': motif,
                'ELBO': elbo,
                'seq': sequence,
                'total_length': length
            }, ignore_index=True)
            
        else:
            print(f'Skipping existing sample: {sample_file}')
        
    results_df.to_csv(result_file, index=False)
            
        
            
if __name__ == "__main__":
    chroma = Chroma()
    parser = create_parser()
    args = parser.parse_args()
    number = args.num
    name = args.pdb
    sample_motif_scaffolding(
        name=name,
        num_samples=number
    )
            