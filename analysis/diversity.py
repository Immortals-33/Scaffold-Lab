import os
import shutil
import numpy as np
import pandas as pd
import subprocess
import random
import shutil
import logging
import typing as T
from typing import Optional, Union, List, Tuple, Dict
from pathlib import Path

"""
Diversity Calculation.

Use Foldseek-Cluster to hierarchically cluster given sets of structures.
Input: A dictionary consisting of different methods and their corresponding paths,
       where structures of different lengths are saved under those paths.
Output: A `pd.DataFrame` object.

See example usage below and Foldseek-cluster (https://github.com/steineggerlab/foldseek?tab=readme-ov-file#cluster)
for further information and customized parameters.

NOTE: Foldseek might sometimes fail to cluster the structures.
(Emprically speaking, this is probably due to the low diversity of the folder)
"""

log = logging.getLogger(__name__)

def foldseek_cluster(
    input: Union[str, Path],
    assist_protein_path: Union[str, Path],
    tmscore_threshold: float = 0.6,
    alignment_type: int = 1,
    output_mode: str = 'FLOAT',
    save_tmp: bool=False,
    foldseek_path: Optional[Union[str, Path]] = None
) -> Union[float, dict]:

    if not os.listdir(input):
        return {"Clusters": 0, "Samples": 0, "Diversity": 0} if output_mode == 'DICT' else 0
    tmp_path = os.path.join(input, 'tmp')
    os.makedirs(tmp_path, exist_ok=True)

    output_prefix = os.path.join(input, 'diversity')

    cmd = f'foldseek easy-cluster \
            {input} \
            {output_prefix} \
            {tmp_path} \
            --alignment-type {alignment_type} \
            --tmscore-threshold {tmscore_threshold} \
            --alignment-mode 2 \
            -v 0'

    if foldseek_path is not None:
        cmd.replace('foldseek', str(foldseek_path))

    assist_num = 0
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        for failed_tmp_file in os.listdir(tmp_path):
            abs_path = os.path.join(tmp_path, failed_tmp_file)
            if os.path.islink(abs_path):
                target_path = os.readlink(abs_path)
                os.unlink(abs_path)
                pass
            else:
                shutil.rmtree(abs_path)

        shutil.copy(assist_protein_path, os.path.join(input, 'assist_protein.pdb'))
        log.info(f'Foldseek-clusters encountered an error. \
        Copied an assistant protein to resume clustering.')
        subprocess.run(cmd, shell=True, check=True)
        assist_num += 1


    result = pd.read_csv(f'{output_prefix}_cluster.tsv', sep='\t', header=None, names=['clusters', 'members'])
    unique_clusters = result['clusters'].nunique() - assist_num
    total_members = len(result) - assist_num
    diversity = round(unique_clusters / total_members, 3)

    if not save_tmp:
        os.remove(f'{output_prefix}_cluster.tsv')
        os.remove(f'{output_prefix}_rep_seq.fasta')
        os.remove(f'{output_prefix}_all_seqs.fasta')

    # Remove assistant protein
    if os.path.exists(os.path.join(input, 'assist_protein.pdb')):
        os.remove(os.path.join(input, 'assist_protein.pdb'))
    shutil.rmtree(tmp_path)

    if output_mode == 'FLOAT':
        return diversity
    elif output_mode == 'DICT':
        return {"Clusters": unique_clusters, "Samples": total_members, "Diversity": diversity}


def process_directories(methods_dict: Dict[str, List[str]]) -> pd.DataFrame:
    rows = []

    for method, dirs in methods_dict.items():
        for dir_path in dirs:
            diversity_result = foldseek_cluster(input=dir_path, output_mode='DICT',tmscore_threshold=0.5)
            length = dir_path.split('/')[-2].split('_')[-1]

            rows.append({
                'backbone_dir': dir_path,
                'method': method,
                'length': length,
                **diversity_result
            })

    results_df = pd.concat([pd.DataFrame([row]) for row in rows], ignore_index=True)

    return results_df

# Example Usage
#results_df = process_directories(methods_dict)
"""
Here, the `methods_dict` is a dictionary like:
{
    "Method_1": [
        "method_1/length_50",
        "method_1/length_100",
        ......
    ],
    "Method_2": [
        "method_2/length_50",
        "method_2/length_100",
        ......
    ],
    ......
}
where it is organized like {method_name}-[List of PDB folders to be clustered]
Then the results will be written into a `pd.DataFrame` object.
Also,
"""
#results_df.to_csv("TM_0.5_diversity_results.csv", index=False)

if __name__ == "__main__":
    results_df = process_directories(methods_dict)
    results_df.to_csv("TM_0.5_diversity_results.csv", index=False)
