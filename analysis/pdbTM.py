import os
import numpy as np
import pandas as pd
import subprocess
import typing as T
from typing import *
from pathlib import Path

def pdbTM(
    input: Union[Path, str],
    foldseek_path: Union[Path, str] = "/dssg/home/acct-clschf/clschf/zzq/methods/SaProt/bin/foldseek",
    save_tmp: bool = False,
    output: Optional[Union[Path, str]] = "./",
) -> Union[float, dict]:
    
    tmp_path = "/dssg/home/acct-clschf/clschf/zzq/se3_diffusion/tmp/"
    pdb100 = "/dssg/home/acct-clschf/clschf/zzq/foldseek/database/pdb100/pdb"
    # Check whether input is a directory or a single file
    if ".pdb" in input:
        cmd = f'{foldseek_path} easy-search \
                {input} \
                {pdb100} \
                {output}/{input}.m8 \
                {tmp_path} \
                --format-mode 4 \
                --format-output query,target,evalue,alntmscore,rmsd \
                --alignment-type 1 \
                --cluster-search 1'
        os.system(cmd)
        # "--cluster-search 1 doesn't work for unexpected reasons"
        # with open(output, "r") as result:
            