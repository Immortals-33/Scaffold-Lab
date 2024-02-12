import os
import numpy as np
import subprocess
import typing as T
import random
import shutil
from typing import *

import mdtraj as md
import MDAnalysis as mda
import biotite
import biotite.structure.io as strucio
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.application.dssp as dssp
from biotite.structure.residues import get_residues
from biotite.structure import get_chains
from biotite.sequence import ProteinSequence


def motif_extract(
    position: str,
    structure_path: Union[str, struc.AtomArray],
    atom_part: Optional[str] = "all-atom",
    split_char: str = "/"
) -> struc.AtomArray:
    """Extracting motif positions from input protein structure.

    Args:
        position (str): Motif region of input protein. DEMO: "A1-7/A28-79" corresponds defines res1-7 and res28-79 in chain A to be motif. 
        structure_path (Union[str, None]): Input protein structure, can either be a path or an AtomArray.
        atom_part (str, optional): _description_. Defaults to "all".
        split_char (str): Spliting character between discontinuous motifs. Defaults to "/".

    Returns:
        motif (biotite.structure.AtomArray): The motif positions extracted by user-specified way (all-atom / CA / backbone)
    """
    
    position = position.split(split_char)
    if isinstance(structure_path, str):
        array = strucio.load_structure(structure_path, model=1)
    else:
        array = structure_path
        
    motif_array = []
    for i in position:
        chain_id = i[0]
        i = i.replace(chain_id, "")
        if "-" not in i: # Single-residue motif
            start = end = int(i)
        else:
            start, end = i.split("-")
            print(start, end)
            start, end = int(start), int(end)
            
        if atom_part == "all-atom":
            motif_array.append(array[(array.chain_id==chain_id) & (array.res_id <= end) & (array.res_id >= start) & (array.hetero==False)])
        elif atom_part == "CA":
            motif_array.append(array[(array.chain_id==chain_id) & (array.res_id <= end) & (array.res_id >= start) & (array.hetero==False) & (array.atom_name=="CA")])
        elif atom_part == "backbone":
            motif_array.append(array[(array.chain_id==chain_id) & (array.res_id <= end) & (array.res_id >= start) & (array.hetero==False) & ((array.atom_name=="N") | (array.atom_name=="CA")| (array.atom_name=="C") | (array.atom_name=="O"))])
            
    motif = motif_array[0]
    for i in range(len(motif_array) - 1):
        motif += motif_array[i + 1]
    return motif

def rmsd(
    reference: Union[str, struc.AtomArray],
    target: Union[str, struc.AtomArray],
    
) -> float:
    
    # Handle input protein type
    if isinstance(reference, str):
        reference = strucio.load_structure(reference, model=1)
    if isinstance(target, str):
        target = strucio.load_structure(target, model=1)
        
    mask_reference = (((reference.atom_name == "N") | (reference.atom_name == "CA") | (reference.atom_name == "C") | (reference.atom_name == "O")) & (biotite.structure.filter_amino_acids(reference)))
    reference = reference[mask_reference]
    mask_target = (((target.atom_name == "N") | (target.atom_name == "CA") | (target.atom_name == "C") | (target.atom_name == "O")) & (biotite.structure.filter_amino_acids(target)) )
    target = target[mask_target]
    superimposed, _ = struc.superimpose(reference, target)
    rms = struc.rmsd(reference, superimposed)
    return rms

def calculate_secondary_structure(
    input: Union[str, struc.AtomArray] = None,
) -> Optional[list]:
    """
    Calculate protein secondary structure.
    """
    array = strucio.load_structure(input)
    app = dssp.DsspApp(array)
    app.start()
    app.join()
    sse = app.get_sse()
    sse_format = "".join(sse)
    
    alpha_composition = (sse_format.count("H") + sse_format.count("G") + sse_format.count("I")) / len(sse_format) if len(sse_format) > 0 else 0
    beta_composition = (sse_format.count("B") + sse_format.count("E")) / len(sse_format) if len(sse_format) > 0 else 0
    coil_composition = (sse_format.count("C") + sse_format.count("S")) / len(sse_format) if len(sse_format) > 0 else 0
    turn_composition = sse_format.count("T") / len(sse_format) if len(sse_format) > 0 else 0

    return [sse_format, alpha_composition, beta_composition, coil_composition, turn_composition]

def radius_of_gyration(
    input: str,
    atom_part: Optional[str] = "all-atom",
    ) -> float:
    t = mda.Universe(input)
    if atom_part == 'all-atom':
        rg = t.atoms.radius_of_gyration()
    if atom_part == 'backbone':
        backbone = t.select_atoms("backbone") 
        rg = backbone.atoms.radius_of_gyration()
    if atom_part == 'CA':
        c_alpha = t.select_atoms("name CA")
        rg = c_alpha.atoms.radius_of_gyration()
    return round(rg, 3)

def randomly_select_and_copy_pdb_files(src_folder, dst_folder, num_files):
    """
    Randomly select and copy a specified number of PDB files 
    from a source folder to a destination folder.

    :param src_folder: Path to the source folder containing PDB files.
    :param dst_folder: Path to the destination folder where files will be copied.
    :param num_files: Number of PDB files to randomly select and copy.
    """
    # Get a list of all PDB files in the source folder
    pdb_files = [file for file in os.listdir(src_folder) if file.endswith('.pdb')]

    # Check if the source folder has enough files
    if len(pdb_files) < num_files:
        raise ValueError(f"Source folder has only {len(pdb_files)} files, but {num_files} files are requested.")

    # Randomly select the specified number of PDB files
    selected_files = random.sample(pdb_files, num_files)

    # Create the destination folder if it doesn't exist
    if not os.path.exists(dst_folder):
        os.makedirs(dst_folder)

    # Copy each selected file to the destination folder
    for file in selected_files:
        shutil.copy(os.path.join(src_folder, file), os.path.join(dst_folder, file))