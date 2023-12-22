import os
import logging
import random
import typing as T
from typing import *
from pathlib import Path

import torch
from chroma import Chroma, Protein, conditioners, api
from chroma.constants.sequence import AA20_3_TO_1


device = 'cuda' if torch.cuda.is_available() else 'cpu'
chroma = Chroma()

ALPHABET = "ABCDEFGHJKLMNOPQRSTUVWXYZ"

class MotifScaffolding():
    """
    Use Chroma Substructure conditioner to do motif-scaffolding task on protein design.
    The class handle the following things:
    1. Read the motif information (1D AA-types along with 3D structure positions)
    2. Use the input contig to initialize a dummy protein structure with manually-favorable length
    3. Fill the motif region with motif information and create a corresponding SubstructureConditioner
    4. Use SubstructureConditioner to generate a new motif-scaffolding protein.
    
    Args:
        `input_protein`: Can either be a path of reference protein or a PDB id.
        `contig`: Use `parse_contig()` to read the contig.
    
    """
    def __init__(self, input_protein:Union[str, Path]) -> None:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.protein = Protein(input_protein, canonicalize=True, device=device)
        self.X, self.C, self.S = self.protein.to_XCS()
        self.alphabet = ALPHABET
        
    def parse_contig(
        self,
        contig: str = None,
        desired_total_length: Union[None, Tuple[int, int]] = None,
        split_char: Optional[str] = "/",
    ) -> [Tuple[List, List, List], int, List, str]:
        components = contig.split(split_char)
        random_sample_length = []
        motif_length = []
        motif_positions = []
        
        # These are designed to be saved and returned
        resolved_contig = contig
        total_motif_positions = []
        
        # Handle total length situation
        if desired_total_length is not None:
            positions = [(part[0], *map(int, part[1:].split('-'))) for part in components if part[0] in ALPHABET]
            motif_total_length = 0
            
            for i, information in enumerate(positions):
                start, end = positions[i][1:]
                single_motif_length = end - start + 1
                motif_total_length += single_motif_length

            min_total_length, max_total_length = desired_total_length
            print(min_total_length, max_total_length, motif_total_length)
            scaffold_length = random.randint(min_total_length - motif_total_length, max_total_length - motif_total_length)
        else:
            scaffold_length = None
            
        # Handle scaffold components
        scaffold_components = [part for part in components if part[0] not in ALPHABET]
        for i, part in enumerate(scaffold_components):
            start, end = map(int, part.split("-"))
            if desired_total_length is not None:
                if i < len(scaffold_components) - 1:
                    max_length = min(end, scaffold_length)
                    sampled_length = random.randint(start, max_length)
                    scaffold_length -= sampled_length
                else:
                    sampled_length = scaffold_length  # Assign the remaining length to the last scaffold component
            else:
                sampled_length = random.randint(start, end)
            random_sample_length.append(sampled_length)
            resolved_contig = resolved_contig.replace(part, f"{sampled_length}-{sampled_length}", 1)
                

        # Handle motif components 
        for part in components:
            if part[0] in ALPHABET:
                chain_id = part[0]
                positions = part[1:]
                start, end = map(int, positions.split("-"))
                motif_length.append(end - start + 1)
                motif_positions.append((chain_id, start, end))
        
        total_length = sum(random_sample_length) + sum(motif_length)
        components = resolved_contig.split(split_char)
        #insert_positions = [random_sample_length[0]]
        #for i in range(1, len(random_sample_length)):
         #   insert_positions.append(insert_positions[i-1] + random_sample_length[i] + motif_length[i-1])
            
        return (random_sample_length, motif_length, motif_positions), total_length, components, resolved_contig
    
    def get_motif_indices(
        self,
        resolved_contig: str = None
    ) -> List:
        cumulative_length = 0
        total_motif_positions = []
        for part in resolved_contig.split("/"):
            if part[0] in ALPHABET:
                # Motif part
                start, end = map(int, part[1:].split("-"))
                motif_start_final = cumulative_length + 1
                motif_end_final = cumulative_length + (end - start + 1)
                total_motif_positions.extend(range(motif_start_final, motif_end_final + 1))
                cumulative_length += (end - start + 1)
            else:
                # Scaffold part
                scaffold_length = int(part.split("-")[0])
                cumulative_length += scaffold_length
        return total_motif_positions
    
    def adjust_motif(
        self,
        chain_value: str,
        index_value_1: str,
        index_value_2: str
    ) -> Tuple[str, str, str]:
        return chain_value, index_value_1, index_value_2
        
            
    def build_scaffold(
        self,
        contig: Tuple[List, List, List],
        total_length: Optional[int],
        components: List,
        index: Tuple[str, str, str]
    ):
        """
        Build a scaffold with motif positions filled, and scaffold positions with dummy variables (0).
        
        Args:
        contig: Contig after (randomly) sampled. It's a string with format like input contig.
        total_length: Total length of protein after length (randomly) sampled.
        components: Different part of the sampled contig.
        index: Used to extract motif information from original protein.
               NOTE: As Chroma protein system use different chain identification from standard PDB files,
                     the selection of motif differs case by case.
                     We suggest users to test the selection of motif index before using them.
                     The current way is very inelegent and is yet to be developed.
        Returns:
        A new protein object with scaffold region filled with dummy variables, 
        and motif region filled with exact coordinates in original protein.
        """
        random_sample_length, motif_length, motif_positions = contig[:]
        
        # Build a full-length empty protein
        X_scaffold = torch.zeros(1, total_length, 4, 3).to(device)
        C_scaffold = torch.full((1, total_length), 1).to(device)
        S_scaffold = torch.full((1, total_length), 0).to(device)
        
        current_position = 0
        motif_index = 0
        scaffold_index = 0
        c_index, m_index_1, m_index_2 = index # Handle different PDB cases
        
        for part in components:
            if part[0] not in self.alphabet:
                length = random_sample_length[scaffold_index]
                scaffold_index += 1
            else:
                chain_id, start, end = motif_positions[motif_index]
                chain_idx = ord(chain_id.upper()) - ord('A') + 1 
                motif_start_idx = eval(f'{c_index}.nonzero(as_tuple=True)[1]{m_index_1}')
                motif_end_idx = eval(f'{c_index}.nonzero(as_tuple=True)[1]{m_index_2}') # Very inelegant but just for handling cases right now
                length = motif_length[motif_index]
                motif_index += 1
                
                X_scaffold[:, current_position:current_position+length, :, :] = self.X[:, motif_start_idx:motif_end_idx, :, :].to(device)
                C_scaffold[:, current_position:current_position+length] = self.C[:, motif_start_idx:motif_end_idx].to(device)
                S_scaffold[:, current_position:current_position+length] = self.S[:, motif_start_idx:motif_end_idx].to(device)
                
            current_position += length
            
        # This is important, otherwise when the motif position will be lost when creating the new protein object
        C_scaffold[:, :] = torch.full((1, total_length), 1).to(device)
        
        protein_new = Protein.from_XCS(X_scaffold, C_scaffold, S_scaffold)
        
        return protein_new
    
    def create_sequence_mask(
        self,
        protein,
        total_length: int,
        index: List = None
    ):
        """
        Define which positions in sequence to be fixed.
        
        Args:
        `protein`: Protein object created by `build_scaffold()` function.
        `total_length`: Total length of the protein.
        `index`: Motif positions extracted by `get_motif_indices` function.
               Used to create sequence mask for fixing motif sequence and structure
               when making conditional design.
        Returns:
        `motif_seq`: Indices of motif positions in the final protein.
                     Fed into the `save_selection()` function in `Protein` class.
        `modified_mask_aa`: Amino acids information for positions need to be fixed.
                            Fed into `chroma.sample()` function during sampling.
        """
        X2, C2, S2 = protein.to_XCS()

        mask_aa = torch.Tensor(total_length * [[1] * 20]).to(device)
        motif_indices = (S2 != 0).nonzero(as_tuple=True)[1]
        
        # REMEMBER THIS, the "-1" is very important in this case
        motif_indices = torch.Tensor([i - 1 for i in index]).long().to(device)
        
        for idx in motif_indices:
            aa_index = S2[0, idx].item()
            mask_aa[idx] = torch.Tensor([0] * 20).to(device)
            mask_aa[idx][aa_index] = 1
            
        motif_seq = motif_indices.tolist()
        modified_mask_aa = mask_aa.unsqueeze(0).to(device)

        return motif_seq, modified_mask_aa
    
    def scaffold_sample(
        self,
        motif_seq,
        modified_mask_aa,
        protein_new,
        random_seed
    ):
        protein_new.sys.save_selection(gti=motif_seq, selname="motif")
        
        
        
        conditioner = conditioners.SubstructureConditioner(
            protein=protein_new,
            backbone_model=chroma.backbone_network,
            #rg=True,
            selection='namesel motif').to(device)
        
        torch.manual_seed(random_seed)
        protein_final = chroma.sample(
            protein_init=protein_new,
            conditioner=conditioner,
            design_selection=modified_mask_aa,
            langevin_factor=4.0,
            langevin_isothermal=True,
            inverse_temperature=8.0,
            sde_func='langevin',
            steps=500,
        )
        
        return protein_final