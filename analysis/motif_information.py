import re
import pandas as pd

"""
This script is written for extracting useful information from RFdiffusion output log.
Some specific parts were added to be printed manually, so it might not be generalized.
*Watch out! Some of the information can't be used directly.

However, the general output log does show some of the information, so I put it here
in case it can be adapted for personal usage.

Patterns:
    `making_design`: This is for extracting the START of a single sample. 
                        Whenever this pattern is caught, a new design is being sampled.
    `sample_mask`: This is for extracting the CONTIG for a specific design.
                    For motif-scaffolding task in RFdiffusion, contigs are sampled randomly each time,
                    so this pattern can catch contig information for each sample.
    `mask_1d_pattern': Also known as DIFFUSION_MASK, where the motifs are fixed and the scaffold 
                        are designed. This is a list of boolean value where True=motif, False=scaffold.
    `sampled_motif_rmsd': The final sampled MOTIF-RMSD value for each specific sample.
    `finished_design`: This is for catching the TIME usage for each sample.
"""
# Regex patterns for extracting information
making_design_pattern = re.compile(r'\[INFO\] - Making design (.*?/([^/]+)/([^_]+)_([^_]+))')
sample_mask_pattern = re.compile(r"'sampled_mask': \['([^']+)'\]")
mask_1d_pattern = re.compile(r"'mask_1d': (\[.*?\])")
sampled_motif_rmsd_pattern = re.compile(r'Sampled motif RMSD: (\d+\.\d+)')
finished_design_pattern = re.compile(r'Finished design in (.+) minutes')

# Process the log file
def process_log_file(file_path):
    results = []
    between_timesteps = False
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        current_design = {}
        for i, line in enumerate(lines):
            # Check for 'Making design' line
            if making_design_match := making_design_pattern.search(line):
                if current_design:
                    results.append(current_design)
                aa, bb, pdb_name, sample_num = making_design_match.groups()
                pdb_name = pdb_name.split('/')[0]
                
                current_design = {'pdb_name': pdb_name, 'sample_num': sample_num}
            
            # Check for 'sampled_mask' line
            if sample_mask_match := sample_mask_pattern.search(line):
                current_design['contig'] = sample_mask_match.group(1)
            
            # Check for 'mask_1d' line
            if mask_1d_match := mask_1d_pattern.search(line):
                mask_1d = eval(mask_1d_match.group(1))
                current_design['mask'] = mask_1d
                current_design['motif_indices'] = [i + 1 for i, val in enumerate(mask_1d) if val]

            # When catching motif-RMSD, set `between_timesteps=True`, 
            # After catching, reset its value to continue to another sample.
            if 'Timestep 3,' in line:
                between_timesteps = True

            # Check for 'Sampled motif RMSD' within the right timestep
            if between_timesteps:
                if sampled_motif_rmsd_match := sampled_motif_rmsd_pattern.search(line):
                    current_design['motif_RMSD'] = sampled_motif_rmsd_match.group(1)

            # Unset flag when 'Timestep 2' is encountered
            if 'Timestep 2,' in line:
                between_timesteps = False

            # Check for 'Finished design' line
            if finished_design_match := finished_design_pattern.search(line):
                current_design['time'] = finished_design_match.group(1)

    if current_design:  # Append the last design
        results.append(current_design)

    return results

# Process the log file and create a DataFrame
log_file_path = './motif_scaffolding_5.out'
data = process_log_file(log_file_path)
df = pd.DataFrame(data)

# Save to CSV
csv_file_path = 'motif_results_5.csv'
df.to_csv(csv_file_path, index=False)
