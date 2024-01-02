import pandas as pd
import os
import glob

def merge_data_with_backbone_path_and_method(root_dir, method_name):
    merged_data = pd.DataFrame()

    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file == 'sc_results.csv':
                csv_path = os.path.join(root, file)
                df = pd.read_csv(csv_path)

                parent_dir = os.path.abspath(os.path.join(root, os.pardir))
                print(parent_dir)
                pdb_files = glob.glob(os.path.join(parent_dir, '*.pdb'))
                
                # Check if there is more than one .pdb file
                if len(pdb_files) > 1:
                    raise RuntimeError(f"More than one .pdb file found in {parent_dir}")
                elif pdb_files:
                    df['backbone_path'] = os.path.abspath(pdb_files[0])
                else:
                    df['backbone_path'] = None

                # Add the method column
                df['method'] = method_name

                merged_data = pd.concat([merged_data, df], ignore_index=True)
    
    return merged_data
            
root_dir = '/dssg/home/acct-clschf/clschf/zzq/Scaffold-Lab/refolding/Framediff'
method_name = 'FrameDiff'  # Replace with your desired method name
merged_data = merge_data_with_backbone_path_and_method(root_dir, method_name)
merged_data.to_csv('../unconditional/framediff.csv', index=False)