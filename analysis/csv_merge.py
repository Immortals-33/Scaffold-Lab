import pandas as pd
import os

merged_data = pd.DataFrame()

for root, dirs, files in os.walk('/dssg/home/acct-clschf/clschf/zzq/methods/chroma/benchmark/motif_scaffolding/'):
    for file in files:
        if file == 'backbone_results.csv':
            csv_path = os.path.join(root, file)
            df = pd.read_csv(csv_path)
            merged_data = pd.concat([merged_data, df], ignore_index=True)
            
merged_data.to_csv('../total_results.csv', index=False)
