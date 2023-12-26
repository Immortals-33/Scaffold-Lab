import pandas as pd
import os

merged_data = pd.DataFrame()

for root, dirs, files in os.walk('/dssg/home/acct-clschf/clschf/zzq/se3_diffusion/inference_outputs/21D_05M_2023Y_16h_28m_08s/'):
    for file in files:
        if file == 'sc_results.csv':
            csv_path = os.path.join(root, file)
            df = pd.read_csv(csv_path)
            merged_data = pd.concat([merged_data, df], ignore_index=True)
            
merged_data.to_csv('../unscore_results.csv', index=False)