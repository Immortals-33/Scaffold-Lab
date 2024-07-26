import sys
import os
import shutil

sys.path.append('./evaluations/pipeline') # This is for correct import of Pipeline()
from evaluations.pipeline.pipeline import Pipeline

inverse_fold_model = 1
fold_model = 2 # dummy models for initializing Pipeline() class

p = Pipeline(inverse_fold_model, fold_model)
# Now suggest the directory storing .crd files are './opt/'
p._preprocess(
    coords_dir='./opt',
    output_dir='./', # This will automatically save output files in os.path.join(output_dir, 'pdbs')
    verbose=True
)

for i, file in enumerate(os.listdir('./pdbs')):
    length = file.split('_')[0]
    num = int(file.split('_')[1].split('.')[0])
    original_path = os.path.join('./pdbs', f'{length}_{num}.pdb')
    new_path = os.path.join(f'./pdbs/length_{length}', f'length_{length}_{num + 1}.pdb')
    if not os.path.exists(f'./pdbs/length_{length}'):
        os.mkdir(f'./pdbs/length_{length}')
    shutil.copy(original_path, new_path)
    print(f'Successfully moved {file}!')