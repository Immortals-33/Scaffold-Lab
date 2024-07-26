import os
import time
import logging
from chroma import Chroma

def sample_chroma(chroma, lengths, num_samples):
    for length in lengths:
        os.makedirs(f"./benchmark/length_{length}", exist_ok=True) 
        start_time = time.time()
        for sample_number in range(1, num_samples + 1):
            filename = f"length_{length}_{sample_number}.pdb"
            sample_file = os.path.join(os.path.join(f"./benchmark/length_{length}/{filename}"))
            if not os.path.exists(sample_file):
                protein = chroma.sample(chain_lengths=[length])
                protein.to(os.path.join(f"./benchmark/length_{length}/{filename}"))
            else:
                print(f'Skipping existing sample: {sample_file}')
            
        elapsed_time = (time.time() - start_time) / num_samples
        logging.info(f'Average sample time for length {length}: {elapsed_time:.2f}.')
            
if __name__ == "__main__":
    chroma = Chroma()
    
    lengths = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    num_samples = 100
    
    sample_chroma(chroma, lengths, num_samples)
