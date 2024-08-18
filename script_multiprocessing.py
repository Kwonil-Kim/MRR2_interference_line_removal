import os
import glob
import numpy as np
import multiprocessing
import subprocess

def process_file(input_file):
    output_dir = "./output/"
    config_file = "config_mrrqc.yaml"
    
    cmd = [
        "python",
        "main_mrr_interference_removal.py",
        config_file,
        input_file,
        output_dir]
    
    subprocess.run(cmd)

if __name__ == "__main__":
    input_dir = "./input/"
    fname_search = "*.nc"
    cnt_cpu = 28
    
    files = np.sort(glob.glob(os.path.join(input_dir, fname_search)))
    
    pool = multiprocessing.Pool(processes=cnt_cpu)
    pool.map(process_file, files)
    pool.close()
    pool.join()
