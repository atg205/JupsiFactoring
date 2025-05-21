import dwave_networkx as dnx
import os
import re
import pandas as pd
from collections import defaultdict

class Reader:
    def __init__(self, mc = False) -> None:
        self.mc =mc
        if mc:
            self.base_path = f'{os.getcwd()[:-4]}/mc/results/'
        else:
            self.base_path = f'{os.getcwd()[:-4]}/results/'
        
    def read_results(self):
        result_dirs = [os.path.join(self.base_path, bit_dir) for bit_dir in os.listdir(self.base_path)]
        result_file_paths = [os.path.join(path,file) for path in result_dirs for file in os.listdir(path)]

        df_dict = {'N': [], 'blocksize':[], 'N_dim':[],'success_rate': [], 'min_energy': []}
        for d in result_file_paths:
            if not re.findall('\.pickle',d):
                continue
            N = int(re.findall('(?<=results\/)\d+',d)[0])
            blocksize = int(re.findall('(?<=blocksize\=)\d+',d)[0])
            with open(d, 'rb') as f:
                df = pd.read_pickle(f)
                df_dict['N'].append(N)
                # get dimensions 
                if not self.mc:
                    p_dim = len([col for col in df.columns if col.startswith('p')])
                    q_dim = len([col for col in df.columns if col.startswith('q')])
                else:
                    p_dim = len([col for col in df.columns if re.match(r'^a\d+$',col)])
                    q_dim = len([col for col in df.columns if re.match(r'^b\d+$',col)])
                df_dict['N_dim'].append(p_dim + q_dim)
                df_dict['blocksize'].append(blocksize)
                df_dict['min_energy'].append(df['energy'].min())
                if df['energy'].min() == 0 and N == 1507:
                    print("hello")
                valid = df[(df['valid'] == True)]
                success_rate = 0 if valid.empty else valid['num_occurrences'].sum()/ df['num_occurrences'].sum() *100
                df_dict['success_rate'].append(success_rate)
    
        self.df=pd.DataFrame(df_dict)
    
    # return the dfs from sampling
    def return_dfs(self, N=0, blocksize=0):
        results = defaultdict(list)
        if N > 0:
            result_dirs = [os.path.join(self.base_path, str(N))]
        else:
            result_dirs = [os.path.join(self.base_path, bit_dir) for bit_dir in os.listdir(self.base_path)]

        result_file_paths = [os.path.join(path,file) for path in result_dirs for file in os.listdir(path)]
        for d in result_file_paths:
            N = int(re.findall('(?<=results\/)\d+',d)[0])
            df_blocksize = int(re.findall('(?<=blocksize\=)\d+',d)[0])
            if blocksize >0 and df_blocksize != blocksize:
                continue

            with open(d, 'rb') as f:
                df = pd.read_pickle(f)
                df["blocksize"] = df_blocksize
                results[N].append(df)
        
        return results
