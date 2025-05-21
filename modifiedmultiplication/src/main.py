from helpers import generator, ProblemCreator
import os
import os.path
import pickle 
import dimod.generators
import pandas as pd
import numpy as np
import re
from dwave.system import DWaveSampler, EmbeddingComposite
import sys


N = 493
def serialize_bqm():

    filepath = os.path.join(os.getcwd(), 'problems',f'{dim[0]}_{dim[1]}_{N}.pickle')


    G = generator.Generator(N=N)

    G.get_bqm_from_file()
    G.convert_bqm_to_spin()
    G.serialize_bqm()
    G.save_output_to_file()
    print(G.bqm)

def checks():
    G = generator.Generator()

    for size in range(10,20): 
        Ns = [sp for sp in G.get_semiprimes(upper_places=size, exact_places=True) if len(bin(sp[0]))-2 > 3 and len(bin(sp[0])) > 3 ]
        spacing = set([int(s) for s in np.linspace(0, len(Ns)-1, num=10)])
        for s in spacing:
            if s >= len(Ns):
                break
            N_i = Ns[s][2]
            dim = (len(bin(Ns[s][0]))-2, len(bin(Ns[s][1]))-2)
            #P = ProblemCreator.Problemcreator(N=N_i)
            #P.generate_full_equation()
            #N_len = P.get_width() + P.get_length()
            N_len = len(bin(N_i))-2
            for block_size in [N_len-1]:
                P = ProblemCreator.Problemcreator(N=N_i, blocksize=block_size, dim=dim)
                P.generate_full_equation()
                P.set_result(Ns[s][0], Ns[s][1])
                result,_ = P.run_get_success_rates(num_reads = 0)
                if result.head(1)['valid'].iloc[0]:
                    print("valid")
                else:
                    print("not valid")

                # equation check
                if P.check_equation(Ns[s][0], Ns[s][1]):
                    print("eq valid")
                else:
                    print("eq not valid")



def generate_Ising():
    for lamda in range(32,35):
        lamda = lamda / 100
        P = ProblemCreator.Problemcreator(N=N, lamda = lamda)
        P.generate_full_equation()
        G = generator.Generator(N = N, bqm = P.return_bqm(), lamda = lamda, width = P.get_width(), length = P.get_length())
        G.convert_bqm_to_spin()
        G.serialize_bqm()
        G.save_output_to_file()

def compare_block_sizes():
    N_arr = []
    block_size_arr = []
    binary_arr = []
    quadratic_arr = []

    for N_i in [N]:
    #for N in [143]:
        P = ProblemCreator.Problemcreator(N=N_i)
        N_len = P.get_width() + P.get_length()

        for i in range(1,N_len):

            P = ProblemCreator.Problemcreator(N=N_i, blocksize=i)
            P.generate_full_equation()
            print(i)
            bqm = P.return_bqm()
            print(bqm.num_variables)
            print(bqm.num_interactions)
            print("-----------------")
            N_arr.append(N_i)
            block_size_arr.append(i)
            binary_arr.append(bqm.num_variables)
            quadratic_arr.append(bqm.num_interactions)

    data = {'N': N_arr, 'block_size': block_size_arr,'num_bin': binary_arr, 'num_quadratic': quadratic_arr}
    return pd.DataFrame(data=data)

def create_histograms():
    dfs = []
    for N_i in [N]:
    #for N in [143]:
        P = ProblemCreator.Problemcreator(N=N_i)
        N_len = P.get_width() + P.get_length()
        for i in range(1,N_len):
            P = ProblemCreator.Problemcreator(N=N_i, blocksize=i)
            P.generate_full_equation()
            bqm = P.return_bqm()
            linears = list(bqm.linear.values())
            quadratics = list(bqm.quadratic.values())
            # make linears equal length
            linears += [None]*(len(quadratics)-len(linears))

            dfs.append(pd.DataFrame({'N':[N]*len(quadratics),'blocksize':[i]*len(quadratics), 'linear_terms': linears, 'quadratic_terms': quadratics}))
    return dfs

def serialize_result_df(df, N_i, blocksize=0, method='mod'):
        mc_path = ''
        if method == 'mc':
            mc_path = 'mc/'
        
        path = os.path.join(os.getcwd(),mc_path, 'results', str(N_i))
        if os.path.isdir(path):
            max_round = list(map(int,re.findall(f'(?<=round\=)\d*(?=_blocksize={blocksize})', ' '.join(os.listdir(path)))))
            max_round = 0 if len(max_round) == 0 else max(max_round)
        else:
            max_round = 0
            os.makedirs(path)
                    
        with open(os.path.join(path, f'round={max_round+1}_blocksize={blocksize}.pickle'),'wb') as f:
            df.to_pickle(f)



def run():
    runs = 3
    blocksize = 17
    for N_i in [253009]:
        P = ProblemCreator.Problemcreator(N=N_i, blocksize=blocksize)
        P.generate_full_equation()
    
        for run_ in range(runs):   
            #result,success_rate = P.run_get_success_rates()
            result, _ = P.run_get_success_rates()
            #result = P.number_violated_constraints(503,503).sort_values(by='energy')
            serialize_result_df(result, N_i, blocksize=blocksize)
            



def run_compare_block_sizes(dim, N_i):
    runs = 5


    N_len = len(bin(N_i))-2
    for block_size in [N_len-1]:
        P = ProblemCreator.Problemcreator(N=N_i, blocksize=block_size, dim=dim)
        P.generate_full_equation()
        for _ in range(runs):        
            result,success_rate = P.run_get_success_rates()
            result = result.sort_values(by='energy')
            serialize_result_df(result, N_i, block_size)

def run_dwave_mc(dim, N, Ta=20):
    sampler = DWaveSampler(region='eu-central-1')  # alternativ um eine QPU in Nordamerika auszuwählen: sampler = DWaveSampler(region='na-west-1')
    solver = EmbeddingComposite(sampler)
    rounds = 10
    bqm = dimod.generators.multiplication_circuit(dim[0], dim[1])
    p_dict = {f'p{i}':0 for i in range(dim[0] + dim[1])}
    for i,bi in enumerate(bin(N)[2:][::-1]):
        p_dict[f'p{i}'] = int(bi)
    bqm.fix_variables(p_dict)
    props =solver.properties
    print(p_dict)
    for _ in range(rounds):
        print(p_dict)

        response = solver.sample(bqm, num_reads = 1000)
        result = response.to_pandas_dataframe().sort_values(by='energy')

        p_columns = [f'a{i}' for i in range(dim[0])]
        q_columns = [f'b{i}' for i in range(dim[1])]

        p_df = result[p_columns]
        q_df = result[q_columns]

        result['a'] = p_df.apply(lambda row: int(''.join(list(row.values.astype(str)))[::-1],2), axis=1) 
        result['b'] = q_df.apply(lambda row: int(''.join(list(row.values.astype(str)))[::-1],2), axis=1)
        
        result['valid'] = result['a'] * result['b'] == int(N)
        serialize_result_df(df=result,N_i=N,blocksize=0,method='mc')
        print(result)



#dim = (int(sys.argv[-3]), int(sys.argv[-2]))
#N_i = (int(sys.argv[-1]))
#dim = (4,4)
#N_i = 9


#dim = (16, 4)

#N_i = 524491
#print(N_i)
#print(dim)

#run_dwave_mc(dim=dim, N=N_i)
#run_compare_block_sizes(dim=dim, N_i=N_i)

for N in [91,143,437,493]:
    generate_Ising()