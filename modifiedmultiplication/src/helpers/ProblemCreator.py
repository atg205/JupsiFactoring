import sympy as sp
import numpy as np
import os
import os.path
import pickle

from sympy import *
import re
from dwave.system import DWaveSampler, EmbeddingComposite
import dimod
import pandas as pd
import pickle
import os
import os.path
import dwave.inspector

import json
from dimod.serialization.json import DimodEncoder, DimodDecoder

class Problemcreator:
    def __init__(self, N, dim=None, lamda = 1, blocksize=0):

        self.N = N
        self.solver = ""

        config = {
            91:((3,4),2),
            143:((4,4),2),

            437: ((5,5),3),

            493: ((5,5),3),

            1591: ((6,6),3),

            59989:((8,8),3),
            
            61937:((8,9),3),
            253009:((9,9),3),
            376289:((10,10),3),

        }
        if dim:
            self.dim = dim
        else:
            self.dim, self.BLOCKSIZE = config[self.N]
            
        if blocksize > 0: # overwrite with parameter
            self.BLOCKSIZE = blocksize

        self.p_dim, self.q_dim = self.dim
        self.lamda = lamda

    def generate_full_equation(self):
        # Read multiplication table
        self.generate_multiplication_table()
        self.calculate_max_blocks_sums()
        self.generate_carry_rows()
        self.derive_cost_function()

        # Generate Ising Model
        self.square_terms()
        self.higher_to_second_order()
        self.generate_bqm()

    def get_width(self):
        return self.dim[0]
    
    def get_length(self):
        return self.dim[1]
    
    def subs(self, term, sub_rule):
        return term if term ==1  or term == 0 else term.subs(sub_rule)

    def split_2d_array(self, array):
        # append zeros so that array can be equally splitted

        if len(array[0])%self.BLOCKSIZE > 0:
        
            missing_columns = self.BLOCKSIZE-(len(array[0])%self.BLOCKSIZE)
        else:
            missing_columns = 0
        array = np.concatenate((np.array(array), np.zeros((array.shape[0],missing_columns))),axis=1 )

        split = len(array[0]) / self.BLOCKSIZE

        return np.split(array[::-1], indices_or_sections=split, axis=1)
    
    def get_block_overflow(self, blocksize, block_sum):
        overflow = len(bin(block_sum))-2 -  blocksize 
        return overflow if overflow > 0 else 0

    # append to already exisiting c_arr
    def append_to_c_array(self, c_arr, number, subs_rule):
        new_cs = sp.symbols(' '.join(f'c{i}' for i in range(len(c_arr), len(c_arr)+number)))    
        #update subs rule
        new_cs = list(new_cs) if number > 1 else [new_cs]
        for c in new_cs:
            subs_rule[c] = 1

        c_arr += new_cs

    def set_result(self, p,q):
        ps = {f'p{i}':bit for i, bit in enumerate(bin(p)[3:-1][::-1])}
        qs = {f'q{i}':bit for i, bit in enumerate(bin(q)[3:-1][::-1])}
        ps.update(qs)
        self.result = pd.DataFrame(ps,index=[0])

              
                                    



    def generate_multiplication_table(self):
        self.q = sp.symbols(' '.join([f'q{i}' for i in range(self.q_dim-2)])) #-2 because we dont consider the leading and trailing ones
        self.p = sp.symbols(' '.join([f'p{i}' for i in range(self.p_dim-2)]))
        
        if self.q_dim == 3:
            self.q=tuple([self.q])
        if self.p_dim == 3:
            self.p=tuple([self.p])
        p_total = [1] + list(self.p)+ [1]

        # pad with zeros to make all rows equal in length
        # so that they can be split into blocks
        self.rows = [[0]* i + [q_sym * p_sym for p_sym in p_total] +
                      [0] * (len(self.q)-i) for i, q_sym in enumerate(self.q)]

        # add first and last table row
        self.rows = [list(self.p)+[1] + [0] * (len(self.q)+1)] + self.rows + [[0] * len(self.q) + p_total]

    def calculate_max_blocks_sums(self):
        self.ones_subs_rule = {sym:1 for sym in list(self.p) + list(self.q)}

        ones_blocks= np.vectorize(self.subs)(np.array(self.rows), self.ones_subs_rule)

        ones_blocks = self.split_2d_array(ones_blocks)


        # int conversion necessary because otherwise zeros are 0.0
        #self.max_block_sums = [[''.join((row.astype(int).astype(str))) for row in block] for block in ones_blocks  ] # don't ask me about this ever

        # convert to binary and sum up rows
        # turn it around for binary conversion
        #self.max_block_sums = [sum(int(row[::-1],2) for row in block) for block in self.max_block_sums] 

    def generate_carry_rows(self):
        self.rows_with_carry = np.array(self.rows.copy())
        self.blocks = self.split_2d_array(self.rows_with_carry)
        self.c = []

        self.N_bin = bin(self.N)

        for i in range(len(self.blocks)):



            ones_blocks= np.vectorize(self.subs)(np.array(self.blocks[i]), self.ones_subs_rule)
            
            # int conversion necessary because otherwise zeros are 0.0
            max_block_sum = [''.join((row.astype(int).astype(str))) for row in ones_blocks]

            # convert to binary and sum up rows,
            # turn it around for binary conversion,
            max_block_sum = sum(int(row[::-1],2) for row in max_block_sum)
            overflow_bits = self.get_block_overflow(blocksize=self.BLOCKSIZE, block_sum=max_block_sum)
            
            if overflow_bits:
                # catch overflow,
                carry_index_start = (i+1)*self.BLOCKSIZE
                carry_index_end = carry_index_start + overflow_bits 
                
                if carry_index_end >= len(self.N_bin)-3: #CHANGED from 4 to 3
                    carry_index_end = len(self.N_bin)-3  #CHANGED
                    overflow_bits = carry_index_end -carry_index_start

                if carry_index_end <= carry_index_start:
                    continue

                column_diff = carry_index_end - len(self.rows_with_carry[0])
                if column_diff > 0:
                    zeros = np.zeros((self.rows_with_carry.shape[0],column_diff))
                    self.rows_with_carry = np.concatenate((self.rows_with_carry, zeros),axis=1)

                # append carries to initial matrix,
                carry_line = np.zeros(len(self.rows_with_carry[0]), dtype=object)
                self.append_to_c_array(self.c, overflow_bits, self.ones_subs_rule)
                carry_line[carry_index_start:carry_index_end] = self.c[-overflow_bits:]
                #immediately substract carries from last column of current block
                carry_line[(i+1)*self.BLOCKSIZE-1] = sp.Add(*[-carry_bit*2**(j+1) for j, carry_bit in enumerate(self.c[-overflow_bits:])])
                
                self.rows_with_carry = np.concatenate((self.rows_with_carry,np.array([carry_line])))

                # IMPROVEMENT: cut the carry for the 1 (last column in general)
                            # update blocks with new carries,
            self.blocks = self.split_2d_array(self.rows_with_carry)
    
    def derive_cost_function(self):
        # derive cost function
        block = self.split_2d_array(self.rows_with_carry)
        N_bin_rev = self.N_bin[2:-1][::-1] # reverse N_bin to because blocks are ordered from LSB to MSB also discard LSB
        # column wise addition, with 2 to the power of column position in block
        # BLOCK not updated in prev function !!
        self.cost_eq = [[sp.Add(*column) *2**i for i, column in enumerate(block.T)]
                         for block in self.blocks ]
        # add the columns of each block together
        self.cost_eq = [sp.Add(*column) for column in self.cost_eq]

        self.cost_eq = [block_equation-int(N_bin_rev[block_inx*self.BLOCKSIZE:(block_inx+1) 
                                                     * self.BLOCKSIZE][::-1],2)
                         for block_inx, block_equation in enumerate(self.cost_eq)]
        self.cost_eq 
    
    def set_lower_order_eq(self, lower_order_expr):
        self.lower_order_expr = lower_order_expr

    def square_terms(self):
        self.c = tuple(self.c)
        self.cost_expr = Add(*[block_expr ** 2  for block_expr in self.cost_eq])
        # expand and simplify with x**2  = x for x = 0,1
        self.cost_expr = expand(self.cost_expr).subs([(term**2, term) for term in self.c + self.p + self.q])

    def higher_to_second_order(self):
        # higher order to second order
        # there are no higher order c terms 
        # therefore we only consider combinations of p and q
        pq_combos = expand(Mul(*[sum(sym_list) for sym_list in [self.p,self.q]])).args



        # generate necessary t symbols for substitution
        self.t = symbols(' '.join([f't{i}' for i in range(len(pq_combos))]))

        # construct subs dict
        self.sub_dict = [[(pq, t_sym)] for pq, t_sym in zip(pq_combos,self.t) ]
       
        # for dynamic penalty weight give coeff 
        def construct_equality_constraint(pqt, coeff = 1):
            return (
                coeff* self.lamda*(pqt[0]*pqt[1]-2*pqt[0]*pqt[2]-2*pqt[1]*pqt[2]+3*pqt[2]) 
            )
        self.lower_order_expr = self.cost_expr.copy()


        constraints = []
        modified = True

        # run until no more modifiable terms
        while modified:
            modified = False
            modified_terms = []
            terms =[_arg for _arg in self.lower_order_expr.args]

            for term in terms:
                if (len(term.free_symbols) >= 3):
                    for sub_rule in self.sub_dict:
                        sub_term = term.subs(sub_rule)
                        # variable replacement starting
                        if sub_term != term:
                            pqt = list(sub_rule[0][0].args) + [sub_rule[0][1]]
                            constraints.append(construct_equality_constraint(pqt, coeff=abs(term.args[0])))
                            term = sub_term
                            modified = True # loop until no more modifications made
                            break
                
                modified_terms.append(term)
            self.lower_order_expr = Add(*modified_terms) 
        self.lower_order_expr += Add(*constraints) 

    def to_qubo(self):
        # construct qubo with named variables s, remove offset 
        lower_order_coeff = self.lower_order_expr.as_coefficients_dict()
        offset = lower_order_coeff[1]
        qubo = {(tuple(str(k).split('*'))*2)[0:2] : v  for k,v in lower_order_coeff.items() if len(str(k).split('*')) < 3 and k != 1}
        return (qubo, offset)

    def generate_bqm(self):
        qubo, offset = self.to_qubo()
        self.bqm = dimod.BinaryQuadraticModel.from_qubo(Q = qubo, offset=offset)

    def return_bqm(self):
        return self.bqm

    
    def init_solving(self):
        self.sampler = DWaveSampler(region='eu-central-1')  # alternativ um eine QPU in Nordamerika auszuwählen: sampler = DWaveSampler(region='na-west-1')
        self.solver = EmbeddingComposite(self.sampler)

        #print(f'Selected QPU {self.sampler.properties["chip_id"]} with {len(self.sampler.nodelist)} qubits')
    
    # Run bqms
    def run_get_success_rates(self, num_reads = 1000):
        if num_reads > 0:
            if not self.solver:
                self.init_solving()


            print(self.N)
            success_rates = {}
            results = {}
            valids = {}
            response = self.solver.sample(self.bqm, num_reads = num_reads, return_embedding=True)
            self.result = response.to_pandas_dataframe().sort_values(by='energy')

        
        p_columns = [f'p{i}' for i in range(len(self.p))]
        q_columns = [f'q{i}' for i in range(len(self.q))]

        p_df = self.result[p_columns]
        q_df = self.result[q_columns]

        self.result['a'] = p_df.apply(lambda row: int(''.join(['1']+list(row.values.astype(str))+['1'])[::-1],2), axis=1) 
        self.result['b'] = q_df.apply(lambda row: int(''.join(['1']+list(row.values.astype(str))+['1'])[::-1],2), axis=1)
        
        self.result['valid'] = self.result['a'] * self.result['b'] == int(self.N)
        valid = self.result[(self.result['valid'] == True)]
        success_rate = 0 if valid.empty  or num_reads == 0 else valid['num_occurrences'].sum()/ num_reads *100

        return (self.result, success_rate)
    
    def check_equation(self,p,q):
        p_bin = bin(p)[3:-1][::-1]
        q_bin = bin(q)[3:-1][::-1]
        sub_dict_t = [(el[0][1], el[0][0]) for el in self.sub_dict]
        sub_dict = [(self.p[i],p_bit) for i,p_bit in enumerate(p_bin)] + [(self.q[i],q_bit) for i,q_bit in enumerate(q_bin)]
        check = self.lower_order_expr.subs(sub_dict_t)
        check = check.subs(sub_dict)


        if len(self.c) == 0:
            return check == 0
        
        max_c = 2**len(self.c)
        for c in range(0, max_c):
            zeros = "0"*(len(self.c)-(len(bin(c))-2))
            c_subs = [(self.c[i],c_bit) for i,c_bit in enumerate(zeros+bin(c)[2:])]
            if check.subs(c_subs) == 0:
                return True
        return False
        print(check)
        print(sub_dict)

    def check_t(self, row):
        #p = row['a']
        #q = row['b']
        p = self.pres
        q = self.qres
        p_bin = bin(p)[3:-1][::-1]
        q_bin = bin(q)[3:-1][::-1]
        #sub_dict_t = [(el[0][1], el[0][0]) for el in self.sub_dict]
        sub_dict = [(self.p[i],p_bit) for i,p_bit in enumerate(p_bin)] + [(self.q[i],q_bit) for i,q_bit in enumerate(q_bin)]
        
        violations = 0
        t_substituted = {el[0][1]:el[0][0].subs(sub_dict) for el in self.sub_dict}


        for t, t_subs in t_substituted.items():
            if row[str(t)] != t_subs:
                violations +=1

        pdf = bin(int(row['a']))
        qdf = bin(int(row['b']))
        violations += sum([pdf[i] != bin(p)[i] for i in range(len(pdf))])
        violations += sum([qdf[i] != bin(q)[i] for i in range(len(qdf))])

        return violations

    def number_violated_constraints(self, pres, qres):
        self.pres = pres
        self.qres = qres
        if len(self.c) != 0: # calc only possible for max blocksize (no carries)
            return -1
         
        #valids = self.result[(self.result['valid'] == True)]
        #if not valids.empty:
        #    valids['violations'] = valids.apply(self.check_t, axis=1)
         #   print(valids)
        
        self.result['violations'] = self.result.apply(self.check_t, axis=1)
        return self.result