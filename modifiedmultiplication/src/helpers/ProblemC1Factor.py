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

class ProblemC1Factor:
    def __init__(self, N):
        self.N = N
        lenN = len(bin(N))-2
        dim = (lenN-1, lenN //2 +1)
    def create_equation(self):
        self.p = sp.symbols([f'p{i}' for i in range(self.dim[0])])
        self.q = sp.symbols([f'q{i}' for i in range(self.dim[1])])
        self.n = sp.symbols([f'n{i}' for i in range(len(bin(N))-2)])
        symbol_arr = [self.p, self.q, self.n]
    
        symbol_arr_factors = [sp.Add(*[sym*2**i for i,sym in enumerate(arr)]) for arr in symbol_arr] 

        self.eq = (symbol_arr_factors[2]-symbol_arr_factors[0]*symbol_arr_factors[1])**2

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
                coeff *(pqt[0]*pqt[1]-2*pqt[0]*pqt[2]-2*pqt[1]*pqt[2]+3*pqt[2]) 
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


P = ProblemC1Factor(N=143)
print(P)
