import os 
from helpers import generator
import numpy as np
import time

G = generator.Generator()


for size in range(20,23):#range(16,20): 
    Ns = [sp for sp in G.get_semiprimes(upper_places=size, exact_places=True) if len(bin(sp[0]))-2 > 3 and len(bin(sp[1]))-2 > 3 ]
    spacing = set([int(s) for s in np.linspace(0, len(Ns)-1, num=10)])
    for s in spacing:
        run_success = False
        error_count = 0
        if s >= len(Ns):
            break
        N_i = Ns[s][2]
        dim = (len(bin(Ns[s][1]))-2, len(bin(Ns[s][0]))-2)
        while not run_success:
            
            return_val =  os.system(f"/usr/bin/env /home/atg205/Documents/modifiedmultiplication/.venv/bin/python /home/atg205/Documents/modifiedmultiplication/src/main.py  -F {dim[0]} {dim[1]} {N_i}")
            if return_val != 0:
                print("-----------------error-----------------")
                time.sleep(2**error_count)
                error_count+=1
                if error_count == 10:
                    raise ValueError("error")
            else:
                run_success = True 