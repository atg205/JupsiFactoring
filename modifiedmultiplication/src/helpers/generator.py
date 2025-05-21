import dimod
import os
import re
import json

class Generator():
    def __init__(self, N=0, lamda = "", bqm = "", width = "", length = ""):
        self.base_path = os.path.join(os.getcwd())
        self.N = N
        self.output = ""
        self.bqm = bqm
        self.width = width
        self.length = length
        self.lamda = lamda


# source : https://www.delftstack.com/howto/python/sieve-of-eratosthenes-in-python/
    def sieveoferatosthenes(self, n):
        max = n + 1
        d = dict()
        for i in range(2, max):
            d[i] = True

        for i in d:
            factors = range(i, max, i)
            for f in factors[1:]:
                d[f] = False
        lst = [i for i in d if d[i] == True]
        return lst

    # return semiprimes either up to a certain number of places
    # or with exact places if variable set to True
    # format (p,q,N)
    def get_semiprimes(self, upper_places = 10, exact_places = False, limit = None):
        upper_boundary = 2 ** (upper_places)
        lower_boundary = upper_boundary / 2
        primes = self.sieveoferatosthenes(int(upper_boundary / 2) +1)
        semiprimes = [(p,q,p*q) for i, p in enumerate(primes) for q in primes[i:] if p * q < upper_boundary and (not exact_places or p*q > lower_boundary) ]  
        return sorted(semiprimes)[:limit]


    def get_multiplication_bqm(self, p):
        p_bin = bin(p[2])[2:]
        #dim = int(len(p_bin)/2)
        #if p_bin[0]:
        #    dim += 1
        #print(p)
        #print(len(bin(p[0]))-2)
        #print(len(bin(p[1]))-2)
        bqm = dimod.generators.multiplication_circuit(len(bin(p[0]))-2,len(bin(p[1]))-2)
        
        for i,place in enumerate(p_bin[::-1]):
            bqm.fix_variable('p' + str(i), int(place))
        return bqm.spin
    
        # def generate xubo readable file from bqm
    def get_bqm_from_file(self):
        bqm_path = os.path.join(self.   base_path, 'bqms')
        files = os.listdir(bqm_path)
        files = re.findall(f'\S+_{self.N}.json', ' '.join(files))
        assert(len(files) == 1)
        self.width = re.findall('\d+', files[0])[0]
        self.length = re.findall('\d+(?=_)', files[0])[1]
        with open(os.path.join(bqm_path, files[0]), 'rb') as f:
            self.bqm = dimod.binary.BinaryQuadraticModel.from_serializable(json.load(f))

    def convert_bqm_to_spin(self):
        self.bqm = self.bqm.spin

    # store linear and quadratic terms h and J
    def serialize_bqm(self):
        # map linear terms
        lin_map = {}
        for i,key in enumerate(self.bqm.linear):
            lin_map[key] = i
        
        self.output = [f'# QUBITS {len(self.bqm.linear)}\n']
        self.output += [f'# offset {self.bqm.offset}\n']
        self.output += [f'# quibitmap {lin_map}\n']
        self.output += [f'# vartype {self.bqm.vartype}\n']
        self.output += [f"{lin_map[k]} {lin_map[k]} {v}\n" for k, v in self.bqm.linear.items()]
        self.output += [f"{lin_map[k[0]]} {lin_map[k[1]]} {v}\n" for k, v in self.bqm.quadratic.items()]
    
    def save_output_to_file(self):
        filename_ext = f'_lamda_{self.lamda}' if self.lamda  else ''
        filepath = os.path.join(os.getcwd(), 'IsingTerms','dynamic',str(self.N), f'{self.width}_{self.length}_{self.N}{filename_ext}.ising')
        with open(filepath, 'w') as f:
            f.writelines(self.output)


    # serialize semiprimes with at most upper places
    def serialize(self, upper = 10, exact_places = False, limit = None):
        output = []
        for semip in self.get_semiprimes(upper_places = upper, exact_places = exact_places, limit = limit):
            bqm = self.get_multiplication_bqm(semip)
            output = self.serialize_bqm(bqm)

            file = open(f'IsingTerms/{semip[0]}x{semip[1]}_{semip[2]}.ising', 'w')
            file.writelines(output)
            file.close()

    # write a provided bqm to file in Ising (xubo) format
    def bqm_to_file(self, bqm):
        self.serialize_bqm(bqm)
        with open(os.path.join(self.base_path, "IsingTerms", f"{self.N}")) as f:
            json.dump(dimod.binary.BinaryQuadraticModel.to_serializable(self.bqm),f)
#for i in range(10,16):
#serialize(upper=10)

#print("..")

