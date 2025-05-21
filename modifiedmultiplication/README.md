# ModifiedMultiplication

This repository contains the code for the implementation of the 'modified multiplication table method' by Jiang et al. (arXiv:1804.02733).

The main logic is in src/helpers/ProblemCreator.py

## Example usage
```python
from helpers import ProblemCreator
P = ProblemCreator.Problemcreator(N=143, blocksize=2, dim=(4,4))
P.generate_full_equation()
result,success_rate = P.run_get_success_rates(num_reads = 1000)
```