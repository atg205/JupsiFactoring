# Comparison of Factoring Algorithms on the D-Wave Quantum Annealer

## Repository Overview

This repository contains data and code related to the bachelor thesis titled **"Comparison of Factoring Algorithms on the D-Wave Quantum Annealer"**, conducted at Forschungszentrum Jülich and FH Aachen by Philipp Hanussek, and published in 2024. The thesis compares various methods for factoring integers using the D-Wave Quantum Annealer.

The thesis can be accessed via the following DOI link: [https://doi.org/10.34734/FZJ-2024-05254](https://doi.org/10.34734/FZJ-2024-05254)

## Factoring Methods

The repository includes data generated using the following methods:

### 1. Multiplication Circuit
- **Description**: This method utilizes a multiplication circuit to reverse the multiplication process to obtain factors of a given number.
- **Code Reference**: The data was generated using code provided by D-Wave, available in their documentation.
- **Documentation Link**: [D-Wave Multiplication Circuit Documentation](https://docs.ocean.dwavesys.com/en/stable/docs_dimod/reference/generated/dimod.generators.multiplication_circuit.html)

### 2. CFA Method (Ding et al.)
- **Description**: The CFA method, introduced by Ding et al., involves using a custom embedding to factorize integers on a quantum annealer.
- **Code Reference**: Data was generated using the code provided in the repository linked below.
- **Repository Link**: [Multiplier Encoder by Ding et al.](https://gitlab.com/jingwen.ding/multiplier-encoder/)
- **Paper Reference**: Ding et al., arXiv:2310.17574

### 3. Modified Multiplication Table (Direct Method) (Jiang et al.)
- **Description**: This method modifies the traditional multiplication table approach by splitting it into blocks for reduced complexity
- **Code Reference**: The data was generated using code available in the `JupsiFactoring/modifiedmultiplication` directory of this repository.
- **Paper Reference**: Jiang et al., "Quantum Annealing for Prime Factorization", arXiv:1804.02733

## Data Availability

All data generated from the aforementioned methods is available in the `data` folder of this repository. This folder is organized as follows:

- **`*_success.pckl` Files**: These files contain success probabilities for each method, grouped by problem size.
- **`ground_states.pckl`**: This file contains the percentage of samples in which the lowest energy state was found.
- **`success_rates_per_semiprime/` Folder**: This folder contains success rates for each semiprime, not grouped by problem size.

