import numpy as np
import hamiltonian_2 as hm
import one_matrix_model as om

def matrix_terms(qubits_per_dim, dimension, lattice_spacing, pow):
    """
    Returns the matrix x^2 terms for test.
    """
    lambdas = []
    for i in range(dimension+1):
        lambdas.append(om.create_lambda_2(i, dimension, qubits_per_dim))
    hamiltonian = 0
    for i in range(dimension+1):
        hamiltonian += (lattice_spacing ** 2) * lambdas[i].power(pow)
    return hamiltonian 


