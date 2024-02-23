# Actions for various papers


import numpy as np
import one_matrix_model as om
import hamiltonian_2 as hm

def lambdas(
  d,
  qubits_per_dim      
):
    """
    Creates the eigenvalues for one matrix 
    """
    lambdas = []
    for i in range(d+1):
        lambdas.append(om.create_lambda_2(i, d, qubits_per_dim))
    return lambdas


def interaction_terms(
  d,
  a,
  lambdas,
  pow
):
    """
    creates interaction(pow=n) / propagator terms(pow=2)
    """
    hamiltonian = 0
    for i in range(d+1):
        hamiltonian += (a ** 2) * lambdas[i].power(pow)
    return hamiltonian


if __name__  == "__main__":
    d = 1
    qubits_per_dim = 2
    a = 1
    
    lambdas = lambdas(d,qubits_per_dim)
    
    action_simple =  - interaction_terms(d, a, lambdas, 2)     # action sum -lambda_1^2 - lambda_2^2
    
    




