import numpy as np
from qiskit.circuit.library import *
import hamiltonian_2 as hm2

def create_lambda_2(
        index,
        dimension,
        qubits_per_dim,
        
):
    if index == qubits_per_dim:
        raise Exception("your index is too high, index starts from 0. So maximal Z insertion is num_qubits-1")
    if index == 0:
        matrix = hm2.position_x(qubits_per_dim)
        for j in range(dimension):
            matrix = matrix ^ hm2.insert_i(qubits_per_dim)
    else:
        matrix = hm2.insert_i(qubits_per_dim)
        for j in range(dimension):
            if index  == j +1:
                matrix = matrix ^ hm2.position_x(qubits_per_dim)
            else:
                matrix = matrix ^ hm2.insert_i(qubits_per_dim)
    return matrix 

def vandermonde_2(
        dimension,
        qubits_per_dim
):
    lambdas = []
    for i in range(dimension+1):
        lambdas.append(create_lambda_2(i, dimension, qubits_per_dim))
    determinant = hm2.insert_i(qubits_per_dim * (dimension+1))
    for i in range(dimension + 1):
        for j in range(i+1 , dimension+1):
            determinant = determinant @ (lambdas[j] - lambdas[i])

    return determinant


if __name__ == "__main__":
    lambdas = []
    for i in range(1+1):
        lambdas.append(create_lambda_2(i, 1, 5))
    print(len(lambdas))

