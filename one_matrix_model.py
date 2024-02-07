from Hamilton import position_x, krId
import numpy as np
from qiskit.circuit.library import *
from create_distribution import create_distribution_varqite
from time import process_time
from qiskit.quantum_info import SparsePauliOp
import hamiltonian_2 as hm2



def create_lambda(pos,d, N):
    matrix = np.array(1)
    for j in range(d+1):
        if pos == j:
            matrix = np.kron(matrix, position_x(N))
        else:
            matrix = np.kron(matrix, krId(N))
    return np.real(matrix)

def create_lambda_2(
        index,
        dimension,
        qubits_per_dim
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



def vandermonde(d, N):
    lambdas = []
    for i in range(d+1):
        lambdas.append(create_lambda(i,d, N))
    determinant = np.identity(2**(N*(d+1)))
    for i in range(d+1):
        for j in range(i + 1, d+1):
            determinant = np.matmul(determinant, lambdas[j] - lambdas[i])

    return determinant

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
    a = vandermonde_2(2,3)
    b= vandermonde(2,3)
    c = create_lambda(2, 3,3)
    d = create_lambda_2(1,2,5)
    #print(position_x(3).shape)
    #print(vandermonde(2, 2)==np.matmul((b-a),np.matmul((c-a),(c-b))))
    print(a.to_matrix()== b)
    print(d)

