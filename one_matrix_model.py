from Hamilton import position_x, krId
import numpy as np
from qiskit.circuit.library import *
from create_distribution import create_distribution_varqite


def create_lambda(pos,d, N):
    matrix = np.array(1)
    for j in range(d+1):
        if pos == j:
            matrix = np.kron(matrix, position_x(N))
        else:
            matrix = np.kron(matrix, krId(N))
    return np.real(matrix)


def vandermonde(d, N):
    lambdas = []
    for i in range(d+1):
        lambdas.append(create_lambda(i,d, N))
    determinant = np.identity(2**(N*(d+1)))
    for i in range(d+1):
        for j in range(i + 1, d+1):
            determinant = np.matmul(determinant, lambdas[j] - lambdas[i])

    return determinant


if __name__ == "__main__":
    #a = create_lambda(0, 2,2)
    #b = create_lambda(1, 2,2)
    #c = create_lambda(2,2,2)
    #print(position_x(3).shape)
    #print(vandermonde(2, 2)==np.matmul((b-a),np.matmul((c-a),(c-b))))

    lambdas = []
    for i in range(d + 1):
        lambdas.append(create_lambda(i, d, N))
    ham = 0
    for i in range(d+1):
        ham += 1 ** 2 * np.linalg.matrix_power(lambdas[i], 2)
    #print(len(lambdas))
    print(lambdas[0])
    print(ham)
    
    varqite = create_distribution_varqite(N, 2, ham, 1)


