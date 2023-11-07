import numpy as np
from qiskit.circuit.library import *
from qiskit.quantum_info import Operator


def Ham(N, a, b):
    def krZ(m, N):
        matrix = np.array(1)
        for j in range(N):
            if m == j:
                matrix = np.kron(matrix, ZGate())
            else:
                matrix = np.kron(matrix, IGate())
        return matrix

    def krId(N):
        matrix = np.array(1)
        for j in range(N):
            matrix = np.kron(matrix, IGate())
        return matrix

    def ktot(N):
        def k_m(N, a, m):
            k = (krId(N) - krZ(m, N)) / 2.0
            return k

        sum = 0
        for h in range(N):
            sum += 2 ** h * k_m(N, 1, h)
        return sum

    mat = a ** 2 * np.linalg.matrix_power((0.5 * np.kron(ZGate(), krId(N - 1)) + np.kron(ZGate(), ktot(N - 1))), b)
    op = Operator(mat)
    return op
