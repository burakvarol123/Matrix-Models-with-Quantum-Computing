import numpy as np
from qiskit.circuit.library import *
from qiskit.extensions import UnitaryGate
from qiskit import *
from qiskit.quantum_info import SparsePauliOp


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

def position_x(N):
    return 0.5 * np.kron(ZGate(), krId(N - 1)) + np.kron(ZGate(), ktot(N - 1))


def Ham(N, a, b):
    x = position_x(N)
    mat = a ** 2 * np.linalg.matrix_power(x, b)
    return mat


def create_Hamilton_2(D): #creates x**2 with pauli strings
    c = ["I" + "I" * D]
    coeffs = [1 / 4]
    for m in range(D):
        c.append("I" + "I" * D)
        coeffs.append(2 ** (m - 1) + 2 ** (2 * m - 1))
        c.append("I" + "I" * m + "Z" + "I" * (D - m - 1))
        coeffs.append(-2 ** (m - 1) - 2 ** (2 * m - 1))
    for m in range(D):
        for n in range(D):
            if m > n:
                c.append("I" + "I" * D)
                coeffs.append(2 ** (m + n - 1))
                c.append("I" + "I" * (m) + "Z" + "I" * (D - m - 1))
                coeffs.append(-2 ** (m + n - 1))
                c.append("I" + "I" * (n) + "Z" + "I" * (D - n - 1))
                coeffs.append(-2 ** (m + n - 1))
                c.append("I" + "I" * (n) + "Z" + "I" * (m - n - 1) + "Z" + "I" * (D - m - 1))
                coeffs.append(2 ** (m + n - 1))
    e= [s + "I"*(D+1) for s in c]

    d = SparsePauliOp(e, coeffs)

    return d

if __name__ == "__main__":
    print(position_x(2))
    print(Ham(2,1,2))

