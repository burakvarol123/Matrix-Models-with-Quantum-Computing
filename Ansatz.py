from qiskit import *
from qiskit.circuit import Parameter
import numpy as np


def ansatz_review_exact(N, depth):  # i only use ry gates now, exactly like in the paper
    circuit = QuantumCircuit(N)
    thetas = []
    for l in range(N * (depth + 1)):
        thetas.append(Parameter('θ_' + str(chr(l))))
    counter = 0
    for i in range(depth):
        for l in range(N):
            if l + counter < len(thetas):
                circuit.ry(thetas[l + counter + N], l )
        counter = counter + N
        for p in range(N):
            circuit.cx( (p+N-1)%N,(p + N) % N)
        circuit.barrier()
    for j in range(N):
        circuit.ry(thetas[j], j)
    for k in range(N):
        if k + int(N / 2) < N:
            circuit.cx(k, (k + int(N / 2)))
    return circuit


def ansatz_varqite_rev(N, depth):  # the order is now like the one in the paper
    if N % 2 != 0:
        raise ValueError('You cant give uneven total qubit number')
    circuit = QuantumCircuit(N)
    thetas = []
    for l in range(N * (depth + 1)):
        thetas.append(Parameter('θ_' + str(l)))
    counter_3 = 0
    for i in range(depth):
        counter_2 = 0

        for l in range(N):
            if l + counter_2 + counter_3 < len(thetas):
                circuit.ry(thetas[l + counter_2 + counter_3 + N], l + counter_2)
        counter_3 = counter_3 + N
        for p in range(N):
            circuit.cx(p, (p + 1) % N)
        circuit.barrier()
    counter = 0
    for j in range(N):
        if j + counter < N:
            circuit.ry(thetas[j + counter], j)
        if j + 1 + counter < N:
            circuit.rx(thetas[j + 1 + counter], j)
        if j + 1 + counter < N:
            circuit.cx(j, (j + int(N / 2)))
        counter = counter + 1

    return circuit


def ansatz_varqite(N, depth):
    if N % 2 != 0:
        raise ValueError('You cant give uneven total qubit number')
    circuit = QuantumCircuit(N)
    thetas = []
    for l in range(N * (depth + 1)):
        thetas.append(Parameter('θ_' + str(l)))
    counter = 0
    for j in range(N):
        if j + counter < N:
            circuit.ry(thetas[j + counter], j)
        if j + 1 + counter < N:
            circuit.rx(thetas[j + 1 + counter], j)
        if j + 1 + counter < N:
            circuit.cx(j, (j + int(N / 2)))
        counter = counter + 1

    counter_3 = 0
    for i in range(depth):
        circuit.barrier()
        counter_2 = 0

        for l in range(N):
            if l + counter_2 + counter_3 < len(thetas):
                circuit.ry(thetas[l + counter_2 + counter_3 + N], l + counter_2)
        counter_3 = counter_3 + N
        for p in range(N):
            circuit.cx(p, (p + 1) % N)
    return circuit



