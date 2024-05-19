"""
Ansatz.py

This module provides functions to generate various quantum circuits (ansatz) 
for use in quantum algorithms, particularly for variational quantum algorithms 
like Quantum Boltzmann Machines and imaginary variational time evolution.

Functions:
- ansatz_review_exact(qubits, depth): Generates an ansatz designed based on the 
  paper "Variational QBMs", with a specified number of qubits and depth.
- ansatz_review_exact_full_cry(qubits, depth): Generates a full entanglement 
  ansatz with Controlled RY gates, optimized for accuracy but time-consuming 
  to compute.
- ansatz_cry_optimized(qubits, depth, drawmode=False): Generates a full 
  entanglement ansatz with Controlled RY gates, optimized for imaginary 
  variational time evolution of diagonal Hamiltonians.

Each function returns a QuantumCircuit object configured according to the 
specified parameters.

Examples:
    # Generate an ansatz based on the paper "Variational QBMs"
    circuit = ansatz_review_exact(qubits=4, depth=2)

    # Generate a full entanglement ansatz with Controlled RY gates
    circuit_full_cry = ansatz_review_exact_full_cry(qubits=4, depth=2)

    # Generate a full entanglement ansatz with Controlled RY gates,
    # optimized for imaginary variational time evolution
    circuit_optimized = ansatz_cry_optimized(qubits=4, depth=2, drawmode=True)

Dependencies:
- qiskit
"""
from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from typing import List


def ansatz_review_exact(qubits: int, depth: int) -> QuantumCircuit:  
    """
    Ansatz designed in the way of the paper Variational QBMs. Close but just
    an approximate.
    Parameters:
        qubits (int): The number of qubits in the quantum circuit.
        depth (int): The depth of the quantum circuit.
    Returns:
        QuantumCircuit: A quantum circuit with the specified number of qubits
        and depth.
    """
    circuit = QuantumCircuit(qubits)
    thetas: List[Parameter] = []

    for l in range(qubits * (depth + 1)):
        thetas.append(Parameter('θ_' + str(chr(l))))
    
    counter = 0
    for i in range(depth):
        for l in range(qubits):
            if l + counter < len(thetas):
                circuit.ry(thetas[l + counter + qubits], l)
        counter += qubits
        
        for p in range(qubits):
            circuit.cx((p + qubits - 1) % qubits, (p + qubits) % qubits)
        
        circuit.barrier()
    
    for j in range(qubits):
        circuit.ry(thetas[j], j)
    
    for k in range(qubits):
        if k + int(qubits / 2) < qubits:
            circuit.cx(k, k + int(qubits / 2))
    
    return circuit


def ansatz_review_exact_full_cry(qubits: int, depth: int) -> QuantumCircuit: 
    """
    Full Entanglement Ansatz with Controlled Ry Gates. The best by far, but
    takes too much time. 
    Can be massively optimized.
    Args:
        qubits (int): Number of qubits.
        depth (int): Depth of the circuit.
    Returns:
        QuantumCircuit: A quantum circuit with the specified number of qubits
        and depth.
    """
    circuit = QuantumCircuit(qubits)
    thetas: List[Parameter] = []
    
    num_thetas = (
        qubits * (depth + 1) 
        + int(depth * qubits * (qubits - 1) / 2) 
        + 100
    )
    for l in range(num_thetas):
        thetas.append(Parameter('θ_' + str(chr(l))))
    
    counter = 0
    for i in range(depth):
        for l in range(qubits):
            if l + counter < len(thetas):
                circuit.h(l)
                circuit.ry(thetas[l + counter + qubits], l)
        
        counter_3 = 0
        for p in range(qubits):
            counter_2 = 0
            for j in range(qubits):
                if p != p + j and p + j <= qubits - 1:
                    circuit.cry(
                        thetas[j + p + counter + 2 * qubits + counter_3 - 1],
                        p,
                        p + j
                    )
                    counter_2 += 1
            counter_3 += counter_2 - 1
            if p == qubits - 1:
                circuit.cx(p, 0)
            circuit.h(p)
        
        counter += 2 * qubits + counter_3
        circuit.barrier()
    
    for j in range(qubits):
        circuit.ry(thetas[j], j)
    
    for k in range(qubits):
        if k + int(qubits / 2) < qubits:
            circuit.cx(k, k + int(qubits / 2))
    
    return circuit


def ansatz_cry_optimized(
        qubits: int, 
        depth: int, 
        drawmode: bool = False) -> QuantumCircuit:
    """
    This Ansatz is currently used to generate the full entanglement ansatz with 
    control R_y gates. Mainly used for imaginary variational time evolution of 
    diagonal Hamiltonians.
    Args:
        qubits (int): Number of qubits.
        depth (int): Depth of the circuit.
        drawmode (bool, optional): If True, parameters theta_i are represented
        with numbers (e.g., theta_3). In the code mode, the chr values of 
        numbers are used since they are sorted alphabetically.
        Defaults to False.
    Returns:
        QuantumCircuit: A quantum circuit with the specified number of qubits 
        and depth.
    """
    circuit = QuantumCircuit(qubits)
    thetas: List[Parameter] = []

    num_thetas = int(qubits * (qubits - 1) / 2) * depth * qubits
    if drawmode:
        thetas = [Parameter('θ_' + str(l)) for l in range(num_thetas)]
    else:
        thetas = [Parameter('θ_' + str(chr(l))) for l in range(num_thetas)]

    counter = 0
    for i in range(depth):
        for l in range(qubits):
            if l + counter < len(thetas):
                if i == 0:
                    circuit.h(l)
                circuit.ry(thetas[l + counter + qubits], l)
        
        counter_3 = 0
        for p in range(qubits):
            counter_2 = 0
            for j in range(qubits):
                if p != p + j and p + j <= qubits - 1:
                    circuit.cry(
                        thetas[j + p + counter + 2 * qubits + counter_3 - 1],
                        p,
                        p + j
                    )
                    counter_2 += 1
            counter_3 += counter_2 - 1
            if p == qubits - 1:
                circuit.cx(p, 0)
        
        counter += (2 * qubits) + counter_3
        circuit.barrier()
    
    return circuit


if __name__ == "__main__":
    ansatz_review_exact(4, 2)
