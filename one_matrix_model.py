"""
one_matrix_model.py

This module provides functions to generate operators related to random 
matrices and their eigenvalues, for use in quantum algorithms and simulations. 
It includes methods for creating lambda operators, the Vandermonde determinant, 
and matrix terms parameterized by eigenvalues.

Functions:
- create_lambda(index, dimension, qubits_per_dim): Creates a discretized 
  eigenvalue operator for a random matrix and tensor products it with identity 
  operators to match the total system dimension.
- vandermonde(dimension, qubits_per_dim): Computes the Vandermonde determinant 
  for a given dimension and number of qubits per dimension.
- matrix_terms(qubits_per_dim, dimension, lattice_spacing, pow): Generates the 
  parameterized matrix terms, M**power, using eigenvalues.

Each function returns an Operator object configured according to the specified 
parameters.

Examples:
    # Create a lambda operator for a specified index, dimension, and qubits 
    # per dimension
    lambda_op = create_lambda(index=0, dimension=2, qubits_per_dim=3)

    # Compute the Vandermonde determinant for a given dimension and qubits per 
    # dimension
    vandermonde_det = vandermonde(dimension=2, qubits_per_dim=3)

    # Generate parameterized matrix terms with eigenvalues
    matrix_term_op = matrix_terms(
        qubits_per_dim=3, 
        dimension=2, 
        lattice_spacing=0.1, 
        pow=2
    )

Dependencies:
- qiskit
- discretisation.py (imported as dc)
"""


import discretisation as dc
from qiskit.quantum_info import Operator
from typing import List


def create_lambda(
    index: int,
    dimension: int,
    qubits_per_dim: int
) -> Operator:
    """
    For a (dimension+1) * (dimension+1) random matrix, creates the discretised 
    eigenvalues and tenor productis it with identity to match the total 
    dimension of the system.

    Args:
        index (int): the index of the eigenvalue.
        dimension (int):Stands for (dimension+1) * (dimension +1) random matrix
        qubits_per_dim (int): How many qubits to be used for a discretised 
        eigenvalue.

    Returns:
        Operator: the eigenvalute operator as pauli string.
    """
    if index == 0:
        matrix = dc.position_x(qubits_per_dim)
        for _ in range(dimension):
            matrix = matrix ^ dc.insert_i(qubits_per_dim)
    else:
        matrix = dc.insert_i(qubits_per_dim)
        for j in range(dimension):
            if index == j + 1:
                matrix = matrix ^ dc.position_x(qubits_per_dim)
            else:
                matrix = matrix ^ dc.insert_i(qubits_per_dim)
    return matrix


def vandermonde(
    dimension: int,
    qubits_per_dim: int
) -> Operator:
    """ The Vandermonde determninant, arises when we parametrisize the random 
    matrix integral with eigenvalues and integrate out the unitary matrices.

    Args:
       dimension (int):Stands for (dimension+1) * (dimension +1) random matrix
        qubits_per_dim (int): How many qubits to be used for a discretised 
        eigenvalue.

    Returns:
        Operator: Vandermonde deterinant as Pauli string.
    """
    lambdas: List[Operator] = []
    for i in range(dimension + 1):
        lambdas.append(create_lambda(i, dimension, qubits_per_dim))
    
    determinant = dc.insert_i(qubits_per_dim * (dimension + 1))
    for i in range(dimension + 1):
        for j in range(i + 1, dimension + 1):
            determinant = determinant @ (lambdas[j] - lambdas[i])

    return determinant


def matrix_terms(
    qubits_per_dim: int,
    dimension: int,
    lattice_spacing: float,
    pow: int
) -> Operator:
    """Creates the eigenvalue parameterised matrix terms, M**power

    Args:
        dimension (int):Stands for (dimension+1) * (dimension +1) random matrix
        qubits_per_dim (int): How many qubits to be used for a discretised 
        eigenvalue.
        lattice_spacing (float): Lattice spacing for discrete integral
        pow (int): _description_

    Returns:
        Operator: Parameterised Random Matrices as Pauli strings.
    """
    lambdas: List[Operator] = []
    for i in range(dimension + 1):
        lambdas.append(create_lambda(i, dimension, qubits_per_dim))
    
    hamiltonian = 0
    for i in range(dimension + 1):
        hamiltonian += (lattice_spacing ** pow) * lambdas[i].power(pow)
    return hamiltonian


if __name__ == "__main__":
    create_lambda(6,0,3)

