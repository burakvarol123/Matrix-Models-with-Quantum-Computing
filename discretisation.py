"""
discretisation.py 

This module provides the operators which will discretisize the real line using 
quantum computing. Details of the discretisation can be found in the report.

Functions:
- insert_z(qubits, index): Inserts a Z operator at the specified index.
- insert_i(qubits): Generates an identity matrix for the given number of qubits.
- k_total(qubits): Generates a non 1/2 pushed away position operator, used in 
  the position_x function.
- position_x(qubits): Generates the position operator for use in Hamiltonians.
- hamiltonian(qubits, lattice_spacing, power): Generates a Hamiltonian for 
  use in varQITE.

Each function returns a SparsePauliOp representing the corresponding operator 
as a Pauli string. These functions handle the construction of complex quantum 
operators required for advanced quantum algorithms and simulations.

Examples:
    # Insert a Z operator at index 3 for a 5-qubit system
    z_operator = insert_z(5, 3)

    # Generate an identity matrix for a 4-qubit system
    identity_operator = insert_i(4)

    # Generate a position operator for a 3-qubit system
    pos_operator = position_x(3)

    # Generate a Hamiltonian with specific lattice spacing and power
    hamiltonian_op = hamiltonian(4, 0.5, 2)

Exceptions:
- ValueError: Raised when an invalid index is provided to the insert_z function.

Dependencies:
- qiskit.quantum_info: The module requires Qiskit's quantum_info library for 
  SparsePauliOp.

"""
from qiskit.quantum_info import SparsePauliOp


def insert_z(
        qubits: int, 
        index: int) -> SparsePauliOp:
    """
    Inserts the Z operator to indexed place. Eg index = 3 => 1*1*1*Z*1..

    Args:
        qubits (int): Number of Qubits that system has.
        index (int): The place where Z has to be inserted.

    Raises:
        ValueError

    Returns:
        SparsePauliOp: The Z operator as Pauli String.
    """
    if index >= qubits:
        raise ValueError(
            "Index is too high; index starts from 0. So the maximal Z" 
            "insertion is num_qubits - 1."
        )
    if index == 0:
        matrix = SparsePauliOp(["Z"], [1])
        for _ in range(qubits - 1):
            matrix = matrix ^ SparsePauliOp(["I"], [1])
    else:
        matrix = SparsePauliOp(["I"], [1])
        for j in range(qubits - 1):
            if index == j + 1:
                matrix = matrix ^ SparsePauliOp(["Z"], [1])
            else:
                matrix = matrix ^ SparsePauliOp(["I"], [1])
    return matrix


def insert_i(
        qubits: int) -> SparsePauliOp:
    """
    ´Generates identity matrix in system size.

    Args:
        qubits (int): Number of qubits in the system

    Returns:
        SparsePauliOp: Returns identity as Pauli string
    """
    matrix = SparsePauliOp(["I"], [1])
    for _ in range(qubits - 1):
        matrix = matrix ^ SparsePauliOp(["I"], [1])
    return matrix


def k_total(
        qubits: int) -> SparsePauliOp:
    """
    The Non 1/2 pushed away position operator. Just used in position_x function

    Args:
        qubits (int): Number of Qubits in the system.

    Returns:
        SparsePauliOp: position operator as Pauli string.
    """
    def k_m(
            qubits: int, 
            index: int) -> SparsePauliOp:

        k = (insert_i(qubits) - insert_z(qubits, index)) / 2.0
        return k
    sum_ = 0
    for i in range(qubits):
        sum_ += 2 ** i * k_m(qubits, i)
    return sum_


def position_x(
        qubits: int) -> SparsePauliOp:
    """Generates the position operator which is used in hamiltonians.

    Args:
        qubits (int): Number of Qubits in the system.

    Returns:
        SparsePauliOp: Returns position operator as Pauli string.
    """
    id_1 = SparsePauliOp(["Z"], [1]) ^ insert_i(qubits - 1)
    id_2 = SparsePauliOp(["Z"], [1]) ^ k_total(qubits - 1)
    return 0.5 * id_1 + id_2


def hamiltonian(
        qubits: int, 
        lattice_spacing: float, 
        power: int) -> SparsePauliOp:
    """
    Hamiltonian which will be used in varqite. 

    Args:
        qubits (int): Number of Qubits in the system.
        lattice_spacing (float): Lattice Spacing a.
        power (int): Power of the term.

    Returns:
        SparsePauliOp: Hamiltonian as a Pauli String.
    """
    x = position_x(qubits)
    mat = (lattice_spacing ** power) * x.power(power)
    return mat

if __name__ == "__main__":
    k_total(3)