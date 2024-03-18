
from qiskit.circuit.library import *
from qiskit import *
from qiskit.quantum_info import SparsePauliOp

def insert_z(number_of_qubits,  index):
    if index == number_of_qubits:
        raise Exception("your index is too high, index starts from 0. So maximal Z insertion is num_qubits-1")
    if index == 0:
        matrix = SparsePauliOp(["Z"] ,[1])
        for j in range(number_of_qubits-1):
            matrix = matrix ^ SparsePauliOp(["I"],[1])
    else:
        matrix = SparsePauliOp(["I"] ,[1])
        for j in range(number_of_qubits-1):
            if index  == j +1:
                matrix = matrix ^ SparsePauliOp(["Z"],[1])
            else:
                matrix = matrix ^ SparsePauliOp(["I"],[1])
    return matrix


def insert_i(
    number_of_qubits
):
    
    matrix = SparsePauliOp(["I"] ,[1])
    for j in range(number_of_qubits-1):
        matrix = matrix ^ SparsePauliOp(["I"] ,[1])
    return matrix


def k_total(
    number_of_qubits
):
    def k_m(
        number_of_qubits,
        index):
        k = (insert_i(number_of_qubits) - insert_z(number_of_qubits, index)) / 2.0
        return k

    sum_ = 0
    for i in range(number_of_qubits):
        sum_ += 2 ** i * k_m(number_of_qubits, i)
    return sum_


def position_x(
        number_of_qubits
        ):
    id = SparsePauliOp(["Z"], ["1"] ) ^ insert_i(number_of_qubits- 1) 
    id_2 = SparsePauliOp(["Z"], ["1"] ) ^ k_total(number_of_qubits - 1)
    return 0.5 * id +id_2

def hamiltonian(
        number_of_qubits,
        lattice_spacing,
        pow
        ):
    x = position_x(number_of_qubits)
    mat = lattice_spacing ** 2 * x.power(pow) * 1/10
    return mat

if __name__ == "__main__":
   matrix = hamiltonian(3,1,2)
   print(matrix.to_matrix())
   
 
   
