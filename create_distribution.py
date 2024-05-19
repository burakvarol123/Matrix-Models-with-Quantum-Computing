"""
create_distribution.py

This module provides functions to perform and compare quantum imaginary time 
evolution using Qiskit. It includes functionalities for creating distributions 
using Variational Quantum Imaginary Time Evolution (varQITE) and SciPy's 
imaginary time evolver, as well as plotting and comparing the results.

Functions:
- create_distribution_varqite(qubits, depth, hamilton, beta, ansatz, 
  num_timesteps): Evolves the given initial state using varQITE and returns 
  the result.
- create_distribution_scipy(qubits, depth, hamilton, beta, ansatz, 
  num_timesteps): Evolves the given initial state using SciPy's imaginary time 
  evolver and returns the result.
- plot_distribution(evolution_result, sol, qubits, depth, a, g, save, 
  filename): Visualizes the distribution for given varQITE and SciPy results.
- quick_compare(evolution_result, sol, qubits, dimension, depth, a, 
  num_timesteps, g, save, filename): Compares the evolution of the energy 
  between varQITE and SciPy evolved distributions.

Each function facilitates different aspects of quantum time evolution 
simulations and their visualization, aiding in the analysis and comparison 
of results from different methods.

Examples:
    # Create a distribution using varQITE
    result_varqite = create_distribution_varqite(
        qubits=4, depth=2, hamilton=hamiltonian, beta=1.0, 
        ansatz=quantum_circuit, num_timesteps=10
    )

    # Create a distribution using SciPy's imaginary time evolver
    result_scipy = create_distribution_scipy(
        qubits=4, depth=2, hamilton=hamiltonian, beta=1.0, 
        ansatz=quantum_circuit, num_timesteps=10
    )

    # Plot the distribution comparison between varQITE and SciPy results
    plot_distribution(result_varqite, result_scipy, qubits=4, depth=2, a=0.5, 
                      g=1.0)

    # Quickly compare the energy evolution between varQITE and SciPy results
    quick_compare(result_varqite, result_scipy, qubits=4, dimension=1, 
                  depth=2, a=0.5, num_timesteps=10, g=1.0)

Dependencies:
- numpy
- os
- matplotlib.pyplot
- pylab
- typing
- qiskit
- utility.py
"""
#!/lustre/fs24/group/cqta/bvarol/package/miniconda3/bin/python3.11
import sys
module_directory_path = "/lustre/fs24/group/cqta/bvarol/workspace/matrixmodels"
sys.path.append(module_directory_path)

import numpy as np
import os
import matplotlib.pyplot as plt
import pylab
from typing import Any, List
from qiskit.algorithms import (SciPyImaginaryEvolver, TimeEvolutionProblem, 
                               VarQITE, TimeEvolutionResult)
from qiskit.algorithms.gradients import ReverseEstimatorGradient, ReverseQGT
from qiskit.algorithms.time_evolvers.variational import (
                            ImaginaryMcLachlanPrinciple, VariationalPrinciple)
from qiskit.circuit import QuantumCircuit
from qiskit.primitives import Estimator
from qiskit.quantum_info import Statevector
from qiskit.opflow import OperatorBase
from utility import generate_combinations


def create_distribution_varqite(
        qubits: int, 
        depth: int, 
        hamilton: OperatorBase,
        beta: float, 
        ansatz: QuantumCircuit,
        num_timesteps: int) -> TimeEvolutionResult: 
    """ 
    Evolves the given initial state according to Variational Quantum 
    Imaginary Time Evolution. Creates rho = exp(-beta * H)

    Args:
        qubits (int): The total numbers of qubits of the evolved system.
        depth (int): depth
        hamilton (OperatorBase): The Hamilton operator which will lead the 
        evolution. 
        beta (float): Inverse Temparature
        ansatz (QuantumCircuit): The initial ansatz for the parameteritzed 
        circuit
        num_timesteps (int): Number of timesteps that the varQITE will use.
        Higher timesteps will result in more pression but will be slower.

    Returns:
        TimeEvolutionResult: This object has everything one needs to reproduce 
        the expectation value and produce expectation values.
        """
    hamiltonian: OperatorBase = hamilton  
    var_principle: VariationalPrinciple = ImaginaryMcLachlanPrinciple(
        qgt=ReverseQGT(),
        gradient=ReverseEstimatorGradient()
    )
    time: float = (1 / 2.0) * beta
    aux_ops: List[OperatorBase] = [hamiltonian]
    init_param_values: List[float] = [0] * (
        depth * qubits + int(depth * qubits * (qubits - 1) / 2)
    )
    
    evolution_problem: Any = TimeEvolutionProblem(
        hamiltonian, 
        time, 
        aux_operators=aux_ops
    )
    var_qite: VarQITE = VarQITE(
        ansatz, 
        init_param_values, 
        var_principle, 
        Estimator(),
        num_timesteps=num_timesteps
    )
    evolution_result: TimeEvolutionResult = var_qite.evolve(evolution_problem)
    return evolution_result


def create_distribution_scipy(qubits: int, depth: int, hamilton: OperatorBase,
                              beta: float, ansatz: QuantumCircuit,
                              num_timesteps: int) -> TimeEvolutionResult:
    """
    Evolves the given initial state by scipy's imaginary time evolver.
    Creates rho = exp(-beta * H)
    Args:
        qubits (int): The total numbers of qubits of the evolved system.
        depth (int): depth
        hamilton (OperatorBase): The Hamilton operator which will lead the 
        evolution. 
        beta (float): Inverse Temparature
        ansatz (QuantumCircuit): The initial ansatz for the parameteritzed 
        circuit
        num_timesteps (int): Number of timesteps that the varQITE will use.
        Higher timesteps will result in more pression but will be slower.

    Returns:
        TimeEvolutionResult: This object has everything one needs to reproduce 
        the expectation value and produce expectation values.
    """
    hamiltonian: OperatorBase = hamilton  
    time: float = beta / 2.0
    aux_ops: List[OperatorBase] = [hamiltonian]
    init_param_values: List[float] = [0] * (
        depth * qubits + int(depth * qubits * (qubits - 1) / 2)
    )
    
    init_state: Statevector = Statevector(
        ansatz.assign_parameters(init_param_values))
    
    evolution_problem: Any = TimeEvolutionProblem(
        hamiltonian, 
        time, 
        initial_state=init_state, 
        aux_operators=aux_ops
    )
    
    exact_evol: SciPyImaginaryEvolver = SciPyImaginaryEvolver(
        num_timesteps=int(beta * num_timesteps)
    )
    sol: TimeEvolutionResult = exact_evol.evolve(evolution_problem)
    return sol


def plot_distribution(evolution_result: TimeEvolutionResult, 
                      sol: TimeEvolutionResult,
                      qubits: int, 
                      depth: int, 
                      a: float, 
                      g: float,
                      save: str = None, 
                      filename: str = None) -> None:
    """
    Visualises the  distribution for given varQITE and scipy results. Only 
    1*1 Matrix Models (just normal integrals) can be properly visualised.

    Args:
        evolution_result (TimeEvolutionResult): The distribution one would get 
        from varQITE.
        sol (TimeEvolutionResult): The distribution one would get from scipy.
        qubits (int): Number of Qubits in the system.
        depth (int): Depth of the circuit.
        a (float): Lattice Spacing.
        g (float): Interaction Constant.
        save (str): The directory to save the visuals. Defaults to None.
        filename (str): The name of the saved visual. Defaults to None.
    """
    
    varqite_result = Statevector(
                        evolution_result.evolved_state).probabilities_dict()
    scipy_result = Statevector(sol.evolved_state).probabilities_dict()
    x = generate_combinations(qubits - 1)
    y = []
    y_2 = []
    for i in np.arange(len(x)):
        y.append(varqite_result[x[i]])
        y_2.append(scipy_result[x[i]])
    plt.bar(x, y, color="red", label="varqite")
    plt.bar(x, y_2, color="blue", label='scipy')
    plt.xticks(ticks=range(len(x)), labels=x, rotation=90)
    plt.xlabel("Position")
    plt.ylabel("Probability")
    plt.title("Scipy vs VarQITE: " + "Qubits= " + str(qubits) + " a= " + str(a) 
              + " g= " + str(g) + " depth= " + str(depth))
    plt.legend()
    
    if save:
        save_path = os.path.join(save, filename + ".png")
        plt.savefig(save_path)
        plt.show()


def quick_compare(evolution_result: TimeEvolutionResult, 
                  sol: TimeEvolutionResult,
                  qubits: int, 
                  dimension: int, 
                  depth: int, 
                  a: float, 
                  num_timesteps: int,
                  g: float, 
                  save: str = None, 
                  filename: str = None) -> None:
    """
    This compares the evolution of the Energy <H> between the 
    varQITE and scipy evolved distributions. The desired varQITE evolution 
    would be close to the scipy one.

    Args:
       evolution_result (TimeEvolutionResult): The distribution one would get 
        from varQITE.
        sol (TimeEvolutionResult): The distribution one would get from scipy.
        qubits (int): Number of Qubits in the system.
        dimension (int): dimension of the matrix -1 . Dimension = 0 => 1*1 
        matrix.
        depth (int): Depth of the circuit.
        a (float): Lattice Spacing.
        g (float): Interaction Constant.
        save (str): The directory to save the visuals. Defaults to None.
        filename (str): The name of the saved visual. Defaults to None.
        """
    qubits_per_dim = int(qubits / (dimension + 1))
    times = evolution_result.times[:(num_timesteps + 1)]
    h_exp_val = np.array([ele[0][0] for ele in evolution_result.observables])
    h_exp_val = h_exp_val[:(num_timesteps + 1)]
    exact_h_exp_val = sol.observables[0][0].real[:(num_timesteps + 1)]
    errors = []
    for i in range(len(h_exp_val)):
        error = h_exp_val[i] - exact_h_exp_val[i]
        errors.append(error)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(times, h_exp_val, label="VarQITE")
    ax1.plot(times, exact_h_exp_val, label="Exact", linestyle='--', color="m")
    ax2.plot(times, errors, linestyle="dotted", label="Error", color="r")
    ax1.set_xlabel('Time')
    ax1.set_ylabel(r"$\langle H \rangle$ (energy)")
    ax2.set_ylabel("Difference")
    plt.title("Qubits per Dimension = " + str(qubits_per_dim) + " Dimension =" 
              + str(dimension) + " Depth =" + str(depth) + " a =" + 
              str(round(a, 3)) + " g =" + str(round(g, 3)))
    ax1.legend(loc=2)
    ax2.legend(loc=1)
    plt.grid()
    if save:
        save_path = os.path.join(save, filename + ".png")
        plt.savefig(save_path)


if __name__ == "__main__":
    """
    two_qubits = create_distribution_varqite(
        3, 1, np.kron(hm.Ham(3,1,2),hm.krId(3)), 1)
    two_qubits_exp = create_distribution_scipy(3,1,hm.Ham(3,1,2),1)
    evolution_result = two_qubits[3]
    h_exp_val = two_qubits[2]
    exact_h_exp_val = two_qubits_exp[0]
    sol = two_qubits_exp[1]
    compare_results(evolution_result, sol, h_exp_val , 
    exact_h_exp_val, qubits = 3 , beta =)

qubits = 3
betas = [0.1, 0.5 , 0.75 , 1.0, 1.25 , 1.5]

for i in betas:
    two_qubits = create_distribution_varqite(qubits, 1, 
    np.kron(hm.Ham(qubits,1,2),hm.krId(3)), i)
    two_qubits_exp = create_distribution_scipy(qubits,1,hm.Ham(qubits,1,2), i)
    evolution_result = two_qubits[3]
    h_exp_val = two_qubits[2]
    exact_h_exp_val = two_qubits_exp[0]
    sol = two_qubits_exp[1]
    compare_results(evolution_result, sol, h_exp_val , exact_h_exp_val, qubits 
    , i)
    """
"""
lambdas = []
for i in range(2):
    lambdas.append(om.create_lambda_2(i, 1, 2))
hamiltonian = 0
for i in range(2):
    hamiltonian += (1 ** 2) * lambdas[i].power(2)

file = "/Users/salsa/MatrixModels/matrixmodels/remote/results/scipy_D_1_N_2_depth_1_evolution "
np.save(file, create_distribution_scipy(4, 1, hamiltonian, 0.7))
"""
