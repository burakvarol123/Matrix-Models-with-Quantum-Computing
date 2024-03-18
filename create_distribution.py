#!/lustre/fs24/group/cqta/bvarol/package/miniconda3/bin/python3.11
import sys
module_directory_path = "/lustre/fs24/group/cqta/bvarol/workspace/matrixmodels"
sys.path.append(module_directory_path)


import hamiltonian_2 as hm2
import numpy as np
import Ansatz as an
from qiskit.primitives import Estimator
from qiskit.algorithms import VarQITE
from qiskit.algorithms.time_evolvers.variational import ImaginaryMcLachlanPrinciple
from qiskit.algorithms import TimeEvolutionProblem
from qiskit.quantum_info import *
import pylab
from qiskit.algorithms import SciPyImaginaryEvolver
from qiskit.quantum_info import Statevector
import matplotlib.pyplot as plt
from utility import generate_combinations
import os 
import one_matrix_model as om
from qiskit.algorithms.gradients import ReverseEstimatorGradient, ReverseQGT



def create_distribution_2_varqite(qubits, depth, hamilton, beta, param):
    hamiltonian = hamilton ^ hm2.insert_i(qubits)
    var_principle = ImaginaryMcLachlanPrinciple(qgt=ReverseQGT(), gradient=ReverseEstimatorGradient())
    time = (1 /2.0)* beta
    aux_ops = [hamiltonian]
    init_param_values = [np.pi / 2 +param] * qubits + [0] * (depth * 2 * qubits  + qubits + int(depth*2*qubits*(2*qubits-1)/2))
    ansatz = an.ansatz_review_exact_full_cry(qubits * 2, depth)
    evolution_problem = TimeEvolutionProblem(hamiltonian, time, aux_operators=aux_ops)
    var_qite = VarQITE(ansatz, init_param_values, var_principle, Estimator())
    evolution_result = var_qite.evolve(evolution_problem)
    return  evolution_result


def create_distribution_scipy(qubits, depth, hamilton, beta, param):
    """
    creates the distribution via scipy, so not "quantumly"
    :returns scipy calculated exp val, and the evolution object
    """
    ansatz = an.ansatz_review_exact(qubits * 2, depth)
    hamiltonian =  hamilton ^ hm2.insert_i(qubits)
    time = beta / 2.0
    aux_ops = [hamiltonian]
    init_param_values = [np.pi / 2 + param] * qubits + [0] * (depth * qubits * 2 + qubits)
    init_state = Statevector(ansatz.assign_parameters(init_param_values))
    evolution_problem = TimeEvolutionProblem(hamiltonian, time, initial_state=init_state, aux_operators=aux_ops)
    exact_evol = SciPyImaginaryEvolver(num_timesteps=int(time * 100))
    sol = exact_evol.evolve(evolution_problem)
    exact_h_exp_val = sol.observables[0][0].real
    return exact_h_exp_val, sol


def compare_results(evolution_result, sol, h_exp_val, exact_h_exp_val, qubits, beta):
    """
    compares the results of varqite with computer
    :param evolution_result: evolution result from varqite
    :param sol: evolution result from scipy
    :param h_exp_val: exp value from varqite
    :param exact_h_exp_val: exp value from scipy
    :return: plots expvals varqite vs scipy, distributions varqite vs scipy
    """
    times = evolution_result.times
    pylab.plot(times, h_exp_val, label="VarQITE")
    pylab.plot(times, exact_h_exp_val, label="Exact", linestyle='--')
    plt.axhline(y=0.25, color='r', linestyle='-')
    pylab.xlabel("Time")
    pylab.ylabel(r"$\langle H \rangle$ (energy)")
    pylab.legend(loc="upper right")
    plt.subplots()
    evolved_state = evolution_result.evolved_state
    varqite_result = Statevector(evolved_state).probabilities_dict(np.arange(qubits))
    scipy_result = Statevector(sol.evolved_state).probabilities_dict(np.arange(qubits))
    
    x = generate_combinations(qubits - 1)
    y = []
    y_2 = []
    for i in range(len(x)):
        y.append(varqite_result[x[i]])
    for i in range(len(x)):
        y_2.append(scipy_result[x[i]])
    
    plt.bar(x, y, width=0.9, color = "red", label = "varqite")
    plt.bar(x, y_2, width=0.9, color = "blue", label = 'scipy')
    plt.xlabel("Position")
    plt.ylabel("Probability")
    plt.title("Comparison Scipy vs VarQITE with" + "" + str(beta))
    plt.legend()
    plt.show()
    save_directory = '/Users/salsa/MatrixModels/matrixmodels/results_gaussian'
    file_name = "beta=" + str(beta) + ".png"
    save_path = os.path.join(save_directory, file_name)
    plt.savefig(save_path)


def compare_results_exp(evolution_result, exact_h_exp_val):
    """
    Compares the evolutions of ground state energies of scipy and varqite
    """
    times = evolution_result.times
    h_exp_val = np.array([ele[0][0] for ele in evolution_result.observables])
    pylab.plot(times, h_exp_val, label="VarQITE")
    pylab.plot(times, exact_h_exp_val, label="Exact", linestyle='--')
    pylab.xlabel("Time")
    pylab.ylabel(r"$\langle H \rangle$ (energy)")
    pylab.legend(loc="upper right")
   

def get_expval(binded, observable, qubits):
    """
    calculates the expextation value of the observable wrt. probability distribution
    :param binded: e**bH that we get from varqite
    :param observable: name of the observable
    :param qubits: num qubits
    :return: expectation value
    """
    estimator = Estimator()
    op = observable ^ hm2.insert_i(qubits)
    expectation_value = estimator.run(binded, op).result().values
    return expectation_value


def get_expval_shot(binded, observable, shots, qubits):
    """
    calculates the expextation value of the observable wrt. probability distribution
    :param binded: e**bH that we get from varqite
    :param observable: name of the observable
    :param qubits: num qubits
    :return: expectation value
    """
    estimator = Estimator()
    op = observable ^ hm2.insert_i(qubits)
    expectation_value = estimator.run(binded, op, shots=shots).result().values
    return expectation_value

if __name__ == "__main__":
    """
    two_qubits = create_distribution_varqite(3, 1, np.kron(hm.Ham(3,1,2),hm.krId(3)), 1)
    two_qubits_exp = create_distribution_scipy(3,1,hm.Ham(3,1,2),1)
    evolution_result = two_qubits[3]
    h_exp_val = two_qubits[2]
    exact_h_exp_val = two_qubits_exp[0]
    sol = two_qubits_exp[1]
    compare_results(evolution_result, sol, h_exp_val , exact_h_exp_val, qubits = 3 , beta =)

qubits = 3
betas = [0.1, 0.5 , 0.75 , 1.0, 1.25 , 1.5]

for i in betas:
    two_qubits = create_distribution_varqite(qubits, 1, np.kron(hm.Ham(qubits,1,2),hm.krId(3)), i)
    two_qubits_exp = create_distribution_scipy(qubits,1,hm.Ham(qubits,1,2), i)
    evolution_result = two_qubits[3]
    h_exp_val = two_qubits[2]
    exact_h_exp_val = two_qubits_exp[0]
    sol = two_qubits_exp[1]
    compare_results(evolution_result, sol, h_exp_val , exact_h_exp_val, qubits , i)
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
