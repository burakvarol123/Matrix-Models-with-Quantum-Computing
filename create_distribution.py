import numpy as np
import Hamilton as hm
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


def create_distribution_varqite(qubits, depth, hamilton, beta):
    """
    creates the distribution e**-beta H.
    :param qubits: number of qubits. Notice that in the algorithm 2*qubits will be used because of varqite reasons
    :param depth: depth of the circuit
    :param hamilton: hamiltonian that we want to have. use the hm.Ham functions here
    :param beta: 1/kT
    :return: the circuit with correct parameters, parameters, expectation value, the object "evolution result"
    """
    hamiltonian = Operator(np.kron(hamilton, hm.krId(qubits)))
    var_principle = ImaginaryMcLachlanPrinciple()
    time = beta / 2.0
    aux_ops = [hamiltonian]
    init_param_values = [np.pi / 2] * qubits + [0] * (depth * qubits * 2 + qubits)
    ansatz = an.ansatz_review_exact(qubits * 2, depth)
    evolution_problem = TimeEvolutionProblem(hamiltonian, time, aux_operators=aux_ops)
    var_qite = VarQITE(ansatz, init_param_values, var_principle, Estimator())
    evolution_result = var_qite.evolve(evolution_problem)
    param = evolution_result.parameter_values[int(time * 100)][:]
    binded = an.ansatz_review_exact(qubits * 2, depth).bind_parameters(param)
    h_exp_val = np.array([ele[0][0] for ele in evolution_result.observables])

    return binded, param, h_exp_val, evolution_result


def create_distribution_scipy(qubits, depth, hamilton, beta):
    """
    creates the distribution via scipy, so not "quantumly"
    :returns scipy calculated exp val, and the evolution object
    """
    ansatz = an.ansatz_review_exact(qubits * 2, depth)
    hamiltonian = Operator(np.kron(hamilton, hm.krId(qubits)))
    time = beta / 2.0
    aux_ops = [hamiltonian]
    init_param_values = [np.pi / 2] * qubits + [0] * (depth * qubits * 2 + qubits)
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


def get_expval(binded, observable, qubits):
    """
    calculates the expextation value of the observable wrt. probability distribution
    :param binded: e**bH that we get from varqite
    :param observable: name of the observable
    :param qubits: num qubits
    :return: expectation value
    """
    estimator = Estimator()
    op = Operator(np.kron(observable, hm.krId(qubits)))
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
    op = Operator(np.kron(observable, hm.krId(qubits)))
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