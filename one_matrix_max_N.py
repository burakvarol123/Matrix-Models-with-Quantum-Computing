#!/lustre/fs24/group/cqta/bvarol/package/miniconda3/bin/python3.11

"""
This calculates the one matrix model for fixed D, a, g and varying N. Used mainly for benchmarking. The module directory part should come BEFORE imports.
"""
import sys
module_directory_path = "/lustre/fs24/group/cqta/bvarol/workspace/matrixmodels"
sys.path.append(module_directory_path)

import numpy as np
from create_distribution import create_distribution_2_varqite, create_distribution_scipy
from copy import copy
import utility as ut
import hamiltonianstryout as tr
import Ansatz as an


if __name__ == "__main__":  
    filenames = sorted(sys.argv[1:])
    parameter_template = {
        'depth': int,
        'D': int,
        "qubits_per_dim_max": int,
        'beta': float,
        "a": float,
        "g": float,
        "power": int,
        "power_interaction": int,
        'out': str

    }
    if len(filenames) == 0:
        print("At least one filename has to be specified")
    for filename in filenames:
        parameters = copy(parameter_template)
        parameters = ut.read_parameters_from_file(filename, parameters)
        print(parameters["power_interaction"])
        varqite_results = []
        scipy_results = []
        for N in np.arange(parameters["qubits_per_dim_max"]-2):
            qubits = (parameters['D']+1)*(N+3)
            ansatz = an.ansatz = an.ansatz_cry_optimized(N+3, parameters["depth"])
            squareterm = tr.matrix_terms(qubits_per_dim=N+3, dimension=parameters['D'], lattice_spacing=parameters['a'] , pow = 2)
            interaction = tr.matrix_terms(N+3, parameters['D'], parameters['a'], parameters['power_interaction']) 
            hamiltonian = squareterm - parameters["g"] * interaction
            evolution_result = create_distribution_2_varqite(qubits, parameters['depth'], hamiltonian, parameters['beta'], ansatz)
            scipy_result = create_distribution_scipy(qubits, parameters['depth'], hamiltonian, parameters['beta'], ansatz)
            varqite_results.append(evolution_result)
            scipy_results.append(scipy_result)
            print("completed cycle, qubits =", N+3)
        results = varqite_results + scipy_results  # first N terms are varqite, others scipy
        ut.write_parameters_to_file(filename, parameters)
        np.save(filename + '_evolution', results)
