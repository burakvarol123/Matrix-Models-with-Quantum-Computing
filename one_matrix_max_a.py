#!/lustre/fs24/group/cqta/bvarol/package/miniconda3/bin/python3.11

"""one_matrix_max_a.py
This script creates both varqite and scipy distributions for given parameter
values. It is designed such that a max lattice spacing is given and iterated 
over every lattice spacings beginning from 0.1 to a_max, in steps of 0.1.
"""
import sys
module_directory_path = "/lustre/fs24/group/cqta/bvarol/workspace/matrixmodels"
sys.path.append(module_directory_path)

import numpy as np
from create_distribution import (create_distribution_varqite, 
                                 create_distribution_scipy)
from copy import copy
import utility as ut
import one_matrix_model as om
import Ansatz as an


if __name__ == "__main__":  
    filenames = sorted(sys.argv[1:])
    parameter_template = {
        'depth': int,
        'D': int,
        "qubits_per_dim": int,
        'beta': float,
        "a_max": float,
        "g": float,
        "power": int,
        "power_interaction": int,
        'out': str,
        "num_timesteps": int

    }
    if len(filenames) == 0:
        print("At least one filename has to be specified")
    for filename in filenames:
        parameters = copy(parameter_template)
        parameters = ut.read_parameters_from_file(filename, parameters)
        print(parameters["power_interaction"])
        varqite_results = []
        scipy_results = []
        for a in np.arange(0.1,parameters["a_max"], 0.05):
            qubits = (parameters["D"]+1)*(parameters["qubits_per_dim"])
            ansatz = an.ansatz_cry_optimized(qubits, parameters["depth"])
            squareterm = om.matrix_terms(
                qubits_per_dim=parameters["qubits_per_dim"], 
                dimension=parameters["D"], 
                lattice_spacing=a, 
                pow=2)
            interaction = om.matrix_terms(
                parameters["qubits_per_dim"], 
                parameters["D"], 
                a, 
                parameters['power_interaction']) 
            hamiltonian = squareterm - parameters["g"] * interaction
            evolution_result = create_distribution_varqite(
                qubits, 
                parameters['depth'], 
                hamiltonian, 
                parameters['beta'], 
                ansatz, 
                parameters["num_timesteps"])
            scipy_result = create_distribution_scipy(
                qubits, 
                parameters['depth'], 
                hamiltonian, 
                parameters['beta'], 
                ansatz, 
                parameters["num_timesteps"])
            varqite_results.append(evolution_result)
            scipy_results.append(scipy_result)
            print("completed cycle, qubits =", qubits)
        results = varqite_results + scipy_results  
        ut.write_parameters_to_file(filename, parameters)
        np.save(filename + '_evolution', results)
