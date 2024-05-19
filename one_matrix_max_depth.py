#!/lustre/fs24/group/cqta/bvarol/package/miniconda3/bin/python3.11

"""one_matrix_max_depth.py
This script creates both varqite and scipy distributions for given parameter
values. It is designed such that a max depth is given and iterated over
every dept beginning from 1 to depth_max
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
        'depth_max': int,
        'D': int,
        "qubits_per_dim": int,
        'beta': float,
        "a": float,
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
        for depth in np.arange(1, parameters["depth_max"]+1):
            qubits = (parameters['D']+1)*parameters["qubits_per_dim"]
            ansatz = an.ansatz = an.ansatz_cry_optimized(qubits, depth)
            squareterm = om.matrix_terms(
                qubits_per_dim=parameters["qubits_per_dim"], 
                dimension=parameters['D'], 
                lattice_spacing=parameters['a'], 
                pow = 2)
            interaction = om.matrix_terms(
                parameters["qubits_per_dim"], 
                parameters['D'], 
                parameters['a'], 
                parameters['power_interaction']) 
            hamiltonian = squareterm - parameters["g"] * interaction
            evolution_result = create_distribution_varqite(
                qubits,
                depth, 
                hamiltonian, 
                parameters['beta'], 
                ansatz, 
                parameters["num_timesteps"])
            scipy_result = create_distribution_scipy(
                qubits, 
                depth, 
                hamiltonian, 
                parameters['beta'], 
                ansatz, 
                parameters["num_timesteps"])
            varqite_results.append(evolution_result)
            scipy_results.append(scipy_result)
            print("completed cycle, qubits =", parameters["qubits_per_dim"])
        results = varqite_results + scipy_results 
        ut.write_parameters_to_file(filename, parameters)
        np.save(filename + '_evolution', results)
