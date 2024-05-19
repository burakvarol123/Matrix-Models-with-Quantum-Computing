#!/lustre/fs24/group/cqta/bvarol/package/miniconda3/bin/python3.11

"""one_matrix_max_depth.py
This script creates both varqite and scipy distributions for given parameter
values. 
"""
import sys
module_directory_path = "/lustre/fs24/group/cqta/bvarol/workspace/matrixmodels"
sys.path.append(module_directory_path)

import numpy as np
from create_distribution import (create_distribution_varqite,
                                 create_distribution_scipy)
from copy import copy
import sys
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
        "a": float,
        "g": float,
        "power": int,
        "power_interaction": int,
        "num_timesteps":int,
        'out': str

    }
    if len(filenames) == 0:
        print("At least one filename has to be specified")
    for filename in filenames:
        parameters = copy(parameter_template)
        parameters = ut.read_parameters_from_file(filename, parameters)
        print(parameters["power_interaction"])
        results = []
        qubits = (parameters['D']+1)*parameters['qubits_per_dim']
        ansatz = an.ansatz_cry_optimized(qubits, parameters["depth"])
        squareterm = om.matrix_terms(
            qubits_per_dim=parameters['qubits_per_dim'], 
            dimension=parameters['D'], 
            lattice_spacing=parameters['a'], 
            pow = 2)
        interaction = om.matrix_terms(
            parameters['qubits_per_dim'], 
            parameters['D'],
            parameters['a'], 
            parameters['power_interaction']) 
        hamiltonian = squareterm - parameters["g"] * interaction
        evolution_result = create_distribution_varqite(
            qubits, 
            parameters['depth'],
            hamiltonian, 
            parameters['beta'],
            ansatz, 
            parameters["num_timesteps"])
        results.append(evolution_result)
        scipy_result = create_distribution_scipy(
            qubits, 
            parameters['depth'],
            hamiltonian, 
            parameters['beta'],
            ansatz, 
            parameters["num_timesteps"])
        results.append(scipy_result)
        ut.write_parameters_to_file(filename, parameters)
        np.save(filename + '_evolution' , results)
