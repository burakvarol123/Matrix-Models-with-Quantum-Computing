import numpy as np
import Hamilton as hm
from create_distribution import create_distribution_varqite
from copy import copy
import sys
from time import process_time
from one_matrix_model import create_lambda


def read_parameters_from_file(path: str, params: dict):
    '''
    Reads arbitrary parameters into a dictionary. The input file
    has to be of the format:

    key1=value1
    key2=value2

    path: The path of the file to read
    params: A key:constructor parameter template dictionary.
            Only variables specified in params  will be read,
            they will be read in by the constructor specified in params
    returns: A dictionary of variables read from path
    ...
    '''
    with open(path, "r") as fd:
        content = fd.read()
        print(content)
        content = content.replace('\n', '')
        lines = content.split(';')
        for line in lines:
            if len(line) == 0:
                continue
            key, value = line.split('=')
            # print(key, value)
            if key in params.keys():
                params[key] = params[key](value)
    return params


def write_parameters_to_file(path: str, params: dict):
    '''
    Writes arbitrary parameters in a dictionary to a file.
    The resulting file will be of the format

    key1=value1
    key2=value2
    ...

    And is therefor parsable by 'read_parameters_from_file'
    path: The path of the file to read
    params: A key:constructor parameter template dictionary.
            Only variables specified in params  will be read,
            they will be read in by the constructor specified in params
    '''
    with open(path, 'w') as fd:
        for key, value in params.items():
            fd.write(f'{key}={value};\n')

if __name__ == "__main__":
    filenames = sorted(sys.argv[1:])
    parameter_template = {
        'depth': int,
        'D': int,
        "qubits_per_dim": int,
        'beta': float,
        "a": float,
        "power": int,
        'out': str

    }
    if len(filenames) == 0:
        print("At least one filename has to be specified")
    for filename in filenames:
        parameters = copy(parameter_template)
        parameters = read_parameters_from_file(filename, parameters)
        qubits = (parameters['D']+1)*parameters['qubits_per_dim']
        lambdas = []
        for i in range(parameters['D']+1):
            lambdas.append(create_lambda(i, parameters['D'], parameters['qubits_per_dim']))
        hamiltonian = 0
        for i in range(parameters['D']+1):
            hamiltonian += 1 ** 2 * np.linalg.matrix_power(lambdas[i], parameters['power'])

        hamiltonian = hm.Ham( qubits, parameters['a'], parameters['power'])
        start = process_time()
        varqite = create_distribution_varqite(qubits, parameters['depth'], hamiltonian, parameters['beta'])
        stop = process_time()
       
        parameters['thetas'] = varqite[1]
        parameters['h_exp_val'] = varqite[2]
        evolution_result = varqite[3]
        parameters['Process time'] = stop - start
        write_parameters_to_file(filename, parameters)
        np.save(filename + '_evolution', evolution_result)
