

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

def generate_combinations(d):
    """
    creates the right order for the visualization of the distrubition

    """
    combinations = []
    combinations_negative = []
    
    def generate_nested_loops(counters, depth):
        if depth == d:
            concatenated_str = "".join(str(counter) for counter in counters)
            concatenated_str_2 = "".join(str(counter) for counter in counters)
            combinations.append("0" + concatenated_str[::-1])
            combinations_negative.append("1" + concatenated_str_2[::-1])
            return

        for i in range(2):
            counters[depth] = i
            generate_nested_loops(counters, depth + 1)

    counters = [0] * d
    generate_nested_loops(counters, 0)
    return combinations_negative[::-1] + combinations


