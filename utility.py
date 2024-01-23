
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


