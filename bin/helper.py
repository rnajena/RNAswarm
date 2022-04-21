def parse_genome(genome_file):
    """
    """
    genome_dict = {}

    header = ""
    seq = ""

    with open(genome_file) as f:  # read genome file
        for line in f:  # parse genome file
            if line.startswith(">"):  # parse header
            #* if there is a header already, we store the current sequence
            #* to this header.
                if header:
                    genome_dict[header] = seq
                    #* then we flush the sequence
                    seq = ""
                #* and parse the new header
                header = line.strip().split(" ")[0][1:]  # remove '>' and ' '
            #elif line != "\n":  # parse sequence and skips newlines
            else:
                #* if no header is detected, we append the new line
                #* to the current sequence
                seq += line.strip()
                #genome_dict[name] = line.strip()  # append sequence to entry
        #* after the last line, we have to store
        #* the last sequence as well. Since no new line with ">" occurs, we
        #* do this manually
        genome_dict[header] = seq
    return genome_dict

def make_combination_array(genome_dict):
    """
    Creates a dictionarry of numpy array of all possible genome segement combinations.
    Use helper.parse_genome() to create genome_dict.
    """
    combination_arrays = {}
    segments = list(genome_dict.keys())

    # segment_combinations = [
    #     segment_combination
    #     for segment_combination in itertools.combinations_with_replacement(segments, 2)
    # ]
    # * while I usually appreciate the usage of list comprehensions, you can directly transform
    # * the iterator to a list. Actually, we also could just put the iterator in the for loop.
    # * should work as well. Is a tad more memory efficient.
    segment_combinations = list(itertools.combinations_with_replacement(segments, 2))

    for segment_combination in segment_combinations:
        # for segment_combination in itertools.combinations_with_replacement(segments,2): # * this should work as well
        combination_arrays[segment_combination] = np.zeros(
            (
                len(genome_dict[segment_combination[0]]),
                len(genome_dict[segment_combination[1]]),
            )
        )
    return combination_arrays