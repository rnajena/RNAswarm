import numpy as np
import itertools

def parse_genome(genome_file):
    genome = {}
    with open(genome_file) as f: # read genome file
        for line in f: # parse genome file
            print(line)
            if line.startswith('>'): # parse header
                name = line.strip().split(' ')[0][1:] # remove '>' and ' '
                genome[name] = '' # create new entry in genome dict
            elif line != '\n': # parse sequence and skips newlines
                genome[name] = line.strip() # append sequence to entry
    return genome

def make_combination_array(genome_dict):
    """"
    Creates a dictionarry of numpy array of all possible genome segement combinations.
    Use parse_genome() to create genome_dict.
    """
    combination_arrays = {}
    segments = list(genome_dict.keys())
    segment_combinations = [segment_combination for segment_combination in itertools.combinations_with_replacement(segments, 2)]
    for segment_combination in segment_combinations:
        combination_arrays[segment_combination] = np.zeros((len(genome_dict[segment_combination[0]]), len(genome_dict[segment_combination[1]])))
    return combination_arrays

def parse_trns_file(trns_file):
    """
    Parse segemehl .trns.txt file and return a dictionary with mappings and values.
    """
    trns_dict = {}
    with open(trns_file) as f:
        for line in f:
            line = line.strip().split('\t')
            line = line[0].split(',')  + line[1].split(',') + line[2].split(',')
            trns_dict[line[-1]] = {
                'mapping01' : {
                    'ref-chr' : line[0],
                    'ref-pos' : line[1],
                    'ref-strand' : line[3],
                    'start-in-read' : line[4],
                    'align-length' : line[5],
                    'algin-edist' : line[6],
                    'score' : line[7]
                },
                'mapping02' : {
                    'ref-chr' : line[0],
                    'ref-pos' : line[1],
                    'ref-strand' : line[3],
                    'start-in-read' : line[4],
                    'align-length' : line[5],
                    'algin-edist' : line[6],
                    'score' : line[7]
                }
            }
    return trns_dict

def fill_combination_array(combination_arrays, genome_dict, trns_dict):
    """
    Fill combination arrays with values from trns_dict.
    """
    for combination in combination_arrays:
        for mapping in trns_dict:
            if trns_dict[mapping]['mapping01']['ref-chr'] == combination[0] and trns_dict[mapping]['mapping01']['ref-chr'] == combination[1]:
                combination_arrays[combination][int(trns_dict[mapping]['mapping01']['ref-pos'])-1][int(trns_dict[mapping]['mapping01']['start-in-read'])-1] = 1
            elif trns_dict[mapping]['mapping02']['ref-chr'] == combination[0] and trns_dict[mapping]['mapping02']['ref-chr'] == combination[1]:
                combination_arrays[combination][int(trns_dict[mapping]['mapping02']['ref-pos'])-1][int(trns_dict[mapping]['mapping02']['start-in-read'])-1] = 1
    return combination_arrays
