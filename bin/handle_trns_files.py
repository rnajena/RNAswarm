import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import itertools
import sys

def parse_genome(genome_file):
    genome_dict = {}
    with open(genome_file) as f:  # read genome file
        for line in f:  # parse genome file
            if line.startswith(">"):  # parse header
                name = line.strip().split(" ")[0][1:]  # remove '>' and ' '
                genome_dict[name] = ""  # create new entry in genome dict
            elif line != "\n":  # parse sequence and skips newlines
                genome_dict[name] = line.strip()  # append sequence to entry
    return genome_dict


def make_combination_array(genome_dict):
    """ "
    Creates a dictionarry of numpy array of all possible genome segement combinations.
    Use parse_genome() to create genome_dict.
    """
    combination_arrays = {}
    segments = list(genome_dict.keys())
    segment_combinations = [
        segment_combination
        for segment_combination in itertools.combinations_with_replacement(segments, 2)
    ]
    for segment_combination in segment_combinations:
        combination_arrays[segment_combination] = np.zeros(
            (
                len(genome_dict[segment_combination[0]]),
                len(genome_dict[segment_combination[1]]),
            )
        )
    return combination_arrays


def parse_trns_file(trns_file):
    """
    Parse segemehl .trns.txt file and return a dictionary with mappings and values.
    """
    trns_dict = {}
    with open(trns_file) as f:
        for line in f:
            line = line.strip().split("\t")
            line = line[0].split(",") + line[1].split(",") + line[2].split(",")
            trns_dict[line[-1]] = {
                "mapping01": {
                    "ref-chr": line[0],
                    "ref-pos": int(line[1]),
                    "ref-strand": line[2],
                    "start-in-read": int(line[3]),
                    "align-length": int(line[4]),
                    "align-edist": int(line[5]),
                    "score": int(line[6]),
                },
                "mapping02": {
                    "ref-chr": line[7],
                    "ref-pos": int(line[8]),
                    "ref-strand": line[9],
                    "start-in-read": int(line[10]),
                    "align-length": int(line[11]),
                    "align-edist": int(line[12]),
                    "score": int(line[13]),
                },
            }
    return trns_dict


def fill_combination_array(combination_arrays, trns_dict):
    """
    Fill combination arrays with values from trns_dict.
    """
    for read_id in trns_dict.keys():
        segment01_segment02 = tuple(
            [
                trns_dict[read_id]["mapping01"]["ref-chr"],
                trns_dict[read_id]["mapping02"]["ref-chr"],
            ]
        )
        segment02_segment01 = tuple(
            [
                trns_dict[read_id]["mapping02"]["ref-chr"],
                trns_dict[read_id]["mapping01"]["ref-chr"],
            ]
        )
        read01_direction = trns_dict[read_id]["mapping01"]["ref-strand"]
        read02_direction = trns_dict[read_id]["mapping02"]["ref-strand"]
        if segment01_segment02 in combination_arrays.keys():
            if read01_direction == "+" and read02_direction == "+":
                combination_arrays[segment01_segment02][
                    trns_dict[read_id]["mapping01"]["ref-pos"] : trns_dict[read_id][
                        "mapping01"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping01"]["align-length"],
                    trns_dict[read_id]["mapping02"]["ref-pos"] : trns_dict[read_id][
                        "mapping02"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping02"]["align-length"],
                ] += 1
            elif read01_direction == "-" and read02_direction == "-":
                combination_arrays[segment01_segment02][
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    - trns_dict[read_id]["mapping01"]["align-length"] : trns_dict[
                        read_id
                    ]["mapping01"]["ref-pos"],
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    - trns_dict[read_id]["mapping02"]["align-length"] : trns_dict[
                        read_id
                    ]["mapping02"]["ref-pos"],
                ] += 1
            elif read01_direction == "-" and read02_direction == "+":
                combination_arrays[segment01_segment02][
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    - trns_dict[read_id]["mapping01"]["align-length"] : trns_dict[
                        read_id
                    ]["mapping01"]["ref-pos"],
                    trns_dict[read_id]["mapping02"]["ref-pos"] : trns_dict[read_id][
                        "mapping02"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping02"]["align-length"],
                ] += 1
            elif read01_direction == "+" and read02_direction == "-":
                combination_arrays[segment01_segment02][
                    trns_dict[read_id]["mapping01"]["ref-pos"] : trns_dict[read_id][
                        "mapping01"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping01"]["align-length"],
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    - trns_dict[read_id]["mapping02"]["align-length"] : trns_dict[
                        read_id
                    ]["mapping02"]["ref-pos"],
                ] += 1
        elif segment02_segment01 in combination_arrays.keys():
            if read01_direction == "+" and read02_direction == "+":
                combination_arrays[segment02_segment01][
                    trns_dict[read_id]["mapping02"]["ref-pos"] : trns_dict[read_id][
                        "mapping02"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping02"]["align-length"],
                    trns_dict[read_id]["mapping01"]["ref-pos"] : trns_dict[read_id][
                        "mapping01"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping01"]["align-length"],
                ] += 1
            elif read01_direction == "-" and read02_direction == "-":
                combination_arrays[segment02_segment01][
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    - trns_dict[read_id]["mapping02"]["align-length"] : trns_dict[
                        read_id
                    ]["mapping02"]["ref-pos"],
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    - trns_dict[read_id]["mapping01"]["align-length"] : trns_dict[
                        read_id
                    ]["mapping01"]["ref-pos"],
                ] += 1
            elif read01_direction == "-" and read02_direction == "+":
                combination_arrays[segment02_segment01][
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    - trns_dict[read_id]["mapping02"]["align-length"] : trns_dict[
                        read_id
                    ]["mapping02"]["ref-pos"],
                    trns_dict[read_id]["mapping01"]["ref-pos"] : trns_dict[read_id][
                        "mapping01"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping01"]["align-length"],
                ] += 1
            elif read01_direction == "+" and read02_direction == "-":
                combination_arrays[segment02_segment01][
                    trns_dict[read_id]["mapping02"]["ref-pos"] : trns_dict[read_id][
                        "mapping02"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping02"]["align-length"],
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    - trns_dict[read_id]["mapping01"]["align-length"] : trns_dict[
                        read_id
                    ]["mapping01"]["ref-pos"],
                ] += 1

def plot_heatmap(array, output_dir, filename):
    heatmap = sns.heatmap(array)
    plt.figure()
    heatmap.figure.savefig(f"{output_dir}/{filename}")
    plt.close("all")

genome_file_path = sys.argv[1]
trns_file_path = sys.argv[2]

print(genome_file_path)
print(trns_file_path)
print(sys.argv[3])

genome_dict = parse_genome(genome_file_path)
print(genome_dict)
trns_dict = parse_trns_file(trns_file_path)
print(trns_dict)

interaction_arrays = make_combination_array(genome_dict)
print(interaction_arrays)

fill_combination_array(interaction_arrays, trns_dict)
print(interaction_arrays)

for combination, array in interaction_arrays.items():
    plot_heatmap(array, sys.argv[3], f'{combination[0]}_{combination[1]}')