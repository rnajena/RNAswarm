import numpy as np
import itertools


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
                    "algin-edist": int(line[5]),
                    "score": int(line[6]),
                },
                "mapping02": {
                    "ref-chr": line[7],
                    "ref-pos": int(line[8]),
                    "ref-strand": line[9],
                    "start-in-read": int(line[10]),
                    "align-length": int(line[11]),
                    "algin-edist": int(line[12]),
                    "score": int(line[13]),
                },
            }
    return trns_dict


def fill_combination_array(combination_arrays, trns_dict):
    """
    Fill combination arrays with values from trns_dict.
    """
    for read_id in trns_dict.keys():
        if (
            tuple([
                trns_dict[read_id]["mapping01"]["ref-chr"],
                trns_dict[read_id]["mapping02"]["ref-chr"],
            ])
            in combination_arrays.keys()
        ):
            if (
                trns_dict[read_id]["mapping01"]["ref-strand"] == "+"
                and trns_dict[read_id]["mapping02"]["ref-strand"] == "+"
            ):
                combination_arrays[
                    tuple([
                        trns_dict[read_id]["mapping01"]["ref-chr"],
                        trns_dict[read_id]["mapping02"]["ref-chr"],
                    ])
                ][
                    trns_dict[read_id]["mapping01"]["ref-pos"] : trns_dict[read_id][
                        "mapping01"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping01"]["align-length"],
                    trns_dict[read_id]["mapping02"]["ref-pos"] : trns_dict[read_id][
                        "mapping02"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping02"]["align-length"],
                ] += 1
            elif (
                trns_dict[read_id]["mapping01"]["ref-strand"] == "-"
                and trns_dict[read_id]["mapping02"]["ref-strand"] == "-"
            ):
                combination_arrays[
                    tuple([
                        trns_dict[read_id]["mapping01"]["ref-chr"],
                        trns_dict[read_id]["mapping02"]["ref-chr"],
                    ])
                ][
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    - trns_dict[read_id]["mapping01"]["align-length"]
                    : trns_dict[read_id]["mapping01"]["ref-pos"],
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    - trns_dict[read_id]["mapping02"]["align-length"]
                    : trns_dict[read_id]["mapping02"]["ref-pos"],
                ] += 1
            elif (
                trns_dict[read_id]["mapping01"]["ref-strand"] == "-"
                and trns_dict[read_id]["mapping02"]["ref-strand"] == "+"
            ):
                combination_arrays[
                    tuple([
                        trns_dict[read_id]["mapping01"]["ref-chr"],
                        trns_dict[read_id]["mapping02"]["ref-chr"],
                    ])
                ][
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    - trns_dict[read_id]["mapping01"]["align-length"]
                    : trns_dict[read_id]["mapping01"]["ref-pos"],
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    : trns_dict[read_id]["mapping02"]["ref-pos"]
                    + trns_dict[read_id]["mapping02"]["align-length"],
                ] += 1
            elif (
                trns_dict[read_id]["mapping01"]["ref-strand"] == "+"
                and trns_dict[read_id]["mapping02"]["ref-strand"] == "-"
            ):
                combination_arrays[
                    tuple([
                        trns_dict[read_id]["mapping01"]["ref-chr"],
                        trns_dict[read_id]["mapping02"]["ref-chr"],
                    ])
                ][
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    : trns_dict[read_id]["mapping01"]["ref-pos"]
                    + trns_dict[read_id]["mapping01"]["align-length"],
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    - trns_dict[read_id]["mapping02"]["align-length"]
                    : trns_dict[read_id]["mapping02"]["ref-pos"],
                ] += 1
        elif (
            tuple([
                trns_dict[read_id]["mapping02"]["ref-chr"],
                trns_dict[read_id]["mapping01"]["ref-chr"],
            ])
            in combination_arrays.keys()
        ):
            if (
                trns_dict[read_id]["mapping01"]["ref-strand"] == "+"
                and trns_dict[read_id]["mapping02"]["ref-strand"] == "+"
            ):
                combination_arrays[
                    tuple([
                        trns_dict[read_id]["mapping02"]["ref-chr"],
                        trns_dict[read_id]["mapping01"]["ref-chr"],
                    ])
                ][
                    trns_dict[read_id]["mapping02"]["ref-pos"] : trns_dict[read_id][
                        "mapping02"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping02"]["align-length"],
                    trns_dict[read_id]["mapping01"]["ref-pos"] : trns_dict[read_id][
                        "mapping01"
                    ]["ref-pos"]
                    + trns_dict[read_id]["mapping01"]["align-length"],
                ] += 1
            elif (
                trns_dict[read_id]["mapping01"]["ref-strand"] == "-"
                and trns_dict[read_id]["mapping02"]["ref-strand"] == "-"
            ):
                combination_arrays[
                    tuple([
                        trns_dict[read_id]["mapping02"]["ref-chr"],
                        trns_dict[read_id]["mapping01"]["ref-chr"],
                    ])
                ][
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    - trns_dict[read_id]["mapping02"]["align-length"]
                    : trns_dict[read_id]["mapping02"]["ref-pos"],
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    - trns_dict[read_id]["mapping01"]["align-length"]
                    : trns_dict[read_id]["mapping01"]["ref-pos"],
                ] += 1
            elif (
                trns_dict[read_id]["mapping01"]["ref-strand"] == "-"
                and trns_dict[read_id]["mapping02"]["ref-strand"] == "+"
            ):
                combination_arrays[
                    tuple([
                        trns_dict[read_id]["mapping02"]["ref-chr"],
                        trns_dict[read_id]["mapping01"]["ref-chr"],
                    ])
                ][
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    - trns_dict[read_id]["mapping02"]["align-length"]
                    : trns_dict[read_id]["mapping02"]["ref-pos"],
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    : trns_dict[read_id]["mapping01"]["ref-pos"]
                    + trns_dict[read_id]["mapping01"]["align-length"],
                ] += 1
            elif (
                trns_dict[read_id]["mapping01"]["ref-strand"] == "+"
                and trns_dict[read_id]["mapping02"]["ref-strand"] == "-"
            ):
                combination_arrays[
                    tuple([
                        trns_dict[read_id]["mapping02"]["ref-chr"],
                        trns_dict[read_id]["mapping01"]["ref-chr"],
                    ])
                ][
                    trns_dict[read_id]["mapping02"]["ref-pos"]
                    : trns_dict[read_id]["mapping02"]["ref-pos"]
                    + trns_dict[read_id]["mapping02"]["align-length"],
                    trns_dict[read_id]["mapping01"]["ref-pos"]
                    - trns_dict[read_id]["mapping01"]["align-length"]
                    : trns_dict[read_id]["mapping01"]["ref-pos"],
                ] += 1
    return combination_arrays

    