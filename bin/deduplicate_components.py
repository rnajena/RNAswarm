#!/usr/bin/env python3

"""deduplicate_components.py

Takes an arbitrary number of trns files, finds and merges interactions, outputing an annotation table and the GMMs for each combination.

Usage:
    deduplicate_components.py -a <annotation_table> -c <count_table> -o <output_file>
    deduplicate_components.py -h | --help

Options:
    -h --help                                    Show this screen.
    -a --annotation_table=<annotation_table>     The annotation table.
    -c --count_table=<count_table>               The count table.
    -o --output=<output_file>                    The output file.

"""

# takes the 'square' with the most reads
# remove all the squares which overlaps with that one from the list of possible interactions
# take the square with most reads, repeat 2/3 until components are done

from docopt import docopt
import pandas as pd

