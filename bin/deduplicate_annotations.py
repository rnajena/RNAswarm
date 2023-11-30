#!/usr/bin/env python3

"""deduplicate_annotations.py

Takes a count table and an annotation table and removes all components which overlap with the component with the most reads.

Usage:
    deduplicate_annotations.py -a <annotation_table> -c <count_table> -o <output_folder>
    deduplicate_annotations.py -h | --help

Options:
    -h --help                                    Show this screen.
    -a --annotation_table=<annotation_table>     The annotation table.
    -c --count_table=<count_table>               The count table.
    -o --output=<output_folder>                  The output folder.

"""

# takes the 'square' with the most reads
# remove all the squares which overlaps with that one from the list of possible interactions
# take the square with most reads, repeat 2/3 until components are done

from docopt import docopt
import pandas as pd
import os

def parse_annotation_table(annotation_table):
    """Parses the annotation table into a pandas dataframe.

    Args:
        annotation_table (str): The annotation table.

    Returns:
        pandas.DataFrame: The annotation table as a dataframe.
    """

    df = pd.read_csv(annotation_table, sep='\t', header=0)
    return df


def parse_count_table(count_table, calculate_means=True, sort_by_mean=True):
    """Parses the count table into a pandas dataframe.

    Args:
        count_table (str): The count table.
        calculate_means (bool): Whether to calculate the means of the count table.
        sort_by_mean (bool): Whether to sort the count table by the mean.

    Returns:
        pandas.DataFrame: The count table as a dataframe.
    """
    # first line is the header, first (unamed) column is the id
    df = pd.read_csv(count_table, sep='\t', header=0, index_col=0)
    if calculate_means:
        # for each line in the count table, calculate the mean of the counts (ignoring the id column)
        df['mean'] = df.mean(axis=1)
        if sort_by_mean:
            df = df.sort_values(by='mean', ascending=False)
    elif not calculate_means and sort_by_mean:
        # cannot sort by mean if it has not been calculated
        raise ValueError('Cannot sort by mean if it has not been calculated.')
    return df


def check_if_overlap(annotation1, annotation2):
    """Checks if two annotations overlap.

    An annotation is a pandas.Series with the following columns:
    - id
    - segment01
    - start01
    - end01
    - segment02
    - start02
    - end02

    Args:
        annotation1 (pandas.Series): The first annotation.
        annotation2 (pandas.Series): The second annotation.

    Returns:
        bool: True if the annotations overlap, False otherwise.


    """
    # check if the segments are the same
    if annotation1['segment01'] != annotation2['segment01'] or annotation1['segment02'] != annotation2['segment02']:
        return False
    # check if the start of the first annotation is between the start and end of the second annotation
    if annotation1['start01'] >= annotation2['start01'] and annotation1['start01'] <= annotation2['end01']:
        return True
    # check if the end of the first annotation is between the start and end of the second annotation
    if annotation1['end01'] >= annotation2['start01'] and annotation1['end01'] <= annotation2['end01']:
        return True
    # check if the start of the second annotation is between the start and end of the first annotation
    if annotation2['start02'] >= annotation1['start02'] and annotation2['start02'] <= annotation1['end02']:
        return True
    # check if the end of the second annotation is between the start and end of the first annotation
    if annotation2['end02'] >= annotation1['start02'] and annotation2['end02'] <= annotation1['end02']:
        return True
    # if none of the above, there is no overlap
    return False


def main():
    """Main function.
    """

    # parse the arguments
    arguments = docopt(__doc__)
    annotation_table = arguments['--annotation_table']
    count_table = arguments['--count_table']
    output_folder = arguments['--output']

    # parse the annotation table
    annotation_table_df = parse_annotation_table(annotation_table)

    # parse the count table
    count_table_df = parse_count_table(count_table)

    # Convert annotation table to a dictionary for faster access
    annotation_dict = annotation_table_df.set_index('id').to_dict('index')

    # List to store rows for the deduplicated DataFrame
    deduplicated_rows = []

    # Iterate over the count table
    for index, row in count_table_df.iterrows():
        # Get the annotation
        annotation = annotation_dict[index]

        # Check if the annotation overlaps with any of the following annotations
        if not any(check_if_overlap(annotation, annotation_dict[idx]) for idx in deduplicated_rows):
            deduplicated_rows.append(index)

    # Create deduplicated DataFrame
    # set index to the id column
    # count_table_deduplicated_df = count_table_df.loc[deduplicated_rows]
    count_table_deduplicated_df = count_table_df.loc[deduplicated_rows]
    annotation_table_deduplicated_df = annotation_table_df.loc[deduplicated_rows]

    # Write the deduplicated count table to a file
    output_files = [os.path.join(output_folder, 'deduplicated_count_table.tsv'),
                   os.path.join(output_folder, 'deduplicated_annotation_table.tsv')]
    count_table_deduplicated_df.to_csv(output_files[0], sep='\t', header=True)
    annotation_table_deduplicated_df.to_csv(output_files[1], sep='\t', header=True)


if __name__ == '__main__':
    main()
