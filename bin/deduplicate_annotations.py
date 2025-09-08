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
    # Check if the segments are the same
    if annotation1['segment01'] != annotation2['segment01'] or annotation1['segment02'] != annotation2['segment02']:
        return False

    # Check for overlap in both segments
    overlap_segment01 = (annotation1['start01'] <= annotation2['end01'] and annotation1['end01'] >= annotation2['start01'])
    overlap_segment02 = (annotation1['start02'] <= annotation2['end02'] and annotation1['end02'] >= annotation2['start02'])

    return overlap_segment01 and overlap_segment02


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
    annotation_table_df['area'] = (annotation_table_df['end01'] - annotation_table_df['start01']) * (annotation_table_df['end02'] - annotation_table_df['start02'])

    # parse the count table
    count_table_df = parse_count_table(count_table)

    # Reorder the annotation table to match the count table
    annotation_table_df = annotation_table_df.loc[count_table_df.index]

    # calculate the mean by area
    count_table_df['mean_by_area'] = count_table_df['mean'] / annotation_table_df['area']
    annotation_table_df['mean_by_area'] = count_table_df['mean_by_area']
    count_table_df = count_table_df.sort_values(by='mean_by_area', ascending=False)
    # Reorder the annotation table to match the count table
    annotation_table_df = annotation_table_df.loc[count_table_df.index]

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

    # Check wich interactions overlap with each of the deduplicated interactions
    main_interactions = {}
    graveyard_interactions = []
    for index in deduplicated_rows:
        annotation = annotation_dict[index]
        main_interactions[index] = []
        for idx, annotation2 in annotation_dict.items():
            if check_if_overlap(annotation, annotation2) and index != idx and idx not in graveyard_interactions:
                main_interactions[index].append(idx)
                graveyard_interactions.append(idx)


    new_column = {}
    for index, row in annotation_table_df.iterrows():
        if index in main_interactions.keys():
            new_column[index] = 'NA'
        elif index in graveyard_interactions:
            # search which interaction it overlaps with
            for main, overlaping in main_interactions.items():
                if index in overlaping:
                    new_column[index] = main


    # add new_column to annotation table
    annotation_table_df['overlaps_with'] = new_column.values()
    # add the mean column from the count table to the annotation table
    annotation_table_df['mean_readcount'] = count_table_df['mean']

    # export the annotation table to a tsv file
    annotation_table_df.to_csv(os.path.join(output_folder, 'annotation_table_extended.tsv'), sep='\t', index=False)

    # Create deduplicated DataFrame
    count_table_deduplicated_df = count_table_df.loc[deduplicated_rows]
    annotation_table_deduplicated_df = annotation_table_df.loc[deduplicated_rows]

    # export the deduplicated count & annotation tables to a tsv file
    count_table_deduplicated_df.to_csv(os.path.join(output_folder, 'count_table_deduplicated.tsv'), sep='\t')
    annotation_table_deduplicated_df.to_csv(os.path.join(output_folder, 'annotation_table_deduplicated.tsv'), sep='\t', index=False)

if __name__ == '__main__':
    main()
