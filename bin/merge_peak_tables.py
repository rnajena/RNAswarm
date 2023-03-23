#!/usr/bin/env python3

import pandas as pd

alias_dict = {
    "WyNAwt": "/home/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20230218/WyNAwt_correct_peak_cells.tsv",
    "WyNA_Udsub": "/home/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20230218/WyNAUdsub_correct_peak_cells.tsv",
}

def merge_peak_tables(alias_dict, output_file):
    """
    Merge peak tables

    Parameters
    ----------
    alias_dict : dict
        A dictionary with the aliases as keys and the filepaths as values
    output_file : str
        The output filepath
    """
    df_dict = {}
    for alias, filepath in alias_dict.items():
        # Read the peak table
        df_dict[alias] = pd.read_csv(filepath, sep="\t")

    # Check if the headers are the same for all peak tables
    for alias, df in df_dict.items():
        if not df.columns.equals(df_dict["WyNAwt"].columns):
            raise ValueError("The headers of the peak tables are not the same!")

    # Check if the first 9 columns (I am counting index as column 0) are the same for all peak tables
    for alias, df in df_dict.items():
        if not df.iloc[:, :9].equals(df_dict["WyNAwt"].iloc[:, :9]):
            raise ValueError("The first 8 columns of the peak tables are not the same!")

    # Create an empty dataframe with the same columns as the peak tables
    merged_df = pd.DataFrame(columns=df_dict["WyNAwt"].columns)

    # For each alias check if "01type" and "02type" are the same as the alias
    for alias, df in df_dict.items():
        # Iterate over the rows
        for index, row in df.iterrows():
            # Check if the "01type" and "02type" are both the same as the alias
            if row["01type"] == alias and row["02type"] == alias:
                # The frame.append method is deprecated and will be removed from pandas in a 
                # future version. Use pandas.concat instead.
                # merged_df.append(row)
                # make sure to keep the index
                merged_df = pd.concat([merged_df, pd.DataFrame(row).T], ignore_index=False)
            elif row["01type"] != row["02type"]:
                raise ValueError("The 01type and 02type are not the same!")

    # Write the sorted merged peak table (sort using the id column)
    merged_df.sort_values(by=["id"]).to_csv(output_file, sep="\t", index=False)

def main():
    merge_peak_tables(alias_dict, "/home/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20230218/PR8-Wy_merged_peak_tables.tsv")

if __name__ == "__main__":
    main()
