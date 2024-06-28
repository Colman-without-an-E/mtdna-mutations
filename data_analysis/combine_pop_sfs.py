""" Combines population and site frequency spectrum data. """

import numpy as np
import pandas as pd
import sys

### Save compiled data

if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.exit("Wrong number of arguments. Usage : compile_data.py [parameters_set] [sim_type]")

    parameters_set = int(sys.argv[1])
    sim_type = sys.argv[2]
    if sim_type not in ["ra", "ssd", "wildtype"]:
        raise ValueError("sim_type must be one of 'ra', 'ssd', or 'wildtype'.")

    pop_df = pd.read_pickle(f"parameters_set{parameters_set}_{sim_type}_pop_master.pkl")
    sfs_df = pd.read_pickle(f"parameters_set{parameters_set}_{sim_type}_sfs_master.pkl")
    combined_df = pd.merge(pop_df, sfs_df, how="outer", on=["sim", "t", "cell"])
    del pop_df, sfs_df
    print(sim_type.upper(), "pop_df:\n", combined_df.head())
    combined_df.to_pickle(f"parameters_set{parameters_set}_{sim_type}_combined_master.pkl")
    print(f"Saved combined dataframe. Dataframe memory usage: {combined_df.memory_usage(deep = True).sum()/(1024**3)}gb.")