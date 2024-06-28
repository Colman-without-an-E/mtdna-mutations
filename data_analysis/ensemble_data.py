""" Save ensemble results """

import numpy as np
import pandas as pd
import sys
import traceback
from utils import *


### Save ensemble site frequency spectrum

if __name__ == "__main__":
    
    # Default values
    h_threshold = 1.0 # record SFS when heteroplasmy threshold is reached

    match len(sys.argv):
        case 3:
            parameters_set = int(sys.argv[1])
            sim_type = sys.argv[2]
        case 4:
            parameters_set = int(sys.argv[1])
            sim_type = sys.argv[2]
            h_threshold = float(sys.argv[3])
        case _:
            sys.exit("Wrong number of arguments. Usage : ensemble_data.py <parameters_set> [h_threshold]")

    # Load data
    try:
        master_df = pd.read_pickle(f"parameters_set{parameters_set}_{sim_type}_combined_master.pkl")
        print(f"Loaded {sim_type.upper()} data.")
    except Exception as error:
        print(f"Error loading {sim_type.upper()} data:", error)
        print(traceback.format_exc())

    # Ensemble and save data
    try:
        ensemble_h_df, ensemble_sfs_df = ensemble_non_wildtype_data(master_df, h_threshold = h_threshold)
        ensemble_h_df.to_pickle(f"parameters_set{parameters_set}_{sim_type}_h_{100*h_threshold:.0f}.pkl")
        ensemble_sfs_df.to_pickle(f"parameters_set{parameters_set}_{sim_type}_sfs_{100*h_threshold:.0f}.pkl")
        print(f"Saved ensemble {sim_type.upper()} heteroplasmy and site frequency spectrum data.")
    except Exception as error:
        print(f"Error saving ensemble {sim_type.upper()} data:", error)
        print(traceback.format_exc())