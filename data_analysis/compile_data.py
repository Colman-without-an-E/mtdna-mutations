""" Compiles data across simulations """

import pandas as pd
from multiprocessing import Pool
import os
import sys
import traceback
from utils import *

### Save compiled data

if __name__ == "__main__":

    if len(sys.argv) != 4:
        sys.exit("Wrong number of arguments. Usage : compile_data.py <parameters_set> <sim_type> <data_type>")

    parameters_set = int(sys.argv[1])
    sim_type = sys.argv[2]
    if sim_type not in ["ra", "ssd", "wildtype"]:
        raise ValueError("sim_type must be one of 'ra', 'ssd', or 'wildtype'.")
    data_type = sys.argv[3]
    if data_type.lower() not in ["pop", "sfs"]:
        raise ValueError("sim_type must be either 'pop' or 'sfs'.")
    
    num_processes = len(os.sched_getaffinity(0)) # Linux
    # num_processes = os.cpu_count() # Windows

    filenames = os.listdir()
    if len(filenames) == 0:
        raise ValueError("No files found.")
    
    # Obtain data path
    try:
        if data_type.lower() == "pop":
            path_list = [fname for fname in filenames if fname.startswith(f"{sim_type}_sim_pop")]
        else:
            path_list = [fname for fname in filenames if fname.startswith(f"{sim_type}_sim_site")]
    except Exception as error:
        print(f"Error obtaining {sim_type.upper()} data path:", error)
        print(traceback.format_exc())

    # Compile data with parallel processing
    try:
        df = pd.DataFrame()
        if data_type.lower() == "pop":
            if sim_type == "wildtype":
                with Pool(num_processes) as pool:
                    df = pd.concat(pool.imap_unordered(get_pop_data_wildtype, path_list), ignore_index = True)
            else:
                with Pool(num_processes) as pool:
                    df = pd.concat(pool.imap_unordered(get_pop_data, path_list), ignore_index = True)
        else:
            with Pool(num_processes) as pool:
                df = pd.concat(pool.imap_unordered(get_sfs_data, path_list), ignore_index = True)
        print(sim_type.upper(), data_type.lower(), "df:\n", df.head())
        df.to_pickle(f"parameters_set{parameters_set}_{sim_type}_{data_type.lower()}_master.pkl")

    except Exception as error:
        print(f"Error loading and compiling {sim_type.upper()} data:", error)
        print(traceback.format_exc())
