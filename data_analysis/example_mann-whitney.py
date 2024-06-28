import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import *

if __name__ == "__main__":
    # Load data
    parameters_set = 1
    h_threshold = 1.0

    ra_sfs_df = pd.read_pickle(f"parameters_set{parameters_set}_ra_sfs_{100*h_threshold:.0f}.pkl")
    ra_sfs = ra_sfs_df[ra_sfs_df["cell"] == -1]["sfs"][0] # Consider only the aggregated data across all cells
    ra_sfs = remove_homoplasmy(ra_sfs)

    ssd_sfs_df = pd.read_pickle(f"parameters_set{parameters_set}_ssd_sfs_{100*h_threshold:.0f}.pkl")
    ssd_sfs = ssd_sfs_df[ssd_sfs_df["cell"] == -1]["sfs"][0] # Consider only the aggregated data across all cells
    ssd_sfs = remove_homoplasmy(ssd_sfs)

    # Performs the Mann-Whitney U test
    rb, U, n1, n2, p = test_statistic(ra_sfs, ssd_sfs, "mann-whitney")
    print(f"rank-biserial correlation coefficient {rb = }, U-statistic {U = }, {n1 = }, {n2 = }, p-value {p = }")