import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import *

if __name__ == "__main__":

    # Load data
    parameters_set = 1
    sim_type = "ra"
    h_threshold = 1.0
    ra_sfs_df = pd.read_pickle(f"parameters_set{parameters_set}_{sim_type}_sfs_{100*h_threshold:.0f}.pkl")
    sfs = ra_sfs_df[ra_sfs_df["cell"] == -1]["sfs"][0] # Consider only the aggregated data across all cells
    sfs = remove_homoplasmy(sfs)

    # Plot the site frequency spectrum
    n_bins = 20 # set number of bins
    fig, ax = plt.subplots(figsize = (3,3))
    plot_sfs(ax, sfs, n_bins, density = True, label = sim_type, type = "line", marker = "o", markersize = 3)
    ax.set_xlabel("h")
    ax.set_ylabel("density")
    plt.tight_layout()
    plt.show()