import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from utils import get_sfs_data, get_populations_data, plot_sfs, merge_sfs

### Ensemble average site frequency spectrum plot

if __name__ == "__main__":

    # Go to correct directory
    os.chdir("../data/parameters_set1")
    
    # Load data
    sfs_df = pd.DataFrame()
    pop_df = pd.DataFrame()
    for sim in np.arange(4,6):
        sfs_path = f"wildtype_sim_site_frequency_spectrum{sim}.txt"
        pop_path = f"wildtype_sim_populations{sim}.txt"
        sfs_df = sfs_df.append(get_sfs_data(sfs_path, sim))
        pop_df = pop_df.append(get_populations_data(pop_path, sim))

    plot_at_t = [100, 2000, 4000, 8000, 10000]
    fig, ax = plt.subplots()
    for t in plot_at_t:
        sfs_list = sfs_df[sfs_df["t"] == t]["sfs"].to_list()
        pop_sum_list = pop_df[pop_df["t"] == t].groupby(["sim", "t"], as_index = False)["wildtype_population"].sum()["wildtype_population"].to_numpy()
        ensemble_sfs = merge_sfs(sfs_list, pop_sum_list, dp = 4)

        # Plot ensemble sfs at time t
        plot_sfs(ax, ensemble_sfs, 1, label = f"t = {t}", density = True, type = "line")
    ax.set_xlim((0,1))
    ax.set_xlabel("heteroplasmy")
    ax.set_ylabel("frequency")
    ax.legend()
    ax.grid()
    plt.show()