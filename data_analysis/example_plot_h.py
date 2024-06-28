import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utils import *

if __name__ == "__main__":

    # Load data
    parameters_set = 1
    sim_type = "ra"
    h_threshold = 1.0
    h_df = pd.read_pickle(f"parameters_set{parameters_set}_{sim_type}_h_{100*h_threshold:.0f}.pkl")
    h_df = h_df[h_df["cell"] == -1] # Consider only the aggregated data across all cells

    # Plot the increase in mutant fraction against time
    colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, ax = plt.subplots(figsize = (3,3))
    ax.plot(h_df["t"], h_df["h_mean"], linewidth = 0.8, color = colours[0])
    ax.fill_between(h_df["t"], h_df["h_mean"]-h_df["h_sem"], h_df["h_mean"]+h_df["h_sem"], color = colours[0], alpha = 0.25)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\langle f \rangle$")
    plt.tight_layout()
    plt.show()