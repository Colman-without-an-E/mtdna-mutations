import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from utils import *

def plot_wildtype(plot_at_t, data_dir = "../data", parameters_set = 1, sims = range(1,11), fig_save_path = None, show = False):
    """ Write """
    # Go to directory of data
    cwd = os.getcwd()
    os.chdir(f"{data_dir}/parameters_set{parameters_set}")
    
    # Load data
    sfs_df = pd.DataFrame()
    pop_df = pd.DataFrame()
    for sim in sims:
        sfs_path = f"wildtype_sim_site_frequency_spectrum{sim}.txt"
        pop_path = f"wildtype_sim_populations{sim}.txt"
        sfs_df = pd.concat([sfs_df, _get_sfs_data(sfs_path, sim)])
        pop_df = pd.concat([pop_df, get_populations_data(pop_path, sim, just_wildtype = True)])
    
    # Return to initial working directory
    os.chdir(cwd)

    # Create figure and ax to plot in
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
    if fig_save_path is not None:
        plt.savefig(fig_save_path)
    if show:
        plt.show()


def plot_ssd(plot_at_t, data_dir = "../data", parameters_set = 1, sims = range(1,11), extinct_t_threshold = None, fig_save_path = None, show = False):
    """ Write """
    # Go to directory of data
    cwd = os.getcwd()
    os.chdir(f"{data_dir}/parameters_set{parameters_set}")
    
    # Load data
    sfs_path_list = [f"ssd_sim_site_frequency_spectrum{sim}.txt" for sim in sims]
    sfs_df = compile_simulations_sfs(sfs_path_list, sims = sims, extinct_t_threshold = extinct_t_threshold, just_wildtype = False)

    pop_path_list = [f"ssd_sim_populations{sim}.txt" for sim in sims]
    pop_df = compile_simulations_populations(pop_path_list, sims = sims, extinct_t_threshold = extinct_t_threshold, just_wildtype = False)

    # Return to initial working directory
    os.chdir(cwd)

    # Create figure and ax to plot in
    fig, ax = plt.subplots()
    for t in plot_at_t:
        sfs_list = sfs_df[sfs_df["t"] == t]["sfs"].to_list()
        pop_sum_list = pop_df[pop_df["t"] == t].groupby(["sim", "t"])[["wildtype_population", "ssd_population"]].sum().sum(axis = 1).to_numpy()
        ensemble_sfs = merge_sfs(sfs_list, pop_sum_list, dp = 4)

        # Plot ensemble sfs at time t
        plot_sfs(ax, ensemble_sfs, 1, label = f"t = {t}", density = True, type = "line")
    ax.set_xlim((0,1))
    ax.set_xlabel("heteroplasmy")
    ax.set_ylabel("frequency")
    ax.legend()
    ax.grid()
    if fig_save_path is not None:
        plt.savefig(fig_save_path)
    if show:
        plt.show()

    # Returns the number of undiscarded simulations
    return len(sfs_df["sim"].value_counts())


def plot_ra(plot_at_t, data_dir = "../data", parameters_set = 1, sims = range(1,11), extinct_t_threshold = None, fig_save_path = None, show = False):
    """ Write """
    # Go to directory of data
    cwd = os.getcwd()
    os.chdir(f"{data_dir}/parameters_set{parameters_set}")
    
    # Load data
    sfs_path_list = [f"ra_sim_site_frequency_spectrum{sim}.txt" for sim in sims]
    sfs_df = compile_simulations_sfs(sfs_path_list, sims = sims, extinct_t_threshold = extinct_t_threshold, just_wildtype = False)

    pop_path_list = [f"ra_sim_populations{sim}.txt" for sim in sims]
    pop_df = compile_simulations_populations(pop_path_list, sims = sims, extinct_t_threshold = extinct_t_threshold, just_wildtype = False)

    # Return to initial working directory
    os.chdir(cwd)

    # Create figure and ax to plot in
    fig, ax = plt.subplots()
    for t in plot_at_t:
        sfs_list = sfs_df[sfs_df["t"] == t]["sfs"].to_list()
        pop_sum_list = pop_df[pop_df["t"] == t].groupby(["sim", "t"])[["wildtype_population", "ra_population"]].sum().sum(axis = 1).to_numpy()
        ensemble_sfs = merge_sfs(sfs_list, pop_sum_list, dp = 4)

        # Plot ensemble sfs at time t
        plot_sfs(ax, ensemble_sfs, 1, label = f"t = {t}", density = True, type = "line")
    ax.set_xlim((0,1))
    ax.set_xlabel("heteroplasmy")
    ax.set_ylabel("frequency")
    ax.legend()
    ax.grid()
    if fig_save_path is not None:
        plt.savefig(fig_save_path)
    if show:
        plt.show()

    # Returns the number of undiscarded simulations
    return len(sfs_df["sim"].value_counts())


### Ensemble average site frequency spectrum plot

if __name__ == "__main__":
    ssd_sims = range(10001, 20001)
    ra_sims = range(20001, 30001)
    plot_at_t = [10000, 12500, 15000, 17500, 20000]
    
    n_undiscarded_ssd_sims = plot_ssd(plot_at_t, ".", 1, ssd_sims, extinct_t_threshold = 11000, fig_save_path = "parameters1_ssd_sfs1.png")
    n_undiscarded_ra_sims = plot_ra(plot_at_t, ".", 1, ra_sims, extinct_t_threshold = 11000, fig_save_path = "parameters1_ra_sfs1.png")
    
    print(f"Undiscarded ssd simulations: {n_undiscarded_ssd_sims}")
    print(f"Undiscarded ra simulations: {n_undiscarded_ra_sims}")