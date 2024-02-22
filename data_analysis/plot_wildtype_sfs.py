import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os

def get_parameters():
    """ Returns dict with value of saved parameters. """
    parameters = dict()
    with open("parameters.txt", "r") as file:

        # skip first line
        file.readline()

        lines = file.readlines()
        for l in lines:
            parameter, value = l.split(",")
            # Parse value accordingly into either float or int
            if "." in value:
                parameters[parameter] = float(value)
            else:
                parameters[parameter] = int(value)
    return parameters

def sfs_data_one_sim(path, sim):
    """ 
    Get site frequency spectrum data from text file of a simulation
     
    Parameters
    ----------
    path : str
        path of text file
    sim : int
        simulation index

    Returns
    -------
    pandas.core.frame.DataFrame
        with columns ["sim", "t", "sfs"]
    """

    df = pd.DataFrame(columns = ["sim", "t", "sfs"])
    with open(path, "r") as file:
        file.readline() # skip header row
        for line in file.readlines():
            _, t, sfs = line.split(",", maxsplit = 2)

            t = float(t)
            sfs = json.loads(sfs[:-1])

            df = df.append(pd.Series([sim, t, sfs], index = df.columns), ignore_index = True)
    return df

def populations_data_one_sim(path, sim):
    """ 
    Get population data from text file of a simulation
     
    Parameters
    ----------
    path : str
        path of text file
    sim : int
        simulation index

    Returns
    -------
    pandas.core.frame.DataFrame
        with columns ["sim", "t", "cell", "wildtype_population"]
    """

    df = pd.read_csv(path, header = None, names = ["sim", "cell", "t", "wildtype_population", "dummy"], skiprows = 1)
    df = df.reindex(columns = ["sim", "t", "cell", "wildtype_population"])
    df["sim"] = sim
    return df

def plot_sfs(ax, counts, population, n_bins = 20, density = False, label = "", type = "hist"):
    """
    Plots the site frequency spectrum
    
    Parameters
    ----------
    ax : matplotlib.pyplot.Axes
        ax to plot on
    counts : dict
        mutant counts
        counts[i] is the frequency of i
    population : int
        population, must be greater or equal to the maximum frequency in counts
    n_bins : int
        number of bins
    density : bool
        controls whether to plot density or frequency
    label : str
        plot label
    type : str
        type of plot, "hist" or "line"
    """

    # Unique heteroplasmy values
    H_unique = np.array(list(counts.keys())).astype(float) / population

    bins = np.arange(0.0, 1.0, 1/n_bins)
    bin_idx = np.floor(H_unique * (n_bins-1) + 1).astype(int)-1
    bin_values = np.zeros(n_bins, dtype = int)
    bin_values[bin_idx] += np.array(list(counts.values()))

    if density:
        bin_values = bin_values / np.sum(bin_values)

    if type == "hist":
        ax.bar(bins, bin_values, width = np.ones_like(bins) / n_bins, align = "edge", label = label)
    elif type == "line":
        ax.plot(bins, bin_values, label = label)

def merge_sfs(sfs_list, pop_list, dp = 4):
    """ Merges site frequency spectra from different simulations. Also normalises site frequency spectrum accordingly.
    
    Parameters
    ----------
    sfs_list : list
        list of dictionaries containing site frequency spectrum data
    pop_list : list
        list of total population across cells for normalisation
    dp : int
        number of decimal places to retain after normalising sfs keys. Controls precision
    
    Returns
    -------
    dict
        ensemble average sfs (normalised)
    """

    ensemble_sfs = dict()
    for sfs, pop in zip(sfs_list, pop_list):
        normalised_sfs_keys = np.array(list(sfs.keys())).astype(int) / pop

        # Normalised keys with 4 decimal places are keys
        normalised_sfs_keys = [np.format_float_positional(h, precision = dp) for h in normalised_sfs_keys]
        normalised_sfs = dict(zip(normalised_sfs_keys, sfs.values()))

        ensemble_sfs = {item: ensemble_sfs.get(item, 0) + normalised_sfs.get(item, 0) for item in set(ensemble_sfs) | set(normalised_sfs)}
    
    return ensemble_sfs

### Save parameter values
os.chdir("../data/parameters_set1")
parameters = get_parameters()


### Simple site frequency spectrum plot from one simulation

# Load data
sim = 4
sfs_path = f"wildtype_sim_site_frequency_spectrum{sim}.txt"
pop_path = f"wildtype_sim_populations{sim}.txt"
sfs_df = sfs_data_one_sim(sfs_path, sim)
pop_df = populations_data_one_sim(pop_path, sim)

# Plot site frequency spectrum at specified t
plot_at_t = [100, 2000, 4000, 8000, 10000]
sfs_at_t = sfs_df[sfs_df["t"].isin(plot_at_t)]["sfs"].to_list()
pop_sum_at_t = pop_df[pop_df["t"].isin(plot_at_t)].groupby(["sim", "t"], as_index = False)["wildtype_population"].sum()["wildtype_population"].to_numpy()
fig, ax = plt.subplots()
for i, t in enumerate(plot_at_t):
    plot_sfs(ax, sfs_at_t[i], pop_sum_at_t[i], label = f"t = {t}", type = "line")
ax.set_xlim((0,1))
ax.set_xlabel("heteroplasmy")
ax.set_ylabel("frequency")
ax.legend()
ax.grid()
plt.show()


### Ensemble average site frequency spectrum plot

# Load data
sfs_df = pd.DataFrame()
pop_df = pd.DataFrame()
for sim in np.arange(4,6):
    sfs_path = f"wildtype_sim_site_frequency_spectrum{sim}.txt"
    pop_path = f"wildtype_sim_populations{sim}.txt"
    sfs_df = sfs_df.append(sfs_data_one_sim(sfs_path, sim))
    pop_df = pop_df.append(populations_data_one_sim(pop_path, sim))

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