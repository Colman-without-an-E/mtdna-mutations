import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import json

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

def get_sfs_data(path):
    """ 
    Get site frequency spectrum data from text file of a simulation
     
    Parameters
    ----------
    path : str
        path of text file

    Returns
    -------
    pandas.core.frame.DataFrame
        with columns ["sim", "t", "cell", "sfs"]
    """

    df_row_list = []
    with open(path, "r") as file:
        file.readline() # skip header row
        for line in file.readlines():

            # Parse data
            t, cell, sfs = line.split(",", maxsplit = 2)

            # Append each row to df_row_list
            df_row = [float(t), int(cell), json.loads(sfs[:-1])]
            df_row_list.append(df_row)

    # Create and return dataframe
    df = pd.DataFrame(df_row_list, columns = ["t", "cell", "sfs"])
    sim = int(re.search(r"\d+", path).group())
    df["sim"] = sim
    
    return df

def _get_sfs_data(path, sim):
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
        with columns ["t", "sfs", "sim"]
    """

    df_row_list = []
    with open(path, "r") as file:
        file.readline() # skip header row
        for line in file.readlines():

            # Parse data
            t, cell, sfs = line.split(",", maxsplit = 2)

            # Append each row to df_row_list
            df_row = [sim, float(t), int(cell), json.loads(sfs[:-1])]
            df_row_list.append(df_row)

    # Create and return dataframe
    df = pd.DataFrame(df_row_list, columns = ["sim", "t", "cell", "sfs"])
    return df

def get_pop_data_wildtype(path):
    """
    Get wildtype population data from text file of a simulation
     
    Parameters
    ----------
    path : str
        path of text file
    sim : int
        simulation index
    just_wildtype : bool
        controls whether to include last column of text file

    Returns
    -------
    pandas.core.frame.DataFrame
        with columns ["t", "cell", "w", "sim"]
    """
    df = pd.read_csv(path)
    df.drop("m", axis = 1, inplace = True)
    sim = int(re.search(r"\d+", path).group())
    df["sim"] = sim
    return df

def get_pop_data(path):
    """
    Get non-wildtype population data from text file of a simulation
     
    Parameters
    ----------
    path : str
        path of text file
    sim : int
        simulation index
    just_wildtype : bool
        controls whether to include last column of text file

    Returns
    -------
    pandas.core.frame.DataFrame
        with columns ["t", "cell", "w", "m", "sim"]
    """
    df = pd.read_csv(path)
    sim = int(re.search(r"\d+", path).group())
    df["sim"] = sim
    return df

def _get_populations_data(path, sim, just_wildtype = False):
    """
    Get population data from text file of a simulation
     
    Parameters
    ----------
    path : str
        path of text file
    sim : int
        simulation index
    just_wildtype : bool
        controls whether to include last column of text file

    Returns
    -------
    pandas.core.frame.DataFrame
        with columns ["t", "cell", "w", "m", "sim"], or
        with columns ["t", "cell", "w", "sim"] if just_wildtype
    """
    df = pd.read_csv(path)
    if just_wildtype:
        df.drop("m", axis = 1, inplace = True)
    df["sim"] = sim
    return df

def plot_sfs(ax, counts, n_bins = 20, density = True, label = "", type = "hist", alpha = 1, **kwargs):
    """
    Plots the site frequency spectrum
    
    Parameters
    ----------
    ax : matplotlib.pyplot.Axes
        ax to plot on
    counts : dict
        mutant counts
        counts[i] is the frequency of i
    n_bins : int
        number of bins
    density : bool
        controls whether to plot density or frequency
    label : str
        plot label
    type : str
        type of plot, "hist" or "line"
    """

    if type == "hist":
        ax.hist(np.array(list(counts.keys()), dtype = float), bins = n_bins, range = (0, 1), weights = counts.values(), density = density, label = label, alpha = alpha)
    
    elif type == "line":
        hist, bins = np.histogram(np.array(list(counts.keys()), dtype = float), bins = n_bins, range = (0, 1), density = density, weights = list(counts.values()))
        bin_centres = (bins[1:] + bins[:-1])/2
        ax.plot(bin_centres, hist, label = label, alpha = alpha, **kwargs)
    
    else:
        raise ValueError("type can only be 'hist' or 'line'.")
    return

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

def get_extinction_time(df):
    """
    Returns the last time of record of a single simulation
    """
    return df["t"].max()

def compile_simulations_sfs(path_list, sims = np.arange(1, 11), extinct_t_threshold = None):
    """
    Compile sfs data from simulations. 

    Parameters
    ----------
    path_list : list
        list containing the paths to simulations of the same type
    sims : iterable
        simulation indices corresponding to the paths
    extinct_t_threshold : NoneType or float
        discard simulations where extinction time <= extinct_t_threshold
    just_wildtype : bool
        controls whether to include last column of text file

    Returns
    -------
    pandas.core.frame.DataFrame
        compiled simulations site frequency spectrum data
    """

    compiled_df = pd.DataFrame()

    if extinct_t_threshold is None:
        for path, sim in zip(path_list, sims):
            single_sim_df = _get_sfs_data(path, sim)
            compiled_df = pd.concat([compiled_df, single_sim_df])
    else:
        for path, sim in zip(path_list, sims):
            single_sim_df = _get_sfs_data(path, sim)

            # Discard simulations where extinction occurs soon
            extinction_t = get_extinction_time(single_sim_df)
            if extinction_t>extinct_t_threshold:
                compiled_df = pd.concat([compiled_df, single_sim_df])
    
    # Assign column names if no simulations pass through extinction threshold filter
    if compiled_df.shape[1] == 0:
        compiled_df = pd.DataFrame(columns = single_sim_df.columns)
    
    return compiled_df

def compile_simulations_populations(path_list, sims = np.arange(1, 11), extinct_t_threshold = None, just_wildtype = False):
    """
    Compile populations data from simulations. 

    Parameters
    ----------
    path_list : list
        list containing the paths to simulations of the same type
    sims : iterable
        simulation indices corresponding to the paths
    extinct_t_threshold : NoneType or float
        discard simulations where extinction time <= extinct_t_threshold
    just_wildtype : bool
        controls whether to include last column of text file

    Returns
    -------
    pandas.core.frame.DataFrame
        compiled simulations populations data
    """

    compiled_df = pd.DataFrame()

    if extinct_t_threshold is None:
        for path, sim in zip(path_list, sims):
            single_sim_df = _get_populations_data(path, sim, just_wildtype)
            compiled_df = pd.concat([compiled_df, single_sim_df])
    else:
        for path, sim in zip(path_list, sims):
            single_sim_df = _get_populations_data(path, sim, just_wildtype)

            # Discard simulations where extinction occurs soon
            extinction_t = get_extinction_time(single_sim_df)
            if extinction_t>extinct_t_threshold:
                compiled_df = pd.concat([compiled_df, single_sim_df])
    
    # Assign column names if no simulations pass through extinction threshold filter
    if compiled_df.shape[1] == 0:
        compiled_df = pd.DataFrame(columns = single_sim_df.columns)
    
    return compiled_df

# Post data ensembling/compilation functions

def remove_zeros(sfs):
    """
    Remove counts of zero heteroplasmy in site frequency spectrum
    
    Parameters
    ----------
    sfs : dict
        sfs[h] is the number of mutations with heteroplasmy h. h is a string.
    """
    try:
        del sfs["0."]
    except KeyError:
        pass
    return sfs

def remove_zeros_sfs_df(df1, *args):
    """
    Removes zero-mutations in the site frequency spectrum column of dataframes.
    
    Parameters
    ----------
    df1 : pandas.core.dataframe.DataFrame
    *args : df2, df3, ...
    """
    for df in list((df1,) + args):
        df["sfs"] = df["sfs"].map(remove_zeros)

def remove_homoplasmy(sfs):
    """
    Returns a copy of site frequency spectrum with homoplasmic counts deleted
    
    Parameters
    ----------
    sfs : dict
        sfs[h] is the number of mutations with heteroplasmy h. h is a string.
    
    Returns
    sfs_copy : dict
    """
    sfs_copy = sfs.copy()

    if "1." in sfs_copy:
        del sfs_copy["1."]
    return sfs_copy

def remove_homoplasmy_sfs_df(df1, *args):
    """
    Removes homoplasmic mutations in the site frequency spectrum column of dataframes.
    
    Parameters
    ----------
    df1 : pandas.core.dataframe.DataFrame
    *args : df2, df3, ...
    """
    for df in list((df1,) + args):
        df["sfs"] = df["sfs"].map(remove_homoplasmy)

def remove_nan(sfs):
    try:
        del sfs["nan"]
    except KeyError:
        pass
    return sfs

def remove_nan_sfs_df(df1, *args):
    for df in list((df1,) + args):
        df["sfs"].map(remove_nan)


# Pandas functions

def normalise(x, dp = 4, sfs_colname = "sfs", total_pop_colname = "total_population"):
    """ Normalises site frequency spectrum to the interval [0,1]. """
    sfs = x[sfs_colname]
    pop = x[total_pop_colname]
    normalised_sfs_keys = np.array(list(sfs.keys())).astype(int) / pop

    # Normalised keys with 4 decimal places are keys
    normalised_sfs_keys = [np.format_float_positional(h, precision = dp) for h in normalised_sfs_keys]
    normalised_sfs = dict(zip(normalised_sfs_keys, sfs.values()))

    return normalised_sfs

def additive_union(x):
    """ Returns the additive union of counters stored as dicts. """
    union = dict()
    for d in x:
        for key, value in d.items():
            union[key] = union.get(key, 0) + value
    return union

def sfs_weight(sfs):
    """
    Computes weight of a site frequency spectrum.
    
    Parameters
    ----------
    sfs : dict
        sfs[h] is the number of mutations with heteroplasmy h. h is a string.
    
    Returns
    -------
    dict
        dict["h"] is the weight of heteroplasmy count with respect to h. h is a string. 
    """

    # Compute weight associated to each heteroplasmy value
    sfs_counts = np.array(list(sfs.values()))
    tot_counts = np.sum(sfs_counts)
    weights = sfs_counts / tot_counts

    # Return weights as dict
    weights_dict = {h: count for h, count in zip(sfs.keys(), weights)}
    return weights_dict

def df_ffill(df):
    """
    Forward fill simulations data frame.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame

    Returns
    -------
    pandas.core.frame.DataFrame
    """

    dfs = []
    for cell in df["cell"].unique():
        sims = df[df["cell"]==cell]["sim"].unique()
        ts = df[df["cell"]==cell]["t"].unique()

        # Template data frame with pre-set indices
        new_index = pd.MultiIndex.from_product([sims, [cell], ts], names=["sim", "cell", "t"])
        template_df = pd.DataFrame(index = new_index)

        # Forward fill
        dfs.append(pd.merge(template_df, df, left_index=True, right_on=["sim", "cell", "t"], how="left").ffill())
        
    filled_df = pd.concat(dfs)
    return filled_df


# Statistical tests

def sfs_median(sfs):
    """
    Returns the median of a site frequency spectrum dictionary.
    
    Parameters
    ----------
    sfs : dict
        site frequency spectrum
    
    Returns
    -------
    float
        median
    """

    total_count = sum(sfs.values())
    target_count = total_count // 2  # Target count for median

    key_list = list(sfs.keys())
    np.argsort

    cumulative_count = 0
    prev_key = None
    for key_i in np.argsort(np.array(key_list, dtype = float)):
        key = key_list[key_i]
        count = sfs[key]
        cumulative_count += count
        if cumulative_count >= target_count:
            if total_count % 2 == 0 and cumulative_count == target_count:
                # If total count is even and at the first middle count
                return (prev_key + float(key)) / 2
            else:
                return float(key)
        prev_key = float(key)

def test_statistic(sfs1, sfs2, statistic = "rb", centre = "median", n_bins = 20):
    """
    Returns the test statistic of site frequency spectrum. 
    
    Parameters
    ----------
    sfs2 : dict
        site frequency spectrum
    sfs2 : dict
        site frequency spectrum
    statistic : str
        type of test statistic to return
        if "rb", rank-biserial correlation coefficient
        if "levene", Levene's test statistic
        if "hist_l2", L2 (Euclidean) distance between histograms
        if "mann-whitney", performs the Mann-Whitney U test
    centre : str
        centre of Levene's test statistic
        either "median" or "mean", default "median"
        
    Returns
    -------
    float or tuple
        test statistic;
        for mann-whitney, returns a tuple
            (rank-biserial correlation coefficient, Mann-Whitney U statistic, n1, n2, two-tailed p-value)
    """

    # Store counts
    n1 = np.sum(list(sfs1.values()))
    n2 = np.sum(list(sfs2.values()))
    N = n1 + n2

    match statistic:
        case "rb":
            # Aggregate site frequency spectrum by taking additive union
            union = additive_union([sfs1, sfs2])

            # Sort aggregated site frequency by heteroplasmy
            tot_h_counts = sorted([(h, count) for h, count in union.items()], key = lambda x: float(x[0]))
            tot_h, tot_counts = zip(*tot_h_counts)
            tot_counts = np.array(tot_counts)

            # Compute corresponding rank of each unique value of heteroplasmy
            ranks = np.append([0], np.cumsum(tot_counts[:-1])) + (tot_counts+1)/2
            
            # Compute the rank means of each site frequency spectrum
            sfs1_rank_counts = np.array([sfs1.get(h, 0) for h in tot_h])
            sfs1_rank_mean = np.sum(ranks * sfs1_rank_counts) / n1
            sfs2_rank_counts = np.array([sfs2.get(h, 0) for h in tot_h])
            sfs2_rank_mean = np.sum(ranks * sfs2_rank_counts) / n2

            # Glass's formula for rank-biserial correlation coefficient
            stat = 2*(sfs1_rank_mean - sfs2_rank_mean) / N
        
        case "mann-whitney":
            from scipy.stats import norm
            # Aggregate site frequency spectrum by taking additive union
            union = additive_union([sfs1, sfs2])

            # Sort aggregated site frequency by heteroplasmy
            tot_h_counts = sorted([(h, count) for h, count in union.items()], key = lambda x: float(x[0]))
            tot_h, tot_counts = zip(*tot_h_counts)
            tot_counts = np.array(tot_counts)

            # Compute corresponding rank of each unique value of heteroplasmy
            ranks = np.append([0], np.cumsum(tot_counts[:-1])) + (tot_counts+1)/2
            
            # Compute the rank means of each site frequency spectrum
            sfs1_rank_counts = np.array([sfs1.get(h, 0) for h in tot_h])
            sfs1_rank_mean = np.sum(ranks * sfs1_rank_counts) / n1
            sfs2_rank_counts = np.array([sfs2.get(h, 0) for h in tot_h])
            sfs2_rank_mean = np.sum(ranks * sfs2_rank_counts) / n2

            # Glass's formula for rank-biserial correlation coefficient
            rb = 2*(sfs1_rank_mean - sfs2_rank_mean) / N

            # Mann-Whitney U test from rank-biserial
            mu = 0.5*n1*n2
            u = mu*(1-np.abs(rb))
            rank_tie_counts = np.array([sfs1_rank_counts, sfs2_rank_counts], dtype = np.uint64).min(axis = 0) # adjustment for tied ranks
            adjustment = np.sum(rank_tie_counts**3 - rank_tie_counts) / (N*(N-1))
            sigma = np.sqrt(mu*(N+1-adjustment)/6)
            z = (u-mu)/sigma
            p_val = norm.sf(abs(z))*2
            
            stat = (rb, u, n1, n2, p_val)

        case "levene":
            h1 = np.array(list(sfs1.keys()), dtype = float)
            h2 = np.array(list(sfs2.keys()), dtype = float)
            counts1 = np.array(list(sfs1.values()), dtype = int)
            counts2 = np.array(list(sfs2.values()), dtype = int)

            match centre:
                case "median":
                    h_centre1 = sfs_median(sfs1)
                    h_centre2 = sfs_median(sfs2)
                case "mean":
                    h_centre1 = np.sum(h1 * counts1) / n1
                    h_centre2 = np.sum(h2 * counts2) / n2
                case _:
                    raise ValueError("Centre must be 'median' or 'mean'.")
        
            # Distance from centre
            D1 = np.abs(h1 - h_centre1)
            D2 = np.abs(h2 - h_centre2)

            d_sum1 = np.sum(D1 * counts1)
            d_sum2 = np.sum(D2 * counts2)
            d_mean1 = d_sum1 / n1
            d_mean2 = d_sum2 / n2
            d_mean_mean = (d_sum1 + d_sum2) / N

            stat = (N-2) * (n1 * (d_mean1-d_mean_mean)**2 + n2 * (d_mean2-d_mean_mean)**2) / (
                np.sum(counts1 * (D1-d_mean1)**2) + np.sum(counts2 * (D2-d_mean2)**2))
        
        case "hist_l2":
            hist1, _ = np.histogram(np.array(list(sfs1.keys()), dtype = float), bins = n_bins, range = (0,1), density = True, weights = list(sfs1.values()))
            hist2, _ = np.histogram(np.array(list(sfs2.keys()), dtype = float), bins = n_bins, range = (0,1), density = True, weights = list(sfs2.values()))
            stat = np.linalg.norm(hist2-hist1)

        case _:
            raise ValueError("Unknown statistic.")

    return stat

def test_statistic_list(sfs_list, fun):
    """
    Computes test statistic between consecutive site frequency spectrum dictionaries.
    
    Parameters
    ----------
    sfs_list : 1-D array-like
        contains site frequency spectrum as dictionaries
    fun : function
        a function which takes in two arguments and returns a float
        to compute the test statistic between two site frequency spectra
    
    Returns
    -------
    1-D numpy.ndarray
    """

    if len(sfs_list) < 2:
        raise ValueError("sfs_list must have length at least 2.")
    
    stats = [fun(sfs1, sfs2) for sfs1, sfs2 in zip(sfs_list[:-1], sfs_list[1:])]
    return np.array(stats)

# Parameter calculations
def kra_match_speed(delta, mu, gamma, Nss):
    """
    Computes the replicative advantage term needed to match SSD wave speed.
    
    Parameters
    ----------
    delta : float
        density, must be between 0 and 1
    mu : float
        degradation rate, must be between 0 and 1
    gamma : float
        diffusion rate, must be between 0 and 1
    Nss : float or int
        target population, must be larger than 0
        
    Returns
    -------
    kra : float
        replicative advantage term
    """
    return (1-delta) * np.sqrt(mu * gamma) / (Nss**(2/3))

def delta_match_speed(kra, mu, gamma, Nss):
    """
    Computes the density needed to match RA wave speed.
    
    Parameters
    ----------
    kra : float
        density, must be between 0 and 1
    mu : float
        degradation rate, must be between 0 and 1
    gamma : float
        diffusion rate, must be between 0 and 1
    Nss : float or int
        target population, must be larger than 0
        
    Returns
    -------
    kra : float
        replicative advantage term
    """
    return 1 - (kra * Nss**(2/3) / np.sqrt(mu * gamma))

### Example: simple site frequency spectrum plot from simulation 4

if __name__ == "__main__":

    # Go to correct directory
    import os
    os.chdir("../data/parameters_set1")

    # Load data
    sim = 3
    sfs_path = f"wildtype_sim_site_frequency_spectrum{sim}.txt"
    pop_path = f"wildtype_sim_populations{sim}.txt"
    sfs_df = _get_sfs_data(sfs_path, sim)
    pop_df = _get_populations_data(pop_path, sim, just_wildtype = True)

    # Plot site frequency spectrum at specified t
    plot_at_t = [100, 2000, 4000, 8000, 10000]
    sfs_at_t = sfs_df[sfs_df["t"].isin(plot_at_t)]["sfs"].to_list()
    pop_sum_at_t = pop_df[pop_df["t"].isin(plot_at_t)].groupby(["sim", "t"], as_index = False)["wildtype_population"].sum()["wildtype_population"].to_numpy()
    fig, h_axes = plt.subplots()
    for i, t in enumerate(plot_at_t):
        plot_sfs(h_axes, sfs_at_t[i], pop_sum_at_t[i], label = f"t = {t}", type = "line")
    h_axes.set_xlim((0,1))
    h_axes.set_xlabel("heteroplasmy")
    h_axes.set_ylabel("frequency")
    h_axes.legend()
    h_axes.grid()
    plt.show()