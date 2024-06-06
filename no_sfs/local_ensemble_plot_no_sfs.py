import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def _ffill(df):
    """
    Forward fill simulations data frame.

    Parameters
    ----------
    df : pandas.core.frame.DataFrame

    Returns
    -------
    pandas.core.frame.DataFrame
    """

    sims = df["sim"].unique()
    ts = sorted(df["t"].unique())

    # Template data frame with pre-set indices
    new_index = pd.MultiIndex.from_product([sims, ts], names=["sim", "t"])
    template_df = pd.DataFrame(index = new_index)

    # Forward fill
    filled_df = pd.merge(template_df, df, left_index=True, right_on=["sim", "t"], how="left")
    filled_df = filled_df.ffill()

    return filled_df


# SSS
fig, ax = plt.subplots(figsize = (5,5))
seed = 1
df = pd.read_csv(f"sss_sim_no_sfs_1unit_{seed}.txt")
# df = _ffill(df)
df["h"] = df["m"] / (df["w"] + df["m"])
df = df.groupby("t").agg(
    w1_mean = ("w", "mean"), w1_sem = ("w", "sem"),
    m1_mean = ("m", "mean"), m1_sem = ("m", "sem"),
    h_mean = ("h", "mean"), h_sem = ("h", "sem"), n_sims = ("h", "count")
).reset_index()
ax.plot(df["t"], df["h_mean"], linewidth = 1)
ax.fill_between(df["t"], df["h_mean"]-df["h_sem"], df["h_mean"]+df["h_sem"], alpha = 0.25)
ax.set_xlabel("t")
ax.set_ylabel(r"$\langle f \rangle$")
ax.grid()
plt.tight_layout()
plt.savefig("../figures/sss_no_sfs_1unit.png")
plt.show()

# RA
fig, ax = plt.subplots(figsize = (5,5))
seed = 1
df = pd.read_csv(f"ra_sim_no_sfs_1unit_{seed}.txt")
# df = _ffill(df)
df["h"] = df["m"] / (df["w"] + df["m"])
df = df.groupby("t").agg(
    w1_mean = ("w", "mean"), w1_sem = ("w", "sem"),
    m1_mean = ("m", "mean"), m1_sem = ("m", "sem"),
    h_mean = ("h", "mean"), h_sem = ("h", "sem"), n_sims = ("h", "count")
).reset_index()
ax.plot(df["t"], df["h_mean"], linewidth = 0.8, label = "stochastic")
ax.fill_between(df["t"], df["h_mean"]-df["h_sem"], df["h_mean"]+df["h_sem"], alpha = 0.25)
# Plot deterministic logistic growth of h
ts = np.linspace(0, 50, 100)
h0 = 0.02
kra = 0.25
ax.plot(ts, h0*np.exp(kra*ts)/(1-h0 + h0*np.exp(kra*ts)), color = "black", linestyle = "--", label = "deterministic")
ax.plot()
ax.set_xlabel("t")
ax.set_ylabel(r"$\langle f \rangle$")
ax.legend()
ax.grid()
plt.tight_layout()
plt.savefig("../figures/ra_no_sfs_1unit.png")
plt.show()


# 1 UNIT
fig, ax = plt.subplots(figsize = (5,5))
seed = 10
df = pd.read_csv(f"ssd_sim_no_sfs_1unit_{seed}.txt")
# df = _ffill(df)
df["h"] = df["m"] / (df["w"] + df["m"])
df = df.groupby("t").agg(
    w1_mean = ("w", "mean"), w1_sem = ("w", "sem"),
    m1_mean = ("m", "mean"), m1_sem = ("m", "sem"),
    h_mean = ("h", "mean"), h_sem = ("h", "sem"), n_sims = ("sim", "count")
).reset_index()
ax.plot(df["t"], df["h_mean"], linewidth = 1)
ax.fill_between(df["t"], df["h_mean"]-df["h_sem"], df["h_mean"]+df["h_sem"], alpha = 0.25)
ax.set_xlabel("t")
ax.set_ylabel(r"$\langle f \rangle$")
ax.grid()
ax.legend()
plt.tight_layout()
plt.savefig("../figures/ssd_no_sfs_1unit.png")
plt.show()

### 2 UNIT

colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
fig, axes = plt.subplots(1, 2, figsize = (10,5))
gammas = [f"{gamma:.2e}" for gamma in np.logspace(-4, 4, 5)]
seeds = [np.arange(1, 6), np.arange(101, 105)]
for i, ax in enumerate(axes):
    for gamma, seed, color in zip(gammas, seeds[i], colours):
        if i==0:
            df = pd.read_csv(f"ssd_sim_no_sfs_2unit_{gamma}_{seed}.txt")
        else:
            df = pd.read_csv(f"ssd_sim_euler_{gamma}_{seed}.txt")
        # df = _ffill(df)
        df["h"] = (df["m1"] + df["m2"]) / (df["w1"] + df["m1"] + df["w2"] + df["m2"])
        df = df.groupby("t").agg(
            w1_mean = ("w1", "mean"), w1_sem = ("w1", "sem"),
            m1_mean = ("m1", "mean"), m1_sem = ("m1", "sem"),
            w2_mean = ("w2", "mean"), w2_sem = ("w2", "sem"),
            m2_mean = ("m2", "mean"), m2_sem = ("m2", "sem"),
            h_mean = ("h", "mean"), h_sem = ("h", "sem"), n_sims = ("h", "count")
        ).reset_index()
        if i==0:
            # m = np.polyfit(df["t"], df["h_mean"], 1)[0]
            m = np.linalg.lstsq(np.reshape(df["t"], (len(df["t"]), 1)), df["h_mean"]-0.5)[0][0]
            ax.plot(df["t"], df["h_mean"], linewidth = 1, label = r"$\gamma=" + f"{float(gamma):.0e}," + f"m = {m:.2e}$")
        else:
            if float(gamma) in [1e-4, 1e-2]:
                m = np.linalg.lstsq(np.reshape(df["t"], (len(df["t"]), 1)), df["h_mean"]-0.5)[0][0]
                ax.plot(df["t"], df["h_mean"], linewidth = 1, label = r"$\gamma=" + f"{float(gamma):.0e}," + f"m = {m:.2e}$")
            else:
                ax.plot(df["t"], df["h_mean"], linewidth = 1, label = r"$\gamma=" + f"{float(gamma):.0e}" + "$")
        ax.fill_between(df["t"], df["h_mean"]-df["h_sem"], df["h_mean"]+df["h_sem"], color = color, alpha = 0.1)
    ax.set_title(["A", "B"][i])
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\langle f \rangle$")
    ax.set_ylim(0.49, 0.53)
    ax.grid()
    ax.legend()
plt.tight_layout()
plt.savefig("../figures/ssd_no_sfs_2unit.png")
plt.show()
