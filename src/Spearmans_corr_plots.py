# This file will load all MCMC chains in this folder.
# Then it will concattenate all the data and compute the spearmans correlation coefficient between each of the parameters.
# The results will be plotted in a 2d histogram heat map.
# The correlation coefficient will be shown on each of the plots.
# The plots will be saved in individual PNG files and saved to the same folder as this file.
# The PNG files will be named "spearmans_corr_{param1}_{param2}_{date}.png" and saved in a folder named "spearmans_corr_plot_{date}"
# Also a csv file will be saved with the same name as the folder inside the folder which contains.
    # Columns will be date, chain file name, number of samples in chain

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import combinations
from scipy.stats import spearmanr, binned_statistic_2d
from datetime import datetime

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
from param_estim_it4_alteration import load_chain

BURNIN = 0
DPI    = 200

date_str   = datetime.now().strftime("%Y_%m_%d")
out_folder = os.path.join(script_dir, f"spearmans_corr_plot_{date_str}")
os.makedirs(out_folder, exist_ok=True)

chain_files = sorted(
    os.path.join(script_dir, f)
    for f in os.listdir(script_dir)
    if f.endswith(".h5")
)

if not chain_files:
    raise FileNotFoundError("No .h5 chain files found in the script directory.")

all_chains = []
all_costs  = []
csv_rows   = []
names      = None

for fpath in chain_files:
    chain, costs, file_names, _ = load_chain(fpath, burnin=BURNIN)
    if names is None:
        names = file_names
    all_chains.append(np.array(chain))
    all_costs.append(np.array(costs))
    csv_rows.append({
        "date":              date_str,
        "chain_file_name":   os.path.basename(fpath),
        "number_of_samples": len(chain),
    })
    print(f"Loaded {os.path.basename(fpath)}: {len(chain)} accepted samples")

combined        = np.concatenate(all_chains, axis=0)
combined_costs  = np.concatenate(all_costs,  axis=0)
cost_min, cost_max = combined_costs.min(), combined_costs.max()
combined_costs_norm = (combined_costs - cost_min) / (cost_max - cost_min)
print(f"\nTotal samples after concatenation: {len(combined)}")

df = pd.DataFrame(csv_rows)
csv_path = os.path.join(out_folder, f"spearmans_corr_plot_{date_str}.csv")
df.to_csv(csv_path, index=False)
print(f"Saved CSV: {os.path.basename(csv_path)}")

for i, j in combinations(range(len(names)), 2):
    name_x = names[i]
    name_y = names[j]
    x      = combined[:, i]
    y      = combined[:, j]

    rho, pval = spearmanr(x, y)

    BINS = 1000
    stat, xedges, yedges, _ = binned_statistic_2d(
        x, y, combined_costs_norm, statistic="mean", bins=BINS
    )

    fig, ax = plt.subplots(figsize=(5, 5))
    mesh = ax.pcolormesh(xedges, yedges, stat.T, cmap="viridis_r", vmin=0, vmax=1)
    plt.colorbar(mesh, ax=ax, label="mean cost (normalised)")
    if name_x == 'C':
        ax.set_xlabel('c trace')
    else:
        ax.set_xlabel(name_x)
    if name_y == 'C':
        ax.set_ylabel('c trace')
    else:
        ax.set_ylabel(name_y)
    #ax.set_title(f"Spearman $\\rho$ = {rho:.4f}")
    ax.text(
        0.05, 0.95, f"$\\rho$ = {rho:.4f}",
        transform=ax.transAxes,
        fontsize=10, verticalalignment="top", color="black",
    )
    plt.tight_layout()

    fname    = f"spearmans_corr_{name_x}_{name_y}_{date_str}.png"
    out_path = os.path.join(out_folder, fname)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fname}")

print("\nDone.")
