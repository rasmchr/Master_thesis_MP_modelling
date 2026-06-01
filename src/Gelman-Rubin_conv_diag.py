import os
import glob
import numpy as np
from param_estim_it4_alteration import load_chain


def gelman_rubin(*chains):
    """
    Compute the Gelman-Rubin R-hat statistic for n chains of
    potentially different lengths. Chains are trimmed to the
    length of the shortest chain.

    Parameters
    ----------
    *chains : variable number of 1D arrays, one per chain
              e.g. gelman_rubin(s1, s2, s3)

    Returns
    -------
    R_hat : float
        Values close to 1.0 indicate convergence.
        Typically R-hat < 1.1 is considered acceptable.
    """
    # Trim all chains to the length of the shortest
    min_len = min(len(c) for c in chains)
    chains  = np.array([c[:min_len] for c in chains])

    m = chains.shape[0]  # number of chains
    n = chains.shape[1]  # number of samples per chain (after trimming)

    # Within-chain means and variances
    chain_means = np.mean(chains, axis=1)        # shape (m,)
    chain_vars  = np.var(chains,  axis=1, ddof=1) # shape (m,)

    # Between-chain variance B
    grand_mean = np.mean(chain_means)
    B = n / (m - 1) * np.sum((chain_means - grand_mean) ** 2)

    # Within-chain variance W
    W = np.mean(chain_vars)

    # Pooled variance estimate
    var_hat = (1 - 1/n) * W + (1/n) * B

    # R-hat
    R_hat = np.sqrt(var_hat / W)
    return R_hat


# ── Load all chains in this folder ──────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
h5_files   = sorted(glob.glob(os.path.join(script_dir, "*.h5")))

if len(h5_files) < 2:
    raise RuntimeError(
        f"Need at least 2 chain files; found {len(h5_files)} in {script_dir}"
    )

chains_data = []
names       = None

for fpath in h5_files:
    chain, costs, param_names, weights = load_chain(fpath)
    chains_data.append(chain)
    names = param_names          # same for every chain
    print(f"Loaded: {os.path.basename(fpath)}  →  {chain.shape[0]} accepted samples")

# ── Gelman-Rubin per parameter ───────────────────────────────────────────────
print(f"\nGelman-Rubin R-hat  ({len(chains_data)} chains)")
print("-" * 35)
for i, name in enumerate(names):
    param_chains = [cd[:, i] for cd in chains_data]
    R_hat = gelman_rubin(*param_chains)
    flag  = "" if R_hat < 1.1 else "  ← not converged"
    print(f"  {name:10s}  R-hat = {R_hat:.4f}{flag}")
