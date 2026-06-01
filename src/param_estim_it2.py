
import numpy as np
import h5py
from scipy.interpolate import interp1d
from datetime import datetime
from it4_anal_solver import B_funky


# ---------------------------------------------------------------------------
# Cost function
# ---------------------------------------------------------------------------

def log_space_cost(model, obs_val, weights=None, eps=1e-6):
    """
    Weighted misfit in log space.

    weights : array-like, same length as obs_val, or None (uniform weight=1).
              Higher weight → larger penalty for misfit at that observation point.
              The array is used as-is (not normalised), so scale it however is
              convenient (e.g. values in [0, 1], or raw inverse-variance).
    """
    if not np.all(np.isfinite(model)):
        return np.inf
    residuals = (np.log(model + eps) - np.log(obs_val + eps)) ** 2
    if weights is not None:
        residuals = residuals * np.asarray(weights)
    return np.sum(residuals)


# ---------------------------------------------------------------------------
# Proposal distributions
# ---------------------------------------------------------------------------

def propose_lognormal(current, sigma, lo, hi):
    """
    Log-normal proposal: exp(Normal(log(current), sigma)).
    Strictly positive by construction; returns None if outside [lo, hi].
    """
    candidate = current * np.exp(sigma * np.random.randn())
    if lo <= candidate <= hi:
        return candidate
    return None


def propose_uniform(current, half_width, lo, hi):
    """
    Symmetric uniform proposal: Uniform(current - half_width, current + half_width).
    half_width plays the same role as sigma in the log-normal case.
    Returns None if outside [lo, hi].
    """
    candidate = current + np.random.uniform(-half_width, half_width)
    if lo <= candidate <= hi:
        return candidate
    return None


def propose(current, cfg):
    """
    Dispatch to the correct proposal based on cfg['proposal'].

    cfg must contain:
      'proposal' : 'lognormal' (default) or 'uniform'
      'sigma'    : step-size parameter (std for log-normal, half-width for uniform)
      'lo', 'hi' : hard bounds
    """
    proposal_type = cfg.get('proposal', 'lognormal')
    lo, hi = cfg['lo'], cfg['hi']
    sigma  = cfg['sigma']

    if proposal_type == 'lognormal':
        return propose_lognormal(current, sigma, lo, hi)
    elif proposal_type == 'uniform':
        return propose_uniform(current, sigma, lo, hi)
    else:
        raise ValueError(f"Unknown proposal type '{proposal_type}'. Use 'lognormal' or 'uniform'.")


# ---------------------------------------------------------------------------
# MCMC
# ---------------------------------------------------------------------------

def run_mcmc(obs_z, obs_val, init_p, param_cfg, num_iterations, output_name: str,
             weights=None):
    """
    Metropolis-Hastings MCMC for PDE parameter estimation.

    Parameters
    ----------
    obs_z          : observed depths
    obs_val        : observed concentrations
    init_p         : dict with all fixed solver inputs
    param_cfg      : dict keyed by param name, each value is a dict with keys:
                       'init'     : float  — starting value
                       'sigma'    : float  — proposal step size
                       'lo'/'hi'  : float  — hard bounds
                       'proposal' : str    — 'lognormal' (default) or 'uniform'
    num_iterations : int
    output_name    : str  — base filename for the .h5 output (timestamp appended)
    weights        : array-like or None
                     Per-observation penalty weights for the cost function.
                     Must match len(obs_val). None = uniform weights.
                     Example: np.linspace(1, 5, len(obs_val)) punishes errors
                     at deeper observations more heavily.

    Returns
    -------
    state       : dict of final parameter values
    output_name : str, the actual file path written (includes timestamp)
    """
    obs_z   = np.asarray(obs_z,   dtype=float)
    obs_val = np.asarray(obs_val, dtype=float)

    stem, _, ext = output_name.rpartition('.')
    timestamp   = datetime.now().strftime("%y_%m_%d_%H_%M")
    output_name = f"{stem}_{timestamp}.{ext}" if stem else f"{output_name}_{timestamp}"

    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        if weights.shape != np.asarray(obs_val).shape:
            raise ValueError("weights must have the same shape as obs_val.")

    z_grid  = init_p['z_grid']
    t       = init_p['t']
    t_array = init_p['t_array']

    def forward(**params):
        _u1    = params.get('u1',    init_p.get('u1'))
        _u2    = params.get('u2',    init_p.get('u2'))
        _alpha = params.get('alpha', init_p.get('alpha'))
        _Da    = params.get('Da',    init_p.get('Da'))
        _kappa = params.get('kappa', init_p.get('kappa'))
        _C     = params.get('C',    init_p.get('C'))

        B = B_funky(z_grid, t, _u1, _u2, _alpha, _Da, _kappa, _C, t_array)
        if not np.all(np.isfinite(B)):
            return np.full(len(obs_z), np.nan)
        return interp1d(z_grid, B, kind="linear", bounds_error=False, fill_value="extrapolate")(obs_z)

    state     = {k: v['init'] for k, v in param_cfg.items()}
    model     = forward(**state)
    cost_curr = log_space_cost(model, obs_val, weights)
    if not np.isfinite(cost_curr):
        raise ValueError(
            f"Initial parameters produce non-finite cost ({cost_curr}). "
            "Check that the solver is stable for the starting parameter values."
        )

    n_params    = len(param_cfg)
    param_names = list(param_cfg.keys())

    # Store the proposal types used so they are recoverable from the file
    proposal_types = [param_cfg[k].get('proposal', 'lognormal') for k in param_names]

    with h5py.File(output_name, "w") as f:
        f.create_dataset("cost",      (num_iterations,),           dtype='f8')
        f.create_dataset("params",    (num_iterations, n_params),  dtype='f8')
        f.create_dataset("accepted",  (num_iterations,),           dtype='i1')
        f.attrs['param_names']    = param_names
        f.attrs['proposal_types'] = proposal_types
        if weights is not None:
            f.create_dataset("weights", data=weights, dtype='f8')

    n_accepted = 0

    with h5py.File(output_name, "r+") as f:
        cost_dset     = f["cost"]
        params_dset   = f["params"]
        accepted_dset = f["accepted"]

        for i in range(num_iterations):
            if i % 10_000 == 0:
                rate = n_accepted / max(i, 1)
                print(f"  iter {i:6d}  cost={cost_curr:.4f}  accept rate={rate:.2f}")

            candidate = dict(state)
            valid = True
            for name in param_names:
                prop = propose(state[name], param_cfg[name])
                if prop is None:
                    valid = False
                    break
                candidate[name] = prop

            if valid:
                model_cand = forward(**candidate)
                cost_cand  = log_space_cost(model_cand, obs_val, weights)

                log_alpha = -(cost_cand - cost_curr)
                accept    = np.log(np.random.rand()) < log_alpha
            else:
                accept = False

            if accept:
                state     = candidate
                cost_curr = cost_cand
                n_accepted += 1

            cost_dset[i]     = cost_curr
            params_dset[i]   = [state[k] for k in param_names]
            accepted_dset[i] = int(accept)

    print(f"\nDone. Acceptance rate: {n_accepted / num_iterations:.3f}  →  saved to {output_name}")
    return state, output_name


# ---------------------------------------------------------------------------
# Chain loading
# ---------------------------------------------------------------------------

def load_chain(path, burnin=20):
    with h5py.File(path, "r") as f:
        params   = f["params"][:]
        accepted = f["accepted"][:]
        costs    = f["cost"][:]
        names    = list(f.attrs["param_names"])
        weights  = f["weights"][:] if "weights" in f else None
    mask  = accepted[burnin:].astype(bool)
    chain = params[burnin:][mask]
    return chain, costs[burnin:][mask], names, weights
