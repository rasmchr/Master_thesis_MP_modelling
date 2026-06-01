import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from it4_anal_solver import B_funky

sys.path.insert(0, r"C:\Users\Træning 3\OneDrive - Roskilde Universitet\12. semester\Speciale - General\Numerik\Scripts\MCMC\it4")
from param_estim_it4_alteration import load_chain

#FILE1 = r"C:\Users"
#FILE2 = r"C:\Users"
File3 = r"C:\Users\"
BASE  = r"C:\Users\"
OUT = Path(BASE) / (Path(File3).stem + "_plots")
OUT.mkdir(exist_ok=True)

alpha = np.mean([1.74, 1.30, 0.87])
u1    = 3.65e5
u2    = 0.01


# N specifies the total number of MCMC iterations pre burn-in (across all chains). This is used to calculate the acceptance rate.
N = 5_000_000

L = 150.0

water_raw = [
    (2.5,   3.5667891682785307),
    (7.5,   2.6969926716060204),
    (12.5,  0.9080153525061379),
    (17.5,  0.5380893150797977),
    (27.5,  0.15000000000000002),
    (32.5,  0.17777777777777778),
    (42.5,  0.07500000000000001),
    (47.5,  0.1),
]
sediment_raw = [
    (0.0125,               634.7414656833154),
    (0.037500000000000006, 478.3568408220815),
    (0.0625,               970.1680149695655),
    (0.08750000000000001,  2343.5408338909124),
    (0.1125,               1420.2905042168202),
    (0.1375,               1237.2282190132125),
    (0.16250000000000003,  1026.581643277907),
    (0.1875,               299.62609220263096),
    (0.21250000000000002,  681.4776292062617),
    (0.2375,               458.20696863375105),
    (0.2625,               422.2530388157952),
    (0.3125,               640.1515738798312),
    (0.3875,               63.55031042845192),
    (0.41250000000000003,  187.55853823063677),
    (0.48750000000000004,  30.846954308574094),
    (0.5375000000000001,   30.516169700062356),
    (0.6125,               17.81769496827207),
]

water_z   = np.array([r[0]        for r in water_raw])
water_val = np.array([r[1]        for r in water_raw])
sed_z     = np.array([r[0] + 50.0 for r in sediment_raw])
sed_val   = np.array([r[1]        for r in sediment_raw])

pad_z   = np.arange(sed_z[-1] + 0.5, L + 0.5, 10)
pad_val = np.zeros(len(pad_z))

obs_z   = np.concatenate([water_z, sed_z, pad_z])
obs_val = np.concatenate([water_val, sed_val, pad_val])

t_final = 70.0
t_array = np.linspace(0, t_final, 2000)
z_grid  = np.linspace(0, 150, 2000)

DPI = 200


def save(fig, path):
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# --- Load chains ---
#chain,  costs,  names,  _ = load_chain(FILE2, burnin=5000)
#chain1, costs1, names1, _ = load_chain(FILE1, burnin=0)
chain3, costs2, names1, _ = load_chain(File3, burnin=50_000)

#combichains = np.concatenate([chain1, chain, chain3])
#combicosts  = np.concatenate([costs1, costs, costs3])
combichains = np.array(chain3)
combicosts  = np.array(costs2)

# --- Cost trace ---
fig, ax = plt.subplots(figsize=(12, 3))
ax.plot(combicosts)
ax.set_ylabel('cost')
ax.set_xlabel('iteration (post-burnin)')
ax.set_title('Cost trace')
plt.tight_layout()
save(fig, OUT / "combi_cost_trace.png")

# --- Parameter traces ---
for i, name in enumerate(names1):
    fig, ax = plt.subplots(figsize=(12, 3))
    ax.plot(combichains[:, i])
    ax.set_ylabel(name)
    ax.set_xlabel('accepted sample')
    if name == 'C':
        ax.set_title('c trace')
    else:
        ax.set_title(f'{name} trace')
    plt.tight_layout()
    save(fig, OUT / f"combi_{name}_trace.png")

# --- Marginal histograms ---
PARAM_COLORS = {'kappa': 'darkorange', 'Da': 'seagreen'}

for i, name in enumerate(names1):
    fig, ax = plt.subplots(figsize=(6, 4))
    color = PARAM_COLORS.get(name, 'tab:blue')
    ax.hist(combichains[:, i], bins=40, density=True, color=color, alpha=0.7)
    mean_val   = np.mean(combichains[:, i])
    median_val = np.median(combichains[:, i])
    lo95, hi95 = np.percentile(combichains[:, i], [2.5, 97.5])
    ax.axvline(mean_val,   color='black',  linestyle='--',  lw=1.4, label=f'mean   = {mean_val:.4g}')
    ax.axvline(median_val, color='red',    linestyle='-.',  lw=1.4, label=f'median = {median_val:.4g}')
    ax.axvline(lo95,       color='purple', linestyle=':',   lw=1.2, label=f'2.5%   = {lo95:.4g}')
    ax.axvline(hi95,       color='purple', linestyle=':',   lw=1.2, label=f'97.5%  = {hi95:.4g}')
    if name == 'C':
        ax.set_xlabel('c trace')
    else:
        ax.set_xlabel(name)
    ax.set_ylabel('density')
    ax.legend(fontsize=8)
    plt.tight_layout()
    save(fig, OUT / f"combi_{name}_hist.png")

# --- Cost histogram ---
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(combicosts, bins=40, density=True, color='dimgray', alpha=0.7)
mean_cost   = np.mean(combicosts)
median_cost = np.median(combicosts)
lo95_cost, hi95_cost = np.percentile(combicosts, [2.5, 97.5])
ax.axvline(mean_cost,   color='black',  linestyle='--',  lw=1.4, label=f'mean   = {mean_cost:.4g}')
ax.axvline(median_cost, color='red',    linestyle='-.',  lw=1.4, label=f'median = {median_cost:.4g}')
ax.axvline(lo95_cost,   color='purple', linestyle=':',   lw=1.2, label=f'2.5%   = {lo95_cost:.4g}')
ax.axvline(hi95_cost,   color='purple', linestyle=':',   lw=1.2, label=f'97.5%  = {hi95_cost:.4g}')
ax.set_xlabel('cost')
ax.set_ylabel('density')
ax.legend(fontsize=8)
plt.tight_layout()
save(fig, OUT / "combi_cost_hist.png")

# --- Pairwise scatter plots ---
from scipy.stats import spearmanr

param_pairs = [('kappa', 'Da')]

for name_x, name_y in param_pairs:
    ix  = names1.index(name_x)
    iy  = names1.index(name_y)
    x   = combichains[:, ix]
    y   = combichains[:, iy]
    rho, _ = spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, s=1, alpha=0.2, color='dimgray', rasterized=True)
    if name_x == 'C':
        ax.set_xlabel('c trace')
    else:
        ax.set_xlabel(name_x)

    if name_y == 'C':
        ax.set_ylabel('c trace')
    else:
        ax.set_ylabel(name_y)

    ax.text(0.05, 0.95, f'ρ = {rho:.3f}', transform=ax.transAxes,
            fontsize=10, verticalalignment='top')
    ax.grid(True, lw=0.4, alpha=0.5)
    plt.tight_layout()
    save(fig, OUT / f"scatter_{name_x}_{name_y}.png")

# --- Pairwise 2D histograms ---
for name_x, name_y in param_pairs:
    ix  = names1.index(name_x)
    iy  = names1.index(name_y)
    x   = combichains[:, ix]
    y   = combichains[:, iy]
    rho, _ = spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.hist2d(x, y, bins=60, density=True, cmap='viridis')
    if name_x == 'C':
        ax.set_xlabel('c trace')
    else:
        ax.set_xlabel(name_x)

    if name_y == 'C':
        ax.set_ylabel('c trace')
    else:
        ax.set_ylabel(name_y)
    ax.text(0.05, 0.95, f'ρ = {rho:.3f}', transform=ax.transAxes,
            fontsize=10, verticalalignment='top', color='white')
    plt.tight_layout()
    save(fig, OUT / f"hist2d_{name_x}_{name_y}.png")

# --- Posterior summaries ---
print("Posterior summaries (95% CI):")
for i, name in enumerate(names1):
    mean = np.mean(combichains[:, i])
    lo, hi = np.percentile(combichains[:, i], [2.5, 97.5])
    print(f"  {name:10s}: mean={mean:.4g}  95% CI=[{lo:.4g}, {hi:.4g}]")
    print(f"  {name:10s}: median={np.median(combichains[:, i]):.4g}")
    print(f"  {name:10s}: standard deviation={np.std(combichains[:, i]):.4f}")
    print(f" Relative CI width (median): {(hi - lo)/np.median(combichains[:, i])*100:.4f}%")
    print(f" Relative CI width (mean): {(hi - lo)/np.mean(combichains[:, i])*100:.4f}%")
    print(f"This is the acceptance rate: {len(combicosts) / N}")
    lo99, hi99 = np.percentile(combichains[:, i], [0.5, 99.5])
    print(f"  {name:10s},  99% CI=[{lo99:.4g}, {hi99:.4g}]")




# --- CI calculation ---
def get_ci_bounds(chain, names, chain1, names1, z_plot, t_final, t_plot,
                  alpha, u1, u2, n_samples=500, ci=95):
    lower_p = (100 - ci) / 2
    upper_p = 100 - lower_p
    if chain1 is not None:

        Da_samples    = np.concatenate([chain[:,  names.index('Da')],
                                        chain1[:, names1.index('Da')]])
        kappa_samples = np.concatenate([chain[:,  names.index('kappa')],
                                        chain1[:, names1.index('kappa')]])
        C_samples = np.concatenate([chain[:,  names.index('C')],
                                        chain1[:, names1.index('C')]])
    else:
        Da_samples    = chain[:, names.index('Da')]
        kappa_samples = chain[:, names.index('kappa')]
        C_samples = chain[:, names.index('C')]

    idx = np.random.choice(len(Da_samples), size=n_samples, replace=False)
    solutions = []
    for i in idx:
        sol = B_funky(z_plot, t_final, u1, u2, alpha,
                      Da_samples[i], kappa_samples[i], C_samples[i], t_plot)
        solutions.append(sol)

    solutions = np.array(solutions)
    lower = np.percentile(solutions, lower_p, axis=0)
    upper = np.percentile(solutions, upper_p, axis=0)
    return lower, upper


# --- Estimated solution ---
alpha_est = alpha
Da_est    = np.median(combichains[:, names1.index('Da')])
kappa_est = np.median(combichains[:, names1.index('kappa')])
C_est     = np.median(combichains[:, names1.index('C')])
u1_est    = u1
u2_est    = u2

print(Da_est, kappa_est)

z_plot = np.linspace(0, 150, 10000)
t_plot = np.linspace(0, t_final, 1000)

total_est = B_funky(z_plot, t_final, u1_est, u2_est, alpha_est,
                    Da_est, kappa_est, C_est, t_plot)

# --- Least squares error (median parameters vs. actual data) ---
data_z_all   = np.concatenate([water_z, sed_z])
data_val_all = np.concatenate([water_val, sed_val])
model_at_data  = np.interp(data_z_all, z_plot, total_est)
model_at_water = np.interp(water_z,    z_plot, total_est)
model_at_sed   = np.interp(sed_z,      z_plot, total_est)

lse       = np.sum((model_at_data  - data_val_all) ** 2)
lse_water = np.sum((model_at_water - water_val)    ** 2)
lse_sed   = np.sum((model_at_sed   - sed_val)      ** 2)
rmse      = np.sqrt(lse / len(data_val_all))
rmse_water = np.sqrt(lse_water / len(water_val))
rmse_sed   = np.sqrt(lse_sed / len(sed_val))

print(f"\nLeast squares error (median parameters):")
print(f"  Total LSE    = {lse:.6g}")
print(f"  Water LSE    = {lse_water:.6g}")
print(f"  Sediment LSE = {lse_sed:.6g}")
print(f"  Total RMSE   = {rmse:.6g}")
print(f"  Water RMSE   = {rmse_water:.6g}")
print(f"  Sediment RMSE = {rmse_sed:.6g}")

ci_lower, ci_upper = get_ci_bounds(chain=combichains, names=names1, z_plot=z_plot, t_final=t_final, t_plot=t_plot,
                                   alpha=alpha_est, u1=u1_est, u2=u2_est,
                                   chain1=None, names1=None,
                                   n_samples=10000, ci=95)

title_str = (
    f'α={alpha_est:.4f},  Da={Da_est:.4f},  κ={kappa_est:.4f},  '
    f'u1={u1_est:.4f},  u2={u2_est:.4f},  c={C_est:.4f}'
)


def make_panel(ax, xscale='linear', ylim=None, xlabel='Concentration'):
    ax.fill_betweenx(z_plot, ci_lower, ci_upper, alpha=0.3, color='gray', label='95% CI')
    ax.plot(total_est, z_plot, label='A + B (estimated)', color='tab:blue')
    ax.plot(ci_lower, z_plot, color='gray', linewidth=0.8, linestyle='--')
    ax.plot(ci_upper, z_plot, color='gray', linewidth=0.8, linestyle='--')
    ax.scatter(obs_val, obs_z, color='k', label='Observations', s=3, zorder=5)
    ax.axhline(y=50, linestyle='--', color='gray', label='Seabed')
    ax.set_xlabel(xlabel)
    ax.set_xscale(xscale)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.invert_yaxis()
    ax.grid()
    ax.legend(fontsize=8)


# Figure 1: Full domain 0-150 m
fig1, axes1 = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
make_panel(axes1[0], xscale='log', xlabel='Concentration (log scale)')
axes1[0].set_ylabel('Depth (m)')
axes1[0].set_title('Log x-axis')
make_panel(axes1[1], xscale='linear', xlabel='Concentration (linear scale)')
axes1[1].set_title('Linear x-axis')
#fig1.suptitle(f'Estimated solution — full domain (0–150 m)\n{title_str}', y=1.01)
plt.tight_layout()
save(fig1, OUT / "solution_full.png")

# Figure 2: Zoomed 0-60 m
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
make_panel(axes2[0], xscale='log', ylim=(0, 60), xlabel='Concentration (log scale)')
axes2[0].set_ylabel('Depth (m)')
axes2[0].set_title('Log x-axis')
make_panel(axes2[1], xscale='linear', ylim=(0, 60), xlabel='Concentration (linear scale)')
axes2[1].set_title('Linear x-axis')
#fig2.suptitle(f'Estimated solution — zoomed (0–60 m)\n{title_str}', y=1.01)
plt.tight_layout()
save(fig2, OUT / "solution_zoom_0_60.png")

# Figure 3: Zoomed 48-52 m
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
make_panel(axes3[0], xscale='log', ylim=(48, 52), xlabel='Concentration (log scale)')
axes3[0].set_ylabel('Depth (m)')
axes3[0].set_title('Log x-axis')
make_panel(axes3[1], xscale='linear', ylim=(48, 52), xlabel='Concentration (linear scale)')
axes3[1].set_title('Linear x-axis')
#fig3.suptitle(f'Estimated solution — zoomed (40–60 m)\n{title_str}', y=1.01)
plt.tight_layout()
save(fig3, OUT / "solution_zoom_48_52.png")

# --- Stitched figures: single depth axis, two x-scales via twiny ---
def make_twiny_figure(ylim, xscale='linear'):
    mask_w = z_plot <= 50.0
    mask_s = z_plot >= 50.0

    fig, ax_s = plt.subplots(figsize=(7, 10))
    ax_w = ax_s.twiny()  # water on top x-axis

    # Sediment compartment — bottom x-axis, orange
    ax_s.fill_betweenx(z_plot[mask_s], ci_lower[mask_s], ci_upper[mask_s],
                       alpha=0.3, color='tab:orange', label='95% CI (sediment)', zorder=1)
    ax_s.plot(total_est[mask_s], z_plot[mask_s], color='tab:orange', label='Solution (sediment)', zorder=3)
    ax_s.plot(ci_lower[mask_s], z_plot[mask_s], color='tab:orange', linewidth=0.8, linestyle='--', zorder=2)
    ax_s.plot(ci_upper[mask_s], z_plot[mask_s], color='tab:orange', linewidth=0.8, linestyle='--', zorder=2)
    ax_s.scatter(sed_val, sed_z, color='tab:orange', s=8, zorder=5, label='Obs (sediment)')

    # Water compartment — top x-axis, blue
    ax_w.fill_betweenx(z_plot[mask_w], ci_lower[mask_w], ci_upper[mask_w],
                       alpha=0.3, color='tab:blue', label='95% CI (water)', zorder=1)
    ax_w.plot(total_est[mask_w], z_plot[mask_w], color='tab:blue', label='Solution (water)', zorder=3)
    ax_w.plot(ci_lower[mask_w], z_plot[mask_w], color='tab:blue', linewidth=0.8, linestyle='--', zorder=2)
    ax_w.plot(ci_upper[mask_w], z_plot[mask_w], color='tab:blue', linewidth=0.8, linestyle='--', zorder=2)
    ax_w.scatter(water_val, water_z, color='tab:blue', s=8, zorder=5, label='Obs (water)')

    ax_s.axhline(y=50, linestyle='--', color='gray', linewidth=0.8, label='Seabed (z=50)')

    ax_s.set_xscale(xscale)
    ax_w.set_xscale(xscale)
    ax_s.set_ylim(ylim)
    ax_s.invert_yaxis()
    ax_s.set_ylabel('Depth (m)')
    ax_s.set_xlabel('Concentration — sediment', color='tab:orange')
    ax_s.tick_params(axis='x', colors='tab:orange')
    ax_w.set_xlabel('Concentration — water', color='tab:blue')
    ax_w.tick_params(axis='x', colors='tab:blue')
    ax_s.grid(alpha=0.3)

    lines_s, labels_s = ax_s.get_legend_handles_labels()
    lines_w, labels_w = ax_w.get_legend_handles_labels()
    ax_s.legend(lines_s + lines_w, labels_s + labels_w, fontsize=8,
                loc='center right', bbox_to_anchor=(1.0, 1/5))

    plt.tight_layout()
    return fig


twiny_configs = [
    {'name': '0_52',  'ylim': (0,  52)},
    {'name': '0_150', 'ylim': (0,  150)},
    {'name': '48_52', 'ylim': (48, 52)},
]

for cfg in twiny_configs:
    for xscale in ('linear', 'log'):
        fig = make_twiny_figure(ylim=cfg['ylim'], xscale=xscale)
        save(fig, OUT / f"twiny_{cfg['name']}_{xscale}.png")


# Source term A
fig_a = plt.figure(figsize=(6, 8))
A = C_est * np.exp(-(alpha_est / Da_est) * z_plot)
plt.plot(A, z_plot, 'k-.', label='A (source)')
plt.gca().invert_yaxis()
plt.xlabel('Concentration')
plt.ylabel('Depth (m)')
plt.legend()
plt.grid()
plt.tight_layout()
save(fig_a, OUT / "source_term_A.png")

print(C_est * np.exp(-(alpha_est / Da_est) * 50))
