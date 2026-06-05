# Vertical Transport of Microplastics in Marine Sediments

A one-dimensional reaction–advection–diffusion model for the vertical transport and
sedimentary accumulation of microplastics (MP), with Bayesian parameter estimation
via Markov-chain Monte Carlo using the Metropolis–Hastings algorithm.

> MSc thesis project — Mathematical Bioscience, Roskilde University (RUC), 2026.
> Authors: Rasmus Willaume Christiansen and Sif Egelund

---

## Overview

This repository contains the model implementation, numerical solvers, the
parameter-estimation pipeline, and the analysis code for a master's thesis on the
vertical transport of microplastics in marine environments.

The aim is to reproduce, with a proposed model, the qualitative concentration profile
observed in field data: an exponential decrease through the water column and a sharp
accumulation just below the water–sediment interface.

The study region covers the Baltic Sea, the inner Danish waters, and the North Sea,
using a compiled dataset drawn from seven independent field studies.

---

## The model(s)

The marine environment is represented as a one-dimensional vertical domain, with the
`z`-axis pointing **downward** from the sea surface. It is split into two coupled
compartments:

- **Water column** (`0 ≤ z ≤ L₁`, with `L₁ = 50 m`)
- **Sediment** (`L₁ ≤ z ≤ L₂`, with `L₂ = 150 m`)


The processes parametrised in the model are:

| # | Process            | Type      | Parameter |
|---|--------------------|-----------|-----------|
| 1 | Buoyancy           | Advection | `S_b`     |
| 2 | Sedimentation rate | Advection | `S_r`     |
| 3 | Biofouling         | Reaction  | `α`       |
| 4 | Bioturbation       | Diffusion | `D_b`     |
| 5 | Bioirrigation      | Diffusion | `D_i`     |
| 6 | Decay              | Reaction  | `κ`       |

<img width="622" height="602" alt="Compartment diagram" src="https://github.com/user-attachments/assets/2f80e6c7-9a9e-4f31-bc32-d285f072eb1d" />

*Schematic overview of the model domain and the fundamental processes governing vertical
MP transport. The water column (blue) spans from $z = 0$ to $z = L_1$, and the sediment
(beige) from $z = L_1$ to $z = L_2$. Particle $A$ (neutrally buoyant) undergoes diffusion
$D_a$ in the water column, indicated by the black bidirectional arrows. The green dashed
arrows indicate the biofouling reaction $\alpha$, converting particle $A$ into the
negatively buoyant particle $B$. The red dashed arrows indicate advective processes:
$S_b$ denotes the settling velocity of particle $B$ through the water column, and $S_r$
denotes the sedimentation rate of non-MP material progressively burying $B$ particles
within the sediment. The dashed black box indicates the bioturbation zone in the
uppermost layer of the sediment, within which the purple bidirectional arrows indicate
diffusion due to bioturbation $D_b$ and bioirrigation $D_i$, acting on $B$ particles. The
blue dashed arrow indicates chemical decay $\kappa$, acting on $B$ particles throughout
the sediment. The arrow to the left of the compartment diagram indicates the direction of
the coordinate system.*

Two model iterations are implemented:

1. **Iteration 1** — a spatially dependent advection formulation that captures the
   interface dynamics directly, including special cases for constant and
   sigmoid-shaped advection. Solved numerically.
2. **Iteration 2** — a two-compartment formulation in which the water column and the
   sediment are modelled separately and coupled through a flux condition at the
   interface. This admits a closed-form analytical solution and is therefore used for
   parameter estimation, where its low computational cost is essential.

The relative importance of advection, diffusion, and reaction is characterised through
the **Péclet** and **Damköhler** numbers; the system is found to be advection-dominated,
with reaction weak in magnitude but necessary to reproduce the observed profile shape.

---

## Parameter estimation

Only the second iteration is used for parameter estimation. Using the method of
characteristics, closed-form analytical expressions can be derived for each compartment,
and for simplicity these expressions are the ones used for estimation.

Three parameters, identified as the least constrained by the literature, are estimated:

- **`c`** — surface microplastic concentration; boundary condition at `z = 0` (particles L⁻¹)
- **`κ`** — effective removal / degradation rate (year⁻¹)
- **`D_a`** — effective diffusion / transport coefficient (m² year⁻¹)

Estimation uses the **Metropolis–Hastings** algorithm with uniform priors. Convergence is
assessed with the **Gelman–Rubin diagnostic** (`R̂`), using multiple chains of 5 million
iterations started from different initial conditions, with a burn-in discard.

Two cost-function weighting schemes are compared:

- **Globally normalised** — normalises by the maximum concentration in the dataset.
  Achieves the lowest RMSE, but leaves `c` and `D_a` poorly identified and produces a
  physically implausible steady-state profile for `A`.
- **Aggressive** — weights every data point relative to the peak of its own compartment.
  Constrains all three parameters and yields a physically realistic profile for `A`, at
  the cost of a higher RMSE.

Parameter interactions and identifiability are examined through pairwise correlation
plots and a **sloppy-model analysis** (eigenvalue / eigenvector decomposition of the
cost-function landscape).

---

## Repository structure

```
.
├── src/                  # model, solvers, MCMC source, MCMC plots
│   ├── B1_ReacAdv_twocomp_solver_implicit_v03.py         # Implicit numerical PDE solver (iteration 2)
│   ├── Gelman-Rubin.py                                   # Gelman Rubin diagnostic chains
│   ├── Spearmans_corr_plots.py                           # Calculates and plots the Spearmans correlations between parameters
│   ├── chain_plot.py                                     # Plots different characteristics about the chains
│   ├── Weights.py                                        # Calculates and plots each of the weighting schemes utilised in the these
│   ├── it4_mcmc.ipynb                                    # Performs the Bayesian Inference on an idealized parameter set
│   ├── it4_mcmc_fixed_C.ipynb                            # Performs the Bayesian Inference on an idealized parameter set while having C fixed.
│   ├── it4_mcmc_nchains.ipynb                            # Performs the Bayesian Inference on an ideliazied parameter set with initial conditions for the estimated parameters being random.      
│   ├── it2_analytcal_twocomp_solver.py                   # Solves the analystical equations of the two compartment model
│   └── param_estim_it2.py                                # performs the parameter estimation of the using the analytical two compartment solver.
├── data/                 # compiled field dataset along with plotting scripts
│   ├── Binned_mean_concentration.csv
│   ├── Data workflow diagram.jpeg
│   ├── Master data_v01.xlsx
│   ├── plot.py
│   └── plot_v1.py       
├── notebooks/            # exploratory analysis of each iteration
│   ├── AR_VS_RAD.ipynb
│   ├── IT1.ipynb
│   ├── IT2.ipynb
│   ├── IT3_fullADV.ipynb
│   ├── IT3_linearADV.ipynb
│   ├── IT3_sigmoidADV.ipynb
│   ├── IT4.ipynb
│   ├── OldANew.ipynb
│   ├── advection_tries_and_characteristics.ipynb
│   ├── it3_sig_damkohler.ipynb       
│   └── plot_diff.ipynb
└── README.md
```

---

## Data

The dataset is compiled from seven independent field studies spanning the Baltic Sea,
the North Sea, and the inner Danish waters. Because the studies use different sampling
methods, units, and spatial resolutions, considerable data harmonisation was required,
which contributes to the observational uncertainty. The master dataset is included as an
`.xlsx` file in the `data/` folder.

---

## Authors & acknowledgements

**Authors:**

- Rasmus Willaume Christiansen — MSc Mathematical Bioscience, Roskilde University (RUC). Willaras@outlook.dk
- Sif Egelund — MSc Mathematical Bioscience, Roskilde University (RUC).

**Supervisor:** Jesper Schmidt, Roskilde University.
