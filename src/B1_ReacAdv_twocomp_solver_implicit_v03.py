import numpy as np
from numba import njit


@njit
def thomas_solve(lower, main, upper, rhs, n):
    """
    Thomas algorithm (tridiagonal matrix algorithm).

    Solves  diag(lower,−1) + diag(main,0) + diag(upper,+1) ) x = rhs.
    Modifies main and rhs in-place; lower and upper are read-only.
    Returns rhs (the solution).
    """
    for i in range(1, n):
        w = lower[i] / main[i - 1]
        main[i] -= w * upper[i - 1]
        rhs[i]  -= w * rhs[i - 1]

    rhs[n - 1] /= main[n - 1]
    for i in range(n - 2, -1, -1):
        rhs[i] = (rhs[i] - upper[i] * rhs[i + 1]) / main[i]

    return rhs


@njit(fastmath=True)
def two_compartment_implicit(
    dz, dt, time_steps, Nx, z_b,
    C, alpha, Da, kappa,
    u1, u2,
    IC_B, M_init
):
    """
    Two-compartment implicit upwind solver with ODE mass transfer at the interface.

    Grid: z[j] = j * dz,  j = 0 .. Nx-1.
    Interface snapped to nearest node: i_b = round(z_b / dz).

    Domain 1 : cells 0 .. i_b       velocity u1  (downward, > 0)
    Domain 2 : cells i_b+1 .. Nx-1  velocity u2  (downward, > 0)

    Interface ODE (backward Euler each step):
        dM/dt = u1 * B[i_b] - u2 * M   →   M_new = (M + dt*u1*B[i_b]) / (1 + dt*u2)

    Boundary conditions:
        B[0]     = 0   Dirichlet inlet, domain 1
        B[i_b+1] = M   Dirichlet inlet, domain 2
        outflow at z = L (upwind stencil is naturally zero-gradient)

    PDE (both domains):
        dB/dt + U dB/dz = alpha * A - kappa * B,   A(z) = C * exp(-sqrt(alpha/Da) * z)

    Parameters
    ----------
    dz, dt       : float
    time_steps   : int
    Nx           : int    – total cells  (L = (Nx-1)*dz)
    z_b          : float  – interface coordinate; snapped to nearest node
    C, alpha, Da : float  – source-profile parameters
    kappa        : float  – first-order decay rate (both domains)
    u1, u2       : float  – advection velocities (both > 0)
    IC_B         : (Nx,)  – initial concentration
    M_init       : float  – initial interface ODE value

    Returns
    -------
    z  : (Nx,)  – spatial grid
    A  : (Nx,)  – source profile
    B  : (Nx,)  – concentration at final time
    M  : float  – interface variable at final time
    """
    i_b = int(z_b / dz + 0.5)

    B = IC_B.copy()
    M = M_init

    z     = np.arange(0, Nx) * dz
    A     = C * np.exp(-np.sqrt(alpha / Da) * z)
    inv_dz = 1.0 / dz

    N1 = i_b + 1        # domain 1: cells 0..i_b
    N2 = Nx - i_b - 1   # domain 2: cells i_b+1..Nx-1

    lo1 = np.zeros(N1); di1 = np.zeros(N1)
    up1 = np.zeros(N1); rh1 = np.zeros(N1)
    lo2 = np.zeros(N2); di2 = np.zeros(N2)
    up2 = np.zeros(N2); rh2 = np.zeros(N2)

    # lo and up are constant throughout — fill once.
    # di and rh are rebuilt each step (thomas_solve overwrites di in-place).
    lo1[0] = 0.0;  up1[0] = 0.0
    for j in range(1, N1):
        lo1[j] = -u1 * inv_dz
        up1[j] = 0.0

    lo2[0] = 0.0;  up2[0] = 0.0
    for j in range(1, N2):
        lo2[j] = -u2 * inv_dz
        up2[j] = 0.0

    diag1 = 1.0 / dt + u1 * inv_dz + kappa
    diag2 = 1.0 / dt + u2 * inv_dz + kappa

    for _n in range(time_steps):

        # --- Domain 1 ---
        di1[0] = 1.0;  rh1[0] = 0.0        # Dirichlet B = 0
        for j in range(1, N1):
            di1[j] = diag1
            rh1[j] = B[j] / dt + alpha * A[j]
        rh1 = thomas_solve(lo1, di1, up1, rh1, N1)
        B[:N1] = rh1

        # --- Interface ODE ---
        M = (M + dt * u1 * B[i_b]) / (1.0 + dt * u2)

        # --- Domain 2 ---
        di2[0] = 1.0;  rh2[0] = M          # Dirichlet B[i_b+1] = M
        for j in range(1, N2):
            di2[j] = diag2
            rh2[j] = B[i_b + 1 + j] / dt + alpha * A[i_b + 1 + j]
        rh2 = thomas_solve(lo2, di2, up2, rh2, N2)
        B[i_b + 1:] = rh2

    return z, A, B, M
